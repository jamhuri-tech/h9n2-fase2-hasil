"""
Kelestarian residu kantong katalitik neuraminidase antara struktur yang dipakai
docking (4K1K, N2 dari H3N2) dan galur H9N2.

Identitas urutan menyeluruh menjawab "seberapa mirip proteinnya". Yang sebenarnya
menentukan sah-tidaknya substitusi struktur untuk docking adalah pertanyaan yang
lebih sempit: apakah residu yang MENYENTUH LIGAN sama persis?

Residu kontak diambil dari pose kristalografi oseltamivir karboksilat pada 4K1K
(<= 4,5 A), lalu dicek satu per satu pada penjajaran terhadap tiap urutan H9N2.
"""

import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
PDB = HERE.parent / "19_docking_targets" / "receptors" / "4K1K.pdb"
LIGAND = "G39"
CHAIN = "A"
CUTOFF = 4.5
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def parse_structure():
    """Urutan teramati (nomor residu -> kode satu huruf) dan koordinat per residu."""
    seq, coords, lig = {}, {}, []
    for line in PDB.read_text().splitlines():
        if line[16] not in (" ", "A"):
            continue
        resn, ch, num = line[17:20].strip(), line[21], line[22:27].strip()
        if line.startswith("ATOM") and ch == CHAIN and resn in THREE2ONE:
            key = int(re.sub(r"[^0-9-]", "", num))
            seq[key] = THREE2ONE[resn]
            coords.setdefault(key, []).append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elif line.startswith("HETATM") and resn == LIGAND and ch == CHAIN:
            lig.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return seq, {k: np.array(v) for k, v in coords.items()}, np.array(lig)


def contact_residues(coords, lig):
    out = []
    for num, xyz in coords.items():
        if np.linalg.norm(xyz[:, None, :] - lig[None, :, :], axis=2).min() <= CUTOFF:
            out.append(num)
    return sorted(out)


def fetch_h9n2(n=40, indonesia_only=False):
    term = ('H9N2[All Fields] AND neuraminidase[Protein Name] '
            'AND "influenza a virus"[Organism] AND 400:500[SLEN]')
    if indonesia_only:
        term += " AND Indonesia[All Fields]"
    url = (f"{EUTILS}/esearch.fcgi?db=protein&term={urllib.parse.quote(term)}"
           f"&retmax={n}&retmode=json")
    ids = json.loads(urllib.request.urlopen(url, timeout=40).read())
    ids = ids["esearchresult"]["idlist"]
    if not ids:
        return []
    url = (f"{EUTILS}/efetch.fcgi?db=protein&id={','.join(ids)}"
           f"&rettype=fasta&retmode=text")
    fasta = urllib.request.urlopen(url, timeout=90).read().decode()
    out = []
    for block in fasta.split(">"):
        if not block.strip():
            continue
        head, *rest = block.splitlines()
        s = "".join(rest).replace(" ", "").upper()
        if 400 <= len(s) <= 500:
            out.append((head.split()[0], head.strip(), s))
    return out


def main():
    seq_map, coords, lig = parse_structure()
    nums = sorted(seq_map)
    ref_seq = "".join(seq_map[n] for n in nums)
    contacts = contact_residues(coords, lig)
    idx_of = {n: i for i, n in enumerate(nums)}

    print(f"struktur 4K1K rantai {CHAIN}: {len(ref_seq)} residu teramati")
    print(f"residu kontak oseltamivir (<= {CUTOFF} A): {len(contacts)}")
    print("  " + ", ".join(f"{seq_map[n]}{n}" for n in contacts) + "\n")

    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"

    seqs = fetch_h9n2(40)
    ind = fetch_h9n2(20, indonesia_only=True)
    known = {a for a, _, _ in seqs}
    seqs += [s for s in ind if s[0] not in known]
    ind_acc = {a for a, _, _ in ind}
    print(f"urutan H9N2 diperiksa: {len(seqs)} (galur Indonesia: {len(ind)})\n")

    per_res = {n: 0 for n in contacts}
    rows = []
    for acc, head, s in seqs:
        aln = al.align(ref_seq, s)[0]
        a, b = aln[0], aln[1]
        # peta indeks referensi -> huruf pasangan
        ri, pair = -1, {}
        for x, y in zip(a, b):
            if x != "-":
                ri += 1
                pair[ri] = y
        same = 0
        for n in contacts:
            if pair.get(idx_of[n]) == seq_map[n]:
                same += 1
                per_res[n] += 1
        rows.append({"accession": acc, "indonesia": acc in ind_acc,
                     "kontak_identik": same, "dari": len(contacts),
                     "persen": round(100 * same / len(contacts), 1)})

    df = pd.DataFrame(rows).sort_values("persen", ascending=False)
    df.to_csv(HERE / "20.2_active_site_conservation.csv", index=False)

    res = pd.DataFrame([{"residu": f"{seq_map[n]}{n}", "lestari_pada": per_res[n],
                         "dari": len(seqs),
                         "persen": round(100 * per_res[n] / len(seqs), 1)}
                        for n in contacts]).sort_values("persen")
    res.to_csv(HERE / "20.3_contact_residue_conservation.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== kelestarian per residu kontak ===")
    print(res.to_string(index=False))
    print(f"\nresidu kontak lestari 100% pada semua galur: "
          f"{(res.persen == 100).sum()} dari {len(res)}")
    print(f"rata-rata kontak identik per galur: {df.persen.mean():.1f}%")
    if df.indonesia.any():
        g = df[df.indonesia]
        print(f"galur Indonesia ({len(g)}): rata-rata {g.persen.mean():.1f}%")
    print(f"\n-> {HERE / '20.2_active_site_conservation.csv'}")
    print(f"-> {HERE / '20.3_contact_residue_conservation.csv'}")


if __name__ == "__main__":
    main()
