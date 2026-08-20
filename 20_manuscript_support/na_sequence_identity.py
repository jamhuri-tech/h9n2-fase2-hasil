"""
Kesamaan urutan neuraminidase: struktur yang dipakai docking (4K1K, N2 dari H3N2)
terhadap neuraminidase galur H9N2 -- termasuk galur Indonesia bila tersedia.

Kenapa ini perlu
----------------
Protein Data Bank tidak memuat satu pun struktur neuraminidase dari galur H9N2,
sehingga docking terpaksa memakai struktur N2 dari galur H3N2. Subtipe
neuraminidasenya tetap N2, tetapi galurnya berbeda -- dan itu titik yang pasti
ditanyakan reviewer. Angka kesamaan urutan mengubah pembelaan kualitatif
("N2 tetap N2") menjadi dasar kuantitatif.

Yang dilaporkan: persen identitas dan persen kemiripan pada penjajaran global
(Needleman-Wunsch, matriks BLOSUM62), terhadap sejumlah galur H9N2 dari NCBI.
"""

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = Path(__file__).resolve().parent
STRUCT_PDB = "4K1K"
N_SEQ = 40                      # jumlah urutan H9N2 yang diambil
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def structure_sequence(pdb_id):
    q = ('{ entry(entry_id:"%s"){ polymer_entities { '
         'entity_poly { pdbx_seq_one_letter_code_can } } } }' % pdb_id)
    req = urllib.request.Request("https://data.rcsb.org/graphql",
                                 data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    d = json.loads(urllib.request.urlopen(req, timeout=40).read())
    ents = d["data"]["entry"]["polymer_entities"]
    seqs = [e["entity_poly"]["pdbx_seq_one_letter_code_can"] for e in ents]
    return max(seqs, key=len)          # rantai neuraminidase = yang terpanjang


def fetch_h9n2_neuraminidase(n=N_SEQ):
    term = ('H9N2[All Fields] AND neuraminidase[Protein Name] '
            'AND "influenza a virus"[Organism] AND 400:500[SLEN]')
    url = (f"{EUTILS}/esearch.fcgi?db=protein&term={urllib.parse.quote(term)}"
           f"&retmax={n}&retmode=json")
    ids = json.loads(urllib.request.urlopen(url, timeout=40).read())
    ids = ids["esearchresult"]["idlist"]
    time.sleep(0.4)
    url = (f"{EUTILS}/efetch.fcgi?db=protein&id={','.join(ids)}"
           f"&rettype=fasta&retmode=text")
    fasta = urllib.request.urlopen(url, timeout=90).read().decode()

    out = []
    for block in fasta.split(">"):
        if not block.strip():
            continue
        head, *rest = block.splitlines()
        seq = "".join(rest).replace(" ", "").upper()
        if 400 <= len(seq) <= 500 and set(seq) <= set("ACDEFGHIKLMNPQRSTVWYX"):
            out.append((head.strip(), seq))
    return out


def make_aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "global"
    return al


def identity(aln, matrix):
    a, b = aln[0], aln[1]
    ident = sim = cols = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        cols += 1
        if x == y:
            ident += 1
            sim += 1
        elif x in matrix.alphabet and y in matrix.alphabet and matrix[x, y] > 0:
            sim += 1
    return 100 * ident / cols, 100 * sim / cols, cols


def strain_of(header):
    m = re.search(r"\(([^()]*H9N2[^()]*)\)", header)
    if m:
        return m.group(1)
    m = re.search(r"A/[A-Za-z0-9_/\-. ]+", header)
    return m.group(0) if m else header[:60]


def main():
    ref = structure_sequence(STRUCT_PDB)
    print(f"struktur {STRUCT_PDB}: {len(ref)} residu\n")

    seqs = fetch_h9n2_neuraminidase()
    print(f"urutan neuraminidase H9N2 diambil: {len(seqs)}\n")

    al = make_aligner()
    matrix = al.substitution_matrix
    rows = []
    for head, seq in seqs:
        aln = al.align(ref, seq)[0]
        pid, psim, cols = identity(aln, matrix)
        acc = head.split()[0]
        rows.append({
            "accession": acc, "galur": strain_of(head), "panjang": len(seq),
            "identitas_%": round(pid, 1), "kemiripan_%": round(psim, 1),
            "kolom_dibandingkan": cols,
            "indonesia": bool(re.search(r"indonesia|jakarta|sulawesi|java|bali|"
                                        r"sumatra|makassar", head, re.I)),
        })

    df = pd.DataFrame(rows).sort_values("identitas_%", ascending=False)
    df.to_csv(HERE / "20.1_na_sequence_identity.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== 10 galur H9N2 dengan kesamaan tertinggi ===")
    print(df.head(10).to_string(index=False))
    print("\n=== ringkasan ===")
    print(f"identitas : {df['identitas_%'].min():.1f}–{df['identitas_%'].max():.1f}% "
          f"(median {df['identitas_%'].median():.1f}%)")
    print(f"kemiripan : {df['kemiripan_%'].min():.1f}–{df['kemiripan_%'].max():.1f}% "
          f"(median {df['kemiripan_%'].median():.1f}%)")
    ind = df[df.indonesia]
    if len(ind):
        print(f"\ngalur Indonesia ({len(ind)}): identitas median "
              f"{ind['identitas_%'].median():.1f}%")
        print(ind.head(5).to_string(index=False))
    else:
        print("\ncatatan: tidak ada galur Indonesia pada cuplikan ini")
    print(f"\n-> {HERE / '20.1_na_sequence_identity.csv'}")


if __name__ == "__main__":
    main()
