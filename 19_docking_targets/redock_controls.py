"""
Kontrol re-docking: memvalidasi protokol pada tiap reseptor.

Ligan ko-kristal dikeluarkan dari struktur, didok ulang ke kotaknya sendiri,
lalu RMSD pose hasil docking dibandingkan dengan pose kristalografinya.
RMSD <= 2,0 A = protokol tervalidasi untuk reseptor tersebut.

Reseptor yang GAGAL kontrol tidak boleh dipakai menilai kandidat: kalau mesin
docking saja tidak mampu menemukan kembali pose yang sudah diketahui benar,
skornya untuk senyawa yang belum diketahui tidak dapat dipercaya.

Kenapa perlu template CCD
-------------------------
Fragmen HETATM dari berkas PDB tidak memuat orde ikatan, sehingga meeko menolak
memprosesnya. Orde ikatan diambil dari berkas "ideal" Chemical Component
Dictionary di RCSB, lalu ditempelkan ke koordinat kristalografi dengan
AssignBondOrdersFromTemplate.
"""

import argparse
import re
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
RECEPTORS = HERE / "receptors"
WORK = HERE / "docking"
CACHE = HERE / "ccd_templates"
CACHE.mkdir(exist_ok=True)
VINA = HERE / "bin" / "vina"
EXHAUSTIVENESS = 32
SEED = 1
SCORING = "vina"      # ditimpa lewat --scoring
RMSD_PASS = 2.0

# heteroatom yang bukan ligan sesungguhnya -- jangan dijadikan kontrol
NOT_A_LIGAND = {"PGO", "BUD", "2HP", "EDO", "GOL", "DMS", "PEG", "SO4", "PO4", "NAG"}


def ccd_template(comp_id):
    path = CACHE / f"{comp_id}_ideal.sdf"
    if not path.exists():
        url = f"https://files.rcsb.org/ligands/download/{comp_id}_ideal.sdf"
        path.write_bytes(urllib.request.urlopen(url, timeout=45).read())
    return Chem.MolFromMolFile(str(path), sanitize=True, removeHs=True)


def native_from_pdb(pdb_path, comp_id, chain):
    """Salinan pertama ligan comp_id pada rantai tertentu, dengan orde ikatan benar."""
    lines, resseq = [], None
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("HETATM") or line[16] not in (" ", "A"):
            continue
        if line[17:20].strip() != comp_id or line[21] != chain:
            continue
        if resseq is None:
            resseq = line[22:27]
        if line[22:27] != resseq:
            continue
        # buang hidrogen: sebagian entri beresolusi tinggi (mis. 8T5W) menyimpan
        # H eksplisit, sedangkan template CCD dibaca tanpa H -- jumlah atomnya
        # harus cocok agar AssignBondOrdersFromTemplate berhasil
        if (line[76:78].strip() or "C").upper() == "H":
            continue
        lines.append(line)
    if not lines:
        return None
    block = "\n".join(lines) + "\nEND\n"
    raw = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
    if raw is None:
        return None
    tpl = ccd_template(comp_id)
    if tpl is None:
        return None
    try:
        return AllChem.AssignBondOrdersFromTemplate(tpl, raw)
    except ValueError:
        return None


def rmsd_to_crystal(native, docked_pdbqt):
    """RMSD sadar-simetri antara pose kristal dan tiap pose hasil docking.

    Molekul hasil docking DIREKONSTRUKSI dari PDBQT lewat meeko, bukan dibaca
    baris demi baris: meeko menyusun ulang urutan atom saat menulis PDBQT,
    sehingga pembandingan indeks-per-indeks membandingkan atom yang berlainan
    dan menghasilkan RMSD palsu 4-8 A untuk pose yang sebenarnya tepat.

    Dikembalikan dua angka:
      rmsd_top  -- pose peringkat 1; menguji FUNGSI SKOR
      rmsd_best -- pose terbaik di antara seluruh pose; menguji PENCARIAN
    Keduanya perlu: pencarian bisa menemukan pose benar tetapi skor gagal
    menaruhnya di peringkat satu.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit.Chem import rdMolAlign

    pm = PDBQTMolecule.from_file(str(docked_pdbqt), skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    if not mols or mols[0] is None:
        return None, None, None
    docked = Chem.RemoveHs(mols[0])
    ref = Chem.RemoveHs(native)
    if docked.GetNumAtoms() != ref.GetNumAtoms():
        return None, None, None

    rms = [rdMolAlign.CalcRMS(docked, ref, prbId=i)
           for i in range(docked.GetNumConformers())]
    best_i = min(range(len(rms)), key=lambda i: rms[i])
    return round(rms[0], 2), round(rms[best_i], 2), best_i + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo", "ad4"])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    tag = "" if args.scoring == "vina" else f"_{args.scoring}"

    boxes = pd.read_csv(HERE / "19.2_receptor_shortlist.csv")
    rows = []
    for box in boxes.itertuples():
        if not isinstance(box.ligand_ids, str):
            print(f"[{box.pdb_id}] {box.target}: tanpa ligan ko-kristal nonpolimer -- "
                  f"kontrol tidak dapat dijalankan")
            rows.append({"pdb_id": box.pdb_id, "target": box.target, "native_ligand": None,
                         "redock_score": None, "rmsd_top_pose_A": None,
                         "rmsd_best_pose_A": None, "control_passed_scoring": None,
                         "control_passed_sampling": None,
                         "note": "situs ditandai glikan (entitas polimer), bukan ligan nonpolimer"})
            continue

        comps = [c for c in box.ligand_ids.split() if c not in NOT_A_LIGAND]
        if not comps:
            continue
        comp = comps[0]
        rdir = WORK / box.pdb_id
        rdir.mkdir(parents=True, exist_ok=True)

        native = native_from_pdb(RECEPTORS / f"{box.pdb_id}.pdb", comp, box.site_chain)
        if native is None:
            print(f"[{box.pdb_id}] {comp}: gagal merekonstruksi orde ikatan")
            continue

        # meeko menolak molekul dengan H implisit -- tambahkan H beserta koordinatnya
        native_h = Chem.AddHs(native, addCoords=True)
        sdf = rdir / f"native_{comp}.sdf"
        w = Chem.SDWriter(str(sdf)); w.write(native_h); w.close()

        nat_q = rdir / f"native_{comp}.pdbqt"
        subprocess.run(["mk_prepare_ligand.py", "-i", str(sdf), "-o", str(nat_q)],
                       capture_output=True, text=True)
        if not nat_q.exists():
            print(f"[{box.pdb_id}] {comp}: meeko gagal menyiapkan ligan")
            continue

        out = rdir / f"native_{comp}_redock{tag}.pdbqt"
        r = None
        if not out.exists():
            r = subprocess.run([str(VINA),
            "--receptor", str(rdir / f"{box.pdb_id}.pdbqt"), "--ligand", str(nat_q),
            "--center_x", str(box.center_x), "--center_y", str(box.center_y),
            "--center_z", str(box.center_z),
            "--size_x", str(box.size_x), "--size_y", str(box.size_y),
            "--size_z", str(box.size_z),
            "--scoring", args.scoring,
            "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", str(SEED),
            "--out", str(out)], capture_output=True, text=True)
        if r is not None:
            scores = [float(m) for m in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", r.stdout, re.M)]
            score = scores[0] if scores else None
        else:
            m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", out.read_text())
            score = float(m.group(1)) if m else None

        rmsd_top, rmsd_best, best_rank = rmsd_to_crystal(native, out)
        passed_top = rmsd_top is not None and rmsd_top <= RMSD_PASS
        passed_sampling = rmsd_best is not None and rmsd_best <= RMSD_PASS
        rows.append({"pdb_id": box.pdb_id, "target": box.target, "native_ligand": comp,
                     "redock_score": score, "rmsd_top_pose_A": rmsd_top,
                     "rmsd_best_pose_A": rmsd_best, "best_pose_rank": best_rank,
                     "control_passed_scoring": passed_top,
                     "control_passed_sampling": passed_sampling,
                     "n_heavy_atoms": Chem.RemoveHs(native).GetNumAtoms(),
                     "scoring": args.scoring})
        print(f"[{box.pdb_id}] {box.target:18s} {comp:5s} skor={score:>8}  "
              f"RMSD pose-1={rmsd_top:>5}  terbaik={rmsd_best:>5} (pose {best_rank})  "
              f"skor:{'LULUS' if passed_top else 'GAGAL'}  "
              f"cari:{'LULUS' if passed_sampling else 'GAGAL'}")

    df = pd.DataFrame(rows)
    dest = HERE / (args.out or f"19.5_redock_controls{tag}.csv")
    df.to_csv(dest, index=False)
    print(f"\n-> {dest}")


if __name__ == "__main__":
    main()
