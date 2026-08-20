"""
Kontrol pengganti untuk MMP2, karena kontrol re-docking sebenarnya tidak mungkin.

Masalahnya
----------
Seluruh PDB hanya memuat dua struktur MMP2 dengan penghambat: 8H78 (L2U, 64 atom
berat -- jauh di luar jangkauan andal Vina) dan 7XJO (lima mutasi rekayasa, dan
"penghambat"-nya ternyata buffer bis-tris propana). Tidak ada ligan ko-kristal
berukuran wajar yang bisa dipakai untuk re-docking.

Penggantinya
------------
ARP101 (kode CCD N73, 28 atom berat) adalah penghambat hidroksamat MMP yang pada
4XCT terikat di situs katalitik MMP9 -- dan MMP9 sudah LULUS kontrol re-docking
(RMSD pose-1 1,21 A). ARP101 didok ke MMP2, lalu diperiksa apakah ia mengunci Zn
katalitik seperti yang seharusnya dilakukan penghambat MMP.

Yang diukur: jarak terpendek atom O/N ligan ke Zn katalitik.
  ~2,0-2,5 A  = koordinasi Zn sejati, setup MMP2 berperilaku benar
  > 4 A       = pose menempel di permukaan, hasil MMP2 tidak dapat dipercaya

Nilai acuan dihitung dari pose kristal ARP101 pada 4XCT.

Ini BUKAN pengganti setara kontrol re-docking: ia tidak menguji ketepatan
geometri pose, hanya apakah situs dan logamnya dikenali dengan benar.
"""

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
VINA = HERE / "bin" / "vina"
RECEPTORS = HERE / "receptors"
WORK = HERE / "docking"


def catalytic_zn(pdb_path, ref_xyz):
    """Zn yang paling dekat ke titik acuan -- yang katalitik, bukan yang struktural."""
    best = None
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "ZN":
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            d = np.linalg.norm(xyz - ref_xyz)
            if best is None or d < best[1]:
                best = (xyz, d)
    return None if best is None else best[0]


def polar_atoms(mol, conf_id=0):
    conf = mol.GetConformer(conf_id)
    return np.array([list(conf.GetAtomPosition(a.GetIdx()))
                     for a in mol.GetAtoms() if a.GetSymbol() in ("O", "N")])


def min_dist_to_zn(mol, zn, conf_id=0):
    p = polar_atoms(mol, conf_id)
    return float(np.linalg.norm(p - zn, axis=1).min()) if len(p) else None


def main():
    boxes = pd.read_csv(HERE / "19.2_receptor_shortlist.csv").set_index("pdb_id")
    rows = []

    # --- acuan: pose kristal ARP101 pada MMP9 ---
    native = Chem.MolFromMolFile(str(WORK / "4XCT" / "native_N73.sdf"), removeHs=True)
    zn9 = catalytic_zn(RECEPTORS / "4XCT.pdb", polar_atoms(native).mean(0))
    rows.append({"reseptor": "MMP9 (4XCT)", "sumber_pose": "kristal",
                 "min_O/N_ke_Zn_A": round(min_dist_to_zn(native, zn9), 2),
                 "skor": None})

    # --- ARP101 didok ulang ke MMP9 (kalibrasi) dan ke MMP2 (uji) ---
    lig = WORK / "4XCT" / "native_N73.pdbqt"
    for pdb_id, label in [("4XCT", "MMP9 (4XCT)"), ("8H78", "MMP2 (8H78)")]:
        box = boxes.loc[pdb_id]
        out = WORK / pdb_id / "surrogate_N73.pdbqt"
        if not out.exists():
            subprocess.run([str(VINA),
                "--receptor", str(WORK / pdb_id / f"{pdb_id}.pdbqt"), "--ligand", str(lig),
                "--center_x", str(box.center_x), "--center_y", str(box.center_y),
                "--center_z", str(box.center_z),
                "--size_x", str(box.size_x), "--size_y", str(box.size_y),
                "--size_z", str(box.size_z),
                "--exhaustiveness", "32", "--seed", "1", "--out", str(out)],
                capture_output=True, text=True)

        from meeko import PDBQTMolecule, RDKitMolCreate
        pm = PDBQTMolecule.from_file(str(out), skip_typing=True)
        docked = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
        score = float(re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)",
                                out.read_text()).group(1))
        zn = catalytic_zn(RECEPTORS / f"{pdb_id}.pdb", polar_atoms(docked).mean(0))
        rows.append({"reseptor": label, "sumber_pose": "docking ARP101",
                     "min_O/N_ke_Zn_A": round(min_dist_to_zn(docked, zn), 2),
                     "skor": score})

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "19.9_mmp2_surrogate_control.csv", index=False)
    print(df.to_string(index=False))
    print(f"\n-> {HERE / '19.9_mmp2_surrogate_control.csv'}")


if __name__ == "__main__":
    main()
