"""
Validasi pose kandidat pada reseptor yang lolos kontrol.

Kontrol re-docking hanya membuktikan mesin docking mampu menemukan kembali pose
ligan PEMBANDING. Ia tidak membuktikan pose KANDIDAT kita masuk akal. Skrip ini
mengujinya dengan dua ukuran yang dapat dihitung, bukan dengan pengamatan mata:

1. TUMPANG-TINDIH RESIDU KONTAK. Residu protein yang bersentuhan (<= 4,0 A) dengan
   pose kandidat dibandingkan dengan residu yang bersentuhan dengan ligan
   ko-kristal. Pose yang menempati kantong yang sama akan berbagi sebagian besar
   residu; pose yang cuma menempel di permukaan tidak.

     recall  = bagian kontak ligan kristal yang berhasil ditiru kandidat
     jaccard = irisan / gabungan, menghukum pose yang melebar ke luar kantong

2. KOORDINASI Zn (khusus MMP9 dan MMP2). Penghambat MMP sejati mengunci Zn
   katalitik lewat atom O atau N pada jarak ~2,0-2,5 A. Pose berskor bagus yang
   Zn-nya 5 A lebih jauh berarti senyawa itu menempel di tepi kantong, bukan
   menghambat.

Fungsi penilai yang dipakai mengikuti hasil kontrol tiap reseptor: Vinardo untuk
neuraminidase, Vina untuk sisanya.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
RECEPTORS = HERE / "receptors"
WORK = HERE / "docking"
CONTACT_CUTOFF = 4.0
ZN_COORD_MAX = 2.6          # batas atas koordinasi Zn sejati

# reseptor tervalidasi -> (ligan kristal, fungsi penilai yang lolos kontrol, cek Zn?)
VALIDATED = {
    "4XCT": ("MMP9", "N73", "vina", True),
    "4A5S": ("DPP4", "N7F", "vina", False),
    "6B1K": ("MIF", "C9G", "vina", False),
    "4K1K": ("Neuraminidase N2", "G39", "vinardo", False),
    "8H78": ("MMP2", "L2U", "vina", True),      # sementara: situs lolos uji pengganti
}


def protein_atoms(pdb_id):
    """(label residu, koordinat) untuk seluruh atom berat protein."""
    labels, xyz = [], []
    for line in (RECEPTORS / f"{pdb_id}.pdb").read_text().splitlines():
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        if (line[76:78].strip() or "C").upper() == "H":
            continue
        labels.append(f"{line[21]}:{line[17:20].strip()}{line[22:27].strip()}")
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(labels), np.array(xyz)


def contacts(lig_xyz, labels, prot_xyz, cutoff=CONTACT_CUTOFF):
    d = np.linalg.norm(prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    return set(labels[(d.min(axis=1) <= cutoff)])


def zinc_sites(pdb_id):
    out = []
    for line in (RECEPTORS / f"{pdb_id}.pdb").read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "ZN":
            out.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return out


def pose_coords(pdbqt, want_polar_only=False):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    conf = mol.GetConformer(0)          # pose peringkat 1
    idx = [a.GetIdx() for a in mol.GetAtoms()
           if (not want_polar_only) or a.GetSymbol() in ("O", "N")]
    return np.array([list(conf.GetAtomPosition(i)) for i in idx])


def main():
    rows = []
    for pdb_id, (target, comp, scoring, check_zn) in VALIDATED.items():
        tag = "" if scoring == "vina" else f"_{scoring}"
        labels, prot = protein_atoms(pdb_id)
        zns = zinc_sites(pdb_id)

        native = Chem.MolFromMolFile(str(WORK / pdb_id / f"native_{comp}.sdf"), removeHs=True)
        nat_xyz = native.GetConformer().GetPositions()
        nat_contacts = contacts(nat_xyz, labels, prot)

        for lig in sorted((HERE / "ligands").glob("*.pdbqt")):
            pose = WORK / pdb_id / f"{lig.stem}_seed1{tag}.pdbqt"
            if not pose.exists():
                continue
            xyz = pose_coords(pose)
            cs = contacts(xyz, labels, prot)
            inter = nat_contacts & cs
            recall = len(inter) / len(nat_contacts) if nat_contacts else None
            jacc = len(inter) / len(nat_contacts | cs) if (nat_contacts | cs) else None

            zn_d = None
            if check_zn and zns:
                polar = pose_coords(pose, want_polar_only=True)
                if len(polar):
                    zn_d = round(float(min(np.linalg.norm(polar - z, axis=1).min()
                                           for z in zns)), 2)

            rows.append({
                "reseptor": target, "pdb_id": pdb_id, "scoring": scoring,
                "ligan": lig.stem.split("_", 1)[1].replace("_", " "),
                "n_kontak_kristal": len(nat_contacts), "n_kontak_pose": len(cs),
                "n_bersama": len(inter),
                "recall": round(recall, 2), "jaccard": round(jacc, 2),
                "min_ON_ke_Zn_A": zn_d,
                "koordinasi_Zn": None if zn_d is None else zn_d <= ZN_COORD_MAX,
            })

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "19.10_pose_validation.csv", index=False)

    pd.set_option("display.width", 250)
    for target in df.reseptor.unique():
        sub = df[df.reseptor == target].sort_values("recall", ascending=False)
        cols = ["ligan", "n_kontak_pose", "n_bersama", "recall", "jaccard"]
        if sub.min_ON_ke_Zn_A.notna().any():
            cols += ["min_ON_ke_Zn_A", "koordinasi_Zn"]
        print(f"\n=== {target} ({sub.pdb_id.iloc[0]}, {sub.scoring.iloc[0]}) — "
              f"{sub.n_kontak_kristal.iloc[0]} residu kontak pada ligan kristal ===")
        print(sub[cols].to_string(index=False))
    print(f"\n-> {HERE / '19.10_pose_validation.csv'}")


if __name__ == "__main__":
    main()
