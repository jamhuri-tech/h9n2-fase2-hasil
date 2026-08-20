"""
Daftar residu kontak bernama untuk kandidat yang lolos validasi pose.

Untuk naskah, tabel pose perlu menyebut residu spesifik ("Zn, Glu402, His401")
bukan sekadar angka recall -- itulah yang dapat dibandingkan pembaca dengan
literatur penghambat yang sudah dikenal.

Dilaporkan tiga hal per pasangan: residu yang juga disentuh ligan ko-kristal
(kontak bersama), residu yang hanya disentuh kandidat, dan residu ligan kristal
yang tidak tersentuh.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DOCK = HERE.parent / "19_docking_targets"
RECEPTORS = DOCK / "receptors"
WORK = DOCK / "docking"
CUTOFF = 4.0

THREE2ONE = {
    "ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys",
    "GLN": "Gln", "GLU": "Glu", "GLY": "Gly", "HIS": "His", "ILE": "Ile",
    "LEU": "Leu", "LYS": "Lys", "MET": "Met", "PHE": "Phe", "PRO": "Pro",
    "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val",
}

TARGETS = {"4XCT": ("MMP9", "N73", ""), "4A5S": ("DPP4", "N7F", ""),
           "6B1K": ("MIF", "C9G", ""), "4K1K": ("Neuraminidase N2", "G39", "_vinardo"),
           "8H78": ("MMP2", "L2U", "")}


def receptor_atoms(pdb_id):
    lab, xyz = [], []
    for line in (RECEPTORS / f"{pdb_id}.pdb").read_text().splitlines():
        if line[16] not in (" ", "A"):
            continue
        resn = line[17:20].strip()
        if line.startswith("ATOM") and resn in THREE2ONE:
            name = f"{THREE2ONE[resn]}{line[22:27].strip()}"
        elif line.startswith("HETATM") and resn in ("ZN", "CA", "MN", "MG"):
            name = resn.capitalize() if len(resn) > 1 else resn
        else:
            continue
        if (line[76:78].strip() or "C").upper() == "H":
            continue
        lab.append(f"{name}({line[21]})")
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(lab), np.array(xyz)


def touched(lig_xyz, lab, xyz):
    d = np.linalg.norm(xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    return sorted(set(lab[d.min(axis=1) <= CUTOFF]))


def pose_xyz(pdbqt):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    return mol.GetConformer(0).GetPositions()


def main():
    verdict = pd.read_csv(DOCK / "19.11_verdict.csv")
    ok = verdict[verdict.pose_ok]
    rows = []
    for pdb_id, (target, comp, tag) in TARGETS.items():
        lab, xyz = receptor_atoms(pdb_id)
        native = Chem.MolFromMolFile(str(WORK / pdb_id / f"native_{comp}.sdf"), removeHs=True)
        nat = set(touched(native.GetConformer().GetPositions(), lab, xyz))

        for r in ok[ok.pdb_id == pdb_id].itertuples():
            stem = [f.stem for f in (DOCK / "ligands").glob("*.pdbqt")
                    if f.stem.split("_", 1)[1].replace("_", " ") == r.ligan]
            if not stem:
                continue
            pose = WORK / pdb_id / f"{stem[0]}_seed1{tag}.pdbqt"
            if not pose.exists():
                continue
            cs = set(touched(pose_xyz(pose), lab, xyz))
            rows.append({
                "reseptor": target, "pdb_id": pdb_id, "ligan": r.ligan,
                "dG": r.affinity_mean, "recall": r.recall,
                "kontak_bersama": ", ".join(sorted(nat & cs)),
                "hanya_kandidat": ", ".join(sorted(cs - nat)),
                "tidak_tersentuh": ", ".join(sorted(nat - cs)),
            })

    df = pd.DataFrame(rows).sort_values(["reseptor", "dG"])
    df.to_csv(HERE / "20.7_contact_residues_named.csv", index=False)

    for r in df.itertuples():
        print(f"\n=== {r.reseptor} — {r.ligan}  (dG {r.dG}, recall {r.recall}) ===")
        print(f"  kontak bersama ligan kristal : {r.kontak_bersama}")
        if r.hanya_kandidat:
            print(f"  hanya disentuh kandidat      : {r.hanya_kandidat}")
        if r.tidak_tersentuh:
            print(f"  kontak kristal yang terlewat : {r.tidak_tersentuh}")
    print(f"\n-> {HERE / '20.7_contact_residues_named.csv'}")


if __name__ == "__main__":
    main()
