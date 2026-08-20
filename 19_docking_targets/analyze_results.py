"""
Menyusun hasil akhir docking: afinitas, efisiensi ligan, dan status kontrol.

Dua koreksi yang wajib ada sebelum hasil dibaca:

1. KONTROL RE-DOCKING. Reseptor yang tidak mampu menemukan kembali pose ligan
   kristalnya sendiri (RMSD > 2 A) tidak dapat dipercaya untuk menilai senyawa
   yang belum diketahui. Statusnya dibawa ke tiap baris hasil.

2. EFISIENSI LIGAN, LE = -dG / jumlah atom berat. Fungsi skor Vina bertambah
   hampir linear terhadap ukuran molekul, sehingga peringkat berdasarkan dG
   mentah otomatis memenangkan molekul besar. LE membuang bias itu dan
   menjawab pertanyaan yang sebenarnya: seberapa produktif tiap atom.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
THRESHOLD = -6.5          # ambang proposal
LE_GOOD = 0.35            # LE di atas ini lazim dianggap titik awal yang sehat

SHORT = {"1JSI": "HA (LSTc)", "1JSH": "HA (LSTa)", "4K1K": "NA-N2",
         "8T5W": "PA-endo", "4DYN": "NP", "4XCT": "MMP9", "8H78": "MMP2",
         "4A5S": "DPP4", "6B1K": "MIF"}


def heavy_atom_counts():
    out = {}
    for f in (HERE / "ligands").glob("*.sdf"):
        m = Chem.MolFromMolFile(str(f))
        if m is not None:
            out[f.stem] = Chem.RemoveHs(m).GetNumAtoms()
    return out


def main():
    d = pd.read_csv(HERE / "19.7_docking_summary.csv")
    ha = heavy_atom_counts()
    d["heavy_atoms"] = d.ligand.map(ha)
    d["ligand_efficiency"] = (-d.affinity_mean / d.heavy_atoms).round(3)
    d["ligan"] = d.ligand.str.replace(r"^\d+_", "", regex=True).str.replace("_", " ")
    d["reseptor"] = d.pdb_id.map(SHORT)

    ctrl_path = HERE / "19.5_redock_controls.csv"
    if ctrl_path.exists():
        c = pd.read_csv(ctrl_path)
        d = d.merge(c[["pdb_id", "native_ligand", "redock_score", "rmsd_top_pose_A",
                       "rmsd_best_pose_A", "control_passed_scoring",
                       "control_passed_sampling"]], on="pdb_id", how="left")

        def verdict(r):
            if pd.isna(r.control_passed_sampling):
                return "tanpa kontrol"
            if r.control_passed_scoring:
                return "tervalidasi"
            if r.control_passed_sampling:
                return "pencarian saja"
            return "tidak tervalidasi"
        d["status_kontrol"] = d.apply(verdict, axis=1)
        # skor relatif terhadap kontrol positif reseptor yang sama
        d["vs_control_kcal"] = (d.affinity_mean - d.redock_score).round(2)

    d["meets_threshold"] = d.affinity_mean <= THRESHOLD
    d["good_efficiency"] = d.ligand_efficiency >= LE_GOOD
    d = d.sort_values(["group", "reseptor", "affinity_mean"])
    d.to_csv(HERE / "19.8_docking_final.csv", index=False)

    pd.set_option("display.width", 250)
    order = [c for c in SHORT.values()]

    print("=" * 78)
    print("AFINITAS rata-rata 3 seed (kcal/mol) -- makin negatif makin kuat")
    print("=" * 78)
    print(d.pivot_table(index="ligan", columns="reseptor",
                        values="affinity_mean")[order].round(2).to_string())

    print("\n" + "=" * 78)
    print("EFISIENSI LIGAN  LE = -dG / atom berat")
    print("=" * 78)
    print(d.pivot_table(index="ligan", columns="reseptor",
                        values="ligand_efficiency")[order].round(3).to_string())

    if "status_kontrol" in d:
        print("\n" + "=" * 78)
        print("KONTROL RE-DOCKING")
        print("=" * 78)
        c = d.drop_duplicates("pdb_id")[["reseptor", "pdb_id", "native_ligand",
                                         "redock_score", "rmsd_top_pose_A",
                                         "rmsd_best_pose_A", "status_kontrol"]]
        print(c.to_string(index=False))

        print("\n" + "=" * 78)
        print("AFINITAS pada reseptor TERVALIDASI saja")
        print("=" * 78)
        ok = d[d.status_kontrol == "tervalidasi"]
        cols = [c for c in order if c in set(ok.reseptor)]
        print(ok.pivot_table(index="ligan", columns="reseptor",
                             values="affinity_mean")[cols].round(2).to_string())
        print("\nefisiensi ligan pada reseptor yang sama:")
        print(ok.pivot_table(index="ligan", columns="reseptor",
                             values="ligand_efficiency")[cols].round(3).to_string())

    print(f"\n-> {HERE / '19.8_docking_final.csv'}")


if __name__ == "__main__":
    main()
