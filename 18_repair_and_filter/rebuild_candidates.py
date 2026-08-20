"""
Membangun ulang peringkat kandidat NDD dengan identitas senyawa yang benar.

Prediksi NDD (probabilitas) sudah benar: fingerprint dihitung dari kolom SMILES
pada berkas 8_FP_*, dan prediksi mengikuti fingerprint tersebut -- terverifikasi
0 inkonsistensi fingerprint di dalam SMILES yang sama. Yang salah hanya label
identitas. Skrip ini membuang label lama dan menempelkan identitas dari
18.1_identity_map_by_canonical_smiles.csv.

Deduplikasi dilakukan ke tingkat struktur: 623 baris konsensus hanya memuat 492
struktur unik, dan daftar 115 "kandidat" hanya memuat 104 struktur unik.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent

CONSENSUS = DATA / "9_ndd_predictions" / "herbal_consensus_across_fingerprints.xlsx"
IDMAP = HERE / "18.1_identity_map_by_canonical_smiles.csv"

PROB_COLS = ["prob_mean_EState", "prob_mean_Klekota-Roth", "prob_mean_PubChem"]
THRESHOLD = 0.9


def main():
    cons = pd.read_excel(CONSENSUS)
    idmap = pd.read_csv(IDMAP)

    # buang seluruh kolom identitas lama yang tidak sebaris dengan struktur
    cons = cons.drop(columns=["IUPAC_Name", "Metabolite", "Organism", "PubChem_ID"])

    # satu baris per struktur; probabilitas identik untuk SMILES yang sama
    per_structure = cons.drop_duplicates(subset="Smiles").copy()
    spread = cons.groupby("Smiles")["prob_mean_across_fingerprints"].agg(lambda x: x.max() - x.min())
    assert spread.max() < 1e-6, f"probabilitas berbeda untuk SMILES yang sama: {spread.max()}"

    # berkas 7 tidak 100% kanonik menurut RDKit (1 baris berbeda), jadi kanonikkan
    # kedua sisi sebelum join
    per_structure["can"] = [
        Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), isomericSmiles=False)
        for s in per_structure["Smiles"]
    ]
    merged = per_structure.merge(idmap, on="can", how="left")
    missing = merged["Metabolite"].isna().sum()
    if missing:
        print(f"PERINGATAN: {missing} struktur tidak ditemukan di peta identitas")

    merged = merged.drop(columns=["can"])
    merged = merged.sort_values("prob_mean_across_fingerprints", ascending=False).reset_index(drop=True)
    merged.insert(0, "rank", range(1, len(merged) + 1))
    merged.to_csv(HERE / "18.2_ndd_ranking_all_structures_corrected.csv", index=False)

    passes = (merged[PROB_COLS] > THRESHOLD).all(axis=1)
    cand = merged[passes].copy()
    cand["rank"] = range(1, len(cand) + 1)
    cand.to_csv(HERE / "18.3_candidates_prob_gt_0.9_corrected.csv", index=False)

    cols = ["rank", "Metabolite", "Organism", "PubChem_ID", "Molecular_formula",
            "prob_mean_across_fingerprints", "identity_ambiguous"]
    print(f"struktur unik dalam konsensus        : {len(merged)}")
    print(f"struktur lolos ketiga FP > {THRESHOLD}      : {len(cand)}")
    print(f"  di antaranya identitas ambigu      : {int(cand.identity_ambiguous.sum())}")
    print()
    print("10 kandidat teratas (identitas TERKOREKSI):")
    print(cand[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
