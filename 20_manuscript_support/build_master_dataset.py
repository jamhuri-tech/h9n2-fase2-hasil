"""
Membangun satu berkas induk terkoreksi yang menggantikan rantai berkas 7-14.

Alasan dibuat satu berkas, bukan meniru rantai lamanya: rantai itu rapuh justru
karena panjang -- identitas berpindah tangan tujuh kali, dan pada satu titik
keselarasannya putus tanpa ketahuan. Berkas induk ini memuat semuanya sekaligus,
dengan kunci tunggal berupa SMILES kanonik.

Isi: identitas terverifikasi (rumus molekul, PubChem), SMILES isomerik untuk
docking, probabilitas ketiga sidik jari, konsensus, penanda jangkauan
keberlakuan model, dan penanda ambiguitas stereokimia.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
REPAIR = DATA / "18_repair_and_filter"
PRED = DATA.parent / "ndd_predictions"

FP_FILES = {"EState": "herbal_predictions_summary_EState.csv",
            "PubChem": "herbal_predictions_summary_PubChem.csv",
            "Klekota-Roth": "herbal_predictions_summary_Klekota_Roth.csv"}


def flat(s):
    m = Chem.MolFromSmiles(str(s))
    return None if m is None else Chem.MolToSmiles(m, isomericSmiles=False)


def main():
    idmap = pd.read_csv(REPAIR / "18.1_identity_map_by_canonical_smiles.csv")
    rank = pd.read_csv(REPAIR / "18.2_ndd_ranking_all_structures_corrected.csv")

    master = rank.merge(
        idmap[["can", "Smiles_isomeric", "Molecular_formula", "Mw",
               "identity_ambiguous", "Metabolite_all", "PubChem_ID_all",
               "Organism_all"]].rename(columns={"Smiles_isomeric": "smiles_iso_map"}),
        left_on="Smiles", right_on="can", how="left", suffixes=("", "_dup"))
    master = master.loc[:, ~master.columns.str.endswith("_dup")]

    # Kunci join HARUS dikanonikkan ulang: satu struktur pada berkas asal tidak
    # persis kanonik menurut RDKit, sehingga pencocokan berbasis string mentah
    # meleset dan memunculkan nilai kosong yang terbaca seolah "di luar jangkauan".
    master["can"] = [flat(x) for x in master["Smiles"]]

    # penanda jangkauan keberlakuan model, per sidik jari
    for fp, fname in FP_FILES.items():
        # pemisah kolom berbeda antar salinan berkas (";" vs ","), jadi dideteksi
        s = pd.read_csv(PRED / fname, sep=None, engine="python")
        s["can"] = [flat(x) for x in s["Smiles"]]
        s = s.drop_duplicates("can")
        cols = {"max_tanimoto_to_train_any": f"tanimoto_{fp}",
                "in_applicability_domain": f"in_AD_{fp}"}
        master = master.merge(s[["can"] + list(cols)].rename(columns=cols),
                              on="can", how="left")

    ad_cols = [c for c in master.columns if c.startswith("in_AD_")]
    master["n_AD_ok"] = master[ad_cols].sum(axis=1)
    master["in_AD_all_fingerprints"] = master["n_AD_ok"] == len(ad_cols)

    # verifikasi PubChem untuk kandidat yang sudah dicek
    ver = HERE / "20.4_pubchem_verification.csv"
    if ver.exists():
        v = pd.read_csv(ver)[["PubChem_ID", "status", "nama_pubchem"]]
        v = v.rename(columns={"status": "pubchem_status", "nama_pubchem": "nama_pubchem"})
        master = master.merge(v, on="PubChem_ID", how="left")

    keep = ["rank", "Metabolite", "Organism", "PubChem_ID", "Molecular_formula", "Mw",
            "Smiles", "smiles_iso_map", "prob_mean_EState", "prob_mean_PubChem",
            "prob_mean_Klekota-Roth", "prob_mean_across_fingerprints",
            "n_fingerprints_active", "consensus_label",
            "tanimoto_EState", "tanimoto_PubChem", "tanimoto_Klekota-Roth",
            "in_AD_EState", "in_AD_PubChem", "in_AD_Klekota-Roth",
            "in_AD_all_fingerprints", "identity_ambiguous",
            "Metabolite_all", "PubChem_ID_all", "Organism_all"]
    if "pubchem_status" in master:
        keep += ["pubchem_status", "nama_pubchem"]
    master = master[[c for c in keep if c in master.columns]]
    master = master.rename(columns={"Smiles": "smiles_canonical_flat",
                                    "smiles_iso_map": "smiles_isomeric"})
    master.to_csv(HERE / "20.5_master_dataset_corrected.csv", index=False)

    cand = master[(master.prob_mean_EState > 0.9) & (master.prob_mean_PubChem > 0.9)
                  & (master["prob_mean_Klekota-Roth"] > 0.9)].copy()
    cand.to_csv(HERE / "20.6_candidates_104_annotated.csv", index=False)

    print(f"berkas induk        : {len(master)} struktur")
    print(f"kandidat > 0,9      : {len(cand)}")
    print(f"di dalam AD ketiga FP: {int(cand.in_AD_all_fingerprints.sum())} / {len(cand)}")
    out = cand[~cand.in_AD_all_fingerprints]
    if len(out):
        print(f"\n=== kandidat DI LUAR jangkauan keberlakuan pada minimal satu FP "
              f"({len(out)}) ===")
        pd.set_option("display.width", 220)
        print(out[["rank", "Metabolite", "Organism", "tanimoto_EState",
                   "tanimoto_PubChem", "tanimoto_Klekota-Roth"]].head(15).to_string(index=False))
        print("\nKandidat ini sebaiknya tidak diangkat sebagai temuan utama: "
              "prediksinya berupa ekstrapolasi di luar ruang kimia data latih.")
    if "pubchem_status" in cand:
        print(f"\nverifikasi PubChem  : "
              f"{int((cand.pubchem_status == 'cocok').sum())} / {len(cand)} cocok")
    print(f"\n-> {HERE / '20.5_master_dataset_corrected.csv'}")
    print(f"-> {HERE / '20.6_candidates_104_annotated.csv'}")


if __name__ == "__main__":
    main()
