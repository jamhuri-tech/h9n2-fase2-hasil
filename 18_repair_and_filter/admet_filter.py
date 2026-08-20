"""
Penyaringan ADMET untuk dua rancangan sediaan: ORAL dan INHALASI.

Sumber ADMET  : 14_df full_admet result to excel.xlsx (Deep-PK, 115 baris).
Kunci join    : kolom SMILES (SATU-SATUNYA kolom pada berkas itu yang dapat
                dipercaya). Kolom IUPAC_Name/Organism/PubChem_ID pada berkas
                tersebut TIDAK sebaris dengan strukturnya -- lihat
                repair_identity_map.py -- sehingga dibuang dan diganti dengan
                identitas hasil perbaikan.

115 baris ADMET hanya memuat 104 struktur unik (11 struktur dikirim ganda ke
Deep-PK). Analisis di sini dilakukan pada tingkat struktur.

Kriteria
--------
ORAL (mengikuti definisi pada dokumen Phase 2):
    Human Oral Bioavailability 50% = Bioavailable
    Blood-Brain Barrier            = Non-Penetrable
    AMES Mutagenesis               = Safe
    Liver Injury II                = Safe

INHALASI (bagian yang sebelumnya "TBD"):
    Bioavailabilitas oral SENGAJA TIDAK dipakai -- rute inhalasi memintas
    absorpsi saluran cerna dan sebagian first-pass hati, sehingga syarat
    itu tidak relevan dan justru membuang kandidat yang sah.

    Inti (core):
        Respiratory Disease   = Safe   <- endpoint paling relevan untuk inhalasi
        AMES Mutagenesis      = Safe
        Liver Injury II       = Safe   <- fraksi yang masuk sirkulasi tetap ke hati
        Blood-Brain Barrier   = Non-Penetrable
                                       <- inhalasi memintas first-pass sehingga
                                          paparan sistemik dapat lebih tinggi
    Gerbang keselamatan penuh (safety gate) = inti + berikut:
        hERG Blockers         = Safe   <- risiko kardiak dari paparan sistemik
        Carcinogenesis        = Safe

    TIDAK dijadikan gerbang, hanya dilaporkan sebagai penanda:
        Skin Sensitisation, Eye irritation
            Deep-PK menandai 8 dari 9 kandidat sebagai "Toxic" pada Skin
            Sensitisation, termasuk asam askorbat (vitamin C). Endpoint ini
            jelas menghasilkan positif palsu untuk polifenol/poliol katekolik,
            sehingga tidak layak dipakai sebagai kriteria eliminasi.
        cLogP / MW / TPSA
            Jendela cLogP 1-5 adalah kaidah untuk obat inhalasi yang bekerja
            SISTEMIK. Untuk terapi antivirus yang bekerja LOKAL di jaringan
            paru, polaritas tinggi justru menguntungkan: retensi di paru naik
            dan paparan sistemik turun. Memakainya sebagai gerbang akan
            membuang seluruh kandidat polifenol tanpa dasar yang benar.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent

ADMET = DATA / "14_df full_admet result to excel.xlsx"
CANDIDATES = HERE / "18.3_candidates_prob_gt_0.9_corrected.csv"

P = "] Predictions"
ORAL_RULES = {
    f"[Absorption/Human Oral Bioavailability 50%{P}": "Bioavailable",
    f"[Distribution/Blood-Brain Barrier{P}": "Non-Penetrable",
    f"[Toxicity/AMES Mutagenesis{P}": "Safe",
    f"[Toxicity/Liver Injury II{P}": "Safe",
}
INHAL_CORE = {
    f"[Toxicity/Respiratory Disease{P}": "Safe",
    f"[Toxicity/AMES Mutagenesis{P}": "Safe",
    f"[Toxicity/Liver Injury II{P}": "Safe",
    f"[Distribution/Blood-Brain Barrier{P}": "Non-Penetrable",
}
INHAL_SAFETY = {
    f"[Toxicity/hERG Blockers{P}": "Safe",
    f"[Toxicity/Carcinogenesis{P}": "Safe",
}
ADVISORY = [
    f"[Toxicity/Skin Sensitisation{P}",
    f"[Toxicity/Eye irritation{P}",
]


def descriptors(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return pd.Series({k: None for k in ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"]})
    return pd.Series(
        {
            "MW": round(Descriptors.MolWt(mol), 2),
            "cLogP": round(Crippen.MolLogP(mol), 2),
            "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "HBD": rdMolDescriptors.CalcNumHBD(mol),
            "HBA": rdMolDescriptors.CalcNumHBA(mol),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(mol),
        }
    )


def apply_rules(df, rules):
    mask = pd.Series(True, index=df.index)
    for col, want in rules.items():
        mask &= df[col] == want
    return mask


def main():
    admet = pd.read_excel(ADMET)
    admet = admet.drop(columns=["IUPAC_Name", "Organism", "PubChem_ID"])
    admet["can"] = [
        Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), isomericSmiles=False) for s in admet["SMILES"]
    ]
    n_rows = len(admet)
    admet = admet.drop_duplicates(subset="can").reset_index(drop=True)
    print(f"baris ADMET {n_rows} -> struktur unik {len(admet)}")

    cand = pd.read_csv(CANDIDATES)
    cand["can"] = [
        Chem.MolToSmiles(Chem.MolFromSmiles(str(s)), isomericSmiles=False) for s in cand["Smiles"]
    ]
    ident = cand[
        ["can", "rank", "Metabolite", "Metabolite_all", "Organism", "Organism_all",
         "PubChem_ID", "PubChem_ID_all", "Molecular_formula", "identity_ambiguous",
         "prob_mean_across_fingerprints"]
    ].rename(columns={"rank": "ndd_rank", "prob_mean_across_fingerprints": "ndd_prob"})

    df = ident.merge(admet, on="can", how="inner")
    print(f"struktur ADMET yang cocok dengan kandidat NDD: {len(df)} / {len(ident)}")
    df = pd.concat([df, df["can"].apply(descriptors)], axis=1)

    df["pass_oral"] = apply_rules(df, ORAL_RULES)
    df["pass_inhalation_core"] = apply_rules(df, INHAL_CORE)
    df["pass_inhalation_safety"] = df["pass_inhalation_core"] & apply_rules(df, INHAL_SAFETY)

    df = df.sort_values("ndd_rank").reset_index(drop=True)
    df.to_csv(HERE / "18.4_admet_annotated_corrected.csv", index=False)

    show = ["ndd_rank", "Metabolite", "Organism", "PubChem_ID", "Molecular_formula",
            "ndd_prob", "MW", "cLogP", "TPSA"]

    print("\n" + "=" * 78)
    print("JALUR ORAL")
    print("=" * 78)
    oral = df[df.pass_oral]
    print(f"lolos: {len(oral)} struktur")
    print(oral[show].to_string(index=False))

    print("\n" + "=" * 78)
    print("JALUR INHALASI -- kriteria inti")
    print("=" * 78)
    core = df[df.pass_inhalation_core]
    print(f"lolos: {len(core)} struktur")
    print(core[show].to_string(index=False))

    print("\n" + "=" * 78)
    print("JALUR INHALASI -- gerbang keselamatan penuh")
    print("=" * 78)
    safe = df[df.pass_inhalation_safety]
    print(f"lolos: {len(safe)} struktur")
    adv = [c.split("/")[-1].replace(P, "") for c in ADVISORY]
    tbl = safe[show + ADVISORY].rename(columns=dict(zip(ADVISORY, adv)))
    print(tbl.to_string(index=False))

    print("\nPenyusutan per kriteria inhalasi (kumulatif, dari %d struktur):" % len(df))
    mask = pd.Series(True, index=df.index)
    for col, want in {**INHAL_CORE, **INHAL_SAFETY}.items():
        mask &= df[col] == want
        label = col.split("/")[-1].replace(P, "")
        print(f"  {label:<24} = {want:<16} -> {mask.sum():>3}")
    print("\nPenanda (dilaporkan, TIDAK menggugurkan):")
    for col in ADVISORY:
        label = col.split("/")[-1].replace(P, "")
        n = (safe[col] != "Safe").sum()
        print(f"  {label:<24} : {n} dari {len(safe)} ditandai Toxic")


if __name__ == "__main__":
    main()
