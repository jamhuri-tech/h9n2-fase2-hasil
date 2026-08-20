"""
Perbaikan alignment identitas senyawa pada pipeline Phase 2.

Masalah
-------
Pada `7_compounds - sambiloto temulawak chanonical smiles.xlsx`, kolom `Smiles`
(hasil kanonikalisasi) tidak lagi sebaris dengan kolom identitas
(IUPAC_Name / Metabolite / Organism / PubChem_ID). Kolom SMILES merupakan
permutasi dari himpunan SMILES yang benar: multiset formula molekulnya identik
dengan sumber `6.1_...smiles.xlsx` (625/625), namun hanya 159/624 baris yang
formula SMILES-nya cocok dengan formula molekul yang dideklarasikan.

Akibatnya seluruh berkas turunan (8_FP_*, 9_ndd_predictions/*, 10_*, 13_*, 14_*)
membawa struktur yang benar tetapi NAMA/CID/ORGANISME yang salah.

Perbaikan
---------
`6.1_...smiles.xlsx` terverifikasi benar (619/625 SMILES cocok dengan formula
molekul yang dideklarasikan). Skrip ini membangun peta

    canonical SMILES (tanpa stereo)  ->  identitas senyawa

dari 6.1, sehingga setiap struktur pada berkas turunan dapat dikembalikan ke
identitas yang benar.

Catatan: kanonikalisasi pada berkas 7 membuang stereokimia, sehingga 44 SMILES
kanonik memetakan ke >1 PubChem CID (stereoisomer). Baris seperti itu ditandai
`identity_ambiguous = True` dan seluruh kandidat identitasnya dicantumkan.
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
OUT = HERE

SRC_TRUSTED = DATA / "6.1_compounds - sambiloto temulawak smiles.xlsx"


def flat_canonical(smiles):
    """SMILES kanonik tanpa stereokimia -- sama dengan yang dipakai berkas 7."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def build_identity_map():
    src = pd.read_excel(SRC_TRUSTED)
    src["can"] = [flat_canonical(s) for s in src["Smiles"]]
    if src["can"].isna().any():
        raise ValueError(f"{src['can'].isna().sum()} SMILES gagal diparsing di {SRC_TRUSTED.name}")

    def agg(group):
        cids = sorted(set(group["PubChem_ID"]))
        names = list(dict.fromkeys(group["Metabolite"]))
        orgs = sorted(set(group["Organism"]))
        return pd.Series(
            {
                "Metabolite": names[0],
                "Metabolite_all": " | ".join(names),
                "PubChem_ID": cids[0],
                "PubChem_ID_all": " | ".join(str(c) for c in cids),
                "Organism": orgs[0],
                "Organism_all": " | ".join(orgs),
                "IUPAC_Name": group["IUPAC_Name"].iloc[0],
                "Smiles_isomeric": group["Smiles"].iloc[0],
                "Molecular_formula": group["Molecular formula"].iloc[0],
                "Mw": group["Mw"].iloc[0],
                "n_source_rows": len(group),
                "n_distinct_cid": len(cids),
                "identity_ambiguous": len(cids) > 1,
            }
        )

    idmap = src.groupby("can", sort=False).apply(agg, include_groups=False).reset_index()
    return src, idmap


if __name__ == "__main__":
    src, idmap = build_identity_map()
    idmap.to_csv(OUT / "18.1_identity_map_by_canonical_smiles.csv", index=False)
    print(f"baris sumber terpercaya : {len(src)}")
    print(f"SMILES kanonik unik     : {len(idmap)}")
    print(f"identitas tunggal       : {(~idmap.identity_ambiguous).sum()}")
    print(f"identitas ambigu (stereo): {idmap.identity_ambiguous.sum()}")
    print(f"-> {OUT / '18.1_identity_map_by_canonical_smiles.csv'}")
