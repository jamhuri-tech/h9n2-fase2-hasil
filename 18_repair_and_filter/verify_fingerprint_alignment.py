"""
Verifikasi independen: fingerprint yang tersimpan mengikuti kolom SMILES,
bukan nama metabolit.

EState fingerprint dihitung ulang dari nol dengan RDKit lalu dibandingkan dengan
bit yang tersimpan pada berkas fingerprint. Dijalankan dua kali:
  (a) dari SMILES yang ada pada baris tersebut;
  (b) dari struktur asli metabolit yang dinamai pada baris tersebut (menurut
      berkas sumber 6.1 yang terverifikasi benar).

Jika (a) cocok hampir sempurna dan (b) tidak, maka model menilai struktur di
kolom SMILES dan label identitasnyalah yang salah -- sehingga prediksi tidak
perlu dihitung ulang, cukup identitasnya yang ditempel ulang.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.EState import Fingerprinter

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
ROOT = DATA.parent

FP_FILE = DATA / "8_Herbal Sambiloto Temulawak_FP_Estate.csv"
TRUSTED = DATA / "6.1_compounds - sambiloto temulawak smiles.xlsx"


def estate_bits(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    counts, _ = Fingerprinter.FingerprintMol(mol)
    return (np.asarray(counts) > 0).astype(int)


def main():
    fp = pd.read_csv(FP_FILE, sep=";")
    fp.columns = [str(c).strip() for c in fp.columns]
    bit_cols = [c for c in fp.columns if c.isdigit()]

    trusted = pd.read_excel(TRUSTED)
    true_smiles = dict(zip(trusted["Metabolite"], trusted["Smiles"]))

    n = match_row_smiles = match_true_structure = 0
    for _, row in fp.iterrows():
        stored = row[bit_cols].astype(int).to_numpy()

        from_row = estate_bits(row["Smiles"])
        if from_row is None:
            continue
        n += 1
        match_row_smiles += bool((from_row == stored).all())

        from_true = estate_bits(true_smiles.get(row["Metabolite"], ""))
        match_true_structure += bool(from_true is not None and (from_true == stored).all())

    print(f"baris diperiksa: {n}")
    print(f"  fingerprint == dihitung dari SMILES pada baris itu   : {match_row_smiles}/{n}")
    print(f"  fingerprint == dihitung dari struktur asli metabolit : {match_true_structure}/{n}")


if __name__ == "__main__":
    main()
