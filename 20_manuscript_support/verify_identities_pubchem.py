"""
Verifikasi 104 identitas kandidat terhadap PubChem -- sumber otoritatif di luar
berkas kita sendiri.

Sampai sekarang perbaikan identitas hanya diverifikasi secara INTERNAL: rumus
molekul yang dihitung dari struktur dibandingkan dengan rumus yang tercatat di
berkas yang sama. Itu menangkap cacat keselarasan, tetapi tidak membuktikan
PubChem ID yang tercantum memang milik senyawa tersebut.

Skrip ini menanyakan langsung ke PubChem: untuk tiap CID, apa rumus molekulnya,
apa strukturnya, dan apa namanya. Lalu dibandingkan dengan yang kita klaim.
"""

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
REPAIR = HERE.parent / "18_repair_and_filter"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BATCH = 50
PROPS = "MolecularFormula,ConnectivitySMILES,SMILES,Title"


def fetch(cids):
    out = {}
    for i in range(0, len(cids), BATCH):
        chunk = [str(c) for c in cids[i:i + BATCH]]
        url = f"{PUG}/compound/cid/{','.join(chunk)}/property/{PROPS}/JSON"
        for attempt in range(3):
            try:
                d = json.loads(urllib.request.urlopen(url, timeout=60).read())
                for p in d["PropertyTable"]["Properties"]:
                    out[int(p["CID"])] = p
                break
            except urllib.error.HTTPError as e:
                if attempt == 2:
                    print(f"  ! gagal untuk {len(chunk)} CID: HTTP {e.code}")
                time.sleep(2)
        time.sleep(0.3)          # PubChem: maksimum 5 permintaan/detik
        print(f"  {min(i + BATCH, len(cids))}/{len(cids)} CID diambil", flush=True)
    return out


def flat(smiles):
    m = Chem.MolFromSmiles(str(smiles))
    return None if m is None else Chem.MolToSmiles(m, isomericSmiles=False)


def main():
    cand = pd.read_csv(REPAIR / "18.3_candidates_prob_gt_0.9_corrected.csv")
    cids = sorted({int(c) for c in cand.PubChem_ID.dropna()})
    print(f"kandidat: {len(cand)} | CID unik: {len(cids)}\n")

    data = fetch(cids)
    print(f"\nditerima dari PubChem: {len(data)}/{len(cids)}\n")

    rows = []
    for r in cand.itertuples():
        cid = int(r.PubChem_ID)
        p = data.get(cid)
        if p is None:
            rows.append({"PubChem_ID": cid, "Metabolite": r.Metabolite,
                         "status": "CID tidak ditemukan"})
            continue
        pub_flat = flat(p.get("ConnectivitySMILES") or p.get("SMILES"))
        ours_flat = flat(r.Smiles)
        formula_ok = str(p["MolecularFormula"]) == str(r.Molecular_formula)
        struct_ok = pub_flat is not None and pub_flat == ours_flat
        rows.append({
            "PubChem_ID": cid, "Metabolite": r.Metabolite,
            "nama_pubchem": p.get("Title", ""),
            "formula_kita": r.Molecular_formula,
            "formula_pubchem": p["MolecularFormula"],
            "formula_cocok": formula_ok,
            "struktur_cocok": struct_ok,
            "identity_ambiguous": r.identity_ambiguous,
            "status": "cocok" if (formula_ok and struct_ok)
                      else ("struktur beda" if formula_ok else "formula beda"),
        })

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "20.4_pubchem_verification.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== hasil ===")
    print(df.status.value_counts().to_string())
    n_ok = int((df.status == "cocok").sum())
    print(f"\ncocok penuh (formula DAN struktur): {n_ok}/{len(df)} "
          f"({100 * n_ok / len(df):.1f}%)")

    bad = df[df.status != "cocok"]
    if len(bad):
        print("\n=== perlu diperiksa manual ===")
        print(bad[["PubChem_ID", "Metabolite", "nama_pubchem", "formula_kita",
                   "formula_pubchem", "status"]].to_string(index=False))
    print(f"\n-> {HERE / '20.4_pubchem_verification.csv'}")


if __name__ == "__main__":
    main()
