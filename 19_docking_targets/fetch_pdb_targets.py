"""
Mengumpulkan kandidat struktur reseptor untuk tahap docking, langsung dari RCSB PDB.

Dua kelompok target:
  JALUR VIRUS  -- protein Influenza A: hemaglutinin H9, neuraminidase N2,
                  PA endonuklease, nukleoprotein, M2. Diutamakan galur H9N2.
  JALUR INANG  -- MMP9, MMP2, DPP4, MIF (hasil intersection
                  SwissTargetPrediction x GeneCards), dicari lewat aksesi UniProt.

Untuk tiap entri dicatat metode, resolusi, organisme, dan ligan ko-kristal.
Ligan ko-kristal penting: posisinya dipakai untuk menempatkan grid box docking,
dan ia menjadi kontrol re-docking (redock RMSD) untuk memvalidasi protokol.
"""

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

# heteroatom yang bukan ligan bermakna: air, ion, buffer, krioprotektan, glikan
JUNK = {
    "HOH", "DOD", "SO4", "PO4", "CL", "NA", "K", "CA", "MG", "MN", "ZN", "FE", "NI",
    "CD", "CU", "HG", "IOD", "BR", "F", "ACT", "GOL", "EDO", "PEG", "PG4", "PGE",
    "1PE", "DMS", "MES", "EPE", "TRS", "IMD", "FMT", "ACY", "CIT", "TLA", "MPD",
    "BME", "DTT", "AZI", "NO3", "SCN", "NH4", "UNX", "UNL",
    "NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "GLC", "BGC", "SIA", "XYP", "A2G",
}

VIRAL_QUERIES = {
    "Hemagglutinin H9": ["H9N2 hemagglutinin", "hemagglutinin H9"],
    "Neuraminidase N2": ["N2 neuraminidase influenza", "H9N2 neuraminidase",
                         "influenza neuraminidase N2 subtype"],
    "PA endonuclease": ["influenza PA endonuclease", "influenza polymerase acidic protein N-terminal"],
    "Nucleoprotein": ["influenza A nucleoprotein"],
    "M2 channel": ["H9N2 M2 proton channel"],
}

HOST_UNIPROT = {
    "MMP9": "P14780",
    "MMP2": "P08253",
    "DPP4": "P27487",
    "MIF": "P14174",
}


def post(url, payload, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(payload).encode(),
                headers={"Content-Type": "application/json"})
            resp = urllib.request.urlopen(req, timeout=45)
            if resp.status == 204:
                return None
            return json.loads(resp.read())
        except (urllib.error.URLError, TimeoutError):
            if i == tries - 1:
                raise
            time.sleep(2)


def search_full_text(value, rows=200):
    body = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": value}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}},
    }
    d = post(SEARCH_URL, body)
    return [h["identifier"] for h in d["result_set"]] if d else []


def search_uniprot(accession, rows=300):
    body = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers"
                         ".reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": accession}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}},
    }
    d = post(SEARCH_URL, body)
    return [h["identifier"] for h in d["result_set"]] if d else []


GQL = """
{ entries(entry_ids: %s) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined deposited_polymer_entity_instance_count }
    rcsb_accession_info { initial_release_date }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      rcsb_entity_source_organism { scientific_name }
    }
    nonpolymer_entities { nonpolymer_comp { chem_comp { id name formula_weight } } }
} }"""


def fetch_details(ids):
    out = []
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        d = post(GRAPHQL_URL, {"query": GQL % json.dumps(chunk)})
        out.extend(d["data"]["entries"])
    return out


def summarise(entry, group, target):
    res = (entry["rcsb_entry_info"] or {}).get("resolution_combined") or []
    ligs, junk = [], []
    for ne in entry.get("nonpolymer_entities") or []:
        cc = ne["nonpolymer_comp"]["chem_comp"]
        (junk if cc["id"] in JUNK else ligs).append(cc)
    chains = entry.get("polymer_entities") or []
    orgs, descs = [], []
    for pe in chains:
        descs.append((pe["rcsb_polymer_entity"] or {}).get("pdbx_description") or "")
        for o in pe.get("rcsb_entity_source_organism") or []:
            if o.get("scientific_name"):
                orgs.append(o["scientific_name"])
    return {
        "group": group,
        "target": target,
        "pdb_id": entry["rcsb_id"],
        "title": (entry["struct"] or {}).get("title", ""),
        "method": (entry["exptl"] or [{}])[0].get("method", ""),
        "resolution_A": res[0] if res else None,
        "released": (entry["rcsb_accession_info"] or {}).get("initial_release_date", "")[:10],
        "organism": " | ".join(dict.fromkeys(orgs)),
        "chains": " | ".join(dict.fromkeys(d for d in descs if d)),
        "n_ligands": len(ligs),
        "ligand_ids": " ".join(c["id"] for c in ligs),
        "ligand_names": " | ".join(c["name"] for c in ligs)[:300],
        "heaviest_ligand_mw": max([c["formula_weight"] or 0 for c in ligs], default=0),
        "other_het": " ".join(sorted({c["id"] for c in junk})),
    }


def main():
    rows, seen = [], set()

    for target, queries in VIRAL_QUERIES.items():
        ids = []
        for q in queries:
            ids.extend(search_full_text(q))
        ids = [i for i in dict.fromkeys(ids)]
        print(f"[virus] {target:22s} {len(ids):3d} entri")
        for e in fetch_details(ids):
            key = (target, e["rcsb_id"])
            if key in seen:
                continue
            seen.add(key)
            rows.append(summarise(e, "virus", target))

    for target, acc in HOST_UNIPROT.items():
        ids = search_uniprot(acc)
        print(f"[inang] {target:22s} {len(ids):3d} entri  (UniProt {acc})")
        for e in fetch_details(ids):
            key = (target, e["rcsb_id"])
            if key in seen:
                continue
            seen.add(key)
            rows.append(summarise(e, "inang", target))

    df = pd.DataFrame(rows)
    df["is_xray"] = df["method"].str.contains("X-RAY", na=False)
    df["is_h9n2"] = df["organism"].str.contains("H9N2", na=False)
    df["has_ligand"] = df["n_ligands"] > 0
    df["druglike_ligand"] = df["heaviest_ligand_mw"] >= 150
    df = df.sort_values(
        ["group", "target", "is_h9n2", "druglike_ligand", "is_xray", "resolution_A"],
        ascending=[True, True, False, False, False, True])
    df.to_csv(HERE / "19.1_pdb_candidates_raw.csv", index=False)
    print(f"\ntotal entri: {len(df)}  ->  {HERE / '19.1_pdb_candidates_raw.csv'}")


if __name__ == "__main__":
    main()
