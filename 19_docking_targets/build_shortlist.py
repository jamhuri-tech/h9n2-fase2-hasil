"""Menggabungkan metadata RCSB dengan grid box menjadi satu daftar reseptor terkurasi."""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
raw = pd.read_csv(HERE / "19.1_pdb_candidates_raw.csv")
box = pd.read_csv(HERE / "19.3_grid_boxes.csv")

meta = raw.drop_duplicates("pdb_id").set_index("pdb_id")
cols = ["title", "method", "resolution_A", "released", "organism", "ligand_ids",
        "ligand_names", "other_het"]
out = box.join(meta[cols], on="pdb_id")
out["box_volume_A3"] = (out.size_x * out.size_y * out.size_z).round(0)
order = ["group", "target", "pdb_id", "resolution_A", "method", "organism", "title",
         "site_res", "site_chain", "n_site_atoms", "ligand_ids", "ligand_names",
         "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
         "box_volume_A3", "metal_x", "metal_y", "metal_z", "other_het", "note"]
out = out[order]
out.to_csv(HERE / "19.2_receptor_shortlist.csv", index=False)

pd.set_option("display.width", 250)
print(out[["group", "target", "pdb_id", "resolution_A", "ligand_ids",
           "site_chain", "box_volume_A3"]].to_string(index=False))
print(f"\n-> {HERE / '19.2_receptor_shortlist.csv'}")
