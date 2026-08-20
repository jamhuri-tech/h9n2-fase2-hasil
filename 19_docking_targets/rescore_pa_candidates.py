"""
Penilaian ulang pose kandidat pada PA endonuklease (8T5W) dengan protokol
dua langkah: pose dicari AutoDock Vina, lalu dinilai ulang Vinardo.

Dasarnya ada di 19.12_cross_rescore.csv: pada kontrol re-docking, pose benar
ligan baloxavir (RMSD 0,30 A) hanya menempati peringkat 4 menurut skor Vina,
tetapi naik ke PERINGKAT 1 ketika pose yang sama dinilai ulang oleh Vinardo.
Protokol dua langkah inilah yang tervalidasi untuk reseptor ini -- bukan Vina
sendirian, bukan pula Vinardo yang mencari sendiri.
"""

import re
import subprocess
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
WORK = HERE / "docking" / "8T5W"
VINA = HERE / "bin" / "vina"
SPLIT = HERE / "_poses"
SPLIT.mkdir(exist_ok=True)
PDB_ID = "8T5W"
SEEDS = [1, 2, 3]


def split_models(pdbqt, stem):
    out, buf, n = [], [], 0
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith("MODEL"):
            buf, n = [], n + 1
            continue
        if line.startswith("ENDMDL"):
            f = SPLIT / f"{stem}_p{n}.pdbqt"
            f.write_text("\n".join(buf) + "\n")
            out.append(f)
            continue
        buf.append(line)
    return out


def score(receptor, ligand, box, scoring="vinardo"):
    r = subprocess.run([str(VINA), "--receptor", str(receptor), "--ligand", str(ligand),
                        "--scoring", scoring, "--score_only",
                        "--center_x", str(box.center_x), "--center_y", str(box.center_y),
                        "--center_z", str(box.center_z),
                        "--size_x", str(box.size_x), "--size_y", str(box.size_y),
                        "--size_z", str(box.size_z)], capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None


def main():
    box = pd.read_csv(HERE / "19.2_receptor_shortlist.csv").set_index("pdb_id").loc[PDB_ID]
    rec = WORK / f"{PDB_ID}.pdbqt"
    rows = []
    for lig in sorted((HERE / "ligands").glob("*.pdbqt")):
        best = []
        for seed in SEEDS:
            src = WORK / f"{lig.stem}_seed{seed}.pdbqt"
            if not src.exists():
                continue
            for pose in split_models(src, f"{PDB_ID}_{lig.stem}_s{seed}"):
                s = score(rec, pose, box)
                if s is not None:
                    best.append(s)
        if not best:
            continue
        rows.append({"pdb_id": PDB_ID, "ligand": lig.stem,
                     "ligan": lig.stem.split("_", 1)[1].replace("_", " "),
                     "n_pose_dinilai": len(best),
                     "vinardo_rescore_best": round(min(best), 2),
                     "vinardo_rescore_mean_top3": round(sum(sorted(best)[:3]) / 3, 2)})
        print(f"  {rows[-1]['ligan']:<20} {len(best):2d} pose -> "
              f"terbaik {rows[-1]['vinardo_rescore_best']:6.2f}", flush=True)

    df = pd.DataFrame(rows).sort_values("vinardo_rescore_best")
    df.to_csv(HERE / "19.13_pa_rescored.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n=== PA endonuklease, protokol Vina cari + Vinardo nilai ulang ===")
    print(df[["ligan", "n_pose_dinilai", "vinardo_rescore_best",
              "vinardo_rescore_mean_top3"]].to_string(index=False))
    print("\nacuan: baloxavir acid (ligan ko-kristal) -7,31 kcal/mol pada skala yang sama")
    print(f"\n-> {HERE / '19.13_pa_rescored.csv'}")


if __name__ == "__main__":
    main()
