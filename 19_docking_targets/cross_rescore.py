"""
Rescoring silang untuk reseptor yang menemukan pose benar tetapi gagal
memeringkatnya: nukleoprotein (3RO5) dan PA endonuklease (8T5W).

Gagasannya: pose dicari oleh satu fungsi penilai, lalu dinilai ulang oleh yang
lain. Kalau pose yang benar naik ke peringkat satu di bawah penilaian silang --
atau kalau kedua fungsi sepakat menaruhnya di atas -- reseptor itu masih dapat
dipakai dengan protokol konsensus.

Catatan teknis: `--score_only` menolak berkas PDBQT berisi banyak pose, jadi
tiap pose dipecah ke berkasnya sendiri lebih dulu.
"""

import re
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
WORK = HERE / "docking"
VINA = HERE / "bin" / "vina"
SPLIT = HERE / "_poses"
SPLIT.mkdir(exist_ok=True)

# reseptor -> (ligan kontrol, fungsi pencari yang menemukan pose benar)
TARGETS = {"3RO5": ("LGH", "vina"), "8T5W": ("E4Z", "vina"),
           "4K1K": ("G39", "vinardo")}      # 4K1K sebagai kalibrasi: sudah lolos
OTHER = {"vina": "vinardo", "vinardo": "vina"}


def split_models(pdbqt, stem):
    """Pecah PDBQT multi-pose menjadi satu berkas per pose."""
    out, buf, n = [], [], 0
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith("MODEL"):
            buf, n = [], n + 1
            continue
        if line.startswith("ENDMDL"):
            f = SPLIT / f"{stem}_pose{n}.pdbqt"
            f.write_text("\n".join(buf) + "\n")
            out.append(f)
            continue
        buf.append(line)
    if not out and buf:
        f = SPLIT / f"{stem}_pose1.pdbqt"
        f.write_text("\n".join(buf) + "\n")
        out.append(f)
    return out


def score_only(receptor, ligand, scoring, box):
    """Vina 1.2.7 tetap menuntut dimensi grid meski hanya menilai satu pose --
    berbeda dari 1.1.2 yang tidak memerlukannya."""
    r = subprocess.run([str(VINA), "--receptor", str(receptor), "--ligand", str(ligand),
                        "--scoring", scoring, "--score_only",
                        "--center_x", str(box.center_x), "--center_y", str(box.center_y),
                        "--center_z", str(box.center_z),
                        "--size_x", str(box.size_x), "--size_y", str(box.size_y),
                        "--size_z", str(box.size_z)],
                       capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    if not m:
        m = re.search(r"Affinity:\s*(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None


def rmsd_per_pose(native_sdf, docked_pdbqt):
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit.Chem import rdMolAlign
    ref = Chem.RemoveHs(Chem.MolFromMolFile(str(native_sdf), removeHs=False))
    pm = PDBQTMolecule.from_file(str(docked_pdbqt), skip_typing=True)
    docked = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
    return [round(rdMolAlign.CalcRMS(docked, ref, prbId=i), 2)
            for i in range(docked.GetNumConformers())]


def main():
    boxes = pd.read_csv(HERE / "19.2_receptor_shortlist.csv").set_index("pdb_id")
    rows = []
    for pdb_id, (comp, finder) in TARGETS.items():
        rdir = WORK / pdb_id
        tag = "" if finder == "vina" else f"_{finder}"
        multi = rdir / f"native_{comp}_redock{tag}.pdbqt"
        if not multi.exists():
            print(f"[{pdb_id}] {multi.name} tidak ada -- dilewati")
            continue

        rmsds = rmsd_per_pose(rdir / f"native_{comp}.sdf", multi)
        poses = split_models(multi, f"{pdb_id}_{comp}")
        rec = rdir / f"{pdb_id}.pdbqt"
        other = OTHER[finder]

        for i, (pose, rmsd) in enumerate(zip(poses, rmsds), start=1):
            rows.append({
                "pdb_id": pdb_id, "ligan_kontrol": comp, "peringkat_asli": i,
                "rmsd_A": rmsd,
                f"skor_{finder}": score_only(rec, pose, finder, boxes.loc[pdb_id]),
                f"skor_{other}": score_only(rec, pose, other, boxes.loc[pdb_id]),
            })
        print(f"[{pdb_id}] {len(poses)} pose dinilai ulang "
              f"({finder} -> {other})", flush=True)

    df = pd.DataFrame(rows)
    df["skor_konsensus"] = df[["skor_vina", "skor_vinardo"]].mean(axis=1)
    df.to_csv(HERE / "19.12_cross_rescore.csv", index=False)

    pd.set_option("display.width", 250)
    for pdb_id in df.pdb_id.unique():
        sub = df[df.pdb_id == pdb_id].copy()
        print(f"\n=== {pdb_id} ({sub.ligan_kontrol.iloc[0]}) ===")
        print(sub[["peringkat_asli", "rmsd_A", "skor_vina", "skor_vinardo",
                   "skor_konsensus"]].round(2).to_string(index=False))
        best = sub.loc[sub.rmsd_A.idxmin()]
        for col in ["skor_vina", "skor_vinardo", "skor_konsensus"]:
            if sub[col].isna().any() or pd.isna(best[col]):
                print(f"  {col.replace('skor_', '')}: skor tidak lengkap, "
                      f"peringkat tidak dihitung")
                continue
            rank = int((sub[col] < best[col]).sum()) + 1
            verdict = "LULUS" if rank == 1 else "gagal"
            print(f"  pose terbaik (RMSD {best.rmsd_A} A) -> peringkat {rank} "
                  f"menurut {col.replace('skor_', '')}  [{verdict}]")
    print(f"\n-> {HERE / '19.12_cross_rescore.csv'}")


if __name__ == "__main__":
    main()
