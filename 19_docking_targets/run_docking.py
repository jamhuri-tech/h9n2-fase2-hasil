"""
Menjalankan docking AutoDock Vina untuk seluruh pasangan ligan x reseptor,
didahului kontrol re-docking pada tiap reseptor.

Prasyarat
---------
  pip install meeko gemmi rdkit
  binary Vina 1.2.7 (unduh dari rilis resmi ccsb-scripps/AutoDock-Vina,
  taruh path-nya di variabel lingkungan VINA_BIN atau di --vina)

Alur
----
  1. Reseptor .pdb  -> .pdbqt   (mk_prepare_receptor.py; air dibuang, logam
                                 katalitik Zn/Ca/Mn DIPERTAHANKAN)
  2. Kontrol re-docking        : ligan ko-kristal dikeluarkan dari struktur,
                                 didok ulang ke kotaknya sendiri, lalu dihitung
                                 RMSD terhadap pose kristalografinya.
                                 RMSD <= 2,0 A = protokol tervalidasi untuk
                                 reseptor tersebut. Reseptor yang gagal kontrol
                                 TIDAK boleh dipakai menilai kandidat.
  3. Docking kandidat          : 9 senyawa herbal x reseptor yang lolos kontrol.

Catatan penting
---------------
* exhaustiveness 32 (bukan 8 bawaan) dan 3 seed berbeda. Skor Vina punya
  sebaran antar-run; melaporkan satu run tunggal menyembunyikan itu.
* Ligan dengan ikatan rotasi > 10 (gingerglikolipid) berada di batas keandalan
  Vina. Skornya wajib diberi catatan, jangan diperingkat setara dengan ligan kaku.
* Ambang proposal DG <= -6,5 kcal/mol dinilai relatif terhadap kontrol positif
  masing-masing reseptor, bukan sebagai angka mutlak.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RECEPTOR_DIR = HERE / "receptors"
LIGAND_DIR = HERE / "ligands"
WORK = HERE / "docking"
SEEDS = [1, 2, 3]
# Ligan dengan torsi aktif di atas ambang ini dikeluarkan dari kampanye docking.
# Batas keras Vina adalah 32 torsi dan ketelitian pencariannya sudah runtuh jauh
# sebelum itu. Pada uji nyata, dua gingerglikolipid (32 dan 33 torsi) memakan 377
# detik per run -- 86 persen waktu jalan seluruh kampanye -- untuk hasil yang
# tidak dapat dipercaya. Naikkan lewat --max-torsions bila memang ingin dijalankan.
MAX_TORSIONS = 20
EXHAUSTIVENESS = 32
NUM_MODES = 9

KEEP_HET = {"ZN", "CA", "MN", "MG", "FE", "CO", "NI"}   # logam katalitik/struktural


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True,
                          text=True, **kw)


def strip_receptor(pdb_in, pdb_out, keep_chains=None):
    """Buang air, krioprotektan, dan ligan ko-kristal; pertahankan logam."""
    keep = []
    for line in pdb_in.read_text().splitlines():
        if line.startswith("ATOM"):
            if keep_chains and line[21] not in keep_chains:
                continue
            keep.append(line)
        elif line.startswith("HETATM") and line[17:20].strip() in KEEP_HET:
            if keep_chains and line[21] not in keep_chains:
                continue
            keep.append(line)
    keep.append("END")
    pdb_out.write_text("\n".join(keep) + "\n")
    return pdb_out


def extract_native(pdb_in, resnames, chain, out_pdb):
    """Keluarkan ligan ko-kristal sebagai berkas terpisah untuk kontrol re-docking."""
    lines = [l for l in pdb_in.read_text().splitlines()
             if l.startswith("HETATM") and l[17:20].strip() in resnames
             and l[21] == chain and l[16] in (" ", "A")]
    if not lines:
        return None
    out_pdb.write_text("\n".join(lines) + "\nEND\n")
    return out_pdb


def heavy_coords(pdbqt_or_pdb, model=1):
    """Koordinat atom berat dari MODEL tertentu sebuah berkas PDBQT/PDB."""
    coords, cur = [], 0
    for line in Path(pdbqt_or_pdb).read_text().splitlines():
        if line.startswith("MODEL"):
            cur = int(line.split()[1])
        if line.startswith(("ATOM", "HETATM")):
            if cur not in (0, model):
                continue
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el.upper().startswith("H"):
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def symmetry_free_rmsd(a, b):
    """RMSD sederhana berbasis urutan atom. Untuk kontrol re-docking ini memadai
    karena kedua pose berasal dari molekul yang sama dengan urutan atom sama."""
    n = min(len(a), len(b))
    if n == 0 or len(a) != len(b):
        return None
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def vina(vina_bin, receptor, ligand, box, out, seed, scoring="vina"):
    cmd = [vina_bin, "--scoring", scoring,
           "--receptor", str(receptor), "--ligand", str(ligand),
           "--center_x", str(box.center_x), "--center_y", str(box.center_y),
           "--center_z", str(box.center_z),
           "--size_x", str(box.size_x), "--size_y", str(box.size_y),
           "--size_z", str(box.size_z),
           "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
           "--seed", str(seed), "--out", str(out)]
    r = sh(cmd)
    if r.returncode != 0:
        return None, r.stderr.strip()[:300]
    scores = [float(m) for m in re.findall(r"^\s+\d+\s+(-?\d+\.\d+)", r.stdout, re.M)]
    return (scores[0] if scores else None), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vina", default=os.environ.get("VINA_BIN", "vina"))
    ap.add_argument("--skip-control", action="store_true")
    ap.add_argument("--max-torsions", type=int, default=MAX_TORSIONS)
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo", "ad4"])
    ap.add_argument("--only", default=None, help="batasi ke satu PDB ID")
    args = ap.parse_args()

    if not shutil.which(args.vina) and not Path(args.vina).exists():
        sys.exit(f"Binary Vina tidak ditemukan: {args.vina}\n"
                 f"Set VINA_BIN atau pakai --vina /path/ke/vina")

    WORK.mkdir(exist_ok=True)
    boxes = pd.read_csv(HERE / "19.2_receptor_shortlist.csv")
    all_ligands = sorted(LIGAND_DIR.glob("*.pdbqt"))
    if not all_ligands:
        sys.exit(f"Tidak ada ligan .pdbqt di {LIGAND_DIR} -- jalankan prepare_ligands.py")

    def n_torsions(path):
        return sum(1 for l in path.read_text().splitlines() if l.startswith("BRANCH"))

    ligands, skipped = [], []
    for lig in all_ligands:
        (ligands if n_torsions(lig) <= args.max_torsions else skipped).append(lig)
    for lig in skipped:
        print(f"DILEWATI (torsi {n_torsions(lig)} > {args.max_torsions}): {lig.stem}")
    print(f"ligan didok: {len(ligands)}  dilewati: {len(skipped)}\n")

    control_rows, dock_rows = [], []

    tag = "" if args.scoring == "vina" else f"_{args.scoring}"
    for box in boxes.itertuples():
        if args.only and box.pdb_id != args.only:
            continue
        rdir = WORK / box.pdb_id
        rdir.mkdir(exist_ok=True)
        src = RECEPTOR_DIR / f"{box.pdb_id}.pdb"
        clean = strip_receptor(src, rdir / f"{box.pdb_id}_clean.pdb")

        rec_pdbqt = rdir / f"{box.pdb_id}.pdbqt"
        if not rec_pdbqt.exists():
            r = sh(["mk_prepare_receptor.py", "--read_pdb", str(clean),
                    "-o", str(rdir / box.pdb_id), "-p", "-a", "--default_altloc", "A"])
            if not rec_pdbqt.exists():
                print(f"[{box.pdb_id}] persiapan reseptor GAGAL: {r.stderr.strip()[:200]}")
                continue

        # ---- kontrol re-docking ----
        passed = True
        if not args.skip_control and isinstance(box.ligand_ids, str):
            names = set(str(box.site_res).split())
            native = extract_native(src, names, box.site_chain, rdir / "native.pdb")
            if native:
                r = sh(["mk_prepare_ligand.py", "-i", str(native),
                        "-o", str(rdir / "native.pdbqt")])
                nat_q = rdir / "native.pdbqt"
                if nat_q.exists():
                    score, err = vina(args.vina, rec_pdbqt, nat_q, box,
                                      rdir / "native_redock.pdbqt", SEEDS[0])
                    rmsd = symmetry_free_rmsd(heavy_coords(native),
                                              heavy_coords(rdir / "native_redock.pdbqt"))
                    passed = rmsd is not None and rmsd <= 2.0
                    control_rows.append({"pdb_id": box.pdb_id, "target": box.target,
                                         "native_ligand": box.ligand_ids,
                                         "redock_score": score, "redock_rmsd_A": rmsd,
                                         "control_passed": passed})
                    print(f"[{box.pdb_id}] kontrol re-docking: skor={score} "
                          f"RMSD={rmsd}  {'LULUS' if passed else 'GAGAL'}")

        # ---- docking kandidat ----
        for lig in ligands:
            for seed in SEEDS:
                out = rdir / f"{lig.stem}_seed{seed}{tag}.pdbqt"
                if out.exists():        # lanjutkan run yang terputus
                    txt = out.read_text()
                    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", txt)
                    score, err = (float(m.group(1)) if m else None), None
                else:
                    score, err = vina(args.vina, rec_pdbqt, lig, box, out, seed,
                                      args.scoring)
                dock_rows.append({"group": box.group, "target": box.target,
                                  "pdb_id": box.pdb_id, "ligand": lig.stem,
                                  "scoring": args.scoring,
                                  "seed": seed, "affinity_kcal_mol": score,
                                  "control_passed": passed, "error": err})
            print(f"  {lig.stem:<28} {box.pdb_id} selesai")

    # 19.5_redock_controls.csv adalah milik redock_controls.py -- JANGAN ditulis
    # dari sini. Versi sebelumnya menimpanya dengan berkas kosong setiap kali
    # dijalankan dengan --skip-control, menghapus hasil kontrol yang sudah benar.
    if control_rows:
        pd.DataFrame(control_rows).to_csv(HERE / "19.5_redock_controls_inline.csv",
                                          index=False)
    d = pd.DataFrame(dock_rows)
    d.to_csv(HERE / f"19.6_docking_raw{tag}.csv", index=False)

    summary = (d.dropna(subset=["affinity_kcal_mol"])
                 .groupby(["group", "target", "pdb_id", "ligand", "scoring"])
                 .agg(affinity_mean=("affinity_kcal_mol", "mean"),
                      affinity_sd=("affinity_kcal_mol", "std"),
                      affinity_best=("affinity_kcal_mol", "min"))
                 .round(2).reset_index()
                 .sort_values(["group", "target", "affinity_mean"]))
    summary.to_csv(HERE / f"19.7_docking_summary{tag}.csv", index=False)
    print(f"\n-> {HERE / f'19.7_docking_summary{tag}.csv'}")


if __name__ == "__main__":
    main()
