"""Melaporkan kemajuan docking + ETA setiap interval, berhenti saat kampanye selesai."""
import os, sys, time, glob, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
TOTAL = 189            # 7 ligan x 9 reseptor x 3 seed
INTERVAL = int(sys.argv[1]) if len(sys.argv) > 1 else 900
SUMMARY = os.path.join(HERE, "19.7_docking_summary.csv")
LOG = os.path.join(HERE, "docking_run.log")

def done():
    return len(glob.glob(os.path.join(HERE, "docking", "*", "*seed*.pdbqt")))

def alive():
    r = subprocess.run(["pgrep", "-f", "run_docking.py"], capture_output=True)
    return r.returncode == 0

# Laju diukur sejak PEMANTAU ini mulai, bukan sejak berkas log dibuat.
# docking_run.log ditulis ulang saat kampanye diluncurkan ulang, sehingga
# st_birthtime-nya masih menunjuk ke run pertama yang jauh lebih lambat --
# memakainya membuat estimasi meleset berkali lipat.
start = time.time()
base = done()                      # run yang sudah ada sebelum pemantau ini

while True:
    time.sleep(INTERVAL)
    n = done()
    el = time.time() - start
    new = max(n - base, 1)
    rate = el / new                 # detik per run, terukur pada jendela pemantauan ini
    left = TOTAL - n
    eta = left * rate
    pct = 100.0 * n / TOTAL
    finish = time.strftime("%H:%M", time.localtime(time.time() + eta))

    if os.path.exists(SUMMARY) or (not alive()):
        state = "SELESAI" if os.path.exists(SUMMARY) else "BERHENTI (proses tidak aktif)"
        print(f"{state} | {n}/{TOTAL} run", flush=True)
        break

    print(f"docking {n}/{TOTAL} ({pct:.0f}%) | {rate:.0f} det/run | "
          f"sisa ~{eta/60:.0f} menit | perkiraan kelar {finish}", flush=True)
