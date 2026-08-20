"""
Mengunduh struktur reseptor terpilih dari RCSB dan menghitung grid box docking
dari koordinat ligan ko-kristalnya.

Pemilihan reseptor (lihat 19.2_receptor_shortlist.csv untuk alasannya) sudah
diverifikasi: bebas mutasi rekayasa, X-ray, dan -- kecuali HA -- membawa ligan
ko-kristal yang mendefinisikan situs pengikatan sekaligus menjadi kontrol
re-docking.
"""

import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PDB_DIR = HERE / "receptors"
PDB_DIR.mkdir(exist_ok=True)

PAD = 5.0  # angstrom, ditambahkan di tiap sisi kotak

# target, pdb, residu HETATM yang menandai situs, catatan
RECEPTORS = [
    # ---------------- jalur virus ----------------
    dict(group="virus", target="Hemagglutinin H9", pdb="1JSI",
         site_res=["SIA", "GAL"], chain_hint=None,
         note="H9N2 asli (A/swine/Hong Kong/9/98); situs pengikatan reseptor "
              "ditandai analog LSTc (tipe manusia). NAG sengaja TIDAK dipakai "
              "sebagai penanda: N-glikan tersebar di seluruh permukaan HA"),
    dict(group="virus", target="Hemagglutinin H9", pdb="1JSH",
         site_res=["SIA", "GAL"], chain_hint=None,
         note="H9N2 asli; analog LSTa (tipe unggas) -- pasangan pembanding 1JSI"),
    dict(group="virus", target="Neuraminidase N2", pdb="4K1K",
         site_res=["G39"], chain_hint=None,
         note="N2 tipe liar 1,6 A dengan oseltamivir karboksilat -- "
              "kontrol positif melekat pada strukturnya"),
    dict(group="virus", target="PA endonuclease", pdb="8T5W",
         site_res=["E4Z"], chain_hint=None,
         note="PA-Nter dengan baloxavir acid; domain ini lestari antar subtipe"),
    dict(group="virus", target="Nucleoprotein", pdb="3RO5",
         site_res=["LGH"], chain_hint=None,
         note="NP dengan analog nukleozin LGH (32 atom berat). Menggantikan 4DYN, "
              "yang ligan kontrolnya 0MR berukuran 41 atom berat dan gagal kontrol "
              "(RMSD terbaik 5,84 A). LGH adalah inti nukleozin dari 0MR tanpa "
              "perpanjangan piridina-karboksamidanya. Catatan: situs NP berada di "
              "ANTARMUKA dua rantai NP -- kedua rantai wajib dipertahankan. "
              "9OTW (2,24 A, lebih baik) tidak dapat dipakai: kode CCD-nya 5 karakter "
              "sehingga identitas ligan hilang saat mmCIF dikonversi ke format PDB"),
    # ---------------- jalur inang ----------------
    dict(group="inang", target="MMP9", pdb="4XCT",
         site_res=["N73"], chain_hint=None,
         note="domain katalitik tipe liar 1,30 A dengan ARP101; Zn katalitik + Ca struktural"),
    dict(group="inang", target="MMP2", pdb="8H78",
         site_res=["L2U"], chain_hint=None, center_on_metal=24.0,
         note="domain katalitik tipe liar; MMP2 memang miskin struktur beresolusi "
              "tinggi. Kotak dipusatkan ke Zn katalitik, BUKAN ke ligan ko-kristal: "
              "L2U berbobot 916 Da sehingga kotaknya menjadi 19.900 A^3 -- terlalu "
              "lebar untuk polifenol 160-680 Da dan menurunkan ketelitian pencarian"),
    dict(group="inang", target="DPP4", pdb="4A5S",
         site_res=["N7F"], chain_hint=None,
         note="tipe liar 1,62 A dengan penghambat heterosiklik"),
    dict(group="inang", target="MIF", pdb="6B1K",
         site_res=["C9G"], chain_hint=None,
         note="tipe liar 1,17 A; situs tautomerase berada di ANTARMUKA TRIMER -- "
              "ketiga rantai wajib dipertahankan"),
]


def download(pdb_id):
    """Ambil struktur. Entri baru/besar tidak lagi punya format PDB warisan,
    jadi mmCIF diunduh lalu dikonversi dengan gemmi."""
    path = PDB_DIR / f"{pdb_id}.pdb"
    if path.exists():
        return path
    try:
        with urllib.request.urlopen(
                f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as r:
            path.write_bytes(r.read())
        return path
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    cif = PDB_DIR / f"{pdb_id}.cif"
    with urllib.request.urlopen(
            f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=90) as r:
        cif.write_bytes(r.read())
    import gemmi
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.write_pdb(str(path))
    print(f"  ({pdb_id}: tanpa format PDB warisan, dikonversi dari mmCIF)")
    return path


METALS = {"ZN", "MN", "MG", "FE", "CO", "NI"}


def parse_het(pdb_path):
    """HETATM -> daftar (resn, chain, resseq, xyz), altLoc pertama saja."""
    out = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("HETATM") or line[16] not in (" ", "A"):
            continue
        out.append((
            line[17:20].strip(), line[21].strip(), line[22:27].strip(),
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        ))
    return out


def cluster(points, cutoff=6.0):
    """Single-linkage clustering; memisahkan salinan ligan yang berjauhan.

    Diperlukan karena beberapa entri memuat lebih dari satu salinan ligan yang
    sama (8T5W punya dua E4Z). Tanpa ini, kotak akan melar melingkupi keduanya.
    """
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(points[i] - points[j]) <= cutoff:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def site_atoms(pdb_path, resnames, chain_hint=None):
    """Koordinat penanda situs: satu klaster saja.

    Jika struktur punya logam katalitik (Zn/Mn), klaster yang dipilih adalah yang
    paling dekat dengan logam tersebut -- pilihan yang benar secara kimia untuk
    metaloproteinase (MMP) dan endonuklease PA yang bergantung Mn.
    """
    het = parse_het(pdb_path)
    cand = [(c, xyz) for resn, c, _, xyz in het
            if resn in resnames and (chain_hint is None or c == chain_hint)]
    if not cand:
        return None, None, None
    metals = np.array([xyz for resn, _, _, xyz in het if resn in METALS])

    pts = np.array([xyz for _, xyz in cand])
    chains = [c for c, _ in cand]
    best = None
    for idx in cluster(pts):
        sub = pts[idx]
        if metals.size:
            score = min(np.linalg.norm(metals - a, axis=1).min() for a in sub)
        else:
            score = -len(idx)          # tanpa logam: ambil klaster terbesar
        if best is None or score < best[0]:
            best = (score, idx)
    idx = best[1]
    chain = max(set(chains[i] for i in idx), key=lambda c: sum(chains[i] == c for i in idx))
    metal_xyz = None
    if metals.size:
        d = [min(np.linalg.norm(pts[idx] - m, axis=1)) for m in metals]
        metal_xyz = metals[int(np.argmin(d))]
    return chain, pts[idx], metal_xyz


def main():
    rows = []
    for rec in RECEPTORS:
        path = download(rec["pdb"])
        chain, xyz, metal = site_atoms(path, set(rec["site_res"]), rec["chain_hint"])
        if xyz is None:
            print(f"  ! {rec['pdb']}: residu situs {rec['site_res']} tidak ditemukan")
            rows.append({**rec, "site_res": " ".join(rec["site_res"]), "site_chain": None})
            continue
        cube = rec.get("center_on_metal")
        if cube and metal is not None:
            center = metal
            size = np.array([cube, cube, cube])
        else:
            center = xyz.mean(0)
            size = (xyz.max(0) - xyz.min(0)) + 2 * PAD
            size = np.maximum(size, 18.0)  # kotak minimum agar ligan besar tetap muat
        rows.append({
            "group": rec["group"], "target": rec["target"], "pdb_id": rec["pdb"],
            "site_res": " ".join(rec["site_res"]), "site_chain": chain,
            "n_site_atoms": len(xyz),
            "center_x": round(float(center[0]), 2),
            "center_y": round(float(center[1]), 2),
            "center_z": round(float(center[2]), 2),
            "size_x": round(float(size[0]), 1),
            "size_y": round(float(size[1]), 1),
            "size_z": round(float(size[2]), 1),
            "metal_x": None if metal is None else round(float(metal[0]), 2),
            "metal_y": None if metal is None else round(float(metal[1]), 2),
            "metal_z": None if metal is None else round(float(metal[2]), 2),
            "note": rec["note"],
        })
        print(f"  {rec['pdb']}  {rec['target']:18s} rantai {chain}  "
              f"center=({center[0]:7.2f},{center[1]:7.2f},{center[2]:7.2f})  "
              f"size=({size[0]:.0f},{size[1]:.0f},{size[2]:.0f})")

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "19.3_grid_boxes.csv", index=False)
    print(f"\n-> {HERE / '19.3_grid_boxes.csv'}")
    print(f"-> berkas PDB di {PDB_DIR}")


if __name__ == "__main__":
    main()
