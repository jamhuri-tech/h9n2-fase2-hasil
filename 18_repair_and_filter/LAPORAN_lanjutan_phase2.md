# Penelitian Phase 2 — koreksi data dan lanjutan

Dokumen ini melanjutkan `Penelitian Phase 2.docx`. Bagian 1 melaporkan cacat data
yang ditemukan pada rantai berkas 7 → 14 beserta perbaikannya; bagian 2–4 adalah
lanjutan analisis (ADMET inhalasi, target protein, PPI) di atas data yang sudah
diperbaiki; bagian 5 adalah rencana tahap berikutnya.

---

## 1. Temuan kritis: kolom SMILES tidak sebaris dengan identitas senyawa

### 1.1 Apa yang terjadi

Pada `7_compounds - sambiloto temulawak chanonical smiles.xlsx`, kolom `Smiles`
hasil kanonikalisasi **tidak lagi sebaris** dengan kolom `IUPAC_Name`,
`Metabolite`, `Organism`, dan `PubChem_ID`. Kolom SMILES tersusun ulang secara
mandiri, sementara kolom identitas tetap pada urutan asal.

Bukti:

| Pemeriksaan | Hasil |
|---|---|
| Formula molekul SMILES == formula yang dideklarasikan, berkas **6.1** (sumber) | **619 / 625** ✔ |
| Formula molekul SMILES == formula yang dideklarasikan, berkas **7** (kanonik) | **159 / 624** ✘ |
| Multiset formula molekul berkas 7 vs berkas 6.1 | **625 / 625 identik** |
| Himpunan SMILES kanonik berkas 7 vs kanonikalisasi berkas 6.1 (RDKit, tanpa stereo) | **identik**, 493 = 493 |

Kesimpulan: tidak ada struktur yang hilang atau bertambah. Ini murni cacat
*alignment* — seluruh struktur yang benar tetap ada, hanya berpindah baris.
Karena itu cacatnya dapat diperbaiki sepenuhnya.

Contoh (baris 1–4 berkas 7):

| Metabolite | CID | Formula dideklarasikan | SMILES yang tertempel | Formula SMILES |
|---|---|---|---|---|
| Farnesyl cyanide | 5365886 | C16H25N | `C=C(C)C(CO)CC=C(C)C` | C10 O1 |
| 1,3,5-Tris(methylene)cycloheptane | 562636 | C10H14 | `C=C(C)C1C(=C)CCCC1O` | C10 O1 |
| alpha-terpineol acetate | 111037 | C12H20O2 | `C=C(C)C1C=CC(C)(O)CC1` | C10 O1 |
| alpha-Santalene | 94164 | C15H24 | `C=C(C)C1CC=C(C)C(O)C1` | C10 O1 |

### 1.2 Seberapa jauh menyebar

Cacat ini merambat ke seluruh berkas turunan: `8_*_FP_*.csv`,
`9_ndd_predictions/*`, `10_*`, `10.1_*`, `12_*`, `13_*`, `14_*`.

Namun **prediksi NDD-nya sendiri tetap sahih**. Fingerprint dihitung dari kolom
`Smiles` di setiap baris, dan prediksi mengikuti fingerprint tersebut — diperiksa
pada berkas PubChem FP: 0 kasus fingerprint berbeda untuk SMILES yang sama,
sedangkan 34 metabolit bernama sama justru menerima fingerprint berbeda. Artinya
model menilai **struktur**, dan strukturnya konsisten. Yang salah hanya
**label identitas** yang menempel pada skor.

Maka: nilai probabilitas, peringkat, dan seleksi kandidat **tidak perlu dihitung
ulang**. Yang perlu diganti hanya nama/CID/organismenya.

### 1.3 Perbaikan

`6.1_compounds - sambiloto temulawak smiles.xlsx` terverifikasi benar dan dipakai
sebagai sumber kebenaran. Peta `SMILES kanonik → identitas` dibangun dengan RDKit
(`repair_identity_map.py`), menghasilkan 493 struktur unik:

- **449** struktur → identitas tunggal, tidak ambigu;
- **44** struktur → lebih dari satu PubChem CID, karena kanonikalisasi pada
  berkas 7 **membuang stereokimia** sehingga stereoisomer melebur. Baris seperti
  ini ditandai `identity_ambiguous = True` dan seluruh kandidat identitasnya
  dicantumkan pada kolom `*_all`.

Catatan: peleburan stereoisomer tidak memengaruhi prediksi, karena ketiga
fingerprint (PubChem, EState, Klekota-Roth) memang berbasis topologi 2D dan buta
terhadap stereokimia. Tetapi stereokimia **wajib dipulihkan** sebelum docking dan
dinamika molekuler — gunakan kolom `Smiles_isomeric` pada peta identitas.

### 1.4 Dua koreksi jumlah

- Daftar "115 senyawa kandidat" sebenarnya berisi **104 struktur unik**; 11
  struktur terkirim ganda ke Deep-PK.
- Berkas konsensus memuat 623 baris / **492 struktur**, bukan 625. Penyebabnya
  dijelaskan pada §1.6.

### 1.5 Verifikasi independen dengan menghitung ulang fingerprint

Klaim "prediksi menilai struktur, bukan nama" pada §1.2 diuji langsung: EState
fingerprint dihitung ULANG dari nol dengan RDKit
(`rdkit.Chem.EState.Fingerprinter`) untuk seluruh 625 baris, lalu dibandingkan
dengan bit yang tersimpan di `Herbal Sambiloto Temulawak_FP_Estate.csv`.

| Fingerprint tersimpan cocok dengan... | Hasil |
|---|---|
| SMILES yang ada pada baris itu | **624 / 625** |
| struktur asli metabolit yang dinamai pada baris itu | **141 / 625** |

Ini bukti langsung — bukan lagi inferensi — bahwa model menilai struktur di kolom
`Smiles`, dan bahwa nama/CID/organisme yang menempel padanyalah yang salah.
Karena itu strategi perbaikan yang dipakai (pertahankan probabilitas, ganti
identitas) memang yang benar; pelatihan ulang tidak diperlukan.

### 1.6 Letak akar masalah, dan dua baris yang hilang diam-diam

`predict_ndd_local.py` **bukan** sumber cacat. Skrip itu mengambil
`base = herbal[meta_cols]` dan `X_herb = herbal[fp_cols]` dari baris yang sama,
jadi ia hanya meneruskan alignment apa adanya dari berkas masukan. Cacatnya sudah
ada sejak berkas 7 → 8, yakni pada tahap kanonikalisasi SMILES sebelum
fingerprinting.

Namun skrip itu memuat satu cacat tersendiri. Pada tahap konsensus:

```python
wide = combined.pivot_table(index=key_cols, columns="dataset", ...)
```

`pivot_table` secara bawaan membuang baris yang punya NaN pada kolom index. Dua
baris memiliki `IUPAC_Name` kosong, sehingga **hilang tanpa peringatan**:

| Metabolite | CID | Organisme | SMILES pada baris | prob EState / PubChem / Klekota-Roth |
|---|---|---|---|---|
| (−)-beta-Sitosterol | 222284 | *Andrographis paniculata* | `CCC(CCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C)C(C)C` | 0,433 / 0,744 / 0,425 |
| Citronellyl acetate | 9017 | *Zingiber officinale* | `CCSCC` | 0,015 / 0,252 / 0,035 |

Keduanya jauh di bawah ambang 0,9, jadi **seleksi kandidat tidak terpengaruh**.
Tetap perlu diperbaiki agar tidak menggigit pada dataset lain — cukup
`pivot_table(..., dropna=False)`, atau isi `IUPAC_Name` yang kosong lebih dulu.
Sebaiknya ditambah juga penjagaan: `assert len(wide) == len(herbal)`.

---

## 2. Peringkat kandidat NDD setelah dikoreksi

Kriteria seleksi tidak diubah: rata-rata probabilitas > 0,9 pada **ketiga**
fingerprint (EState, PubChem, Klekota-Roth), masing-masing sudah merupakan
rata-rata tiga strategi uji (stratified, scaffold, similarity cluster).

**104 struktur** lolos. Sebarannya: *Zingiber officinale* 32, *Andrographis
paniculata* 28, *Melissa officinalis* 19, *Phyllanthus niruri* 15, *Forsythia
suspensa* 10 — **tidak ada** dari *Curcuma xanthorrhiza*.

Lima peringkat teratas (bandingkan dengan tabel lama yang identitasnya salah):

| No | Metabolit | Organisme | PubChem CID | Probability score |
|---|---|---|---|---|
| 1 | Glucogallin | *Phyllanthus niruri* | 124375 | 0,99479 |
| 2 | 5-O-Caffeoylshikimic acid | *Zingiber officinale* | 5281762 | 0,99403 |
| 3 | Adoxosidic acid | *Forsythia suspensa* | 13892717 | 0,99371 |
| 4 | Labiatenic acid | *Melissa officinalis* | 5281792 | 0,99334 |
| 5 | Melitric acid A | *Melissa officinalis* | 10459878 | 0,99302 |

Skor probabilitasnya sama persis dengan tabel lama — hanya nama, organisme, dan
CID-nya yang berubah, sesuai penjelasan §1.2.

Peringkat teratas kini didominasi asam fenolik dan polifenol (glukogalin,
turunan kafeoil, asam salvianolat, asam litospermat). Ini jauh lebih koheren
secara kimia dengan struktur yang benar-benar diskor daripada daftar lama yang
menyebut seskuiterpen dan flavon termetilasi.

---

## 3. Penyaringan ADMET: jalur oral dan jalur inhalasi

Berkas ADMET (`14_*.xlsx`) di-*join* ulang lewat kolom SMILES — satu-satunya
kolom yang dapat dipercaya pada berkas itu. Kolom identitasnya dibuang dan
diganti dengan hasil perbaikan.

### 3.1 Jalur oral (mengulang kriteria dokumen lama)

Kriteria: Human Oral Bioavailability 50% = Bioavailable; Blood-Brain Barrier =
Non-Penetrable; AMES Mutagenesis = Safe; Liver Injury II = Safe.

**2 struktur lolos** — jumlahnya sama dengan laporan lama, tetapi identitasnya
berbeda:

| Metabolit | Organisme | PubChem CID | CID versi lama (salah) | MW | cLogP | NDD prob |
|---|---|---|---|---|---|---|
| p-Coumaric acid | *Melissa officinalis* | **637542** | 636822 | 164,16 | 1,49 | 0,98195 |
| Caffeic acid | *Andrographis paniculata*, *Melissa officinalis* | **689043** | 643820 | 180,16 | 1,20 | 0,97830 |

Strukturnya memang benar sejak awal (asam p-kumarat dan asam kafeat); yang salah
adalah CID dan organisme yang dilekatkan padanya. Dokumen lama menyebut keduanya
berasal dari *Zingiber officinale* — itu keliru.

### 3.2 Jalur inhalasi (bagian yang sebelumnya "TBD")

Bioavailabilitas oral **sengaja tidak dipakai**: rute inhalasi memintas absorpsi
saluran cerna dan sebagian metabolisme lintas-pertama di hati, sehingga syarat
itu tidak relevan dan hanya akan membuang kandidat yang sah.

Gerbang keselamatan yang dipakai:

| Kriteria | Alasan | Sisa kandidat (kumulatif, dari 104) |
|---|---|---|
| Respiratory Disease = Safe | endpoint paling langsung untuk sediaan hirup | 49 |
| AMES Mutagenesis = Safe | genotoksisitas | 37 |
| Liver Injury II = Safe | fraksi yang masuk sirkulasi tetap melewati hati | 19 |
| Blood-Brain Barrier = Non-Penetrable | inhalasi memintas *first-pass*, paparan sistemik dapat lebih tinggi | 9 |
| hERG Blockers = Safe | risiko kardiak dari paparan sistemik | 9 |
| Carcinogenesis = Safe | keamanan jangka panjang | **9** |

**9 struktur lolos:**

| NDD rank | Metabolit | Organisme | PubChem CID | Formula | MW | cLogP | TPSA |
|---|---|---|---|---|---|---|---|
| 1 | Glucogallin | *Phyllanthus niruri* | 124375 | C13H16O10 | 332,26 | −1,63 | 166,1 |
| 3 | Adoxosidic acid | *Forsythia suspensa* | 13892717 | C16H24O10 | 376,36 | −2,24 | 166,1 |
| 23 | Gingerglycolipid C | *Zingiber officinale* | 10259020 | C33H60O14 | 680,83 | 0,57 | 225,1 |
| 24 | Gingerglycolipid B | *Zingiber officinale* | 10009754 | C33H58O14 | 678,81 | 0,34 | 225,1 |
| 29 | p-Coumaric acid | *Melissa officinalis* | 637542 | C9H8O3 | 164,16 | 1,49 | 57,5 |
| 35 | Caffeic acid | *Andrographis paniculata* | 689043 | C9H8O4 | 180,16 | 1,20 | 77,8 |
| 41 | (+)-Ascorbic acid | *Zingiber officinale* | 54670067 | C6H8O6 | 176,12 | −1,41 | 107,2 |
| 55 | Salvianic acid A | *Melissa officinalis* | 439435 | C9H10O5 | 198,17 | 0,09 | 98,0 |
| 85 | Forsythoside D | *Forsythia suspensa* | 5317383 | C20H30O13 | 478,45 | −3,20 | 219,0 |

Kedua kandidat oral termasuk di dalam sembilan ini, jadi keduanya layak untuk
kedua rancangan sediaan.

**Dua endpoint sengaja dilaporkan sebagai penanda, bukan sebagai penggugur:**

1. *Skin Sensitisation.* Deep-PK menandai 8 dari 9 kandidat sebagai "Toxic" —
   termasuk **asam askorbat (vitamin C)**. Untuk polifenol dan poliol katekolik,
   endpoint ini jelas menghasilkan positif palsu, sehingga tidak layak dijadikan
   kriteria eliminasi. (3 dari 9 juga ditandai pada Eye irritation:
   p-kumarat, kafeat, salvianat A.)
2. *Jendela fisikokimia cLogP 1–5 dan MW ≤ 500.* Kaidah ini berlaku untuk obat
   inhalasi yang bekerja **sistemik**. Untuk antivirus yang diharapkan bekerja
   **lokal di jaringan paru**, polaritas tinggi justru menguntungkan: retensi di
   paru naik dan paparan sistemik turun. Menerapkannya sebagai gerbang akan
   menggugurkan seluruh sembilan kandidat tanpa dasar yang benar. Nilainya tetap
   dicantumkan di tabel agar dapat dipertimbangkan pada tahap formulasi
   (dua gingerglikolipid, MW ≈ 680, paling berisiko dari sisi disolusi).

---

## 4. Target protein dan jaringan PPI

### 4.1 Hasil intersection — tervalidasi

Perhitungan ulang dari berkas mentah mereproduksi hasil dokumen lama:

- SwissTargetPrediction, probability > 0,5 → **28** gen unik (dokumen lama
  menulis 27; selisih 1, perlu dicek ambang yang dipakai).
- GeneCards "Influenza A Virus" → 2824 baris, **2820** *protein coding*
  (dokumen lama menulis 2824), relevance score > 2,1 → **741** gen.
- Intersection → **DPP4, MMP9, MMP2, MIF** ✔ persis sama.

Keluaran SwissTargetPrediction (anhidrase karbonat, ESR2) konsisten dengan asam
fenolik kecil, sehingga struktur yang disubmit hampir pasti sudah benar — yang
salah hanya penamaan berkasnya (`636822`, `643820`). **Tetap perlu dijalankan
ulang** dengan CID yang benar (637542, 689043) agar jejak audit bersih.

### 4.2 Analisis hub PPI (lanjutan)

Jaringan STRING (4 seed + 1 lapis ketetanggaan) berisi 9 simpul dan 18 sisi,
densitas 0,50. Peringkat sentralitas:

| Protein | Degree | Weighted degree | Betweenness | Closeness | Seed |
|---|---|---|---|---|---|
| **MMP9** | 8 | 6,36 | 0,393 | 1,000 | ✔ |
| **MMP2** | 6 | 5,07 | 0,125 | 0,800 | ✔ |
| TIMP1 | 6 | 4,12 | 0,125 | 0,800 | |
| TIMP4 | 4 | 2,89 | 0,000 | 0,667 | |
| RECK | 4 | 2,53 | 0,000 | 0,667 | |
| SCUBE3 | 2 | 1,68 | 0,000 | 0,571 | |
| DMP1 | 2 | 1,57 | 0,000 | 0,571 | |
| **DPP4** | 2 | 1,07 | 0,000 | 0,571 | ✔ |
| **MIF** | 2 | 0,90 | 0,000 | 0,571 | ✔ |

Urutan prioritas seed untuk docking: **MMP9 > MMP2 > DPP4 > MIF**. MMP9 adalah
satu-satunya simpul yang terhubung ke seluruh jaringan (closeness = 1,0) dan
memegang betweenness tertinggi, sekaligus relevance score GeneCards tertinggi
(6,67).

### 4.3 Catatan strategi: target inang vs target virus

Keempat protein tersebut — DPP4, MMP9, MMP2, MIF — semuanya protein **inang**,
bukan protein virus. Yang sedang dibangun adalah terapi *host-directed*
(pengendalian remodeling matriks dan inflamasi paru), bukan penghambatan
langsung siklus replikasi virus.

Ini sah dan menarik, tetapi **menyimpang dari usulan proposal**, yang menyatakan
"Data protein target Avian Influenza A/H9N2 diperoleh dari Protein Data Bank
(PDB)" dan menargetkan skor docking ΔG ≤ −6,5 kcal/mol pada protein target
H9N2. Rekomendasi: jalankan **kedua** jalur pada tahap docking —

- **jalur virus**: neuraminidase dan hemaglutinin H9N2, ditambah PA endonuklease
  dan NP, diambil dari PDB (sesuai proposal);
- **jalur inang**: MMP9 dan MMP2 sebagai hub utama, DPP4 dan MIF sebagai
  pendukung,

lalu laporkan keduanya. Dengan begitu indikator capaian proposal tetap terpenuhi
sementara temuan *host-directed* tetap dapat dipublikasikan.

---

## 5. Tahap berikutnya

1. **Jalankan ulang SwissTargetPrediction** untuk CID 637542 dan 689043, dan
   perluas ke tujuh kandidat inhalasi lainnya — saat ini target hanya berasal
   dari dua senyawa oral.
2. **Pulihkan stereokimia** seluruh kandidat terpilih dari kolom
   `Smiles_isomeric` sebelum pembuatan struktur 3D. Untuk 44 struktur ambigu,
   pilih stereoisomer secara eksplisit berdasarkan CID sumbernya.
3. **Siapkan reseptor**: unduh struktur PDB untuk jalur virus dan jalur inang,
   bersihkan air dan ligan bawaan, tentukan grid box.
4. **Docking** (AutoDock Vina / PyRx) untuk 9 kandidat inhalasi + 2 kandidat oral
   terhadap kedua kelompok target, ambang ΔG ≤ −6,5 kcal/mol.
5. **Optimasi metaheuristik** (GA / HOA / BA) untuk kandidat yang belum memenuhi
   ambang.
6. **Simulasi MD** (GROMACS, 50–100 ns) untuk minimal tiga kompleks terbaik,
   dievaluasi dengan RMSD, RMSF, dan MM/PBSA.
7. **Perbaiki `predict_ndd_local.py`** sesuai §1.6 (`dropna=False` pada
   `pivot_table` + penjagaan jumlah baris).
8. **Perbaiki berkas 7–14** atau tandai sebagai *deprecated*, dan bangun ulang
   turunannya dari berkas 6.1 agar tidak ada yang memakai identitas lama.

---

## Berkas yang dihasilkan

| Berkas | Isi |
|---|---|
| `verify_fingerprint_alignment.py` | Menghitung ulang EState FP dengan RDKit untuk membuktikan §1.5 |
| `repair_identity_map.py` | Membangun peta SMILES kanonik → identitas dari berkas 6.1 |
| `18.1_identity_map_by_canonical_smiles.csv` | 493 struktur + identitas + penanda ambiguitas + SMILES isomerik |
| `rebuild_candidates.py` | Menempelkan ulang identitas ke hasil prediksi NDD |
| `18.2_ndd_ranking_all_structures_corrected.csv` | 492 struktur, terurut, identitas benar |
| `18.3_candidates_prob_gt_0.9_corrected.csv` | 104 kandidat lolos ambang 0,9 |
| `admet_filter.py` | Penyaring ADMET jalur oral dan jalur inhalasi |
| `18.4_admet_annotated_corrected.csv` | 104 struktur + seluruh endpoint ADMET + deskriptor + kolom lolos/tidak |
| `ppi_hub_analysis.py` | Sentralitas jaringan STRING |
| `18.5_ppi_hub_ranking.csv` | Peringkat hub 9 protein |
