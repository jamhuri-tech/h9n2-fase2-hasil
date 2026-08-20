# Hasil komputasi Fase 2 — antivirus herbal H9N2

> **Status:** hasil antara dari penelitian yang sedang berjalan; naskahnya belum
> terbit. Angka-angka di sini masih dapat berubah sebelum publikasi. Bila Anda
> memakainya, mohon hubungi tim penelitinya lebih dulu.

Seluruh berkas hasil dari tahap perbaikan data, penyaringan ADMET, penentuan
target, docking, dan verifikasi lanjutan. Setiap angka yang dikutip di laporan
tim dapat ditelusuri ke salah satu berkas di sini.

Skrip disertakan bersama datanya, jadi setiap tahap dapat dijalankan ulang dari
nol — bukan sekadar hasil akhir yang harus dipercaya begitu saja.

## Kalau hanya ingin melihat yang penting

| Ingin tahu | Buka berkas ini |
|---|---|
| Daftar 104 kandidat, identitas sudah benar & terverifikasi PubChem | `20_manuscript_support/20.6_candidates_104_annotated.csv` |
| Kandidat mana yang lolos ADMET oral / inhalasi | `18_repair_and_filter/18.4_admet_annotated_corrected.csv` |
| Hasil docking + efisiensi ligan + status validasi | `19_docking_targets/19.8_docking_final.csv` |
| Pasangan mana yang lolos pemeriksaan pose | `19_docking_targets/19.11_verdict.csv` |
| Reseptor mana yang boleh dipercaya | `19_docking_targets/19.5_redock_controls.csv` (+ `_vinardo`) |
| Residu kontak bernama untuk naskah | `20_manuscript_support/20.7_contact_residues_named.csv` |

## Isi per folder

### `18_repair_and_filter/` — perbaikan identitas dan penyaringan ADMET

| Berkas | Isi |
|---|---|
| `18.1_identity_map_by_canonical_smiles.csv` | Peta SMILES kanonik → identitas benar, 493 struktur |
| `18.2_ndd_ranking_all_structures_corrected.csv` | 492 struktur, terurut, identitas benar |
| `18.3_candidates_prob_gt_0.9_corrected.csv` | 104 kandidat lolos ambang 0,9 pada ketiga sidik jari |
| `18.4_admet_annotated_corrected.csv` | Seluruh endpoint ADMET + deskriptor + kolom lolos/tidak |
| `18.5_ppi_hub_ranking.csv` | Sentralitas jaringan STRING, 9 protein |
| `verify_fingerprint_alignment.py` | Bukti cacat identitas: hitung ulang sidik jari dari nol |
| `repair_identity_map.py`, `rebuild_candidates.py` | Perbaikan identitas |
| `admet_filter.py`, `ppi_hub_analysis.py` | Penyaringan ADMET dan analisis jaringan |

### `19_docking_targets/` — pemilihan reseptor, docking, validasi

| Berkas | Isi |
|---|---|
| `19.1_pdb_candidates_raw.csv` | 708 struktur PDB yang ditimbang untuk kedua jalur |
| `19.2_receptor_shortlist.csv` | 9 reseptor terpilih + grid box + alasan |
| `19.4_ligands_prepared.csv` | 9 ligan: stereopusat, muatan formal, torsi |
| `19.5_redock_controls.csv` / `_vinardo.csv` | Kontrol re-docking, dua fungsi penilai |
| `19.6_docking_raw.csv` | 231 run mentah, per seed |
| `19.8_docking_final.csv` | Afinitas + efisiensi ligan + status kontrol |
| `19.9_mmp2_surrogate_control.csv` | Uji koordinasi Zn sebagai kontrol pengganti MMP2 |
| `19.10_pose_validation.csv` | Recall, Jaccard, jarak ke Zn tiap pasangan |
| `19.11_verdict.csv` | Gabungan skor + LE + validasi pose, kolom lolos/gugur |
| `19.12_cross_rescore.csv` | Rescoring silang untuk NP dan PA endonuklease |
| `19.13_pa_rescored.csv` | PA endonuklease, protokol Vina-cari + Vinardo-nilai |
| `receptors/` | Struktur PDB yang dipakai, apa adanya dari RCSB |
| `ligands/` | Ligan siap docking (SDF + PDBQT), stereokimia dipulihkan |
| `docking/` | Seluruh pose hasil docking dan kontrol |

### `20_manuscript_support/` — bahan pendukung naskah

| Berkas | Isi |
|---|---|
| `20.1_na_sequence_identity.csv` | Kesamaan urutan N2 terhadap 40 galur H9N2 |
| `20.2_active_site_conservation.csv` | Kelestarian kantong per galur, 60 galur |
| `20.3_contact_residue_conservation.csv` | Kelestarian per residu kontak |
| `20.4_pubchem_verification.csv` | Verifikasi 104 identitas ke PubChem |
| `20.5_master_dataset_corrected.csv` | Berkas induk terkoreksi, menggantikan berkas 7–14 |
| `20.6_candidates_104_annotated.csv` | 104 kandidat + jangkauan model + status PubChem |
| `20.7_contact_residues_named.csv` | Residu kontak bernama per pasangan |
| `20.8_pains_screen.csv` | Penyaringan PAINS dan Brenk |

## Yang TIDAK ada di sini

Berkas 7–14 yang lama **sengaja tidak disertakan**. Identitas senyawa di dalamnya
salah pada sekitar 75% baris, dan menyimpannya di sini hanya akan membuat orang
mengutipnya kembali. Penggantinya `20.5_master_dataset_corrected.csv`.

Binary AutoDock Vina dan berkas pose sementara juga tidak disertakan — yang
pertama dapat diunduh dari rilis resminya, yang kedua dapat dibuat ulang.

## Menjalankan ulang

Perlu Python dengan `rdkit`, `pandas`, `networkx`, `biopython`, `meeko`, dan
binary AutoDock Vina 1.2.7. Setiap skrip berdiri sendiri dan mencetak apa yang
dikerjakannya.
