"""
Menyiapkan ligan untuk docking: 9 kandidat inhalasi (dua kandidat oral termasuk
di dalamnya) plus ligan ko-kristal tiap reseptor sebagai kontrol re-docking.

Tiga hal yang dikerjakan di sini dan tidak boleh dilewati:

1. STEREOKIMIA DIPULIHKAN. Berkas hasil skrining memakai SMILES kanonik tanpa
   stereo -- itu tidak masalah untuk fingerprint 2D, tetapi fatal untuk docking.
   SMILES isomerik diambil kembali dari peta identitas 18.1.

2. PROTONASI pH 7,4. Gugus karboksilat dideprotonasi (-COOH -> -COO^-). Seluruh
   kandidat teratas adalah asam fenolik, sehingga pada pH fisiologis mereka
   bermuatan negatif; mendok mereka dalam bentuk netral akan salah menilai
   interaksi elektrostatik, terutama pada situs neuraminidase yang kaya arginin
   dan pada Zn katalitik MMP. Fenol (pKa ~9-10) dibiarkan terprotonasi.

3. KONFORMER 3D dicari dengan ETKDGv3 (banyak konformer, dioptimasi MMFF94s,
   diambil yang berenergi terendah).
"""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
REPAIR = HERE.parent / "18_repair_and_filter"
OUT = HERE / "ligands"
OUT.mkdir(exist_ok=True)

N_CONF = 120   # batas atas; jumlah sebenarnya diskalakan ke ikatan rotasi
VINA_TORSION_LIMIT = 32
FLEX_WARN = 15
SEED = 0xF00D

# asam karboksilat -> karboksilat (pH 7,4)
ACID = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
# enol yang terkonjugasi ke karbonil lakton -- asam vinilog seperti askorbat
# (pKa1 4,17), jadi pada pH 7,4 ia berbentuk anion juga
ENOL_LACTONE = Chem.MolFromSmarts("[OX2H1]-[CX3]=[CX3]-[CX3](=[OX1])-[OX2]")


def deprotonate_acids(mol):
    """Deprotonasi gugus -COOH menjadi -COO^- di tempatnya.

    Dikerjakan dengan menyunting atom, BUKAN dengan RunReactants: reaksi SMARTS
    hanya mengembalikan atom yang dipetakan, sehingga sisa molekul ikut terbuang
    dan yang tersisa hanya fragmen [OH-].
    """
    rw = Chem.RWMol(mol)
    targets = [m[2] for m in mol.GetSubstructMatches(ACID)]
    targets += [m[0] for m in mol.GetSubstructMatches(ENOL_LACTONE)]
    for idx in set(targets):
        o = rw.GetAtomWithIdx(idx)
        o.SetFormalCharge(-1)
        o.SetNumExplicitHs(0)
        o.SetNoImplicit(True)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def embed_best(mol, n_conf=N_CONF):
    """Jumlah konformer diskalakan ke jumlah ikatan rotasi -- molekul kaku tidak
    perlu ratusan konformer, sedangkan gingerglikolipid (rantai panjang) perlu."""
    mol = Chem.AddHs(mol)
    from rdkit.Chem import rdMolDescriptors
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    # dibatasi 120: molekul dengan puluhan ikatan rotasi tidak akan terwakili
    # cukup oleh berapa pun jumlah konformer yang masih masuk akal waktunya,
    # dan ligan seperti itu memang di luar batas keandalan Vina (lihat
    # kolom dockable_vina pada 19.4)
    n = int(min(n_conf, 120, max(30, 12 * rot)))
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.useSmallRingTorsions = True
    params.pruneRmsThresh = 0.5
    params.numThreads = 0
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n, params=params)
    if not cids:
        return None, None
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s",
                                            maxIters=1000, numThreads=0)
    energies = [e for _, e in res]
    best = min(range(len(energies)), key=lambda i: energies[i])
    return mol, (cids[best], energies[best])


def write_pdbqt(mol, conf_id, path):
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        return False
    prep = MoleculePreparation()
    setups = prep.prepare(mol, conformer_id=conf_id)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        print(f"    ! meeko gagal: {err}")
        return False
    path.write_text(pdbqt)
    return True


def main():
    admet = pd.read_csv(REPAIR / "18.4_admet_annotated_corrected.csv")
    idmap = pd.read_csv(REPAIR / "18.1_identity_map_by_canonical_smiles.csv")

    sel = admet[admet.pass_inhalation_safety].copy()
    sel = sel.merge(idmap[["can", "Smiles_isomeric"]], on="can", how="left")
    print(f"kandidat disiapkan: {len(sel)}\n")

    rows = []
    for r in sel.sort_values("ndd_rank").itertuples():
        name = str(r.Metabolite).replace(" ", "_").replace("/", "-").replace("(", "").replace(")", "")
        mol = Chem.MolFromSmiles(str(r.Smiles_isomeric))
        if mol is None:
            print(f"  ! {r.Metabolite}: SMILES isomerik gagal diparsing")
            continue

        n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
        charged = deprotonate_acids(mol)
        formal_charge = Chem.GetFormalCharge(charged)

        mol3d, best = embed_best(charged)
        if mol3d is None:
            print(f"  ! {r.Metabolite}: gagal embedding 3D")
            continue
        conf_id, energy = best

        stem = f"{int(r.ndd_rank):03d}_{name}"
        sdf_path = OUT / f"{stem}.sdf"
        w = Chem.SDWriter(str(sdf_path))
        mol3d.SetProp("_Name", str(r.Metabolite))
        mol3d.SetProp("PubChem_ID", str(r.PubChem_ID))
        mol3d.SetProp("Organism", str(r.Organism))
        mol3d.SetProp("NDD_prob", f"{r.ndd_prob:.5f}")
        w.write(mol3d, confId=conf_id)
        w.close()

        pdbqt_ok = write_pdbqt(mol3d, conf_id, OUT / f"{stem}.pdbqt")

        n_rot = int(r.RotB)
        rows.append({
            "dockable_vina": n_rot <= VINA_TORSION_LIMIT,
            "flexibility_warning": n_rot > FLEX_WARN,
            "ndd_rank": r.ndd_rank, "Metabolite": r.Metabolite, "Organism": r.Organism,
            "PubChem_ID": r.PubChem_ID, "ndd_prob": round(r.ndd_prob, 5),
            "smiles_isomeric": r.Smiles_isomeric,
            "n_stereocenters": n_stereo, "formal_charge_pH7.4": formal_charge,
            "mmff_energy": round(energy, 2), "n_rotatable": r.RotB,
            "MW": r.MW, "cLogP": r.cLogP, "TPSA": r.TPSA,
            "sdf": sdf_path.name, "pdbqt": f"{stem}.pdbqt" if pdbqt_ok else "",
        })
        print(f"  {r.ndd_rank:>3}  {r.Metabolite:<22} stereocenter={n_stereo}  "
              f"muatan={formal_charge:+d}  rotB={int(r.RotB)}  E={energy:8.2f}")

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "19.4_ligands_prepared.csv", index=False)
    print(f"\n-> {HERE / '19.4_ligands_prepared.csv'}")
    print(f"-> SDF + PDBQT di {OUT}")


if __name__ == "__main__":
    main()
