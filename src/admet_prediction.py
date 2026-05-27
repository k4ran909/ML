import os
import urllib.parse
import requests
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import FilterCatalog

# Paths
HIT_PREDICTIONS_CSV = os.path.join(os.path.dirname(__file__), "..", "results", "hit_predictions_results.csv")
OUTPUT_ADMET_CSV = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "admet_predictions.csv")

# Active controls for comparison
CONTROLS = [
    {"name": "Paclitaxel", "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"},
    {"name": "Docetaxel", "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O"},
    {"name": "Cabazitaxel", "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)OC)C)OC"}
]

def calculate_esol_logs(mol):
    """
    Calculate ESOL LogS (Aqueous Solubility) based on Delaney's model:
    LogS = 0.16 - 0.63*LogP - 0.0062*MW - 0.0062*RotatableBonds + 0.05*AromaticProportion
    """
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rtb = Descriptors.NumRotatableBonds(mol)
    
    # Calculate aromatic proportion
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    aromatic_prop = len(aromatic_atoms) / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0
    
    logs = 0.16 - (0.63 * logp) - (0.0062 * mw) - (0.0062 * rtb) + (0.05 * aromatic_prop)
    return logs

def check_pains_filter(mol):
    """Check if the molecule triggers any PAINS toxicity alerts."""
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    
    matches = catalog.GetMatches(mol)
    if matches:
        # Return the description of the first matching rule
        return [m.GetDescription() for m in matches]
    return []

def main():
    if not os.path.exists(HIT_PREDICTIONS_CSV):
        print(f"Error: {HIT_PREDICTIONS_CSV} not found. Run resolve_and_predict.py first.")
        return
        
    hits_df = pd.read_csv(HIT_PREDICTIONS_CSV)
    
    all_compounds = []
    
    # Add hits
    for idx, row in hits_df.iterrows():
        all_compounds.append({
            "name": row["CompoundName"],
            "smiles": row["SMILES"],
            "type": "Hit Candidate",
            "ml_active_prob": row["Consensus_Active_Prob"]
        })
        
    # Add controls
    for ctrl in CONTROLS:
        all_compounds.append({
            "name": ctrl["name"],
            "smiles": ctrl["smiles"],
            "type": "Control Drug (Active)",
            "ml_active_prob": 1.0  # reference
        })
        
    admet_data = []
    
    # Initialize PAINS catalog
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    
    for comp in all_compounds:
        name = comp["name"]
        smiles = comp["smiles"]
        comp_type = comp["type"]
        prob = comp["ml_active_prob"]
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Invalid SMILES for '{name}' ignored.")
            continue
            
        # Lipinski Parameters
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rtb = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Calculate Lipinski Violations
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        # Solubility
        logs = calculate_esol_logs(mol)
        solubility_mg_l = (10 ** logs) * mw * 1000.0 # convert mol/L to mg/L
        
        # PAINS Alert
        pains_matches = catalog.GetMatches(mol)
        pains_alert = "Clean" if not pains_matches else ", ".join([m.GetDescription() for m in pains_matches])
        
        admet_data.append({
            "CompoundName": name,
            "Type": comp_type,
            "ML_Active_Prob": prob,
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "RotatableBonds": rtb,
            "TPSA": round(tpsa, 2),
            "ESOL_LogS": round(logs, 2),
            "Est_Solubility_mg_L": round(solubility_mg_l, 2),
            "Lipinski_Violations": violations,
            "PAINS_Alert": pains_alert
        })
        
    df = pd.DataFrame(admet_data)
    df.to_csv(OUTPUT_ADMET_CSV, index=False)
    print(f"\nADMET predictions completed and saved to: {OUTPUT_ADMET_CSV}")
    
    # Print formatted output comparing hits with active drugs
    print("\n" + "="*95)
    print("                              ADMET PREDICTIONS REPORT")
    print("="*95)
    cols = ["CompoundName", "Type", "MW", "LogP", "Lipinski_Violations", "ESOL_LogS", "TPSA", "PAINS_Alert"]
    print(df[cols].to_string(index=False))
    print("="*95)

if __name__ == "__main__":
    main()
