import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths
MODEL_SAVE_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "trained_models.pkl")
OUTPUT_PLOT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "feature_importance.png")

# Selected human-readable explanations for common RDKit 2D descriptors to make output scientifically meaningful
DESCRIPTOR_EXPLANATIONS = {
    "MolLogP": "Octanol-Water Partition Coefficient (Lipophilicity)",
    "MolWt": "Molecular Weight",
    "HeavyAtomMolWt": "Heavy Atom Molecular Weight",
    "TPSA": "Topological Polar Surface Area (Polarity)",
    "NumHDonors": "Number of Hydrogen Bond Donors",
    "NumHAcceptors": "Number of Hydrogen Bond Acceptors",
    "NumRotatableBonds": "Number of Rotatable Bonds (Flexibility)",
    "FractionCSP3": "Fraction of C atoms that are SP3 hybridized (Complexity/3D shape)",
    "MolMR": "Molecular Refractivity (Polarizability)",
    "LabuteASA": "Labute Approximate Surface Area",
    "NOCount": "Number of Nitrogen and Oxygen atoms",
    "NHOHCount": "Number of NH and OH groups",
    "NumAliphaticRings": "Number of Aliphatic Rings",
    "NumAromaticRings": "Number of Aromatic Rings",
    "NumSaturatedRings": "Number of Saturated Rings",
    "BertzCT": "Bertz Complexity Index"
}

# Standard biological descriptions of selected MACCS keys related to tubulin activity
MACCS_EXPLANATIONS = {
    136: "Presence of Amide Group (O=C-N)",
    139: "Presence of Hydroxyl Group (OH)",
    120: "Presence of Heteroatom Ring (carbon-heteroatom ring)",
    84: "Presence of Carbonyl Group (C=O)",
    72: "Presence of Oxygen Atom bound to Ring",
    65: "Presence of Ester/Ether linkage",
    125: "Presence of Aromatic Ring",
    143: "Presence of Halogen Atom"
}

def get_feature_meaning(feat_name):
    """Return a scientific explanation for a feature name."""
    if feat_name in DESCRIPTOR_EXPLANATIONS:
        return DESCRIPTOR_EXPLANATIONS[feat_name]
    if feat_name.startswith("MACCS_"):
        try:
            bit = int(feat_name.split("_")[1])
            return MACCS_EXPLANATIONS.get(bit, f"MACCS Key structural fragment (Bit {bit})")
        except ValueError:
            pass
    return "Chemical/Structural property"

def main():
    if not os.path.exists(MODEL_SAVE_PATH):
        print(f"Error: {MODEL_SAVE_PATH} not found. Run train_ml.py first.")
        return
        
    with open(MODEL_SAVE_PATH, "rb") as f:
        saved_data = pickle.load(f)
        
    models = saved_data["models"]
    descriptor_names = saved_data["descriptor_names"]
    
    # Create the complete list of feature names
    maccs_names = [f"MACCS_{i}" for i in range(1, 167)]
    feature_names = maccs_names + descriptor_names
    
    # 1. Extract Random Forest Feature Importances
    rf_model = models.get("Random Forest")
    rf_importances = rf_model.feature_importances_
    
    rf_df = pd.DataFrame({
        "Feature": feature_names,
        "Importance": rf_importances
    }).sort_values(by="Importance", ascending=False)
    
    # 2. Extract Logistic Regression Coefficients (L2 Regularized weights)
    lr_model = models.get("Logistic Regression")
    lr_coefs = np.abs(lr_model.coef_[0]) # take absolute value for magnitude of impact
    
    lr_df = pd.DataFrame({
        "Feature": feature_names,
        "Importance": lr_coefs
    }).sort_values(by="Importance", ascending=False)
    
    # Select Top 12 features for plotting
    top_rf = rf_df.head(12).copy()
    top_lr = lr_df.head(12).copy()
    
    # Create publication-ready plot
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Left Plot: Random Forest
    y_pos1 = np.arange(len(top_rf))
    bars1 = ax1.barh(y_pos1, top_rf["Importance"], align='center', color='#2b5c8f', edgecolor='black', alpha=0.9)
    ax1.set_yticks(y_pos1)
    ax1.set_yticklabels(top_rf["Feature"], fontsize=11, fontweight='bold')
    ax1.invert_yaxis()  # top-down
    ax1.set_xlabel('Relative Feature Importance Score', fontsize=12, fontweight='bold')
    ax1.set_title('Random Forest - Top Descriptors', fontsize=14, fontweight='bold', pad=15)
    
    # Right Plot: Logistic Regression
    y_pos2 = np.arange(len(top_lr))
    bars2 = ax2.barh(y_pos2, top_lr["Importance"], align='center', color='#c2593f', edgecolor='black', alpha=0.9)
    ax2.set_yticks(y_pos2)
    ax2.set_yticklabels(top_lr["Feature"], fontsize=11, fontweight='bold')
    ax2.invert_yaxis()  # top-down
    ax2.set_xlabel('Absolute Coefficient Weight (L2-Norm Impact)', fontsize=12, fontweight='bold')
    ax2.set_title('Logistic Regression - Top Descriptors', fontsize=14, fontweight='bold', pad=15)
    
    # Add annotations and grid refinement
    for ax in [ax1, ax2]:
        ax.grid(True, linestyle='--', alpha=0.5, axis='x')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_PATH, dpi=300)
    print(f"\nFeature importance plot successfully generated and saved to: {OUTPUT_PLOT_PATH}")
    
    # Output a structured scientific table of top features
    print("\n" + "="*95)
    print("                    TOP SIGNIFICANT DESCRIPTORS (RANDOM FOREST)")
    print("="*95)
    print(f"{'Feature':<25} | {'Importance':<12} | {'Scientific Definition':<50}")
    print("-"*95)
    for idx, row in top_rf.iterrows():
        feat = row["Feature"]
        val = f"{row['Importance']:.4f}"
        meaning = get_feature_meaning(feat)
        print(f"{feat:<25} | {val:<12} | {meaning:<50}")
    print("="*95)

if __name__ == "__main__":
    main()
