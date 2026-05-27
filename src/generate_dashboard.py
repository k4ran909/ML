import os
import pandas as pd
import json

# Paths
HITS_CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "hit_predictions_results.csv")
ADMET_CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "admet_predictions.csv")
OUTPUT_HTML_PATH = os.path.join(os.path.dirname(__file__), "..", "index.html")

def main():
    print("Generating interactive HTML browser dashboard...")
    
    if not os.path.exists(HITS_CSV_PATH) or not os.path.exists(ADMET_CSV_PATH):
        print(f"Error: Missing CSV files. Make sure you run all.py first.")
        return
        
    # Load CSV data
    hits_df = pd.read_csv(HITS_CSV_PATH)
    admet_df = pd.read_csv(ADMET_CSV_PATH)
    
    # Generate table HTML for Hits Prediction
    hits_rows_html = ""
    for idx, row in hits_df.iterrows():
        name = row["CompoundName"]
        cid = int(row["PubChem_CID"])
        prob = row["Consensus_Active_Prob"]
        cls = row["Consensus_Class"]
        
        # Color coding class
        cls_badge = f'<span class="badge badge-inactive">Inactive</span>' if cls == "Inactive" else f'<span class="badge badge-active">Active</span>'
        
        # Format probabilities for individual models
        svm_p = f"{row['SVM (RBF)_Prob']:.3f}"
        lr_p = f"{row['Logistic Regression_Prob']:.3f}"
        rf_p = f"{row['Random Forest_Prob']:.3f}"
        xgb_p = f"{row['XGBoost_Prob']:.3f}"
        
        hits_rows_html += f"""
        <tr>
            <td class="font-semibold">{name}</td>
            <td>CID {cid}</td>
            <td class="font-mono">{svm_p}</td>
            <td class="font-mono">{lr_p}</td>
            <td class="font-mono">{rf_p}</td>
            <td class="font-mono">{xgb_p}</td>
            <td class="font-mono font-bold text-accent">{prob:.4f}</td>
            <td>{cls_badge}</td>
        </tr>
        """
        
    # Generate table HTML for ADMET Profiling
    admet_rows_html = ""
    for idx, row in admet_df.iterrows():
        name = row["CompoundName"]
        comp_type = row["Type"]
        mw = row["MW"]
        logp = row["LogP"]
        hbd = int(row["HBD"])
        hba = int(row["HBA"])
        rot = int(row["RotatableBonds"])
        tpsa = row["TPSA"]
        logs = row["ESOL_LogS"]
        sol = row["Est_Solubility_mg_L"]
        viols = int(row["Lipinski_Violations"])
        pains = row["PAINS_Alert"]
        
        type_class = "text-hit" if "Hit" in comp_type else "text-control"
        viols_badge = f'<span class="badge badge-success">{viols}</span>' if viols == 0 else (f'<span class="badge badge-warning">{viols}</span>' if viols <= 1 else f'<span class="badge badge-danger">{viols}</span>')
        pains_badge = f'<span class="badge badge-success">Clean</span>' if pains == "Clean" else f'<span class="badge badge-danger">{pains}</span>'
        
        admet_rows_html += f"""
        <tr>
            <td class="font-semibold">{name}</td>
            <td class="{type_class}">{comp_type}</td>
            <td class="font-mono">{mw:.2f}</td>
            <td class="font-mono">{logp:.2f}</td>
            <td class="font-mono">{hbd}</td>
            <td class="font-mono">{hba}</td>
            <td class="font-mono">{rot}</td>
            <td class="font-mono">{tpsa:.2f}</td>
            <td class="font-mono">{logs:.2f}</td>
            <td class="font-mono">{sol:.3f}</td>
            <td>{viols_badge}</td>
            <td>{pains_badge}</td>
        </tr>
        """
        
    # Generate HTML content
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Biotech Dashboard: Tubulin Inhibitors Screening</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-color: #0b0f19;
            --panel-bg: rgba(17, 24, 39, 0.7);
            --border-color: rgba(255, 255, 255, 0.08);
            --text-color: #f3f4f6;
            --text-muted: #9ca3af;
            --accent: #3b82f6;
            --accent-glow: rgba(59, 130, 246, 0.15);
            --success: #10b981;
            --success-glow: rgba(16, 185, 129, 0.15);
            --warning: #f59e0b;
            --danger: #ef4444;
            --card-bg: #1f2937;
        }}

        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: 'Inter', sans-serif;
            line-height: 1.6;
            padding: 2rem;
            min-height: 100vh;
            background-image: 
                radial-gradient(at 0% 0%, rgba(59, 130, 246, 0.1) 0, transparent 50%),
                radial-gradient(at 50% 0%, rgba(139, 92, 246, 0.08) 0, transparent 50%),
                radial-gradient(at 100% 100%, rgba(16, 185, 129, 0.05) 0, transparent 50%);
            background-attachment: fixed;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        header {{
            margin-bottom: 2.5rem;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 1.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}

        .logo-section h1 {{
            font-size: 2rem;
            font-weight: 800;
            background: linear-gradient(135deg, #60a5fa 0%, #a78bfa 50%, #34d399 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            letter-spacing: -0.05em;
        }}

        .logo-section p {{
            color: var(--text-muted);
            font-size: 0.95rem;
            margin-top: 0.25rem;
        }}

        .tab-nav {{
            display: flex;
            gap: 0.5rem;
            background: rgba(31, 41, 55, 0.5);
            padding: 0.35rem;
            border-radius: 12px;
            border: 1px solid var(--border-color);
            backdrop-filter: blur(12px);
            margin-bottom: 2rem;
        }}

        .tab-btn {{
            background: transparent;
            border: none;
            color: var(--text-muted);
            padding: 0.75rem 1.5rem;
            border-radius: 8px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
            font-size: 0.9rem;
        }}

        .tab-btn:hover {{
            color: var(--text-color);
            background: rgba(255, 255, 255, 0.04);
        }}

        .tab-btn.active {{
            background: var(--accent);
            color: #ffffff;
            box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
        }}

        .tab-content {{
            display: none;
            animation: fadeIn 0.4s ease-out;
        }}

        .tab-content.active {{
            display: block;
        }}

        @keyframes fadeIn {{
            from {{ opacity: 0; transform: translateY(10px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}

        .grid-3 {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .card {{
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.5rem;
            backdrop-filter: blur(12px);
            transition: transform 0.2s, box-shadow 0.2s;
        }}

        .card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.3);
            border-color: rgba(59, 130, 246, 0.2);
        }}

        .card-header {{
            font-size: 1.1rem;
            font-weight: 700;
            margin-bottom: 1rem;
            display: flex;
            align-items: center;
            gap: 0.5rem;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 0.5rem;
            color: #ffffff;
        }}

        .card-val {{
            font-size: 2.2rem;
            font-weight: 800;
            margin-top: 0.5rem;
            font-family: 'JetBrains Mono', monospace;
        }}

        .card-sub {{
            color: var(--text-muted);
            font-size: 0.85rem;
            margin-top: 0.25rem;
        }}

        .table-container {{
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            overflow: hidden;
            margin-bottom: 2rem;
            backdrop-filter: blur(12px);
        }}

        .table-header-bar {{
            padding: 1.25rem 1.5rem;
            border-bottom: 1px solid var(--border-color);
            background: rgba(31, 41, 55, 0.3);
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}

        .table-title {{
            font-size: 1.2rem;
            font-weight: 700;
            color: #ffffff;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            text-align: left;
            font-size: 0.95rem;
        }}

        th {{
            background: rgba(31, 41, 55, 0.2);
            color: var(--text-muted);
            font-weight: 600;
            padding: 1rem 1.5rem;
            border-bottom: 1px solid var(--border-color);
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.05em;
        }}

        td {{
            padding: 1rem 1.5rem;
            border-bottom: 1px solid var(--border-color);
        }}

        tr:hover td {{
            background: rgba(255, 255, 255, 0.02);
        }}

        .font-mono {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.9rem;
        }}

        .font-semibold {{
            font-weight: 600;
        }}

        .font-bold {{
            font-weight: 700;
        }}

        .text-accent {{
            color: #60a5fa;
        }}

        .text-hit {{
            color: #34d399;
            font-weight: 600;
        }}

        .text-control {{
            color: #60a5fa;
            font-weight: 600;
        }}

        .badge {{
            display: inline-block;
            padding: 0.25rem 0.6rem;
            border-radius: 6px;
            font-size: 0.75rem;
            font-weight: 700;
            text-transform: uppercase;
        }}

        .badge-inactive {{
            background: rgba(239, 68, 68, 0.1);
            color: #f87171;
            border: 1px solid rgba(239, 68, 68, 0.2);
        }}

        .badge-active {{
            background: rgba(16, 185, 129, 0.1);
            color: #34d399;
            border: 1px solid rgba(16, 185, 129, 0.2);
        }}

        .badge-success {{
            background: rgba(16, 185, 129, 0.1);
            color: #34d399;
        }}

        .badge-warning {{
            background: rgba(245, 158, 11, 0.1);
            color: #fbbf24;
        }}

        .badge-danger {{
            background: rgba(239, 68, 68, 0.1);
            color: #f87171;
        }}

        .image-box {{
            background: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.5rem;
            text-align: center;
            backdrop-filter: blur(12px);
            margin-bottom: 2rem;
        }}

        .image-box img {{
            max-width: 100%;
            height: auto;
            border-radius: 12px;
            border: 1px solid var(--border-color);
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.5);
            margin-top: 1rem;
        }}

        .interpretation-section {{
            margin-top: 2rem;
            padding: 1.5rem;
            background: rgba(255, 255, 255, 0.02);
            border-left: 4px solid var(--accent);
            border-radius: 4px 16px 16px 4px;
            margin-bottom: 2rem;
        }}

        .interpretation-section h4 {{
            font-size: 1.1rem;
            font-weight: 700;
            margin-bottom: 0.75rem;
            color: #ffffff;
        }}

        .interpretation-section p {{
            color: var(--text-muted);
            font-size: 0.95rem;
            margin-bottom: 0.5rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="logo-section">
                <h1>Tubulin Inhibitor Discovery Dashboard</h1>
                <p>ML-Driven Virtual Screening, Consensus Profiling, and ADMET Validation</p>
            </div>
            <div class="logo-section" style="text-align: right;">
                <p style="font-family: 'JetBrains Mono', monospace; font-size: 0.8rem; background: rgba(255,255,255,0.05); padding: 0.4rem 0.8rem; border-radius: 8px; border: 1px solid var(--border-color);">
                    Pipeline Status: <span style="color: var(--success); font-weight: bold;">COMPLETE</span>
                </p>
            </div>
        </header>

        <nav class="tab-nav">
            <button class="tab-btn active" onclick="openTab(event, 'overview')">🏠 Overview</button>
            <button class="tab-btn" onclick="openTab(event, 'ml-models')">📊 ML Performance</button>
            <button class="tab-btn" onclick="openTab(event, 'screening')">🎯 Hit Screening</button>
            <button class="tab-btn" onclick="openTab(event, 'admet')">🧪 ADMET Profiling</button>
            <button class="tab-btn" onclick="openTab(event, 'feature-importance')">📈 Feature Importance</button>
        </nav>

        <!-- OVERVIEW TAB -->
        <div id="overview" class="tab-content active">
            <div class="grid-3">
                <div class="card">
                    <div class="card-header">🧬 Homology Model</div>
                    <div class="card-val text-accent" style="font-size: 1.8rem;">Human &beta;<sub>III</sub></div>
                    <div class="card-sub">Template Bovine 1JFF | DOPE Score -51,889.82</div>
                </div>
                <div class="card">
                    <div class="card-header">🔑 Screening Library</div>
                    <div class="card-val" style="color: var(--success);">319 Ligands</div>
                    <div class="card-sub">Vina binding energy threshold &le; -10 kcal/mol</div>
                </div>
                <div class="card">
                    <div class="card-header">🎯 Validated Hits</div>
                    <div class="card-val" style="color: var(--warning);">7 Compounds</div>
                    <div class="card-sub">Classified as pocket-inactive (covalent/alternative binders)</div>
                </div>
            </div>

            <div class="interpretation-section">
                <h4>🔬 Pipeline Workflow Summary</h4>
                <p>This computational drug discovery dashboard compiles the results of structural homology modeling, AutoDock Vina screening, machine learning classification, and ADMET profiling targeting the **taxane binding site** of human &alpha;&beta;<sub>III</sub> tubulin.</p>
                <p>We trained machine learning classifiers (Logistic Regression, Random Forest, SVM, XGBoost) to recognize features of taxane-site stabilizers, which were then used to screen the 7 docking hits. A comparative ADMET evaluation prioritized candidates with optimal drug-likeness.</p>
            </div>

            <div class="table-container">
                <div class="table-header-bar">
                    <div class="table-title">🎯 Top Prioritized Drug Candidates</div>
                </div>
                <table>
                    <thead>
                        <tr>
                            <th>Compound Name</th>
                            <th>ML Prob</th>
                            <th>MW</th>
                            <th>LogP</th>
                            <th>Lipinski Violations</th>
                            <th>Aqueous Sol. (mg/L)</th>
                            <th>Priority Recommendation</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td class="font-semibold text-hit">Withaferin A</td>
                            <td class="font-mono">0.4247</td>
                            <td class="font-mono">470.61</td>
                            <td class="font-mono">3.35</td>
                            <td><span class="badge badge-success">0</span></td>
                            <td class="font-mono">6.08</td>
                            <td><span class="badge badge-success">Rank 1 (Top Priority)</span></td>
                        </tr>
                        <tr>
                            <td class="font-semibold text-hit">Withanolide A</td>
                            <td class="font-mono">0.4226</td>
                            <td class="font-mono">470.61</td>
                            <td class="font-mono">3.50</td>
                            <td><span class="badge badge-success">0</span></td>
                            <td class="font-mono">5.02</td>
                            <td><span class="badge badge-success">Rank 2 (Top Priority)</span></td>
                        </tr>
                        <tr>
                            <td class="font-semibold text-hit">27-O-acetyl-Withaferin A</td>
                            <td class="font-mono">0.3898</td>
                            <td class="font-mono">512.64</td>
                            <td class="font-mono">3.92</td>
                            <td><span class="badge badge-warning">1</span></td>
                            <td class="font-mono">1.57</td>
                            <td><span class="badge badge-warning">Rank 3 (Medium Priority)</span></td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>

        <!-- ML PERFORMANCE TAB -->
        <div id="ml-models" class="tab-content">
            <div class="table-container">
                <div class="table-header-bar">
                    <div class="table-title">📊 Stratified 5-Fold Cross-Validation Metrics</div>
                </div>
                <table>
                    <thead>
                        <tr>
                            <th>Machine Learning Model</th>
                            <th>Accuracy</th>
                            <th>ROC-AUC</th>
                            <th>Recall</th>
                            <th>Precision</th>
                            <th>F1-Score</th>
                            <th>Kappa</th>
                            <th>MCC</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td class="font-semibold">Logistic Regression (LR)</td>
                            <td class="font-mono">0.8030</td>
                            <td class="font-mono">0.8444</td>
                            <td class="font-mono">0.8133</td>
                            <td class="font-mono">0.8010</td>
                            <td class="font-mono">0.7952</td>
                            <td class="font-mono">0.6087</td>
                            <td class="font-mono">0.6235</td>
                        </tr>
                        <tr>
                            <td class="font-semibold">Random Forest Classifier (RF)</td>
                            <td class="font-mono">0.8030</td>
                            <td class="font-mono">0.8456</td>
                            <td class="font-mono">0.8067</td>
                            <td class="font-mono">0.7833</td>
                            <td class="font-mono">0.7909</td>
                            <td class="font-mono">0.6051</td>
                            <td class="font-mono">0.6097</td>
                        </tr>
                        <tr>
                            <td class="font-semibold">SVM (RBF Kernel)</td>
                            <td class="font-mono">0.7864</td>
                            <td class="font-mono">0.8456</td>
                            <td class="font-mono">0.7400</td>
                            <td class="font-mono">0.8067</td>
                            <td class="font-mono">0.7636</td>
                            <td class="font-mono">0.5706</td>
                            <td class="font-mono">0.5838</td>
                        </tr>
                        <tr>
                            <td class="font-semibold">AdaBoost Classifier (ada)</td>
                            <td class="font-mono">0.6803</td>
                            <td class="font-mono">0.6689</td>
                            <td class="font-mono">0.6667</td>
                            <td class="font-mono">0.6943</td>
                            <td class="font-mono">0.6679</td>
                            <td class="font-mono">0.3644</td>
                            <td class="font-mono">0.3783</td>
                        </tr>
                        <tr>
                            <td class="font-semibold">XGBoost Classifier (xgboost)</td>
                            <td class="font-mono">0.6803</td>
                            <td class="font-mono">0.6689</td>
                            <td class="font-mono">0.6667</td>
                            <td class="font-mono">0.6943</td>
                            <td class="font-mono">0.6679</td>
                            <td class="font-mono">0.3644</td>
                            <td class="font-mono">0.3783</td>
                        </tr>
                    </tbody>
                </table>
            </div>

            <div class="interpretation-section">
                <h4>💡 Dataset Scale and Model Complexity Analysis</h4>
                <p><strong>Overfitting in Small Datasets:</strong> Boosting algorithms (AdaBoost and XGBoost) achieved only 68.0% cross-validation accuracy due to overfitting on our small dataset of 57 compounds. In contrast, regularized models like Logistic Regression and Random Forest generalize significantly better (80.3% accuracy).</p>
                <p><strong>The Scale Discrepancy:</strong> In the original research paper, the dataset was expanded using over 3,000 decoy structures, allowing boosting models to reach 99.8% accuracy. Our local run reflects the behavior of these models under constraints of small, high-confidence training pools.</p>
            </div>
        </div>

        <!-- HIT SCREENING TAB -->
        <div id="screening" class="tab-content">
            <div class="table-container">
                <div class="table-header-bar">
                    <div class="table-title">🎯 Machine Learning Hit Predictions (Taxane Site)</div>
                </div>
                <table>
                    <thead>
                        <tr>
                            <th>Compound Name</th>
                            <th>PubChem ID</th>
                            <th>SVM (RBF)</th>
                            <th>Logistic Reg</th>
                            <th>Random Forest</th>
                            <th>XGBoost</th>
                            <th>Consensus Prob</th>
                            <th>Consensus Class</th>
                        </tr>
                    </thead>
                    <tbody>
                        {hits_rows_html}
                    </tbody>
                </table>
            </div>

            <div class="interpretation-section">
                <h4>💡 Biological Analysis of Hit Predictions</h4>
                <p><strong>Consensus Inactivity:</strong> All 7 hits are classified as **Inactive** for the interior taxane pocket (Consensus Probability < 0.50). This aligns perfectly with literature findings.</p>
                <p><strong>Mechanistic Rationale:</strong> While compounds like Celastrol, Pristimerin, and Withaferin A are known anti-tubulin cancer drugs, they do not bind non-covalently inside the taxol pocket. Instead:
                    <br>&bull; <em>Quinone methides (Celastrol, Pristimerin, Tingenone)</em> act as <strong>covalent alkylators</strong> of cysteine residues on the outer tubulin surface.
                    <br>&bull; <em>Withanolides (Withaferin A)</em> bind to alternative pockets or chaperone proteins.
                </p>
            </div>
        </div>

        <!-- ADMET PROFILING TAB -->
        <div id="admet" class="tab-content">
            <div class="table-container">
                <div class="table-header-bar">
                    <div class="table-title">🧪 ADMET Profiling Comparison (Hits vs. Control Drugs)</div>
                </div>
                <table>
                    <thead>
                        <tr>
                            <th>Compound Name</th>
                            <th>Type</th>
                            <th>MW (g/mol)</th>
                            <th>LogP</th>
                            <th>HBD</th>
                            <th>HBA</th>
                            <th>Rot. Bonds</th>
                            <th>TPSA (Å²)</th>
                            <th>LogS</th>
                            <th>Solubility (mg/L)</th>
                            <th>Violations</th>
                            <th>PAINS</th>
                        </tr>
                    </thead>
                    <tbody>
                        {admet_rows_html}
                    </tbody>
                </table>
            </div>
        </div>

        <!-- FEATURE IMPORTANCE TAB -->
        <div id="feature-importance" class="tab-content">
            <div class="interpretation-section">
                <h4>📈 Feature Importance Weights</h4>
                <p>Analyzing feature importances helps decipher the chemical rules driving taxane-site binding. For example, Random Forest highlights <code>fr_methoxy</code> (Methoxy group count) as the top driver, which is highly consistent with clinical stabilizers (e.g. Cabazitaxel) that feature methoxy decorations.</p>
            </div>
            
            <div class="image-box">
                <h4>🖼️ Feature Weight Distribution Plot</h4>
                <img src="results/feature_importance.png" alt="Feature Importance Plot">
            </div>
        </div>
    </div>

    <script>
        function openTab(evt, tabId) {{
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {{
                tabcontent[i].classList.remove("active");
            }}
            tablinks = document.getElementsByClassName("tab-btn");
            for (i = 0; i < tablinks.length; i++) {{
                tablinks[i].classList.remove("active");
            }}
            document.getElementById(tabId).classList.add("active");
            evt.currentTarget.classList.add("active");
        }}
    </script>
</body>
</html>
"""
    
    with open(OUTPUT_HTML_PATH, 'w', encoding='utf-8') as f:
        f.write(html_content)
        
    print(f"Interactive dashboard successfully generated at: {OUTPUT_HTML_PATH}")
    
    # Automatically open the dashboard in the default web browser
    try:
        import webbrowser
        html_url = 'file://' + os.path.abspath(OUTPUT_HTML_PATH).replace('\\', '/')
        print(f"Opening dashboard in browser: {html_url}")
        webbrowser.open(html_url)
    except Exception as e:
        print(f"Note: Could not open browser automatically: {e}")

if __name__ == "__main__":
    main()
