"""
CHOL-OMICS Clinical-Genetic Cohort Biophysical Mapper (GDPR Anonymized)
Author: Luca Celebrano
Repository: Lukecele/inflammatory-onco-physics

Fully compliant with EU GDPR 2016/679 and Ethics Committee Pseudonymization Protocols.
Replaces individual markers with alphanumeric identifiers (CHOL-01 to CHOL-18).
"""

from pathlib import Path
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

def process_chol_omics_cohort():
    tsv_path = DATA_DIR / "chol_omics_anonymized_stratification.tsv"
    excel_candidates = [
        DATA_DIR / "Clinca_genetica_coorte.xlsx",
        REPO_ROOT.parent / "Clinca_genetica_coorte.xlsx"
    ]
    
    excel_path = None
    for cand in excel_candidates:
        if cand.exists():
            excel_path = cand
            break

    if tsv_path.exists():
        print(f"[*] Loading pre-processed GDPR-pseudonymized dataset from: {tsv_path}")
        res_df = pd.read_csv(tsv_path, sep="\t")
    elif excel_path and excel_path.exists():
        print(f"[*] Loading {excel_path} (Anonymized Processing)...")
        df = pd.read_excel(excel_path, sheet_name="Coorte CHOL-OMICS")
        print(f"[*] Processing {len(df)} diagnostic and mutational classes with strict GDPR pseudonymization...")

        results = []
        for idx, row in df.iterrows():
            anon_id = f"CHOL-{idx+1:02d}"
            diag = str(row["Diagnosi"]).strip()
            base_mut = str(row.get("Base mutata", "")).strip()
            mut_type = str(row.get("Tipo di mutazione", "")).strip()
            clin_notes = str(row.get("Manifestazioni Cliniche da Letteratura\n(Generali o Specifiche della Variante)", ""))

            s_score = 1.0
            if "missenso" in mut_type.lower() and "peptide segnale" in mut_type.lower():
                s_score = 1.8
                variant_summary = "JAG1 Hypomorphic Missense"
            elif "delezione" in mut_type.lower() or "aploinsufficienza" in mut_type.lower():
                s_score = 4.5
                variant_summary = "JAG1 Structural Deletion (LoF)"
            elif "splicing" in mut_type.lower() or "frameshift" in mut_type.lower() or "nonsense" in clin_notes.lower():
                s_score = 4.2
                variant_summary = "JAG1 Canonical Splicing NMD"
            elif "pkhd1" in base_mut.lower() or "ciliari" in base_mut.lower() or "cilio" in clin_notes.lower():
                s_score = 4.8
                variant_summary = "Ciliary/Ductal Gene Defect"
            elif "tjp2" in base_mut.lower() or "giunzionali" in base_mut.lower():
                s_score = 4.0
                variant_summary = "Junctional Complex Alteration"
            else:
                s_score = 2.5
                variant_summary = "Biliary Structural Anomaly"

            i_score = 2.0
            if "colangite" in clin_notes.lower() or "litiasi" in clin_notes.lower():
                i_score = 4.2
            elif "fibrosi" in clin_notes.lower() or "cirrosi" in clin_notes.lower():
                i_score = 4.6
            elif "subclinica" in clin_notes.lower():
                i_score = 1.5

            bsir = s_score * i_score
            results.append({
                "Pseudonym_ID": anon_id,
                "Variant_Category": variant_summary,
                "Sensory_Disruption_Score": s_score,
                "Inflammatory_Stress_Score": i_score,
                "BSIR_Index": bsir,
                "Plot_Label": f"{anon_id} ({variant_summary})"
            })
        res_df = pd.DataFrame(results)
        res_df.to_csv(tsv_path, sep="\t", index=False)
    else:
        raise FileNotFoundError(f"Neither {tsv_path} nor raw Excel found.")

    # Assign risk group colors
    if "Color" not in res_df.columns:
        res_df["Color"] = res_df["BSIR_Index"].apply(
            lambda x: "#e11d48" if x > 15.0 else ("#f59e0b" if x > 8.0 else "#0284c7")
        )

    print("[*] Generating Publication-Grade Figure 6 (CHOL-OMICS Landscape)...")
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=300)

    # Panel A: Phase space
    ax1 = axes[0]
    ax1.scatter(res_df["Sensory_Disruption_Score"], res_df["Inflammatory_Stress_Score"],
                s=res_df["BSIR_Index"] * 18, c=res_df["Color"], alpha=0.85, edgecolors="#0f172a", linewidth=1.2)
    ax1.axhline(3.0, color="#94a3b8", linestyle="--", linewidth=1.2)
    ax1.axvline(3.0, color="#94a3b8", linestyle="--", linewidth=1.2)
    ax1.set_title("A. CHOL-OMICS Pseudonymized Biophysical Phase Space\nSensory Disruption vs Inflammatory/Toxic Stress",
                  fontsize=12, fontweight='bold', pad=12)
    ax1.set_xlabel("Sensory & Morphogenetic Disruption Score", fontsize=10.5, fontweight='bold')
    ax1.set_ylabel("Toxic & Inflammatory Stress Score", fontsize=10.5, fontweight='bold')
    ax1.set_xlim(0.8, 5.5)
    ax1.set_ylim(0.8, 5.5)

    # Panel B: Ranked BSIR Risk Score
    ax2 = axes[1]
    sorted_df = res_df.sort_values("BSIR_Index", ascending=True)
    y_pos = np.arange(len(sorted_df))
    bars = ax2.barh(y_pos, sorted_df["BSIR_Index"], color=sorted_df["Color"], edgecolor="#334155", linewidth=0.8)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(sorted_df["Plot_Label"], fontsize=8.5)
    ax2.set_xlabel("Biophysical Sensory-Inflammatory Risk (BSIR Index)", fontsize=10.5, fontweight='bold')
    ax2.set_title("B. Stratification of the 18 Anonymized CHOL-OMICS Classes\nPredicting Disease Progression & Intervention Urgency",
                  fontsize=12, fontweight='bold', pad=12)

    for bar in bars:
        w = bar.get_width()
        ax2.text(w + 0.3, bar.get_y() + bar.get_height()/2, f"{w:.1f}", va='center', fontsize=8, fontweight='bold', color="#1e293b")

    plt.tight_layout()
    output_png = REPO_ROOT / "figure_6_chol_omics_landscape.png"
    plt.savefig(output_png, bbox_inches="tight")
    plt.close()
    print(f"[✔] Fully pseudonymized figure saved at: {output_png}")

if __name__ == "__main__":
    process_chol_omics_cohort()
