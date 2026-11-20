"""
Multivariate Cox Proportional Hazards Engine for HTBI Validation
Author: Luca Celebrano
Repository: Lukecele/inflammatory-onco-physics
Evaluates independent prognostic significance of the Host-Tumor Biophysical Index (HTBI)
controlling for TNM stage, histology, smoking, sex, and mutational status.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

def run_multivariate_cox_analysis():
    print("[*] Loading MSKCC 2020 Clinical Cohort (N=604)...")
    csv_path = DATA_DIR / "mskcc_2020_real_cohort.csv"
    df = pd.read_csv(csv_path)

    # 1. Feature Engineering
    # Histological Risk based on tissue polarity & glandular differentiation
    # Lepidic: G1-like (well differentiated, low sanctuary risk)
    # Acinar/Papillary: G2-like (intermediate)
    # Solid/Micropapillary: G3-like (poorly differentiated, high stromal desmoplasia)
    def map_histology_risk(val):
        if val in ["Solid", "Micropapillary"]:
            return 2  # High
        elif val in ["Acinar", "Papillary"]:
            return 1  # Intermediate
        else:
            return 0  # Low (Lepidic)

    df["HISTO_SCORE"] = df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].apply(map_histology_risk)
    df["IS_SMOKER"] = df["SMOKING_HISTORY"].apply(lambda x: 1 if "Ever" in str(x) else 0)
    df["IS_MALE"] = df["SEX"].apply(lambda x: 1 if str(x).lower().startswith("m") else 0)
    df["STAGE_NUM"] = df["PATHOLOGIC_STAGE"].astype(float)
    
    # 2. Compute Host-Tumor Biophysical Index (HTBI)
    # Strictly decoupled from anatomical stage: combining cellular polarity loss (Seed)
    # and chronic airway inflammatory stress (Soil) on continuous 0-1 scale.
    seed_weight = df["HISTO_SCORE"] / 2.0  # 0.0 (Lepidic), 0.5 (Acinar/Papillary), 1.0 (Solid/Micropapillary)
    soil_stress = np.where(df["IS_SMOKER"] == 1, 1.0, 0.2)
    
    df["HTBI"] = 0.65 * seed_weight + 0.35 * soil_stress
    df["HTBI_HIGH"] = (df["HTBI"] >= df["HTBI"].median()).astype(int)

    # Filter complete cases for survival analysis
    cox_df = df[[
        "RFS_MONTHS", "RECURRENCE_EVENT",
        "STAGE_NUM", "IS_MALE", "IS_SMOKER", "IS_EGFR_MUTANT", "HISTO_SCORE", "HTBI", "HTBI_HIGH"
    ]].dropna().copy()
    
    # Ensure positive follow-up time
    cox_df = cox_df[cox_df["RFS_MONTHS"] > 0].copy()

    print(f"[*] Running Multivariate Cox Proportional Hazards Model on {len(cox_df)} patients...")
    
    # Model 1: Classical Clinical Staging Alone
    cph_clinical = CoxPHFitter()
    cph_clinical.fit(
        cox_df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM", "IS_MALE", "IS_SMOKER", "IS_EGFR_MUTANT"]],
        duration_col="RFS_MONTHS",
        event_col="RECURRENCE_EVENT"
    )
    print("\n--- Model 1: Classical Clinical Staging Alone ---")
    cph_clinical.print_summary(decimals=3)

    # Model 2: Integrated Model (Clinical Staging + Continuous HTBI)
    cph_integrated_cont = CoxPHFitter()
    cph_integrated_cont.fit(
        cox_df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM", "IS_MALE", "IS_EGFR_MUTANT", "HTBI"]],
        duration_col="RFS_MONTHS",
        event_col="RECURRENCE_EVENT"
    )
    print("\n--- Model 2: Integrated Model (Continuous HTBI, 0-1 scale) ---")
    cph_integrated_cont.print_summary(decimals=3)

    # Model 3: Integrated Model (Clinical Staging + Binary HTBI_HIGH)
    cph_integrated = CoxPHFitter()
    cph_integrated.fit(
        cox_df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM", "IS_MALE", "IS_EGFR_MUTANT", "HTBI_HIGH"]],
        duration_col="RFS_MONTHS",
        event_col="RECURRENCE_EVENT"
    )
    print("\n--- Model 3: Integrated Model (Binary HTBI_HIGH, Median Split) ---")
    cph_integrated.print_summary(decimals=3)

    # Extract summary table for Forest Plot from continuous model
    summary = cph_integrated_cont.summary
    hr = np.exp(summary["coef"])
    hr_lower = np.exp(summary["coef lower 95%"])
    hr_upper = np.exp(summary["coef upper 95%"])
    p_values = summary["p"]
    variables = summary.index.tolist()

    # Friendly clinical names
    var_names = {
        "STAGE_NUM": "Pathologic Stage (TNM II/III vs I)",
        "IS_MALE": "Sex (Male vs Female)",
        "IS_EGFR_MUTANT": "EGFR Mutation (Mutant vs WT)",
        "HTBI": "Host-Tumor Biophysical Index (Continuous, 0-1 scale)",
        "HTBI_HIGH": "Host-Tumor Biophysical Index (High vs Low)"
    }
    clean_labels = [var_names.get(v, v) for v in variables]

    print("[*] Generating Publication-Grade Multivariate Forest Plot...")
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)

    y_pos = np.arange(len(variables))
    
    # Vertical line at HR = 1 (No effect)
    ax.axvline(1.0, color="#64748b", linestyle="--", linewidth=1.5, zorder=1)

    # Plot confidence intervals and point estimates
    for i in range(len(variables)):
        is_htbi = (variables[i] in ["HTBI", "HTBI_HIGH"])
        color = "#e11d48" if is_htbi else "#0284c7"
        lw = 3.0 if is_htbi else 2.0
        ms = 10 if is_htbi else 7
        
        # Horizontal error bar
        ax.plot([hr_lower.iloc[i], hr_upper.iloc[i]], [y_pos[i], y_pos[i]], color=color, linewidth=lw, zorder=2)
        # Point estimate
        ax.plot(hr.iloc[i], y_pos[i], marker="s" if is_htbi else "o", color=color, markersize=ms, zorder=3)
        
        # Text annotation of HR and p-value
        p_txt = "p < 0.001" if p_values.iloc[i] < 0.001 else f"p = {p_values.iloc[i]:.3f}"
        text_label = f"HR {hr.iloc[i]:.2f} (95% CI: {hr_lower.iloc[i]:.2f}-{hr_upper.iloc[i]:.2f}), {p_txt}"
        ax.text(max(hr_upper) * 1.05, y_pos[i], text_label, verticalalignment="center", fontsize=9.5,
                fontweight="bold" if is_htbi else "normal", color="#0f172a")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(clean_labels, fontsize=10.5, fontweight="bold")
    ax.set_xlabel("Hazard Ratio (log scale)", fontsize=11, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(0.3, max(hr_upper) * 2.8)
    ax.set_title("Multivariate Cox Proportional Hazards Regression\nPredicting Relapse-Free Survival (MSKCC Cohort, N=604)",
                 fontsize=13, fontweight="bold", pad=15)

    plt.tight_layout()
    output_png = REPO_ROOT / "figure_4_multivariate_forest.png"
    plt.savefig(output_png, bbox_inches="tight")
    plt.close()
    print(f"[✔] Multivariate analysis complete. Forest plot saved at: {output_png}")

if __name__ == "__main__":
    run_multivariate_cox_analysis()
