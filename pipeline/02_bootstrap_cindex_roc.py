"""
Bootstrap Concordance Index & Time-Dependent ROC Engine
Author: Luca Celebrano
Repository: Lukecele/inflammatory-onco-physics

Performs bootstrap resamples to prove incremental prognostic value of HTBI,
and computes time-dependent ROC performance at 12, 24, 36, and 60 months.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from sklearn.metrics import roc_auc_score, roc_curve

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

def run_bootstrap_and_roc():
    print("[*] Loading Cohort Data for Bootstrap Iterations...")
    csv_path = DATA_DIR / "mskcc_2020_real_cohort.csv"
    df = pd.read_csv(csv_path)

    # Feature engineering: Decoupled from stage
    def map_histology_risk(val):
        if val in ["Solid", "Micropapillary"]:
            return 2
        elif val in ["Acinar", "Papillary"]:
            return 1
        else:
            return 0

    df["HISTO_SCORE"] = df["PREDOMINANT_HISTOLOGIC_SUBTYPE"].apply(map_histology_risk)
    df["IS_SMOKER"] = df["SMOKING_HISTORY"].apply(lambda x: 1 if "Ever" in str(x) else 0)
    df["IS_MALE"] = df["SEX"].apply(lambda x: 1 if str(x).lower().startswith("m") else 0)
    df["STAGE_NUM"] = df["PATHOLOGIC_STAGE"].astype(float)
    
    # Pure biological HTBI: Cellular polarity (Seed) + Airway inflammation (Soil)
    seed_weight = df["HISTO_SCORE"] / 2.0
    soil_stress = np.where(df["IS_SMOKER"] == 1, 1.0, 0.2)
    df["HTBI"] = 0.65 * seed_weight + 0.35 * soil_stress
    df["HTBI_HIGH"] = (df["HTBI"] >= df["HTBI"].median()).astype(int)

    data = df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM", "IS_MALE", "IS_SMOKER", "IS_EGFR_MUTANT", "HTBI_HIGH", "HTBI"]].dropna()
    data = data[data["RFS_MONTHS"] > 0].reset_index(drop=True)

    n_samples = len(data)
    n_bootstraps = 500
    print(f"[*] Running {n_bootstraps} Bootstrap Resamples on N={n_samples} patients...")

    c_clinical_boots = []
    c_integrated_boots = []
    delta_c_boots = []

    np.random.seed(42)

    for b in range(n_bootstraps):
        boot_idx = np.random.choice(n_samples, size=n_samples, replace=True)
        boot_df = data.iloc[boot_idx].reset_index(drop=True)

        # Fit Model 1: Clinical Stage alone
        cph1 = CoxPHFitter()
        cph1.fit(boot_df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM"]],
                  duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        c1 = cph1.concordance_index_

        # Fit Model 2: Integrated (Stage + HTBI)
        cph2 = CoxPHFitter()
        cph2.fit(boot_df[["RFS_MONTHS", "RECURRENCE_EVENT", "STAGE_NUM", "HTBI"]],
                  duration_col="RFS_MONTHS", event_col="RECURRENCE_EVENT")
        c2 = cph2.concordance_index_

        c_clinical_boots.append(c1)
        c_integrated_boots.append(c2)
        delta_c_boots.append(c2 - c1)

    mean_c1 = np.mean(c_clinical_boots)
    mean_c2 = np.mean(c_integrated_boots)
    mean_delta = np.mean(delta_c_boots)
    ci_low = np.percentile(delta_c_boots, 2.5)
    ci_high = np.percentile(delta_c_boots, 97.5)
    p_boot = np.mean(np.array(delta_c_boots) <= 0)

    print(f"\n[✔] Bootstrap Results ({n_bootstraps} iterations):")
    print(f" • TNM Staging C-Index:   {mean_c1:.3f} (SD = {np.std(c_clinical_boots):.3f})")
    print(f" • Integrated C-Index:    {mean_c2:.3f} (SD = {np.std(c_integrated_boots):.3f})")
    print(f" • Incremental Delta C:   +{mean_delta:.3f} (95% CI: [{ci_low:.3f}, {ci_high:.3f}])")
    print(f" • Bootstrap p-value:     p = {p_boot:.4f}")

    # Time-dependent ROC at landmark intervals
    roc_results = {}
    for horizon in [12, 24, 36, 60]:
        y_true = ((data["RFS_MONTHS"] <= horizon) & (data["RECURRENCE_EVENT"] == 1)).astype(int)
        
        # Risk score from Stage
        auc_stage = roc_auc_score(y_true, data["STAGE_NUM"])
        fpr_s, tpr_s, _ = roc_curve(y_true, data["STAGE_NUM"])

        # Combined risk score
        combined_score = 0.6 * data["STAGE_NUM"] + 0.4 * data["HTBI"]
        auc_comb = roc_auc_score(y_true, combined_score)
        fpr_c, tpr_c, _ = roc_curve(y_true, combined_score)

        roc_results[horizon] = {
            "auc_stage": auc_stage,
            "auc_combined": auc_comb,
            "fpr_stage": fpr_s,
            "tpr_stage": tpr_s,
            "fpr_comb": fpr_c,
            "tpr_comb": tpr_c,
            "n_events": sum(y_true),
            "n_total": len(y_true)
        }

    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Panel A: Bootstrap C-Index Distribution
    ax1 = axes[0]
    ax1.hist(c_clinical_boots, bins=35, density=True, alpha=0.6, color="#0284c7", label=f"TNM Stage Alone (Mean C = {mean_c1:.3f})")
    ax1.hist(c_integrated_boots, bins=35, density=True, alpha=0.6, color="#e11d48", label=f"TNM + HTBI Integrated (Mean C = {mean_c2:.3f})")
    ax1.axvline(mean_c1, color="#0369a1", linestyle="--", linewidth=2.0)
    ax1.axvline(mean_c2, color="#be123c", linestyle="--", linewidth=2.0)

    ax1.set_title(f"A. 500-Fold Bootstrap C-Index Validation\n(Δ C-index = +{mean_delta:.3f}, 95% CI: [{ci_low:.3f}, {ci_high:.3f}], p = {p_boot:.4f})",
                  fontsize=11.5, fontweight="bold", pad=12)
    ax1.set_xlabel("Harrell's Concordance Index (C-Index)", fontsize=10.5, fontweight="bold")
    ax1.set_ylabel("Bootstrap Probability Density", fontsize=10.5, fontweight="bold")
    ax1.legend(loc="upper left", fontsize=9, frameon=True)

    # Panel B: Time-Dependent ROC at 24 Months
    ax2 = axes[1]
    res_24 = roc_results[24]
    ax2.plot(res_24["fpr_comb"], res_24["tpr_comb"], color="#e11d48", linewidth=2.8,
             label=f"TNM + HTBI Integrated (AUC = {res_24['auc_combined']:.3f})")
    ax2.plot(res_24["fpr_stage"], res_24["tpr_stage"], color="#0284c7", linewidth=2.0, linestyle="--",
             label=f"TNM Stage Alone (AUC = {res_24['auc_stage']:.3f})")
    ax2.plot([0, 1], [0, 1], color="#94a3b8", linestyle=":", linewidth=1.5, label="Random Chance (AUC = 0.500)")

    ax2.set_title(f"B. 2-Year (24 Months) Time-Dependent ROC Curve\n(Δ AUC = +{res_24['auc_combined'] - res_24['auc_stage']:.3f}, n={res_24['n_events']} relapses / {res_24['n_total']})",
                  fontsize=11.5, fontweight="bold", pad=12)
    ax2.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=10.5, fontweight="bold")
    ax2.set_ylabel("Sensitivity (True Positive Rate)", fontsize=10.5, fontweight="bold")
    ax2.set_xlim(-0.02, 1.02)
    ax2.set_ylim(-0.02, 1.02)
    ax2.legend(loc="lower right", fontsize=9.5, frameon=True)

    plt.tight_layout()
    output_png = REPO_ROOT / "figure_5_bootstrap_roc.png"
    plt.savefig(output_png, bbox_inches="tight")
    plt.close()
    print(f"[✔] Bootstrap & ROC analysis complete. Figure saved at: {output_png}")

if __name__ == "__main__":
    run_bootstrap_and_roc()
