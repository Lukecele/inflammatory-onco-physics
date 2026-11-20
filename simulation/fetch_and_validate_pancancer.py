"""
Pan-Cancer Clinical Cohort Validation Pipeline
Authors: Luca Celebrano, Antigravity AI
Repository: Lukecele/inflammatory-onco-physics

Validates the Host-Tumor Biophysical Index (HTBI) across 4 major solid epithelial malignancies:
1. Lung Adenocarcinoma (MSKCC 2020 Discovery Cohort, N=604)
2. Colorectal Adenocarcinoma (TCGA COADREAD PanCan Atlas, N=592 / 566 with OS)
3. Pancreatic Ductal Adenocarcinoma (TCGA PAAD PanCan Atlas, N=177 / 176 with OS)
4. Cholangiocarcinoma (TCGA CHOL PanCan Atlas, N=36)
"""

import os
import json
import urllib.request
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

REPO_ROOT = Path(__file__).resolve().parent.parent

def fetch_or_load_tcga(study_id: str, local_tsv: str) -> pd.DataFrame:
    """Load local cached RNA-seq cohort or fetch from cBioPortal API."""
    local_path = REPO_ROOT / "data" / local_tsv
    if local_path.exists():
        print(f"[*] Loading cached transcriptomic cohort from {local_path.name}...")
        return pd.read_csv(local_path, sep="\t")
    
    print(f"[*] Fetching {study_id} from cBioPortal API...")
    url_clin = f"https://www.cbioportal.org/api/studies/{study_id}/clinical-data?clinicalDataType=PATIENT&pageSize=100000"
    req_p = urllib.request.Request(url_clin, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req_p, timeout=30) as r:
        p_raw = json.loads(r.read())
    
    p_dict = {}
    for entry in p_raw:
        pid = entry["patientId"]
        attr = entry["clinicalAttributeId"]
        val = entry["value"]
        if pid not in p_dict:
            p_dict[pid] = {"patientId": pid}
        p_dict[pid][attr] = val
    df_patient = pd.DataFrame(list(p_dict.values()))

    # Fetch expression
    gene_map = {"CDH1": 999, "TJP2": 7082, "EPCAM": 4072, "COL1A1": 1277, "ACTA2": 59, "IL1B": 3553}
    url_expr = f"https://www.cbioportal.org/api/molecular-profiles/{study_id}_rna_seq_v2_mrna_median_Zscores/molecular-data/fetch"
    payload = {"entrezGeneIds": list(gene_map.values()), "sampleListId": f"{study_id}_rna_seq_v2_mrna"}
    req_e = urllib.request.Request(url_expr, data=json.dumps(payload).encode("utf-8"),
                                   headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"}, method="POST")
    with urllib.request.urlopen(req_e, timeout=30) as r:
        e_raw = json.loads(r.read())
    
    df_e = pd.DataFrame(e_raw)
    id_to_sym = {v: k for k, v in gene_map.items()}
    df_e["symbol"] = df_e["entrezGeneId"].map(id_to_sym)
    pivot = df_e.pivot_table(index="sampleId", columns="symbol", values="value", aggfunc="mean")
    pivot["patientId"] = pivot.index.str.slice(0, 12)
    pivot = pivot.groupby("patientId").mean().reset_index()

    df = pd.merge(df_patient, pivot, on="patientId", how="inner")
    df.to_csv(local_path, sep="\t", index=False)
    return df


def run_pancancer_analysis():
    print("=" * 80)
    print(" PAN-CANCER CLINICAL COHORT VALIDATION PIPELINE")
    print("=" * 80)

    # ---------------------------------------------------------
    # 1. Cohort 1: Lung Adenocarcinoma (MSKCC 2020, N=604)
    # ---------------------------------------------------------
    print("\n--- [1/4] Processing Lung Adenocarcinoma (MSKCC 2020) ---")
    lung_path = REPO_ROOT / "data" / "mskcc_2020_real_cohort.csv"
    if not lung_path.exists():
        raise FileNotFoundError(f"Missing clinical cohort at: {lung_path}")
    df_lung = pd.read_csv(lung_path)

    def map_lung_histo(val):
        if val in ["Solid", "Micropapillary"]: return 2
        elif val in ["Acinar", "Papillary"]: return 1
        else: return 0

    df_lung["SEED"] = df_lung["PREDOMINANT_HISTOLOGIC_SUBTYPE"].apply(map_lung_histo) / 2.0
    df_lung["SOIL"] = df_lung["SMOKING_HISTORY"].apply(lambda x: 1.0 if "Ever" in str(x) else 0.2)
    # Biophysical decoupled HTBI: Seed polarity loss + Chronic airway inflammatory stress
    df_lung["HTBI"] = 0.65 * df_lung["SEED"] + 0.35 * df_lung["SOIL"]
    df_lung["TIME"] = pd.to_numeric(df_lung["RFS_MONTHS"], errors="coerce")
    df_lung["EVENT"] = pd.to_numeric(df_lung["RECURRENCE_EVENT"], errors="coerce")
    df_lung = df_lung.dropna(subset=["TIME", "EVENT"])
    df_lung = df_lung[df_lung["TIME"] > 0].copy()

    median_htbi_lung = df_lung["HTBI"].median()
    df_lung["RISK_GROUP"] = np.where(df_lung["HTBI"] >= median_htbi_lung, "High HTBI (Adverse Soil)", "Low HTBI (Favorable Soil)")

    # ---------------------------------------------------------
    # 2. Cohort 2: Colorectal Adenocarcinoma (TCGA COADREAD)
    # ---------------------------------------------------------
    print("\n--- [2/4] Processing Colorectal Cancer (TCGA COADREAD) ---")
    df_crc = fetch_or_load_tcga("coadread_tcga_pan_can_atlas_2018", "tcga_coadread_rnaseq.tsv")
    df_crc["TIME"] = pd.to_numeric(df_crc["OS_MONTHS"], errors="coerce")
    df_crc["EVENT"] = df_crc["OS_STATUS"].apply(lambda x: 1 if ("DECEASED" in str(x) or "1" in str(x)) else 0)
    df_crc = df_crc.dropna(subset=["TIME", "EVENT"])
    df_crc = df_crc[df_crc["TIME"] > 0].copy()

    # Biophysical mucosal barrier loss (seed loss) + stromal desmoplasia
    seed_loss = - (df_crc["CDH1"] + df_crc["TJP2"]) / 2.0
    stromal_mech = (df_crc["COL1A1"] + df_crc["ACTA2"]) / 2.0
    df_crc["HTBI"] = 0.50 * seed_loss + 0.50 * stromal_mech
    median_crc = df_crc["HTBI"].median()
    df_crc["RISK_GROUP"] = np.where(df_crc["HTBI"] >= median_crc, "High HTBI (Adverse Soil)", "Low HTBI (Favorable Soil)")

    # ---------------------------------------------------------
    # 3. Cohort 3: Pancreatic Adenocarcinoma (TCGA PAAD)
    # ---------------------------------------------------------
    print("\n--- [3/4] Processing Pancreatic Adenocarcinoma (TCGA PAAD) ---")
    df_paad = fetch_or_load_tcga("paad_tcga_pan_can_atlas_2018", "tcga_paad_rnaseq.tsv")
    df_paad["TIME"] = pd.to_numeric(df_paad["OS_MONTHS"], errors="coerce")
    df_paad["EVENT"] = df_paad["OS_STATUS"].apply(lambda x: 1 if ("DECEASED" in str(x) or "1" in str(x)) else 0)
    df_paad = df_paad.dropna(subset=["TIME", "EVENT"])
    df_paad = df_paad[df_paad["TIME"] > 0].copy()

    # Pancreatic biophysical barrier: Hypoxia transport resistance + tight junction disruption
    hyp = pd.to_numeric(df_paad.get("BUFFA_HYPOXIA_SCORE", 0), errors="coerce").fillna(0)
    hyp_norm = (hyp - hyp.mean()) / (hyp.std() + 1e-6)
    tjp2_norm = df_paad["TJP2"].fillna(0)
    df_paad["HTBI"] = 0.60 * hyp_norm + 0.40 * tjp2_norm
    median_paad = df_paad["HTBI"].median()
    df_paad["RISK_GROUP"] = np.where(df_paad["HTBI"] >= median_paad, "High HTBI (Adverse Soil)", "Low HTBI (Favorable Soil)")

    # ---------------------------------------------------------
    # 4. Cohort 4: Cholangiocarcinoma (TCGA CHOL)
    # ---------------------------------------------------------
    print("\n--- [4/4] Processing Cholangiocarcinoma (TCGA CHOL) ---")
    df_chol = fetch_or_load_tcga("chol_tcga_pan_can_atlas_2018", "tcga_chol_rnaseq.tsv")
    df_chol["TIME"] = pd.to_numeric(df_chol["OS_MONTHS"], errors="coerce")
    df_chol["EVENT"] = df_chol["OS_STATUS"].apply(lambda x: 1 if ("DECEASED" in str(x) or "1" in str(x)) else 0)
    df_chol = df_chol.dropna(subset=["TIME", "EVENT"])
    df_chol = df_chol[df_chol["TIME"] > 0].copy()

    # Biliary fibrogenesis: Desmoplasia (COL1A1/ACTA2) + loss of junctional polarity (-CDH1)
    desm = (df_chol["COL1A1"] + df_chol["ACTA2"]) / 2.0
    pol_loss = - df_chol["CDH1"]
    df_chol["HTBI"] = 0.60 * desm + 0.40 * pol_loss
    median_chol = df_chol["HTBI"].median()
    df_chol["RISK_GROUP"] = np.where(df_chol["HTBI"] >= median_chol, "High HTBI (Adverse Soil)", "Low HTBI (Favorable Soil)")

    # ---------------------------------------------------------
    # Plotting Pan-Cancer 4-Panel Kaplan-Meier Validation
    # ---------------------------------------------------------
    print("\n[*] Plotting Pan-Cancer Validation Grid...")
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=300)

    cohorts_data = [
        ("A. Lung Adenocarcinoma (MSKCC Discovery, N=604)", df_lung, axes[0, 0], "Relapse-Free Survival (Months)", 60),
        (f"B. Colorectal Adenocarcinoma (TCGA COADREAD, N={len(df_crc)})", df_crc, axes[0, 1], "Overall Survival (Months)", 60),
        (f"C. Pancreatic Adenocarcinoma (TCGA PAAD, N={len(df_paad)})", df_paad, axes[1, 0], "Overall Survival (Months)", 48),
        (f"D. Cholangiocarcinoma (TCGA CHOL, N={len(df_chol)})", df_chol, axes[1, 1], "Overall Survival (Months)", 48)
    ]

    for title, df_c, ax, xlabel, max_t in cohorts_data:
        df_low = df_c[df_c["RISK_GROUP"].str.contains("Low")].copy()
        df_high = df_c[df_c["RISK_GROUP"].str.contains("High")].copy()

        # Compute log-rank test
        lr_res = logrank_test(
            df_low["TIME"], df_high["TIME"],
            event_observed_A=df_low["EVENT"], event_observed_B=df_high["EVENT"]
        )
        p_val = lr_res.p_value

        km_low = KaplanMeierFitter()
        km_high = KaplanMeierFitter()

        km_low.fit(df_low["TIME"], df_low["EVENT"], label=f"Low HTBI / Intact Polarity (n={len(df_low)})")
        km_high.fit(df_high["TIME"], df_high["EVENT"], label=f"High HTBI / Adverse Microenvironment (n={len(df_high)})")

        km_low.plot_survival_function(ax=ax, color="#0284c7", linewidth=2.5)
        km_high.plot_survival_function(ax=ax, color="#e11d48", linewidth=2.5)

        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
        ax.set_ylabel("Survival Probability", fontsize=10, fontweight='bold')
        ax.set_xlim(0, max_t)
        ax.set_ylim(-0.02, 1.02)
        ax.legend(loc="lower left", fontsize=8.5, frameon=True)

        # Survival separation annotation
        p_str = f"p = {p_val:.4f}" if p_val >= 0.0001 else f"p = {p_val:.2e}"
        ax.text(0.95, 0.92, f"Log-Rank: {p_str}\nBiophysical Decoupled Model", transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8fafc', edgecolor='#cbd5e1', lw=1.2),
                fontsize=8.5, fontweight='bold', color='#1e293b')

    plt.tight_layout()
    output_png = REPO_ROOT / "figure_3_pancancer_validation.png"
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()
    print(f"\n[✔] Pan-Cancer validation complete! Figure saved at: {output_png}")



if __name__ == "__main__":
    run_pancancer_analysis()
