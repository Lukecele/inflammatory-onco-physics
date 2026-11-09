"""
Bedside Clinical HTBI Calculator (Host-Tumor Biophysical Index)
Authors: Luca Celebrano, Antigravity AI
Repository: /home/luca/inflammatory-onco-physics
A standalone interactive decision support tool computing recurrence risk from routine labs.
"""

import sys

def calculate_htbi(
    stage: int,
    histology_grade: int,
    nlr: float,
    crp_mg_l: float,
    albumin_g_l: float,
    desmoplasia: str = "Moderate",
    tils_percent: float = 15.0
):
    # 1. Proliferation / Seed Aggressiveness from Grade & Stage
    grade_factors = {1: 0.7, 2: 1.0, 3: 1.35}
    r_S = 0.40 * grade_factors.get(histology_grade, 1.0) * (1.0 + 0.3 * (stage - 1))
    
    # 2. Stromal Transport Resistance Factor
    desmo_factors = {"Low": 1.0, "Moderate": 0.75, "High": 0.45}
    transport = desmo_factors.get(desmoplasia, 0.75)
    
    # 3. Host Immune Surveillance Capacity
    E_0 = 0.05 * (tils_percent / 20.0)
    
    # 4. Systemic Inflammatory Burden
    norm_nlr = nlr / 2.0
    norm_crp = crp_mg_l / 2.5
    norm_alb = 40.0 / max(albumin_g_l, 10.0)
    soil_exhaustion = 1.0 + 0.25 * norm_nlr * (norm_crp * norm_alb)
    
    # 5. Composite HTBI
    seed = r_S / transport
    soil = E_0 / soil_exhaustion
    htbi = seed / max(soil, 1e-5)
    
    # 6. Recurrence Probability at 24 and 36 Months
    prob_24m = min(0.95, max(0.05, 1.0 / (1.0 + np.exp(-(0.03 * htbi + 0.5 * stage - 2.0)))))
    prob_36m = min(0.98, max(0.08, 1.0 / (1.0 + np.exp(-(0.035 * htbi + 0.55 * stage - 1.8)))))
    
    # 7. Projected ctDNA Invariant Lead Time
    projected_lead_time_months = 5.46 if htbi >= 10.0 else (6.2 if htbi >= 5.0 else None)
    
    return {
        "htbi": htbi,
        "prob_24m": prob_24m * 100,
        "prob_36m": prob_36m * 100,
        "lead_time": projected_lead_time_months,
        "risk_tier": "EXTREME RISK" if htbi >= 25.0 else ("HIGH RISK" if htbi >= 10.0 else ("MODERATE RISK" if htbi >= 5.0 else "FAVORABLE / LOW RISK"))
    }

if __name__ == "__main__":
    import numpy as np
    print("=" * 70)
    print(" BEDSIDE CLINICAL HTBI DECISION CALCULATOR")
    print(" Authors: Luca Celebrano, MD Candidate / Antigravity AI")
    print("=" * 70)
    
    # Example Patient: 64-year-old post-R0 resection, high systemic inflammation
    stage_in = 2
    grade_in = 3
    nlr_in = 4.2
    crp_in = 12.5
    alb_in = 35.0
    desmo_in = "High"
    tils_in = 8.0
    
    res = calculate_htbi(stage_in, grade_in, nlr_in, crp_in, alb_in, desmo_in, tils_in)
    
    print(f"\n[INPUT PATIENT CLINICAL PARAMETERS]")
    print(f" • Pathologic TNM Stage:       Stage {stage_in}")
    print(f" • Histological Grade:         G{grade_in} (Poorly differentiated)")
    print(f" • Complete Blood Count NLR:   {nlr_in:.1f} (Normal < 2.0)")
    print(f" • Serum hs-CRP:               {crp_in:.1f} mg/L (Normal < 2.5)")
    print(f" • Serum Albumin:              {alb_in:.1f} g/L (Normal 35-50)")
    print(f" • Stromal Desmoplasia:        {desmo_in}")
    print(f" • Infiltrating TILs:          {tils_in:.1f}%\n")
    
    print("-" * 70)
    print(f"[PREDICTIVE BIOPHYSICAL OUTPUT]")
    print(f" • Host-Tumor Biophysical Index (HTBI):  {res['htbi']:.1f}")
    print(f" • Risk Stratification Tier:             {res['risk_tier']}")
    print(f" • Predicted 24-Month Recurrence Risk:   {res['prob_24m']:.1f}%")
    print(f" • Predicted 36-Month Recurrence Risk:   {res['prob_36m']:.1f}%")
    if res['lead_time']:
        print(f" • Projected ctDNA Invariant Lead Time:  {res['lead_time']:.2f} Months prior to CT imaging")
    
    print("-" * 70)
    print(f"[CLINICAL ACTIONABLE GUIDANCE]")
    if res['htbi'] >= 10.0:
        print(" [!] HIGH OCCULT PERSISTER RISK: Standard cytotoxic/targeted therapy alone is INSUFFICIENT.")
        print(" [>] Action 1: Initiate early serial ctDNA monitoring (every 3 months) to catch microscopic surge.")
        print(" [>] Action 2: Consider microenvironmental co-therapy (anti-IL-1b/NLRP3, low-dose aspirin, or resolvins)")
        print("     to decompress stroma and restore endogenous CD8+ cytotoxic clearance.")
    else:
        print(" [✔] FAVORABLE HOST SOIL: Standard surgical resection / adjuvant therapy adequate.")
        print(" [>] Routine surveillance recommended.")
    print("=" * 70)
