"""
Comparative Simulation: Low vs High Chronic Inflammatory Stress in Target Therapy Recurrence
Authors: Luca Celebrano, Antigravity AI
Repository: /home/luca/inflammatory-onco-physics
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from model import ModelParameters, InflammatoryOncoContinuum

def run_comparative_study():
    print("[*] Initializing Biophysical Onco-Continuum Simulation...")
    
    # Cohort A: Low Inflammatory Stress (Clean Homeostasis, NLR < 2, Low CRP)
    params_low = ModelParameters(
        sigma_base=0.04,   # Low baseline systemic stress
        alpha=0.03,        # Low persister transition
        gamma=0.008        # Slow immune exhaustion
    )
    sim_low = InflammatoryOncoContinuum(params_low)
    print(" -> Simulating Cohort A (Low Systemic Inflammation)...")
    res_low = sim_low.simulate(t_span_months=60.0, t_therapy_stop=36.0)

    # Cohort B: High Inflammatory Stress (Chronic Inflammation, NLR > 4, Elevated CRP)
    params_high = ModelParameters(
        sigma_base=0.22,   # Elevated chronic systemic stress
        alpha=0.06,        # High persister induction
        gamma=0.016        # Accelerated immune exhaustion
    )
    sim_high = InflammatoryOncoContinuum(params_high)
    print(" -> Simulating Cohort B (High Chronic Inflammation)...")
    res_high = sim_high.simulate(t_span_months=60.0, t_therapy_stop=36.0)

    print("[*] Generating Publication-Grade Figures...")
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)

    time = res_low["time"]

    # -------------------------------------------------------------
    # Panel 1: Total Tumor Burden (Clinically detectable recurrence)
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1.axvspan(0, 36, color='#e0f2fe', alpha=0.5, label='Active Adjuvant Therapy Window (36m)')
    ax1.axvline(36, color='#dc2626', linestyle='--', linewidth=1.5, label='Therapy Cessation (Month 36)')
    
    ax1.plot(time, res_low["N_tot"], color='#0284c7', linewidth=2.5, label='Low Inflammation (Cure / Dormant Stability)')
    ax1.plot(time, res_high["N_tot"], color='#e11d48', linewidth=2.5, label='High Chronic Inflammation (Catastrophic Relapse)')
    
    ax1.set_title('A. Total Tumor Burden & Post-Cessation Dynamics', fontsize=12, fontweight='bold', pad=10)
    ax1.set_xlabel('Time (Months)', fontsize=10)
    ax1.set_ylabel('Mean Tumor Density (arb. units)', fontsize=10)
    ax1.set_ylim(-0.02, 0.45)
    ax1.legend(loc='upper right', frameon=True, fontsize=9)

    # -------------------------------------------------------------
    # Panel 2: Persister Reservoir (DTPs) Dynamics
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.axvspan(0, 36, color='#e0f2fe', alpha=0.5)
    ax2.axvline(36, color='#dc2626', linestyle='--', linewidth=1.5)

    ax2.plot(time, res_low["N_R"], color='#0284c7', linestyle='-.', linewidth=2.2, label='DTP Reservoir (Low Inflammation)')
    ax2.plot(time, res_high["N_R"], color='#e11d48', linestyle='-.', linewidth=2.2, label='DTP Reservoir (High Inflammation)')

    ax2.set_title('B. Plastic Drug-Tolerant Persister (DTP) Subpopulation', fontsize=12, fontweight='bold', pad=10)
    ax2.set_xlabel('Time (Months)', fontsize=10)
    ax2.set_ylabel('DTP Density (arb. units)', fontsize=10)
    ax2.legend(loc='upper left', frameon=True, fontsize=9)

    # -------------------------------------------------------------
    # Panel 3: History-Dependent Immune Capacity Exhaustion
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3.axvspan(0, 36, color='#e0f2fe', alpha=0.5)
    ax3.axvline(36, color='#dc2626', linestyle='--', linewidth=1.5)

    # Theoretical minimal clearance threshold (r_R)
    ax3.axhline(params_low.r_R, color='#64748b', linestyle=':', linewidth=1.5, label=r'Critical Recurrence Threshold ($r_R$)')

    ax3.plot(time, res_low["E_imm"], color='#0284c7', linewidth=2.2, label=r'Immune Competence $\mathcal{E}_{\mathrm{imm}}$ (Low Stress)')
    ax3.plot(time, res_high["E_imm"], color='#e11d48', linewidth=2.2, label=r'Immune Exhaustion $\mathcal{E}_{\mathrm{imm}} \to 0$ (High Stress)')

    ax3.set_title(r'C. Cumulative Epigenetic Immune Exhaustion ($\mathcal{E}_{\mathrm{imm}}$)', fontsize=12, fontweight='bold', pad=10)
    ax3.set_xlabel('Time (Months)', fontsize=10)
    ax3.set_ylabel('Cytotoxic Clearance Rate (1/month)', fontsize=10)
    ax3.legend(loc='lower left', frameon=True, fontsize=9)

    # -------------------------------------------------------------
    # Panel 4: Spatial Microenvironmental Profile at Cessation (t=36m)
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    snap_low = res_low["snapshots"].get("month_36")
    snap_high = res_high["snapshots"].get("month_36")

    if snap_low and snap_high:
        x = snap_low["x"]
        ax4.plot(x, snap_low["N_R"], color='#0284c7', linewidth=2.2, label='Persisters (Low Inf)')
        ax4.plot(x, snap_high["N_R"], color='#e11d48', linewidth=2.2, label='Persisters (High Inf)')
        
        # Twin axis for drug concentration
        ax4_twin = ax4.twinx()
        ax4_twin.plot(x, snap_low["C"], color='#6b7280', linestyle='--', linewidth=1.8, label='Drug Gradient C(x)')
        ax4_twin.set_ylabel('Local Drug Penetration C(x)', color='#6b7280', fontsize=10)
        ax4_twin.tick_params(axis='y', labelcolor='#6b7280')
        ax4_twin.grid(False)

    ax4.set_title('D. Spatial Sanctuary Profile at Month 36 (Cessation)', fontsize=12, fontweight='bold', pad=10)
    ax4.set_xlabel(r'Distance from Capillary $x$ ($\mu$m)', fontsize=10)
    ax4.set_ylabel('DTP Density at Month 36', fontsize=10)
    ax4.legend(loc='upper left', frameon=True, fontsize=9)

    plt.tight_layout()
    output_png = '/home/luca/inflammatory-onco-physics/figure_1_bifurcation.png'
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()
    print(f"[✔] Simulation complete. Figure saved at: {output_png}")

if __name__ == '__main__':
    run_comparative_study()
