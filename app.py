"""
Inflammatory Onco-Physics: Interactive Streamlit Web Application & Digital Twin
Author: Luca Celebrano, MD Candidate (Università Federico II / TIGEM)
Repository: https://github.com/Lukecele/inflammatory-onco-physics
DOI: https://doi.org/10.5281/zenodo.24118328
"""

import streamlit as st
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="Inflammatory Onco-Physics",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling
st.markdown("""
<style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 800;
        color: #0f172a;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.1rem;
        color: #475569;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 16px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .badge-bar {
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
        margin-bottom: 1.5rem;
    }
</style>
""", unsafe_allow_html=True)

# Top Bar
st.markdown('<div class="main-title">🔬 Reaction-Diffusion Theory of Microvascular Pharmacological Sanctuaries</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Interactive Biophysical Simulation & Bedside HTBI Risk Engine • Luca Celebrano (N = 3,461)</div>', unsafe_allow_html=True)

st.markdown("""
<div class="badge-bar">
    <a href="https://doi.org/10.5281/zenodo.24118328" target="_blank"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.24118328.svg" alt="DOI"></a>
    <a href="https://github.com/Lukecele/inflammatory-onco-physics" target="_blank"><img src="https://img.shields.io/badge/GitHub-Lukecele%2Finflammatory--onco--physics-24292e.svg?logo=github" alt="GitHub"></a>
    <a href="https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"></a>
    <a href="https://lucacelebrano.me/" target="_blank"><img src="https://img.shields.io/badge/%F0%9F%8F%B0%20Luca's%20Cave-Research%20Hub-60d9c5.svg" alt="Luca's Cave"></a>
    <img src="https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg" alt="CC BY 4.0">
    <img src="https://img.shields.io/badge/REMARK-Compliant-success.svg" alt="REMARK">
</div>
""", unsafe_allow_html=True)

# Sidebar Parameters
st.sidebar.header("⚙️ Biophysical Parameters")

r_v = st.sidebar.slider("Capillary Vessel Radius (r_v, μm)", min_value=5.0, max_value=25.0, value=10.0, step=1.0)
r_ox = st.sidebar.slider("Outer Normoxic Boundary (r_ox, μm)", min_value=120.0, max_value=300.0, value=200.0, step=10.0)
C_0 = st.sidebar.slider("Luminal Free Drug (C_0, nM)", min_value=20.0, max_value=300.0, value=100.0, step=10.0)
IC_50 = st.sidebar.slider("Cytocidal Threshold (IC_50, nM)", min_value=2.0, max_value=30.0, value=10.0, step=1.0)
D_eff = st.sidebar.slider("Effective Diffusion (D_eff, μm²/s)", min_value=200.0, max_value=1200.0, value=650.0, step=50.0)
k_deg = st.sidebar.slider("Metabolic Clearance (k × 10⁻⁴ s⁻¹)", min_value=0.5, max_value=3.0, value=1.51, step=0.1) * 1e-4

st.sidebar.markdown("---")
st.sidebar.header("🩸 Clinical Biomarkers (Bedside)")
stage = st.sidebar.selectbox("Pathological TNM Stage", ["Stage IB", "Stage IIA/B", "Stage IIIA"])
histology_grade = st.sidebar.selectbox("Histologic Differentiation", ["G1 (Well-differentiated / Lepidic)", "G2 (Acinar / Papillary)", "G3 (Solid / Micropapillary)"])
nlr = st.sidebar.number_input("Pre-Op Neutrophil-to-Lymphocyte Ratio (NLR)", min_value=0.5, max_value=25.0, value=2.5, step=0.1)
crp = st.sidebar.number_input("Serum hs-CRP (mg/L)", min_value=0.1, max_value=150.0, value=6.5, step=0.5)
albumin = st.sidebar.number_input("Serum Albumin (g/L)", min_value=15.0, max_value=55.0, value=40.0, step=1.0)

# Solver Computation
lam_C = np.sqrt(D_eff / k_deg) # um
r_vals = np.linspace(r_v, r_ox, 1000)
C_r = C_0 * sp.kn(0, r_vals / lam_C) / sp.kn(0, r_v / lam_C)

# Critical boundary
crit_mask = C_r <= IC_50
if np.any(crit_mask):
    idx_crit = np.where(crit_mask)[0][0]
    r_crit = r_vals[idx_crit]
else:
    r_crit = r_ox

V_total = np.pi * (r_ox**2 - r_v**2)
V_sanctuary = np.pi * (r_ox**2 - r_crit**2)
frac_sanctuary = max(0.0, min(100.0, (V_sanctuary / V_total) * 100))

# Tabs
tab1, tab2, tab3 = st.tabs([
    "🔬 1. Cylindrical Bessel K₀ Solver & Spatial Sanctuary",
    "⏱️ 2. Double-Barrier ctDNA Shedding & Lead Time",
    "🏥 3. Bedside HTBI Risk Engine & Survival"
])

with tab1:
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Diffusion Length (λ)", f"{lam_C:.1f} μm")
    col2.metric("Critical Boundary (r_crit)", f"{r_crit:.1f} μm")
    col3.metric("Sanctuary Sleeve Width", f"{r_ox - r_crit:.1f} μm")
    col4.metric("DTP Sanctuary Volume", f"{frac_sanctuary:.1f}%")

    col_plot, col_2d = st.columns([1.2, 0.8])
    with col_plot:
        fig1, ax1 = plt.subplots(figsize=(7, 4), dpi=130)
        ax1.plot(r_vals, C_r, color="#0284c7", lw=2.5, label="Osimertinib $C(r)$")
        ax1.axhline(IC_50, color="#dc2626", ls="--", lw=1.8, label=f"Cytocidal $\\mathrm{{IC}}_{{50}} = {IC_50}$ nM")
        ax1.axvline(r_crit, color="#d97706", ls=":", lw=1.8, label=f"$r_{{\\mathrm{{crit}}}} = {r_crit:.1f}$ μm")
        ax1.fill_between(r_vals[r_vals >= r_crit], 0, C_r[r_vals >= r_crit], color="#fef3c7", alpha=0.85, label=f"DTP Sanctuary ({frac_sanctuary:.1f}% Vol)")
        ax1.set_xlabel("Radial Distance from Capillary (μm)", fontweight="bold")
        ax1.set_ylabel("Free Drug Concentration (nM)", fontweight="bold")
        ax1.set_title("Cylindrical Radial Reaction-Diffusion Profile", fontweight="bold", fontsize=11)
        ax1.set_xlim(r_v, r_ox)
        ax1.set_ylim(0, max(C_0 * 1.05, 110))
        ax1.legend(loc="upper right", fontsize=8, framealpha=0.95)
        ax1.grid(True, linestyle="--", alpha=0.5)
        st.pyplot(fig1)

    with col_2d:
        # 2D cross-section contour map
        grid_dim = 250
        x = np.linspace(-r_ox, r_ox, grid_dim)
        y = np.linspace(-r_ox, r_ox, grid_dim)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        Z = np.zeros_like(R)
        valid = (R >= r_v) & (R <= r_ox)
        Z[valid] = C_0 * sp.kn(0, R[valid] / lam_C) / sp.kn(0, r_v / lam_C)
        Z[R < r_v] = C_0

        fig2, ax2 = plt.subplots(figsize=(5, 4), dpi=130)
        im = ax2.contourf(X, Y, Z, levels=25, cmap="magma")
        # Circle boundaries
        vessel_circ = plt.Circle((0, 0), r_v, color="#38bdf8", fill=True, label="Capillary Vessel")
        crit_circ = plt.Circle((0, 0), r_crit, color="#fbbf24", fill=False, ls="--", lw=1.5, label="r_crit Threshold")
        outer_circ = plt.Circle((0, 0), r_ox, color="#e2e8f0", fill=False, ls="-", lw=1.2)
        ax2.add_patch(vessel_circ)
        ax2.add_patch(crit_circ)
        ax2.add_patch(outer_circ)
        ax2.set_xlim(-r_ox, r_ox)
        ax2.set_ylim(-r_ox, r_ox)
        ax2.set_aspect("equal")
        ax2.set_title("2D Tumor Cord Sanctuary Cross-Section", fontweight="bold", fontsize=11)
        ax2.set_xlabel("x (μm)")
        ax2.set_ylabel("y (μm)")
        plt.colorbar(im, ax=ax2, label="Drug (nM)", fraction=0.046, pad=0.04)
        st.pyplot(fig2)

with tab2:
    st.subheader("Steric Entrapment & Pre-Radiological ctDNA Lead Time")
    st.markdown("Because early persister cells are arrested in $G_0/G_1$ without significant apoptosis, and dense collagen fibers ($\le 11$ nm pore size) trap nucleosomes, plasma ctDNA remains false-negative during initial dormancy.")
    
    t_days = np.linspace(0, 500, 500)
    T_d = 24.6
    N_t = 1e4 * (2.0 ** (t_days / T_d))
    h_shed = 1e-4 * (1.0 + np.tanh((np.log10(N_t) - 6.0) * 1.5))
    ctDNA_vaf = h_shed * (N_t / 1e9) * 100 # %

    idx_mrd = np.where(ctDNA_vaf >= 0.02)[0][0]
    idx_rad = np.where(N_t >= 1e8)[0][0]
    t_mrd = t_days[idx_mrd]
    t_rad = t_days[idx_rad]
    lead_time = (t_rad - t_mrd) / 30.4375

    col_lt1, col_lt2, col_lt3 = st.columns(3)
    col_lt1.metric("ctDNA MRD+ Detection", f"Day {t_mrd:.1f} ({t_mrd/30.4375:.1f} mo)")
    col_lt2.metric("CT Imaging Recurrence", f"Day {t_rad:.1f} ({t_rad/30.4375:.1f} mo)")
    col_lt3.metric("Pre-Radiological Lead Time", f"{lead_time:.2f} Months", delta="Intervention Window")

    fig3, ax3 = plt.subplots(figsize=(8, 3.8), dpi=130)
    ax3.plot(t_days, ctDNA_vaf, color="#0284c7", lw=2.5, label="Plasma ctDNA VAF (%)")
    ax3.axhline(0.02, color="#059669", ls="--", label="Digital NGS MRD Limit (0.02%)")
    ax3.axvline(t_mrd, color="#059669", ls=":")
    ax3.axvline(t_rad, color="#dc2626", ls=":", label="CT Macro-Recurrence")
    ax3.axvspan(t_mrd, t_rad, color="#10b981", alpha=0.18, label=f"Window of Clinical Opportunity ({lead_time:.1f} Mo)")
    ax3.set_xlabel("Days Following Adjuvant Cessation", fontweight="bold")
    ax3.set_ylabel("Plasma ctDNA VAF (%)", fontweight="bold")
    ax3.set_yscale("log")
    ax3.set_ylim(0.001, 10.0)
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(loc="upper left", fontsize=9)
    st.pyplot(fig3)

with tab3:
    st.subheader("Bedside Host-Tumor Biophysical Index (HTBI)")
    
    # Calculate GPS
    gps = 0
    if crp > 10.0 and albumin < 35.0:
        gps = 2
    elif crp > 10.0 or albumin < 35.0:
        gps = 1

    # Grade
    grade_val = 0.50
    if "G1" in histology_grade:
        grade_val = 0.10
    elif "G3" in histology_grade:
        grade_val = 0.90

    # HTBI formula
    htbi_score = 0.40 * grade_val + 0.35 * min(1.0, nlr / 5.0) + 0.25 * (gps / 2.0)
    
    # Risk category
    if htbi_score < 0.35:
        risk_cat = "Low Host-Tissue Risk"
        risk_color = "#10b981"
        rec_3yr = "12 - 18%"
        monitoring = "Standard Follow-up (ctDNA q3m, CT q6m)"
    elif htbi_score < 0.65:
        risk_cat = "Intermediate Host-Tissue Risk"
        risk_color = "#f59e0b"
        rec_3yr = "28 - 38%"
        monitoring = "Intensive Surveillance (ctDNA q8w, anti-inflammatory pre-conditioning)"
    else:
        risk_cat = "High Host-Tissue Risk (Sanctuary Escape)"
        risk_color = "#ef4444"
        rec_3yr = "52 - 68%"
        monitoring = "High-Risk Protocol (ctDNA q6w, consider prolonged adjuvant / trial enrollment)"

    st.markdown(f"""
    <div style="background-color: #f1f5f9; padding: 20px; border-radius: 10px; border-left: 6px solid {risk_color}; margin-bottom: 20px;">
        <h3 style="margin-top:0; color: #0f172a;">Composite HTBI Score: <strong>{htbi_score:.3f}</strong> ({risk_cat})</h3>
        <p><strong>3-Year Post-Adjuvant Recurrence Risk:</strong> {rec_3yr} | <strong>Glasgow Prognostic Score:</strong> GPS {gps}</p>
        <p style="margin-bottom:0;"><strong>Recommended Strategy:</strong> {monitoring}</p>
    </div>
    """, unsafe_allow_html=True)
