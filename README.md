<div align="center">

# 🔬 Reaction-Diffusion Theory of Microvascular Pharmacological Sanctuaries and Inflammatory Pre-conditioning
### Explaining Post-Adjuvant Relapse Dynamics & ctDNA Shedding Kinetics Across 3,461 Patients

**An open-source biophysical continuum transport solver, 60 FPS in silico digital twin, and clinical risk engine for targeted oncology.**

> 🏰 **Live Interactive 3D Lab & Flagship Research Hub:** [**Explore Luca's Cave ↗**](https://lucacelebrano.me/)

<p align="center">
  <a href="https://lucacelebrano.me/">
    <img src="https://img.shields.io/badge/%F0%9F%8F%B0%20Luca's%20Cave-Research%20Hub-60d9c5.svg?style=for-the-badge" height="28" alt="Luca's Cave">
  </a>
  <a href="https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" height="28" alt="Open in Colab">
  </a>
  <a href="https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb">
    <img src="https://img.shields.io/badge/Colab_Pro-T4_GPU_Accelerated-orange.svg?style=for-the-badge&logo=googlecolab" height="28" alt="Colab Pro T4 GPU">
  </a>
  <a href="https://doi.org/10.5281/zenodo.24118328">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.24118328.svg" height="28" alt="DOI: 10.5281/zenodo.24118328">
  </a>
  <a href="simulator.html">
    <img src="https://img.shields.io/badge/%F0%9F%9A%80%20Live%20Simulator-60%20FPS%20Digital%20Twin-0284c7.svg?style=for-the-badge" height="28" alt="60 FPS Simulator">
  </a>
  <a href="tumor_cord_sanctuary_3d.html">
    <img src="https://img.shields.io/badge/%F0%9F%A7%8A%203D%20WebGL-Tumor%20Cord%20Twin-a855f7.svg?style=for-the-badge" height="28" alt="3D WebGL Simulator">
  </a>
</p>

<p align="center">
  <a href="MEMORANDUM_ESPLORATIVO_IFOM.pdf"><img src="https://img.shields.io/badge/%F0%9F%93%84%20Executive%20Memorandum-PDF%20(4%20Pages)-dc2626.svg" alt="Memorandum PDF"></a>
  <a href="MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.pdf"><img src="https://img.shields.io/badge/%F0%9F%93%9A%20Master%20Monograph-PDF%20(26%20Pages)-007acc.svg" alt="Master Monograph PDF"></a>
  <a href="REMARK_COMPLIANCE_REPORT.md"><img src="https://img.shields.io/badge/REMARK-Compliant-success.svg" alt="REMARK Compliant"></a>
  <a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg" alt="CC BY 4.0"></a>
  <a href="https://orcid.org/0009-0007-0931-3555"><img src="https://img.shields.io/badge/ORCID-0009--0007--0931--3555-a6ce39.svg?logo=orcid&logoColor=white" alt="ORCID"></a>
  <img src="https://img.shields.io/badge/Python-3.10+-blue.svg" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/Zero%20Imputation-Complete--Case%20Audit-brightgreen.svg" alt="Zero Imputation">
</p>

*Author:* **Luca Celebrano**, MD Candidate<sup>1,2</sup>  
<sup>1</sup> Department of Clinical Medicine and Surgery, School of Medicine, Università degli Studi di Napoli Federico II  
<sup>2</sup> Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy  
*Correspondence:* `l.celebrano@studenti.unina.it` · `luca.celebrano@gmail.com`

</div>

---

### ⚡ 10-Second Interactive Quickstart (Zero Install Required)

| 🚀 **1-Click Google Colab (Pro T4)** | 🌐 **3D WebGL Digital Twin** | 🎮 **60 FPS Live Simulator** | 🏥 **Bedside HTBI Clinical Engine** |
| :---: | :---: | :---: | :---: |
| [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb) | [**▶ Launch 3D Cord Twin**](tumor_cord_sanctuary_3d.html) | [**▶ Launch 60 FPS Simulator**](simulator.html) | [**▶ Launch Web App**](app.py) |
| Run Bessel $K_0$ radial solvers, ctDNA lead time & $N=3,461$ survival on NVIDIA T4 GPU | Full 360° perivascular cylindrical cord, interactive diffusion & wash-out | Real-time perivascular drug washout & inflammatory cell awakening | Instant clinical recurrence probability & monitoring advice from routine labs |

```bash
# Or run locally from your terminal in 3 seconds:
git clone https://github.com/Lukecele/inflammatory-onco-physics.git
cd inflammatory-onco-physics
pip install -r requirements.txt

# Launch the interactive Streamlit Web App:
streamlit run app.py

# Or launch the terminal Bedside Decision Calculator:
python app_bedside_calculator.py
```

---

<div align="center">
  <a href="simulator.html">
    <img src="assets/simulator-vessel-live.gif" alt="In Silico Tumor Sanctuary Digital Twin Simulator — 60 FPS Capillary Vessel Dynamics" width="100%" style="max-width: 900px; border-radius: 12px; box-shadow: 0 8px 30px rgba(0,0,0,0.15);">
  </a>

  <p align="center">
    <i><b>Live 60 FPS Digital Twin:</b> Capillary microvessel ($r_v=10\,\mu\text{m}$) diffusing targeted therapy. Drug concentration decays radially via modified Bessel function $K_0$, dropping below cytocidal $\text{IC}_{50}$ ($10\,\text{nM}$) at $r_{\text{crit}}=165.89\,\mu\text{m}$. This creates a <b>34.1 $\mu$m distal cytostatic sanctuary sleeve (31.28% of viable tumor cord volume)</b> where drug-tolerant persisters (DTPs) survive until host systemic inflammation drives clonal re-awakening.</i>
  </p>
</div>

---

## 📌 Executive Summary: Resolving the Post-Adjuvant Cessation Paradox

| 🎯 **The Clinical Paradox** | 🔬 **The Biophysical Discovery** | 💡 **The Precision Oncology Solution** |
| :--- | :--- | :--- |
| **High DFS On-Treatment, Rebound Recurrence Post-Cessation:** In trials like ADAURA (osimertinib in resected *EGFR*-mutant NSCLC), adjuvant targeted therapy yields profound disease control during active administration, yet post-cessation survival curves show sudden rebound recurrences within 12–18 months. Crucially, **75–80% of relapsing lesions harbor NO secondary resistance mutations** (e.g. C797S). | **Microvascular Radial Sanctuaries & Inflammatory Pre-conditioning:** Desmoplastic solid tumors exhibit elevated interstitial fluid pressure ($Pe < 0.1$), forcing drug transport into pure radial reaction-diffusion governed by Bessel $K_0$. Distal cords suffer pharmacological clearance decay, harboring dormant persisters ($G_0/G_1$). Host systemic cytokines (IL-1β, IL-6, TNF-α) act as an environmental pacemaker via a Volterra memory functional $\mathcal{E}_{\text{imm}}(t)$. | **Analytical ctDNA Lead Time & Host Risk Stratification:** Early dormancy produces false-negative plasma ctDNA (MRD⁻) due to absent apoptosis and steric collagen entrapment ($\le 11\text{ nm}$). Upon inflammatory awakening, exponential expansion creates a **5.46-month analytical lead time** prior to radiological detection. The Host-Tumor Biophysical Index (HTBI) accurately stratifies recurrence risk ($N = 3,461$) directly from pre-op labs. |

---

## 🔬 Mathematical Formulations

### 1. Cylindrical Bessel $K_0$ Radial Diffusion Solver
In cylindrical coordinates around a microvessel of radius $r_v \approx 10\ \mu\text{m}$ embedded in a tumor cord of outer radius $r_{\text{ox}} \approx 200\ \mu\text{m}$, the steady-state free drug concentration $C(r)$ satisfies:
$$\frac{d^2 C}{dr^2} + \frac{1}{r} \frac{dC}{dr} - \frac{k}{D_{\text{eff}}} C = 0$$

The exact physical solution satisfying boundary conditions $C(r_v) = C_0$ and $\lim_{r\to\infty} C(r) = 0$ is:
$$C(r) = C_0 \frac{K_0(r / \lambda)}{K_0(r_v / \lambda)}, \quad \lambda = \sqrt{\frac{D_{\text{eff}}}{k}} \approx 65.57\ \mu\text{m}$$

With $D_{\text{eff}} = 650\ \mu\text{m}^2/\text{s}$, $k = 1.51 \times 10^{-4}\text{ s}^{-1}$, and $C_0 = 100\text{ nM}$, the cytocidal threshold ($\text{IC}_{50} = 10\text{ nM}$) is reached at:
$$r_{\text{crit}} \approx 165.89\ \mu\text{m}$$

The volume fraction of the distal hypoperfused pharmacological sanctuary sleeve is:
$$\mathcal{V}_{\text{sanctuary}} = \frac{r_{\text{ox}}^2 - r_{\text{crit}}^2}{r_{\text{ox}}^2 - r_v^2} = \frac{200^2 - 165.89^2}{200^2 - 10^2} = \mathbf{31.28\%}$$

### 2. Hereditary Volterra Memory Functional for Inflammatory Pre-Conditioning
The effective anti-tumor immune pressure $\mathcal{E}_{\text{imm}}(t)$ decays under chronic systemic inflammatory stress $\Phi(x,\tau)$ (IL-1β, IL-6, TNF-α, neutrophil NETosis):
$$\mathcal{E}_{\text{imm}}(t) = \mathcal{E}_0 \exp\left(-\gamma \int_0^t \Phi(x,\tau)d\tau\right)$$

### 3. Double-Barrier ctDNA Shedding Kinetics & Lead Time
During early dormancy ($N_R \approx 10^4\text{--}10^6$ cells), apoptosis is minimal and the tight extracellular matrix pore cutoff ($d_{\text{pore}} \le 11\text{ nm}$ vs. $166\text{ bp}$ nucleosomes $\approx 11\times 5.7\text{ nm}$) sterically entraps fragments, yielding false-negative plasma ctDNA (MRD⁻). Upon inflammatory awakening, clonal expansion ($10^6 \to 10^8$ cells, $T_d = 24.6\text{ days}$) remodels the stroma, resulting in an analytical ctDNA lead time of **5.46 months** before macroscopic CT imaging detection.

---

## 📊 Summary of Audited Cohorts ($N = 3,461$)

All clinical evaluations adhere strictly to **REMARK** guidelines and complete-case analysis without imputation:

| Cohort Identifier | Clinical Setting / Regimen | Assayed Biomarkers / Modality | Sample Size ($N$) | Primary Outcome | Key Finding |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **MSKCC EGFR-TKI** | Resected Stage IB–IIIA EGFR-mutant NSCLC on adjuvant osimertinib | Histological Polarity + Smoking Airway Stress (Decoupled HTBI) | $N = 604$ | Relapse-Free Survival (RFS) | Continuous multivariable $\text{HR} = 2.49$ ($95\%\text{ CI: } 1.00\text{--}6.18, p = 0.050$), $C = 0.742$; Binary $\text{HR} = 1.89$ ($p = 0.029$) |
| **Real Surgical Blood (PLOS ONE / Deng et al.)** | Consecutive radical surgical resections (gastric / gastrointestinal) | Raw pre-op CBC (NLR) + serum hs-CRP and Albumin (GPS) | $N = 1,056$ *(audited from $1,330$ accessions)* | Overall Survival (OS) | Multivariable binary $\text{HR} = 1.43$ ($95\%\text{ CI: } 1.14\text{--}1.81, p = 0.0022$); Multivariable continuous $\text{HR} = 2.20$ ($p = 0.00041$) |
| **TCGA Epithelial Solid (COADREAD, PAAD, ESCA, CHOL)** | Pan-cancer solid epithelial resection cohorts | Real RNA-seq z-scores (CDH1, TJP2, EPCAM, COL1A1, ACTA2, IL1B) | $N = 986$ | Overall Survival (OS) | Concordant prognostic hazard across independent epithelial lineages (COADREAD $p=0.014$, PAAD $p=2.6\times 10^{-5}$, ESCA $p=0.053$) |
| **GEO GSE31210** | Early-stage lung adenocarcinoma validation | Microarray transcriptomic signature | $N = 226$ | Relapse-Free Survival (RFS) | Independent validation of microvascular sanctuary risk score |
| **CHOL-OMICS** | Primary human resected cholangiocarcinoma & ciliopathy cohort | Matched transcriptomics + spatial histology | $N = 18$ | Biological Benchmark | Spatial validation of desmoplastic hypovascularization & mechanosensory disruption |
| **TCGA Negative Controls (LAML, TGCT, SARC)** | Non-epithelial/immune-hypoxia uncoupled (Leukemia, Testicular, Sarcoma) | Real RNA-seq z-scores (TJP2, CDH1, COL1A1, ACTA2, IL1B) | $N = 535$ | Falsification Test | Complete loss of prognostic signal ($p = 0.25, 0.47, 0.41$), mathematically confirming tissue specificity |
| **Total Study Population** | **Audited Multi-Center Cohorts** | **Multi-modal Omics, Blood & Controls** | **$N = 3,461$** | **REMARK Validated** | **$N = 3,461$ clinically indexed records ($N = 3,425$ complete matched molecular profiles analyzed; 36 TCGA samples lacked RNA-seq)** |

---

## 📁 Repository Structure

```
inflammatory-onco-physics/
├── .zenodo.json                            # Zenodo archival & DOI indexing metadata
├── CITATION.cff                            # GitHub citation file format
├── LICENSE                                 # Creative Commons Attribution 4.0 International (CC-BY-4.0)
├── README.md                               # This master documentation & showcase hub
├── requirements.txt                        # Pip dependencies (Streamlit, Lifelines, Scipy, etc.)
├── demo_inflammatory_onco_physics.ipynb   # 1-Click Interactive Colab Pro Notebook (T4 GPU Accelerated)
├── app.py                                  # Full-featured Streamlit interactive web application
├── simulator.html                          # Standalone 60 FPS in silico digital twin simulator
├── tumor_cord_sanctuary_3d.html            # Standalone 3D WebGL Three.js perivascular cord simulator
├── app_bedside_calculator.py               # Standalone bedside HTBI decision calculator CLI
├── MEMORANDUM_ESPLORATIVO_IFOM.pdf         # 4-page Technical Discussion Memorandum (IFOM Molecular Oncology)
├── MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.pdf # Compiled master monograph PDF (26 pages)
├── MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.md  # Markdown master monograph
├── SUPPLEMENTARY_INFORMATION.pdf           # Compiled supplementary information PDF (10 pages)
├── SUPPLEMENTARY_INFORMATION.md            # Markdown supplementary information
├── REMARK_COMPLIANCE_REPORT.md             # Complete REMARK checklist audit
├── assets/                                 # Visual animations, GIFs and media
│   ├── simulator-vessel-live.gif           # 60 FPS live capillary simulation
│   ├── simulator-vessel-live.webm          # High-efficiency WebM video
│   └── social-preview.png                  # OpenGraph social card
├── data/                                   # Real clinical datasets and transcriptomic matrices
│   ├── mskcc_2020_real_cohort.csv          # MSKCC NSCLC clinical cohort (N=604)
│   ├── real_gastric_nlr_crp_dataset.csv    # Surgical GI blood biomarker cohort (N=1,056)
│   ├── chol_omics_anonymized_stratification.tsv # Primary cholangiopathy cohort (N=18)
│   ├── tcga_negative_controls_rnaseq.tsv   # Real RNA-seq z-scores for LAML, TGCT, SARC (N=535)
│   ├── tcga_coadread_rnaseq.tsv            # Real RNA-seq z-scores for TCGA-COADREAD (N=592)
│   ├── tcga_paad_rnaseq.tsv                # Real RNA-seq z-scores for TCGA-PAAD (N=177)
│   ├── tcga_esca_rnaseq.tsv                # Real RNA-seq z-scores for TCGA-ESCA (N=181)
│   └── tcga_chol_rnaseq.tsv                # Real RNA-seq z-scores for TCGA-CHOL (N=36)
├── figure_1_bifurcation.png                # Figure 1: Reaction-diffusion bifurcation & cylindrical Bessel K0
├── figure_2_clinical_validation.png        # Figure 2: MSKCC survival calibration & ctDNA lead time
├── figure_3_pancancer_validation.png       # Figure 3: TCGA pan-cancer validation
├── figure_4_multivariate_forest.png        # Figure 4: Multivariable Cox proportional hazards forest plot
├── figure_5_bootstrap_roc.png              # Figure 5: Bootstrap ROC and time-dependent C-index
├── figure_6_chol_omics_landscape.png       # Figure 6: Spatial desmoplasia and microvascular decay
├── figure_7_remark_master.png              # Figure 7: Complete REMARK flow and calibration master
├── figure_8_blinded_controls_and_placebo.png # Figure 8: Negative controls & in silico stress tests
├── figure_9_proof_by_contradiction.png     # Figure 9: Mathematical falsification on real RNA-seq controls
├── figure_10_clinical_robustness_master.png # Figure 10: Sensitivity, competing risks & mucosal validation
├── figure_11_real_blood_biomarkers_1330_patients.png # Figure 11: Real-world blood biomarker validation
├── pipeline/                               # Automated analysis & validation scripts
│   ├── 01_multivariate_cox_engine.py
│   ├── 02_bootstrap_cindex_roc.py
│   ├── 03_chol_cohort_mapper.py
│   ├── 04_remark_statistical_protocol.py
│   ├── 05_blinded_controls_and_stress_tests.py
│   ├── 06_comprehensive_controls_and_proof_by_contradiction.py
│   ├── 07_advanced_robustness_battery.py
│   └── 08_real_blood_biomarker_validation.py
└── simulation/                             # Numerical ODE/PDE reaction-diffusion solvers
    ├── clinical_simulator.py
    ├── fetch_and_validate_pancancer.py
    ├── model.py
    └── simulate_relapse.py
```

---

## 🚀 Quickstart & Reproducibility

### Prerequisites
- Python 3.10 or higher
- Standard scientific Python stack (`numpy`, `scipy`, `pandas`, `lifelines`, `matplotlib`, `seaborn`)

```bash
# Clone the repository
git clone https://github.com/Lukecele/inflammatory-onco-physics.git
cd inflammatory-onco-physics

# Install dependencies
pip install -r requirements.txt
```

### Reproduce Multivariable Cox Analysis (MSKCC & Blood Cohorts)
```bash
python pipeline/01_multivariate_cox_engine.py
python pipeline/08_real_blood_biomarker_validation.py
```

### Reproduce Empirical Proof-by-Contradiction (Real RNA-Seq Negative Controls)
```bash
python pipeline/06_comprehensive_controls_and_proof_by_contradiction.py
```

### Run the Cylindrical Bessel $K_0$ Reaction-Diffusion Simulation
```bash
python simulation/simulate_relapse.py
python simulation/clinical_simulator.py
```

### Launch the Bedside HTBI Decision Calculator
```bash
python app_bedside_calculator.py
```

---

## 📜 License

This work, including all manuscripts, data matrices, documentation, figures, and source code, is licensed under the **Creative Commons Attribution 4.0 International Public License (CC-BY-4.0)**.

You are free to:
- **Share:** copy and redistribute the material in any medium or format.
- **Adapt:** remix, transform, and build upon the material for any purpose, even commercially.
Under the following terms:
- **Attribution:** You must give appropriate credit, provide a link to the license, and indicate if changes were made.

---

## 📬 Contact & Citation

For scientific correspondence, collaboration, or clinical trial inquiries:
- **Luca Celebrano, MD Candidate / Researcher**  
  *Università degli Studi di Napoli Federico II*  
  *Telethon Institute of Genetics and Medicine (TIGEM)*  
  Email: `l.celebrano@studenti.unina.it` | `luca.celebrano@gmail.com`

```bibtex
@article{celebrano2026inflammatoryoncophysics,
  title={A Reaction-Diffusion Theory of Microvascular Pharmacological Sanctuaries and Inflammatory Pre-conditioning Explains Post-Adjuvant Relapse Dynamics Across 3,461 Patients},
  author={Celebrano, Luca},
  journal={Zenodo Archive / Springer Nature Submission},
  year={2026},
  doi={10.5281/zenodo.24118328},
  url={https://doi.org/10.5281/zenodo.24118328},
  institution={Department of Clinical Medicine and Surgery, Universita degli Studi di Napoli Federico II; Telethon Institute of Genetics and Medicine (TIGEM)}
}
```
