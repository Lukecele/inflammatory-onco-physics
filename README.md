# A Reaction-Diffusion Theory of Microvascular Pharmacological Sanctuaries and Inflammatory Pre-conditioning Explains Post-Adjuvant Relapse Dynamics Across 3,461 Patients

[![Research Square Preprint](https://img.shields.io/badge/Research_Square-10.21203%2Frs.3.rs--11695210%2Fv1-blue.svg)](https://doi.org/10.21203/rs.3.rs-11695210/v1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24118328.svg)](https://doi.org/10.5281/zenodo.24118328)
[![YouTube Short Part 1](https://img.shields.io/badge/YouTube_Short-Part_1:_Sanctuaries-red.svg?logo=youtube&logoColor=white)](https://youtube.com/shorts/2IxHC94_Z18)
[![YouTube Short Part 2](https://img.shields.io/badge/YouTube_Short-Part_2:_ctDNA_Breach-red.svg?logo=youtube&logoColor=white)](https://youtube.com/shorts/m5TXgQOB2WA)
[![YouTube Short Part 3](https://img.shields.io/badge/YouTube_Short-Part_3:_Desmoplastic_Shield-red.svg?logo=youtube&logoColor=white)](https://youtube.com/shorts/o4tKdbQV9cU)
[![YouTube Short Part 4](https://img.shields.io/badge/YouTube_Short-Part_4:_Antibody_Strike-red.svg?logo=youtube&logoColor=white)](https://youtube.com/shorts/9EMsP73WSho)
[![YouTube Short Part 5](https://img.shields.io/badge/YouTube_Short-Part_5:_5.5M_Lead_Time-red.svg?logo=youtube&logoColor=white)](https://youtube.com/shorts/0HPGny83MRU)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![REMARK Compliant](https://img.shields.io/badge/REMARK-Compliant-success.svg)](REMARK_COMPLIANCE_REPORT.md)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Author:** Luca Celebrano<sup>1,2</sup>  
<sup>1</sup> Department of Clinical Medicine and Surgery, School of Medicine and Surgery, Università degli Studi di Napoli Federico II, Naples, Italy  
<sup>2</sup> Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy  
*Correspondence:* [l.celebrano@studenti.unina.it](mailto:l.celebrano@studenti.unina.it) / [luca.celebrano@gmail.com](mailto:luca.celebrano@gmail.com) | ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
*Research Hub:* [Luca's Cave (lucacelebrano.me)](https://lucacelebrano.me/)

---

## 📌 Overview & Scientific Scope

This repository provides the complete, open-access computational framework, biophysical continuous reaction-diffusion solvers, biostatistical survival pipelines, and clinical cohort matrices supporting the monograph:

> **Celebrano, L. (2026).** *A Reaction-Diffusion Theory of Microvascular Pharmacological Sanctuaries and Inflammatory Pre-conditioning Explains Post-Adjuvant Relapse Dynamics Across 3,461 Patients.* Springer Nature In Review / Research Square Preprint [DOI: 10.21203/rs.3.rs-11695210/v1](https://doi.org/10.21203/rs.3.rs-11695210/v1) & Zenodo Archive [DOI: 10.5281/zenodo.24118328](https://doi.org/10.5281/zenodo.24118328).
>
> 📄 **Executive Memorandum:** A concise 4-page technical discussion memorandum is available directly in this repository: [`MEMORANDUM_ESPLORATIVO_IFOM.pdf`](MEMORANDUM_ESPLORATIVO_IFOM.pdf).  
> 📚 **Full Monograph:** The complete 26-page monograph is available as [`MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.pdf`](MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.pdf) with [`SUPPLEMENTARY_INFORMATION.pdf`](SUPPLEMENTARY_INFORMATION.pdf).

### The Biological Problem
Adjuvant targeted therapies (e.g., third-generation EGFR-TKIs like osimertinib in the ADAURA trial) demonstrate profound disease-free survival (DFS) benefits during active administration, yet post-cessation survival trajectories exhibit a characteristic rebound recurrence within 12–18 months. Crucially, 75–80% of relapsing lesions harbor **no secondary resistance mutations** (e.g., EGFR C797S).

We reconcile this clinical paradox through a **continuum biophysical reaction-diffusion model** coupled to **hereditary Volterra immune exhaustion**:
1. **Microvascular Pharmacological Sanctuaries:** Due to cylindrical radial diffusion decay ($K_0$ Bessel solution) and metabolic clearance ($k \approx 1.51 \times 10^{-4}\text{ s}^{-1}$), free drug concentrations drop below therapeutic thresholds ($\text{IC}_{50} \approx 10\text{ nM}$) at a radial distance $r_{\text{crit}} \approx 165.89\ \mu\text{m}$ from the microvessel wall ($r_v \approx 10\ \mu\text{m}$). This creates a distal hypoperfused sanctuary sleeve ($165.89 < r \le 200\ \mu\text{m}$) harboring drug-tolerant persister (DTP) cells, accounting for **31.28%** of viable tumor cord volume.
2. **Inflammatory Pre-conditioning:** Persisters maintain dormant viability under metabolic restriction ($G_0/G_1$ arrest) but undergo epigenetic and phenotypic re-awakening under the influence of systemic inflammatory cytokines (IL-1β, IL-6, TNF-α) and NETosis, modeled via a memory functional $\mathcal{E}_{\text{imm}}(t) = \mathcal{E}_0 \exp(-\gamma \int_0^t \Phi(x,\tau)d\tau)$.
3. **Double-Barrier ctDNA Shedding Kinetics:** In early dormancy ($N_R \approx 10^4\text{--}10^6$ cells), negligible apoptosis and tight steric entrapment in dense extracellular matrix (pore size $\le 11\text{ nm}$ vs. $166\text{ bp}$ nucleosomes) lead to false-negative plasma ctDNA (MRD$^{-}$). Upon inflammatory awakening and rapid clonal expansion ($10^6 \to 10^8$ cells, $T_d = 24.6\text{ days}$), the model predicts an analytical ctDNA lead time of **5.46 months** prior to macroscopic imaging detection.

---

## ⚡ Interactive Google Colab Notebook & 3D Simulation

### 🎬 Official 3D Simulation Series (YouTube Shorts Quintilogy)

| 🔴 **Part 1: Microvascular Sanctuaries** | 🔴 **Part 2: The Matrix Breach & ctDNA** | 🔴 **Part 3: The Desmoplastic Shield** | 🔴 **Part 4: Targeted Antibody Strike** | 🔴 **Part 5: 5.5-Month Lead Time** |
| :---: | :---: | :---: | :---: | :---: |
| [![Part 1: Why Cancers Relapse](https://img.shields.io/badge/▶_Part_1-Sanctuaries-red.svg?style=for-the-badge&logo=youtube&logoColor=white)](https://youtube.com/shorts/2IxHC94_Z18) | [![Part 2: The Matrix Breach](https://img.shields.io/badge/▶_Part_2-Matrix_Breach-red.svg?style=for-the-badge&logo=youtube&logoColor=white)](https://youtube.com/shorts/m5TXgQOB2WA) | [![Part 3: Desmoplastic Shield](https://img.shields.io/badge/▶_Part_3-Desmoplastic_Shield-red.svg?style=for-the-badge&logo=youtube&logoColor=white)](https://youtube.com/shorts/o4tKdbQV9cU) | [![Part 4: Antibody Strike](https://img.shields.io/badge/▶_Part_4-Antibody_Strike-red.svg?style=for-the-badge&logo=youtube&logoColor=white)](https://youtube.com/shorts/9EMsP73WSho) | [![Part 5: 5.5M Lead Time](https://img.shields.io/badge/▶_Part_5-5.5M_Lead_Time-red.svg?style=for-the-badge&logo=youtube&logoColor=white)](https://youtube.com/shorts/0HPGny83MRU) |
| *Bessel $K_0$ drug decay & DTP sanctuary sleeve* | *Collagen pore cutoff ($\le 11$ nm) & helical ctDNA entry* | *CAF microvessel squeeze & arrested convection* | *Monoclonal antibodies penetrating matrix to bind tumor cells* | *ctDNA MRD⁺ flags relapse 5.46 months before CT scan* |

### 🚀 Google Colab (Pro T4 GPU)
Run the complete biophysical solvers, Bessel $K_0$ radial simulations, and clinical survival regressions directly in your browser with zero local installation:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lukecele/inflammatory-onco-physics/blob/main/demo_inflammatory_onco_physics.ipynb)

* Notebook file: [`demo_inflammatory_onco_physics.ipynb`](demo_inflammatory_onco_physics.ipynb)
* Hardware acceleration: automatically detects and utilizes NVIDIA T4 GPU (Google Colab Pro) or runs efficiently on standard CPU.

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
├── README.md                               # Master study documentation
├── requirements.txt                        # Python dependencies
├── demo_inflammatory_onco_physics.ipynb   # 1-Click Interactive Colab Notebook (T4 GPU Accelerated)
├── MEMORANDUM_ESPLORATIVO_IFOM.pdf         # 4-page Technical Discussion Memorandum (IFOM Molecular Oncology)
├── MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.pdf # Compiled master monograph PDF (26 pages)
├── MANOSCRITTO_SPRINGER_NATURE_26_PAGINE.md  # Markdown master monograph
├── SUPPLEMENTARY_INFORMATION.pdf           # Compiled supplementary information PDF (10 pages)
├── SUPPLEMENTARY_INFORMATION.md            # Markdown supplementary information
├── REMARK_COMPLIANCE_REPORT.md             # Complete REMARK checklist audit
├── app_bedside_calculator.py               # Standalone bedside HTBI decision calculator CLI
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
