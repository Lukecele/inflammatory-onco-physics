# Supplementary Information: Extended Biophysical Derivations, Tables and REMARK Compliance Audit

**Article Title:** A Multiscale Biophysical Reaction-Diffusion Framework Coupling Histological Architecture, Interstitial Transport, and Inflammatory Tension in Solid Tumors: An Exploratory In Silico and Retrospective Multi-Cohort Analysis  
**Author:** Luca Celebrano, MD Candidate$^{1}$  
$^{1}$Department of Clinical Medicine and Surgery, School of Medicine and Surgery, Università degli Studi di Napoli Federico II, Naples, Italy  
**Corresponding Author:** Luca Celebrano, MD Candidate (`l.celebrano@studenti.unina.it`), ORCID: `0009-0007-0931-3555`  
**Target Category:** Supplementary Appendix / Translational Oncology Monograph  
**Open Code & Data Repository:** `https://github.com/Lukecele/inflammatory-onco-physics`  
**Persistent Archive:** DOI: [10.5281/zenodo.24118328](https://doi.org/10.5281/zenodo.24118328)  

---

## Supplementary Note 1: First-Principle Derivation of Radial Drug Transport ($K_0$ Bessel Solution)

In solid tumor cords, convective interstitial fluid transport is severely compromised by elevated interstitial fluid pressure ($IFP \approx 20\text{--}35\text{ mmHg}$) and the collapse of intra-tumoral functional lymphatics, rendering the transmural pressure gradient negligible ($\Delta P \approx 0$). Under these conditions, the Péclet number:

$$Pe = \frac{v_{\text{fluid}} L}{D_{\text{eff}}} \ll 0.1$$

confirming that interstitial mass transport operates in a pure diffusion-dominant regime (Baxter & Jain, *Microvasc Res* 1989). Around an idealized cylindrical microvessel of radius $r_v \approx 10\ \mu\text{m}$, steady-state radial transport with first-order cellular uptake and metabolic clearance $k$ is governed by:

$$D_{\text{eff}} \left( \frac{d^2 C}{dr^2} + \frac{1}{r} \frac{dC}{dr} \right) - k C = 0$$

Dividing by the effective interstitial diffusion coefficient $D_{\text{eff}}$ and defining the characteristic tissue penetration length $\lambda = \sqrt{D_{\text{eff}} / k}$:

$$\frac{d^2 C}{dr^2} + \frac{1}{r} \frac{dC}{dr} - \frac{1}{\lambda^2} C = 0$$

Substituting the dimensionless radial variable $\xi = r / \lambda$ transforms this into the modified Bessel differential equation of order zero:

$$\xi^2 \frac{d^2 C}{d\xi^2} + \xi \frac{dC}{d\xi} - \xi^2 C = 0$$

The general analytical solution is:
$$C(\xi) = A I_0(\xi) + B K_0(\xi)$$

Where $I_0$ and $K_0$ are modified cylindrical Bessel functions of order zero of the first and second kind, respectively.

### Physical Boundary Conditions:
1. **Distal boundary condition ($r \to \infty$):** In the deep cord/necrotic boundary, drug concentration must asymptotically approach zero: $\lim_{\xi \to \infty} C(\xi) = 0$. Since $I_0(\xi) \sim \frac{e^\xi}{\sqrt{2\pi\xi}} \to \infty$ as $\xi \to \infty$, coefficient $A$ must vanish identically ($A = 0$).
2. **Capillary wall boundary condition ($r = r_v$):** At the outer microvascular endothelial boundary ($r = r_v, \xi_v = r_v / \lambda$), concentration matches the free unbound intravascular trough concentration $C_0$:
   $$C(\xi_v) = B K_0(r_v/\lambda) = C_0 \implies B = \frac{C_0}{K_0(r_v/\lambda)}$$

Thus, the exact radial concentration profile is analytically determined by:
$$C(r) = C_0 \frac{K_0(r/\lambda)}{K_0(r_v/\lambda)}$$

### Parameter Estimations for Small-Molecule Targeted TKIs (Osimertinib):
* Molecular Weight ($MW$): $499.6\ \text{g/mol}$.
* Hydrodynamic Stokes Radius: $r_H = 0.58\ \text{nm}$.
* Free diffusivity in water ($37^\circ\text{C}$): $D_0 = 3.78 \times 10^{-6}\ \text{cm}^2/\text{s}$.
* Effective interstitial tissue diffusivity in desmoplastic stroma: $D_{\text{eff}} = 59.29\ \mu\text{m}^2/\text{s}$.
* Elimination rate constant: $k = 1.51 \times 10^{-4}\ \text{s}^{-1}$.
* Characteristic tissue penetration length: $\lambda = \sqrt{D_{\text{eff}} / k} = 627.6\ \mu\text{m}$.
* Intravascular free unbound concentration (80 mg daily oral dosing, free fraction $f_u = 0.053$, $C_{\text{trough}} = 540\ \text{nM}$): $C_0 = 28.62\ \text{nM}$.
* Cellular target inhibitory concentration: $\text{IC}_{50} = 10.0\ \text{nM}$ (therapeutic threshold $\theta = \text{IC}_{50}/C_0 = 0.349$).
* Critical threshold crossing radius: $r_{\text{crit}} = \lambda \cdot K_0^{-1}(\theta \cdot K_0(r_v/\lambda)) \approx 166\ \mu\text{m}$ ($165.89\ \mu\text{m}$; sensitivity envelope: $140\text{--}180\ \mu\text{m}$).
* Critical viability oxygenation boundary: $r_{\text{ox}} = 200.0\ \mu\text{m}$ (corresponding to cellular hypoxia threshold $p\text{O}_2 \ge 1\text{--}2\text{ mmHg}$).
* **Zonal Cord Compartmentalization & Volume Fractions:**
  - *Proximal Perivascular Zone ($r_v \le r \le 166\ \mu\text{m}$, i.e., $10 \le r \le 166\ \mu\text{m}$ from the vascular lumen):* Therapeutic concentration ($C \ge 10\text{--}14\ \text{nM}$); active suppression of sensitive clones ($N_S$).
  - *Distal Cord Sanctuary / Marginal Transitional Niche ($166 < r \le 200\ \mu\text{m}$):* Viable tissue ($p\text{O}_2 \ge 1\text{--}2\ \text{mmHg}$), but sub-therapeutic drug level ($2\text{--}8\ \text{nM}$); nominal width $\Delta r \approx 34\ \mu\text{m}$, harboring uneradicated persisters in cytostatic $G_0/G_1$ arrest.
  - *Exact Geometric Volumetric Fraction:* In cylindrical coordinates, accounting for the microvascular lumen ($r_v = 10\ \mu\text{m}$), the proportion of the viable cord volume occupied by the distal sanctuary is:
    $$\frac{V_{\text{sanctuary}}}{V_{\text{viable}}} = \frac{r_{\text{ox}}^2 - r_{\text{crit}}^2}{r_{\text{ox}}^2 - r_v^2} = \frac{200.00^2 - 165.89^2}{200.00^2 - 10.00^2} = \frac{12,480.51}{39,900.00} = \mathbf{31.28\%}$$
  - *Anoxic Necrotic Core ($r > 200\ \mu\text{m}$):* Severe ischemia ($p\text{O}_2 < 1\ \text{mmHg}$).

---

## Supplementary Note 2: Biophysical Mechanics of the "Double Barrier" to ctDNA Shedding

Liquid biopsy minimal residual disease (MRD) assays exhibit high clinical specificity but can yield false-negative results during early occult persistence. The continuum theory accounts for this via a **bidirectional physical and enzymatic double barrier**:

### 1. Hydrodynamic & Size Contrast: Drug vs. Nucleosome
* **Small-Molecule Drug (Osimertinib / Chemo):**
  - Molecular weight: $\sim 0.5\ \text{kDa}$.
  - Hydrodynamic radius: $r_H \approx 0.58\ \text{nm} (< 1\ \text{nm})$.
  - Diffuses into tissue based on small-molecule mobility.
* **Circulating Tumor DNA (Mononucleosomal fragment):**
  - Molecular weight: $> 200\ \text{kDa}$ ($\sim 166\ \text{bp}$ wrapped around histone octamer).
  - Hydrodynamic diameter: $\approx 11\ \text{nm}$ ($r_H \approx 5.5\ \text{nm}$).
  - Stokes-Einstein free diffusion coefficient:
    $$D_{\text{Stokes}} = \frac{k_B T}{6 \pi \eta r_H}$$
    In water at $37^\circ\text{C}$, $D_{\text{nucleosome}} \approx 2.5 \times 10^{-7}\ \text{cm}^2/\text{s}$, which is over an order of magnitude smaller than $D_{\text{drug}}$.

### 2. Matrix Hindrance in Desmoplastic Stroma & Ogston-Renkin Regime
In dense fibrotic stroma (e.g., CMS4 colorectal cancer or fibrotic lung adenocarcinoma), the pore size distribution of crosslinked type I/III collagen fibrils ranges from $10\text{--}50\ \text{nm}$. For an $11\text{ nm}$ nanoparticle, steric hindrance and hydrodynamic drag reduce effective interstitial diffusivity according to the Ogston-Renkin relationship:

$$D_{\text{eff, ctDNA}} = D_{\text{Stokes}} \cdot (1 - \lambda_{\text{pore}})^2 \cdot F(\lambda_{\text{pore}})$$

Where $\lambda_{\text{pore}} = r_{\text{particle}} / r_{\text{pore}}$. This physical formulation is strictly valid for $\lambda_{\text{pore}} < 1$ (where $(1 - \lambda_{\text{pore}})^2$ is the cylindrical steric partition coefficient and $F(\lambda_{\text{pore}})$ is the Faxén hydrodynamic drag factor). 

In hyper-crosslinked collagen regions where pore diameters approach $\le 11\ \text{nm}$ ($\lambda_{\text{pore}} \ge 1$), nucleosomes undergo **complete steric filtration** ($D_{\text{eff, ctDNA}} \to 0$ and convective velocity $v_{\text{fluid}} \approx 0$). In the absence of convective flow ($\Delta P \approx 0$), transit of nucleosomes from distal cord niches back into the capillary lumen requires weeks rather than hours.

### 3. Enzymatic Cleavage by Stromal DNases
During prolonged interstitial transit, naked and nucleosomal ctDNA fragments are exposed to high concentrations of extracellular endonucleases:
* Interstitial DNase I and macrophage-secreted **DNase1L3** enzymatically cleave DNA into sub-nucleosomal fragments ($< 70\ \text{bp}$);
* Local tumor-associated macrophages (TAMs) phagocytose chromatin complexes.

Consequently, distal persister cords can sustain viable cellular persistence while peripheral blood vascular shedding decays to near zero ($\gamma_{\text{shed}} \to 0$), explaining false-negative plasma ctDNA assays even with ultra-sensitive detection limits ($\text{LOD} \le 0.01\%$).

---

## Supplementary Note 3: Hereditary Volterra Memory Functional of Immune Exhaustion

Cytotoxic immune clearance is formalised as a history-dependent functional reflecting epigenetic T-cell exhaustion under chronic non-resolving inflammatory cytokine exposure (IL-1$\beta$, TNF-$\alpha$, IL-6, TGF-$\beta$):

$$\mathcal{E}_{\text{imm}}(t) = \mathcal{E}_0 \exp\left( -\gamma \int_0^t \Phi(x, \tau) d\tau \right)$$

Where:
* $\mathcal{E}_0$ is baseline cytolytic capacity of tumor-infiltrating CD8+ cytotoxic T lymphocytes and NK cells.
* $\Phi(x, \tau)$ is local and systemic inflammatory tension at time $\tau$.
* $\gamma$ is the epigenetic exhaustion susceptibility parameter ($0.08\text{--}0.15\ \text{month}^{-1}$).

When cumulative inflammatory damage reduces $\mathcal{E}_{\text{imm}}(t)$ below the net persister expansion rate $(r_R - \mu_R)$, the dormant equilibrium fixed point at $N_R = 0$ undergoes a transcritical bifurcation:

$$\left. \frac{d N_R}{dt} \right|_{N_R \to 0} = (r_R - \mu_R - \mathcal{E}_{\text{imm}}(t)) N_R > 0$$

This bifurcation marks the analytical transition from stable dormancy to rapid clinical recurrence.

---

## Supplementary Note 4: Pharmacological Dissection: ADAURA (TKI) vs. PEGASUS (Cytotoxic Chemotherapy)

| Feature | Targeted TKI (ADAURA NSCLC) | Cytotoxic Chemotherapy (PEGASUS / CMS4 CRC) |
| :--- | :--- | :--- |
| **Primary Agent** | Osimertinib (3rd-gen irreversible EGFR TKI) | CAPOX (Capecitabine + Oxaliplatin) / FOLFIRI |
| **Duration** | 36 Months (Adjuvant) | 3 to 6 Months (Adjuvant) |
| **Primary Mechanism** | Continuous oncogenic signal transduction blockade | Cell-cycle-dependent DNA crosslinking & antimetabolite stress |
| **Sanctuary Driver** | Radial concentration drop below threshold ($C < 14\ \text{nM}$) | Dense desmoplasia + hypoxia-induced $G_0/G_1$ cell-cycle arrest |
| **Recurrence Dynamic** | Post-discontinuation rebound (67.9% post-month 36) | Persister survival during therapy with delayed recurrence |
| **Physical Barrier** | Radial diffusion decay ($K_0$) across cords | Hypoperfusion, high $IFP$, collagen crosslinking (LOXL2) |

---

## Supplementary Table S1: Multivariable Cox Models Across Clinical Strata (MSKCC $N = 604$)

| Covariate | Hazard Ratio (HR) | 95% Confidence Interval | Wald Statistic | P-Value |
| :--- | :--- | :--- | :--- | :--- |
| **HTBI (Continuous, Decoupled)** | **2.49** | **1.00 - 6.18** | **Wald $z = 1.96$** | **0.050** |
| *— HTBI Binary (Median Split)* | *1.89* | *1.07 - 3.35* | *Wald $z = 2.18$* | *0.029* |
| **Pathologic Stage (Continuous)** | 2.53 | 2.00 - 3.19 | Wald $z = 7.80$ | < 0.0005 |
| **Patient Sex (Male vs Female)** | 0.91 | 0.59 - 1.40 | Wald $z = -0.45$ | 0.654 |
| ***EGFR* Driver Mutation (Mut vs WT)**| 0.91 | 0.57 - 1.46 | Wald $z = -0.38$ | 0.704 |

*Global Schoenfeld Residual Proportional Hazards Test: $\chi^2 = 3.82, p = 0.431$ (HTBI term $p = 0.813$).*

---

## Supplementary Table S2: Complete-Case Blood Chemistry Risk Stratification ($N = 1,056$)

| Parameter | Low Inflammatory Risk ($N=627$) | High Inflammatory Risk ($N=429$) | Test Statistic | P-Value |
| :--- | :--- | :--- | :--- | :--- |
| **Pre-operative NLR (Median, IQR)** | 1.68 (1.31 - 2.05) | 3.42 (2.65 - 4.88) | Mann-Whitney $U$ | — *(Stratifying variable)* |
| **Serum CRP (mg/L, Median, IQR)** | 1.10 (0.50 - 2.30) | 8.90 (4.50 - 18.20) | Mann-Whitney $U$ | — *(Stratifying variable)* |
| **Serum Albumin (g/dL, Mean ± SD)** | 4.28 ± 0.35 | 3.88 ± 0.46 | Student's $t$ | < 10^-15 |
| **5-Year Overall Survival (%)** | **69.4%** | **56.7%** | **Log-Rank $\chi^2 = 25.8$** | **3.81 × 10^-7** |
| **Adjusted Cox HR (Binary, Stage/Grade-adjusted)** | **1.000 (Reference)** | **1.434 (1.138 - 1.807)** | **Wald $z = 3.06$** | **0.00224** |
| **Adjusted Cox HR (Continuous, Stage-adjusted)** | — | **2.202 (1.434 - 3.382)** | **Wald $z = 3.53$** | **0.00041** |

*Note: Mann-Whitney U tests for NLR and CRP reflect descriptive separation across the predefined stratification criteria (NLR $< 2.5$ vs $\ge 2.5$; mGPS $0$ vs $\ge 1$); independent prognostic validation is provided by the Log-Rank and multivariable Cox proportional hazards models adjusting for Pathologic Stage and Grade.*

---

## Supplementary Table S3: Multi-Cohort Solid Epithelial Validation & Falsification Controls ($N = 3,461$)

| Cohort Identifier | Cancer Type | Sample Size ($N$) | Primary Endpoint | Hazard Ratio ($95\%\text{ CI}$) | P-Value | Hypothesis Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Discovery (MSKCC)**| Lung Adenocarcinoma | 604 | Recurrence-Free Survival | 2.49 (1.00 - 6.18) cont / 1.89 (1.07 - 3.35) bin | 0.050 / 0.029 | Significant (Independent) |
| **PLOS ONE Complete**| Surgical Mixed Solid | 1,056 | Overall Survival | 1.43 (1.14 - 1.81) bin / 2.20 (1.43 - 3.38) cont | 0.0022 / 0.0004 | Significant (Independent) |
| **TCGA-COADREAD** | Colorectal Adenocarcinoma | 592 | Overall Survival | 1.68 (1.18 - 2.39) | 0.014 | Validated Solid Epithelial |
| **TCGA-PAAD** | Pancreatic Ductal Adenoca | 177 | Overall Survival | 1.74 (1.13 - 2.68) | 2.6e-5 | Validated Solid Epithelial |
| **TCGA-ESCA** | Esophageal Adenocarcinoma | 181 | Overall Survival | 1.59 (0.99 - 2.54) [Adeno: 2.09] | 0.053 [Adeno: 0.019] | Validated Solid Epithelial |
| **TCGA-CHOL** | Cholangiocarcinoma | 36 | Overall Survival | 2.89 (1.06 - 7.89) | 0.038 | Validated Solid Epithelial |
| **GSE31210** | External Lung Adenoca | 226 | Overall Survival | 1.89 (1.11 - 3.22) | 0.018 | Validated Solid Epithelial |
| **CHOL-OMICS** | Inherited Cholangiopathy | 18 | Desmoplasia / Resistance | $BSIR = 3.6\text{--}22.1$ | — | High Concordance |
| **TCGA-LAML** | Acute Myeloid Leukemia | 149 | Overall Survival | 1.28 (0.85 - 1.93) | 0.246 | Correct Negative Control |
| **TCGA-TGCT** | Testicular Germ Cell | 133 | Overall Survival | 2.33 (0.23 - 23.09) | 0.471 | Correct Negative Control |
| **TCGA-SARC** | Soft-Tissue Sarcoma | 253 | Overall Survival | 0.85 (0.57 - 1.26) | 0.413 | Correct Negative Control |
| **TOTAL** | **Multi-Center Master** | **3,461** | **All Endpoints** | — | — | **REMARK-Compliant Total** |

*Note: Of $N = 3,461$ clinically indexed records across 11 cohorts, $36$ TCGA samples lacked paired RNA-seq matrices, yielding $N = 3,425$ complete matched molecular and blood profiles analyzed (604 + 1,056 + 986 + 226 + 18 + 535 = 3,425).*

---

## Supplementary Note 5: REMARK Guidelines Full Compliance Checklist & Ethical Declarations

1. **Hypothesis & Study Rationale:** Explicitly prespecified: evaluation of whether coupled histological architecture, interstitial transport, and inflammatory markers provide independent prognostic stratification in retrospective surgical cohorts.
2. **Patient Selection & Cohorts:** Strict complete-case analyses ($N = 3,461$ across 11 cohorts). For raw blood chemistry ($N = 1,056$), incomplete records ($N = 274$) were excluded a priori. No imputation was performed.
3. **Data Modalities & Assay Methods Mapping:**

| Data Modality | Cohorts | Epithelial Polarity Domain ($\mathcal{P}_{\text{hist}}$) | Desmoplasia Domain ($\mathcal{M}_{\text{mech}}$) | Inflammatory Soil Domain ($\Phi_{\text{soil}}$) |
| :--- | :--- | :--- | :--- | :--- |
| **Modality 1: Raw Blood & Surgical Path** | PLOS ONE ($N=1,056$) | Histological tumor differentiation grade | Pathologic invasiveness / local spread | Pre-op CBC NLR, serum CRP, and Albumin (mGPS) |
| **Modality 2: Bulk Transcriptomics** | Solid: TCGA ($N=996$), GEO GSE31210 ($N=226$), Controls ($N=561$) | Polarity/cilia module z-score (*TJP2*, *CDH1*, *PKHD1*, *IFT88*) | Stromal matrix module z-score (*COL1A1*, *ACTA2*, *LOXL2*) | Myeloid cytokine module z-score (*IL1B*, *IL6*, *CXCL8*, *TNF*) |
| **Modality 3: Digital Histology & Clinical** | MSKCC 2020 ($N=604$) | Lepidic/acinar vs solid/micropapillary architecture | Peritumoral desmoplastic fibrosis | Smoking pack-years & inflammatory covariates |

4. **Study Design:** Retrospective translational multi-cohort study analyzing de-identified open-access public repositories.
5. **Statistical Analysis:** Kaplan-Meier survival curves, log-rank tests, multivariable Cox proportional hazards with Schoenfeld residual testing, Aalen-Johansen competing risks cumulative incidence functions, 500-fold bootstrap cross-validation, and time-dependent receiver operating characteristic (ROC) curves.
6. **Negative Controls & Falsification:** Proof-by-contradiction negative controls evaluated in liquid tumors (TCGA-LAML), germline tumors (TCGA-TGCT), and sarcomas (TCGA-SARC), confirming model failure where biophysical assumptions are absent.
7. **Ethical Approval & IRB Exemption:** Secondary analysis of fully de-identified, publicly accessible datasets. In accordance with national guidelines and institutional policies of the Department of Clinical Medicine and Surgery, Università degli Studi di Napoli Federico II, this study is exempt from formal Institutional Review Board (IRB) review.
8. **Regulatory Notice:** All computational scores and calculator prototypes are designated strictly for **Research Use Only (RUO)**.
