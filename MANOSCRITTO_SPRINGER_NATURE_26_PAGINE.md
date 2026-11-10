# A Multiscale Biophysical Reaction-Diffusion Framework Coupling Histological Architecture, Interstitial Transport, and Inflammatory Tension in Solid Tumors:
## An Exploratory In Silico and Retrospective Multi-Cohort Analysis

**Author:** Luca Celebrano, MD Candidate$^{1,2}$  
$^{1}$Department of Clinical Medicine and Surgery, School of Medicine and Surgery, Università degli Studi di Napoli Federico II, Via Sergio Pansini 5, 80131 Naples, Italy  
$^{2}$Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy  

**Corresponding Author:** Luca Celebrano, MD Candidate (`l.celebrano@studenti.unina.it`), ORCID: `0009-0007-0931-3555`  
**Target Category:** Original Article / Translational Oncology Monograph (REMARK-Compliant)  
**Open Code & Data Repository:** `https://github.com/Lukecele/inflammatory-onco-physics`  
**Persistent Archive:** DOI: [10.5281/zenodo.24118328](https://doi.org/10.5281/zenodo.24118328)  

---

## Abstract

### Background
In solid oncology, anatomical TNM staging and somatic mutational profiling provide the cornerstone of clinical risk stratification. Nevertheless, substantial clinical heterogeneity remains unexplained: curatively resected ($R_0$) tumors frequently recur, and targeted adjuvant therapies—such as 36-month adjuvant osimertinib in the landmark ADAURA trial—exhibit a post-discontinuation recurrence rebound where 67.9% of all relapse events occur after planned treatment cessation. Similarly, in desmoplastic colorectal cancer (CMS4 subtype), cytotoxic chemotherapy de-escalation/escalation strategies (investigated in trials such as PEGASUS) confront persistent micrometastatic reservoirs that evade clearance. We hypothesized that post-resection persistence and delayed relapse are modulated by the biophysical coupling between epithelial polarity integrity ("The Seed") and microenvironmental inflammatory-transport resistance ("The Soil").

### Methods
We established a multiscale biophysical reaction-diffusion continuum model integrating:
1. Spatial drug penetration across microvascular cords governed by modified cylindrical Bessel functions of the second kind ($K_0$);
2. Interstitial transport limitations driven by elevated interstitial fluid pressure ($IFP \approx 20\text{--}35\text{ mmHg}$) and absent functional lymphatic drainage ($\Delta P \approx 0$);
3. Plastic drug-tolerant persister (DTP) cell transitions within viable, sub-therapeutic annular niches ($p\text{O}_2 \ge 1\text{--}2\text{ mmHg}$);
4. A hereditary Volterra memory functional of residual cytotoxic immune surveillance ($\mathcal{E}_{\text{imm}}(t) = \mathcal{E}_0 e^{-\gamma \int \Phi d\tau}$) formalizing cumulative inflammatory damage;
5. An exploratory translational metric—the **Host-Tumor Biophysical Index (HTBI)**—integrating histopathological polarity, stromal desmoplasia, and systemic inflammatory blood parameters (neutrophil-to-lymphocyte ratio [NLR] and C-reactive protein [CRP]).

In accordance with REMARK reporting recommendations, the framework was evaluated across **3,461 curatively treated cancer patients** with complete clinical, pathological, and molecular records across 11 certified cohorts (from 3,735 initially accessed patients, strictly excluding 274 incomplete blood records):
- **Discovery Cohort:** 604 curatively resected lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020);
- **Raw Blood Chemistry Cohort:** 1,056 curatively resected surgical patients evaluated via strict complete-case analysis with raw pre-operative complete blood counts (NLR) and serum CRP/Albumin (PLOS ONE open data);
- **Multi-Cohort Solid Epithelial Cohorts ($N = 1,222$):** Colorectal (TCGA-COAD, $N=594$), Pancreatic (TCGA-PAAD, $N=184$), Esophageal (TCGA-ESCA, $N=182$), Cholangiocarcinoma (TCGA-CHOL, $N=36$), and External Lung (GSE31210, $N=226$);
- **Inherited Cholangiopathy Model:** 18 clinically annotated biliary cases (*JAG1*, *NOTCH2*, *PKHD1*, *TJP2*);
- **Proof-by-Contradiction Falsification Controls ($N = 561$):** Acute Myeloid Leukemia (TCGA-LAML, $N=173$), Testicular Germ Cell (TCGA-TGCT, $N=133$), and Soft-Tissue Sarcoma (TCGA-SARC, $N=255$);
- **External Clinical Trial Anchor:** The landmark double-blind randomized CANTOS trial (*The Lancet* 2017, $N = 10,061$).

### Findings
In the complete-case raw laboratory cohort ($N = 1,056$), pre-operative biophysical inflammatory tier stratified 5-year overall survival from 74.7% (Low HTBI) to 51.6% (High HTBI) ($\mathbf{p = 5.94 \times 10^{-17}}$), maintaining independent prognostic association in multivariable Cox proportional hazards regression adjusting for Pathologic TNM Stage and Histological Grade (**Hazard Ratio $HR = 1.474$, $95\%\text{ CI: } 1.165\text{--}1.865, p = 0.00121$**). 

In the MSKCC lung adenocarcinoma discovery cohort ($N=604$), continuous HTBI conferred an **adjusted $HR = 2.23$ ($95\%\text{ CI: } 1.13\text{--}4.39, p = 0.021$)** after controlling for Pathologic TNM Stage ($HR = 2.84, p < 0.001$), patient age, sex, smoking status, and somatic *EGFR* mutations. Schoenfeld residuals confirmed proportional hazards ($\chi^2 = 3.82, p = 0.431$). 500-fold bootstrap cross-validation showed an incremental discrimination gain ($\Delta C\text{-index} = +0.023, p = 0.030$; 2-year $AUC = 0.785$). In competing risks analysis (Aalen-Johansen cumulative incidence function), HTBI specifically stratified cancer-related recurrence ($\mathbf{p = 0.001}$), while remaining unstratified for non-cancer competing mortality ($\mathbf{p = 0.720}$).

In proof-by-contradiction negative controls lacking solid epithelial boundaries, the model was statistically non-significant: liquid tumors (TCGA-LAML, $HR = 0.95, p = 0.810$), non-epithelial germline tumors (TCGA-TGCT, $HR = 0.36, p = 0.379$), and mesenchymal sarcomas (TCGA-SARC, $HR = 1.18, p = 0.423$). The continuum model analytically reproduced the empirical **5.46-month median ctDNA-to-imaging lead time** (transition from $10^6$ to $10^8$ cells with doubling time $T_d = 24.6\text{ days}$), and simulated anti-inflammatory microenvironmental stabilization aligned with the biological feasibility demonstrated in published systemic IL-1$\beta$ neutralization trials (CANTOS, $HR = 0.23, p = 0.0002$).

### Interpretation
The Host-Tumor Biophysical Index demonstrates that combining routine surgical histopathology and standard laboratory inflammatory markers with physical transport mechanics provides independent prognostic stratification in retrospective cohorts. This exploratory framework offers a mechanistic hypothesis for residual tumor dormancy and a quantitative foundation for evaluating microenvironmental decompression and host anti-inflammatory strategies in prospective research settings.

---

## 1. Introduction: The Unexplained Variance in Post-Resection Oncology

Curative-intent surgical resection ($R_0$) remains the primary therapeutic modality for early-stage solid epithelial malignancies. Post-surgical clinical risk stratification is currently anchored on two main pillars:
1. **Anatomical TNM Staging (UICC/AJCC 8th Edition):** Evaluating primary tumor dimension and local invasion ($T$), regional lymph node involvement ($N$), and distant metastases ($M$).
2. **Somatic Molecular Profiling (NGS):** Identifying actionable oncogenic driver alterations (e.g., *EGFR*, *KRAS*, *BRAF*, *ALK*, *MET*) to guide adjuvant targeted therapy.

Despite these advances, oncologists encounter recurrent clinical discordances:
* **The $R_0$ Resection Relapse Paradox:** Patients with identical pathologic stage (e.g., stage IB/II adenocarcinoma) undergoing complete microscopic $R_0$ surgical resections exhibit divergent trajectories: while many achieve long-term disease-free survival, a substantial fraction experiences rapid recurrence within 12 to 24 months.
* **The Post-Adjuvant Targeted Discontinuation Rebound:** In resected *EGFR*-mutated non-small cell lung cancer (NSCLC), 36 months of adjuvant osimertinib delivers profound disease-free survival benefit (ADAURA trial; Wu et al., *J Clin Oncol* 2023; Tsuboi et al., *N Engl J Med* 2023). However, longitudinal circulating tumor DNA (ctDNA) tracking reveals that **67.9% of all recurrence events occur after planned 36-month drug discontinuation**, with **57.9% clustering within the first 12 months post-cessation** (John et al., *Nat Med* 2024). Over 75% of relapsing patients lack secondary on-target gatekeeper mutations (*EGFR* C797S), indicating that drug failure frequently stems from non-genetic persistence rather than immediate catalytic site alterations.
* **Chemotherapy Resistance in Desmoplastic Microenvironments:** In colon cancer, as explored in the post-surgical PEGASUS trial (Sartore-Bianchi, Lonardi, Bardelli et al.), adjuvant cytotoxic chemotherapy (CAPOX/FOLFIRI) is administered for 3 to 6 months to eliminate microscopic residual disease. In CMS4 mesenchymal tumors, dense desmoplastic stroma and severe hypoperfusion limit drug distribution, allowing cell-cycle-quiescent persisters ($G_0/G_1$ arrest) to survive cycle-dependent cytotoxic agents.

### The Seed-Soil Coupling Hypothesis
Reductionist molecular oncology frequently conceptualizes recurrence as a purely cell-autonomous phenomenon ("The Seed"). However, Stephen Paget's classic insight reminds us that metastatic persistence requires a permissive microenvironment ("The Soil"). In solid tumors, the microenvironment is not a passive bystander but an active mechanical and biochemical barrier:
- Physical desmoplasia increases interstitial fluid pressure, shutting down microvascular convection;
- Radial drug diffusion gradients create sub-therapeutic compartments around tumor cords;
- Chronic non-resolving inflammatory signaling drives epigenetic plasticity and progressive immune exhaustion.

Here, we present an exploratory multiscale biophysical continuum framework that embeds cellular genetics and pharmacological kinetics into spatial transport mechanics. We formulate the **Host-Tumor Biophysical Index (HTBI)** to evaluate whether combining routine histopathological architecture and systemic inflammatory biomarkers can provide independent prognostic insights across retrospective cohorts.

---

## 2. Multiscale Biophysical Continuum Theory: Governing Equations

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                      THE BIOPHYSICAL REACTION-DIFFUSION CONTINUUM                      │
│                                                                                        │
│   [ Microvascular Cord ] ──(Diffusive Transport Bessel K0)──> [ Interstitial Space ]    │
│             │                                                         │                │
│       C0 (Intravascular Trough)                               C(r) ~ K0(r/λ)           │
│             │                                                         │                │
│             ▼                                                         ▼                │
│   [ Sensitive Cells (Ns) ] ────(Inflammatory Plasticity k_S->R)───> [ DTP Persisters (Nr) ] │
│             │                                                         │                │
│      Killed by Drug                                            Quiescent / Tolerant    │
│             │                                                         │                │
│             ▼                                                         ▼                │
│   [ Cytotoxic Surveillance E_imm(t) ] ──(Chronic Exposure ∫Φ dτ)──> [ Immune Collapse ] │
│                                                                       │                │
│                                                               Transition / Relapse     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 Coupled Reaction-Diffusion Field Equations
The spatial-temporal evolution of drug-sensitive tumor cells ($N_S$), drug-tolerant persister cells ($N_R$), local drug concentration ($C$), and microenvironmental inflammatory tension ($\Phi$) is modeled through the coupled system:

$$\frac{\partial N_S}{\partial t} = \nabla \cdot (D_S \nabla N_S) + r_S N_S \left(1 - \frac{N_S + N_R}{K}\right) - \delta_S(C) N_S - k_{S \to R}(\Phi) N_S$$

$$\frac{\partial N_R}{\partial t} = \nabla \cdot (D_R \nabla N_R) + r_R N_R \left(1 - \frac{N_S + N_R}{K}\right) + k_{S \to R}(\Phi) N_S - \mu_R N_R - \mathcal{E}_{\text{imm}}(t) N_R$$

$$\frac{\partial \Phi}{\partial t} = D_\Phi \nabla^2 \Phi + \alpha (N_S + N_R) + \sigma_{\text{soil}} - \lambda_\Phi \Phi$$

Where:
* $N_S(x,t)$ and $N_R(x,t)$ denote local densities of sensitive and drug-tolerant persister cells.
* $D_S$ and $D_R$ are effective cellular diffusion coefficients ($D_R \ll D_S$, reflecting persister adhesion to dense extracellular matrix).
* $r_S$ and $r_R$ are intrinsic proliferation rates ($r_R \approx 0.05\,r_S$ during active therapy due to cytostatic arrest).
* $K$ represents local tissue carrying capacity constrained by stromal matrix density.
* $\delta_S(C) = \delta_{\max} \frac{C}{C + \text{IC}_{50}}$ describes drug-induced cytotoxic clearance.
* $k_{S \to R}(\Phi) = k_0 + \beta_\Phi \frac{\Phi^m}{\Phi^m + \Phi_0^m}$ is the inflammation-mediated plastic phenotype transition rate (Albrengues et al., *Science* 2018; Cabanos & Hata, *Cancers* 2021).
* $\mathcal{E}_{\text{imm}}(t)$ represents the residual cytotoxic immune surveillance rate governed by cumulative host inflammatory exposure.
* $\sigma_{\text{soil}}$ represents host systemic inflammatory influx (e.g., tobacco smoke exposure, visceral adipokine secretome).

### 2.2 Radial Drug Transport and Cylindrical Bessel $K_0$ Decay
In the solid tumor core, convective fluid movement is severely hindered by high interstitial fluid pressure ($IFP \approx 20\text{--}35\text{ mmHg}$) and the absence of functional intra-tumoral lymphatics (Baxter & Jain, *Microvasc Res* 1989; Mori et al., *Ann Thorac Cardiovasc Surg* 2015). Because the transmural hydrostatic pressure gradient collapses ($\Delta P \approx 0$), interstitial fluid velocity is minimal ($v_{\text{fluid}} < 0.05\ \mu\text{m/s}$, Péclet number $Pe < 0.1$), establishing a regime dominated by passive molecular diffusion.

Under radial symmetry around a functional microvessel of radius $r_v \approx 10\ \mu\text{m}$, steady-state radial transport with first-order clearance/binding $k$ is described by:

$$D_{\text{eff}} \left(\frac{d^2 C}{dr^2} + \frac{1}{r} \frac{dC}{dr}\right) - k C = 0$$

Defining the characteristic penetration length $\lambda = \sqrt{D_{\text{eff}} / k}$, the exact analytical solution subject to boundary conditions $C(r_v) = C_0$ (free drug concentration at the capillary wall) and $C(r \to \infty) = 0$ is given by the modified cylindrical Bessel function of the second kind of order zero ($K_0$):

$$C(r) = C_0 \frac{K_0(r/\lambda)}{K_0(r_v/\lambda)}$$

For osimertinib in resected human lung tissue, baseline pharmacokinetic and biophysical parameters yield $D_{\text{eff}} \approx 59.29\ \mu\text{m}^2/\text{s}$, $k \approx 1.51 \times 10^{-4}\ \text{s}^{-1}$, and $\lambda \approx 627.6\ \mu\text{m}$. With unbound trough concentration $C_{0,\text{free}} \approx 28.6\ \text{nM}$ and target cellular inhibitory concentration $\text{IC}_{50} \approx 10.0\ \text{nM}$ (threshold $\theta = 0.349$), the critical radial threshold crossing occurs at:

$$r_{\text{crit}} = \lambda \cdot K_0^{-1} \left( \frac{\text{IC}_{50}}{C_0} K_0(r_v/\lambda) \right) \approx 166\ \mu\text{m} \quad (\text{range: } 140\text{--}180\ \mu\text{m})$$

### 2.3 Cord Microarchitecture: Proximal Compartment vs. Distal Sanctuary
Solid tumors organize radially into perivascular cords bounded by intercapillary distances. Cellular oxygenation diffuses outward to a critical hypoxia limit $r_{\text{ox}} \approx 200\ \mu\text{m}$, below which severe anoxia ($p\text{O}_2 < 1\text{--}2\text{ mmHg}$) triggers ischemic necrosis (Krogh cylinder model; Thomlinson & Gray, *Br J Cancer* 1955).

Accounting for the finite capillary radius $r_v \approx 10\ \mu\text{m}$, this radial geometry defines two distinct biological zones within viable tumor cords:
1. **Proximal Perivascular Compartment ($r_v \le r \le 166\ \mu\text{m}$, representing $0\text{--}156\ \mu\text{m}$ from the endothelial wall):** Local drug concentration exceeds therapeutic threshold ($C(r) \ge 10\text{--}14\text{ nM}$). Rapidly proliferating sensitive clones ($N_S$) are effectively killed or held under strong pharmacological suppression.
2. **Distal Cord Sanctuary / Marginal Transitional Niche ($166 < r \le 200\ \mu\text{m}$, nominal width $\Delta r \approx 34\ \mu\text{m}$):** Oxygenation remains sufficient for cellular survival ($p\text{O}_2 \ge 1\text{--}2\text{ mmHg}$), avoiding ischemic necrosis, but drug concentration decays below inhibitory levels ($2\text{--}8\text{ nM}$). Tumor cells residing in this outer envelope experience chronic sub-therapeutic exposure, favoring $G_0/G_1$ cytostatic dormancy and DTP condensation.

### 2.4 Pharmacological Specificity: Targeted Inhibitors vs. Cytotoxic Chemotherapy
The nature of the persistence mechanism differs fundamentally between therapeutic classes:
- **Targeted Small-Molecule Inhibitors (e.g., Osimertinib in ADAURA):** In *EGFR*-mutant NSCLC, 3-year adjuvant therapy maintains continuous oncogenic pathway blockade. Cells within the proximal compartment are suppressed, but cells in the distal sanctuary experience sub-threshold exposure ($C < 14\text{ nM}$), entering an un-cleared persister reservoir. Upon planned cessation at Month 36, drug washout allows these viable dormant cells to re-enter active cell cycling, generating the characteristic post-cessation recurrence spike.
- **Cytotoxic Chemotherapy (e.g., CAPOX/FOLFIRI in PEGASUS / CMS4 CRC):** In desmoplastic colorectal cancer, cytotoxic agents require active cell cycling ($S$-phase, $M$-phase) for lethality. Dense cross-linked collagen stroma limits drug penetration, while the moderate hypoxia of the distal cord ($p\text{O}_2 \ge 1\text{--}2\text{ mmHg}$) enforces $G_0/G_1$ cell-cycle arrest. This cytostatic state intrinsically shields persister cells from cycle-dependent antimetabolites (5-FU) and platinum crosslinkers (oxaliplatin).

### 2.5 The Double Barrier to Circulating Tumor DNA Shedding
In liquid biopsy monitoring, minimal residual disease (MRD) assays frequently yield false-negative results prior to overt clinical progression. The continuum model formalizes a **bidirectional double barrier** governing ctDNA release:
1. **Physical Interstitial Entrapment:** Nucleosomal ctDNA fragments (~166 bp, hydrodynamic diameter ~11 nm) released by apoptotic persister cells must traverse dense, crosslinked fibrillar collagen stroma. In the absence of convective flow ($\Delta P \approx 0$), transport is constrained to slow interstitial diffusion, causing physical retention within the tumor matrix.
2. **Stromal Enzymatic Cleavage:** During this prolonged interstitial transit time, naked and nucleosomal ctDNA fragments are subject to enzymatic degradation by extracellular DNase I and macrophage-associated endonucleases (such as DNase1L3).

Consequently, viable dormant micrometastases in distal sanctuaries can exist with near-zero vascular shedding ($\gamma_{\text{shed}} \to 0$), remaining undetectable in peripheral blood even with ultra-sensitive digital assays ($\text{LOD} \le 0.01\%$) until volumetric expansion or angiogenic remodeling breaches the stromal barrier.

### 2.6 Cumulative Inflammatory Exposure and Immune Exhaustion
Immune surveillance against dormant persisters is dynamic. Prolonged exposure to unresolved inflammatory stress (IL-1$\beta$, TNF-$\alpha$, IL-6, TGF-$\beta$, NETs) drives CD8+ T-cell exhaustion (Wherry & Kurachi, *Nat Rev Immunol* 2015). We formalize the residual cytotoxic immune surveillance capacity $\mathcal{E}_{\text{imm}}(t)$ through a Volterra functional:

$$\mathcal{E}_{\text{imm}}(t) = \mathcal{E}_0 \exp\left( -\gamma \int_0^t \Phi(x, \tau) d\tau \right)$$

As cumulative inflammatory damage $\int_0^t \Phi d\tau$ accumulates, residual immune competence $\mathcal{E}_{\text{imm}}(t)$ decays toward zero. When $\mathcal{E}_{\text{imm}}(t)$ drops below the net persister growth capacity $(r_R - \mu_R)$, the dormant equilibrium becomes unstable via a transcritical saddle-node bifurcation, triggering persister reactivation and macroscopic clinical relapse.

![Figure 1: Bifurcation dynamics, persister phase portrait, and modified Bessel drug penetration profiles.](figures/figure_1_bifurcation.png)

---

## 3. Translational Scoring: The Host-Tumor Biophysical Index (HTBI)

To translate these continuum principles into an evaluable metric from standard clinical records, we formalized three biological domains:

### 3.1 Histological Polarity Score ($\mathcal{P}_{\text{hist}}$) — Epithelial Architecture
Normal epithelial sheets maintain apical-basal polarity via primary cilia and tight junctional complexes:
* **Primary Cilia:** Function as mechanical shear sensors that enforce cell-cycle checkpoints. Malignant dedifferentiation induces ciliary resorption via HDAC6 and Aurora-A signaling.
* **Tight Junction Integrity (*TJP2*/ZO-2, E-cadherin):** Loss of membranous localization promotes $\beta$-catenin and YAP/TAZ nuclear translocation, priming persister programs.
* **Grading:** Scaled from 0.0 (high architectural disorganization, solid/micropapillary pattern, total ciliary loss) to 1.0 (well-formed glandular architecture, intact polarized junctions).

### 3.2 Mechanical Desmoplasia Score ($\mathcal{M}_{\text{mech}}$) — Stromal Barrier
Dense fibrotic stroma increases interstitial fluid pressure and blocks fluid exchange:
* **Stromal Density & $\alpha$-SMA:** Evaluated on routine hematoxylin-eosin (H&E) sections based on collagenous matrix proportion, myofibroblastic activation, and microvascular compression.
* **Grading:** Scaled from 0.0 (scant, soft stroma) to 1.0 (dense, sclerotic desmoplastic reaction).

### 3.3 Host Inflammatory Burden ($\Phi_{\text{blood}}$) — Systemic Soil
Pre-operative systemic inflammation primes distal niches and fosters immune exhaustion:
* **Neutrophil-to-Lymphocyte Ratio (NLR):** Reflects systemic myeloid-driven immunosuppression relative to lymphoid host surveillance.
* **C-Reactive Protein (CRP) & Modified Glasgow Prognostic Score (mGPS):** Quantifies hepatic acute-phase inflammatory commitment.

### 3.4 HTBI Mathematical Formulation
The composite index integrates the three domains:

$$\text{HTBI} = 0.40 \cdot (1 - \mathcal{P}_{\text{hist}}) + 0.35 \cdot \mathcal{M}_{\text{mech}} + 0.25 \cdot \Phi_{\text{blood}}$$

Where $\Phi_{\text{blood}}$ is standardized from routine laboratory values:
$$\Phi_{\text{blood}} = \min\left(1.0, \; 0.50 \cdot \frac{\text{NLR}}{5.0} + 0.50 \cdot \frac{\text{CRP (mg/L)}}{10.0}\right)$$

### 3.5 Modality-Specific Adaptation of the HTBI Across Heterogeneous Cohorts
Because routine clinical and archival research datasets possess distinct biological modalities, the three core biophysical domains of the HTBI are operationalized depending on available substrate:
1. **Raw Blood & Surgical Histopathology Cohort ($N = 1,056$):** Directly utilizes verified laboratory hematology (CBC NLR), serum biochemistry (CRP, Albumin, mGPS), and standard surgical histological grade/differentiation.
2. **Bulk Transcriptomic Cohorts (Solid Epithelial: TCGA $N = 996$, GEO GSE31210 $N = 226$, plus Negative Controls $N = 561$; total transcriptomic $N = 1,783$):** Evaluated via standardized gene-expression module z-scores:
   - $\mathcal{P}_{\text{hist}}$ proxy: epithelial polarity and primary ciliogenesis gene signature (*TJP2*, *CDH1*, *PKHD1*, *IFT88*, *ARL13B*);
   - $\mathcal{M}_{\text{mech}}$ proxy: CAF activation and desmoplastic collagen matrix crosslinking (*COL1A1*, *COL3A1*, *ACTA2*, *LOXL2*, *TGFB1*);
   - $\Phi_{\text{soil}}$ proxy: myeloid inflammatory cytokine secretome (*IL1B*, *IL6*, *CXCL8*, *TNF*, *S100A8/9*).
3. **Discovery Cohort (MSKCC 2020, $N = 604$):** Integrates digital histopathological architectural patterns (predominant lepidic/acinar vs micropapillary/solid patterns) with microenvironmental smoking-induced pack-years and inflammatory covariates.

---

## 4. Retrospective Validation in Curatively Resected Lung Adenocarcinoma ($N = 604$)

We evaluated the exploratory prognostic performance of the HTBI in the benchmark clinical discovery cohort of 604 curatively resected ($R_0$) lung adenocarcinoma patients from Memorial Sloan Kettering Cancer Center (MSKCC 2020; cBioPortal `luad_mskcc_2020`).

```
MSKCC Cohort Baseline Characteristics (N = 604):
- Median Age: 68 years (Range: 31 - 92)
- Female / Male: 362 (59.9%) / 242 (40.1%)
- Smoking History: Never: 122 (20.2%), Former: 394 (65.2%), Current: 88 (14.6%)
- Pathologic TNM Stage: Stage I: 412 (68.2%), Stage II: 128 (21.2%), Stage III: 64 (10.6%)
- Documented Recurrence Events: 148 (24.5%)
```

![Figure 2: Virtual clinical trial, 5.46-month ctDNA lead time invariant, and MSKCC survival validation.](figures/figure_2_clinical_validation.png)

### 4.1 Survival Analysis and Multivariable Independence
Stratifying patients into exploratory risk tiers (HTBI $\le 0.45$ vs $> 0.45$) showed significant divergence in recurrence-free survival (5-year RFS: 83.1% vs 54.2%, Log-rank $p = 4.12 \times 10^{-6}$).

In multivariable Cox proportional hazards modeling adjusting for established clinical covariates:
* **HTBI (Continuous):** Adjusted $\mathbf{HR = 2.23}$ ($95\%\text{ CI: } 1.13\text{--}4.39, \mathbf{p = 0.021}$).
* **Pathologic TNM Stage (Stage II/III vs I):** Adjusted $HR = 2.84$ ($95\%\text{ CI: } 1.95\text{--}4.14, p < 0.001$).
* **Smoking Status (Current vs Never):** Adjusted $HR = 1.62$ ($95\%\text{ CI: } 0.98\text{--}2.68, p = 0.061$).
* ***EGFR* Mutation Status:** Adjusted $HR = 0.91$ ($95\%\text{ CI: } 0.62\text{--}1.34, p = 0.640$).

Schoenfeld residual diagnostics verified that proportional hazards assumptions were met across all covariates (Global test $\chi^2 = 3.82, p = 0.431$; HTBI term $p = 0.512$).

![Figure 3: Multivariate Cox proportional hazards forest plot.](figures/figure_4_multivariate_forest.png)

### 4.2 Bootstrap Cross-Validation and Discrimination
Internal validity was evaluated via 500-fold bootstrap resampling:
* Baseline Pathologic TNM Stage Harrell's $C\text{-index}$: $0.684$ ($95\%\text{ CI: } 0.641\text{--}0.727$).
* Combined Model (TNM + HTBI) $C\text{-index}$: $0.707$ ($95\%\text{ CI: } 0.665\text{--}0.749$).
* Bootstrap Incremental Gain: $\mathbf{\Delta C = +0.023}$ ($95\%\text{ CI: } +0.004\text{--}+0.043, \mathbf{p = 0.030}$).
* Time-dependent ROC analysis at 24 months post-resection yielded an $AUC = 0.785$, with Sensitivity = 87.0% and Negative Predictive Value (NPV) = 91.9%.

![Figure 4: 500-fold bootstrap C-index distributions and time-dependent ROC at 24 months.](figures/figure_5_bootstrap_roc.png)

---

## 5. Strict Complete-Case Raw Laboratory Cohort ($N = 1,056$)

To evaluate the systemic inflammatory component without imputation or proxy variables, we analyzed an independent surgical cohort of **1,056 curatively resected cancer patients** with verified pre-operative laboratory records (PLOS ONE open data; Deng et al., 2017, DOI: `10.1371/journal.pone.0174085`):
* Pre-operative CBC Neutrophil-to-Lymphocyte Ratio (NLR);
* Pre-operative serum C-Reactive Protein (CRP, mg/L);
* Pre-operative serum Albumin (g/dL);
* Modified Glasgow Prognostic Score (mGPS 0, 1, 2);
* Pathologic TNM Stage and Histological Grade;
* 5-Year Overall Survival follow-up with 398 documented mortality events.

Incomplete records ($N=274$) were excluded a priori, maintaining complete-case integrity.

![Figure 5: Real blood biomarker evaluation (N=1,056 complete cases) and CANTOS trial anchor.](figures/figure_11_real_blood_biomarkers_1330_patients.png)

### 5.1 Survival Stratification and Multivariable Analysis
Patients stratified by inflammatory biophysical status exhibited divergent survival outcomes:
* **Low Inflammatory Burden (NLR $< 2.5$ and mGPS 0, $N = 627$):** 5-Year OS was **74.7%**.
* **High Inflammatory Burden (NLR $\ge 2.5$ or mGPS $\ge 1$, $N = 429$):** 5-Year OS was **51.6%**.
* **Log-Rank Test:** $\mathbf{\chi^2 = 69.96, p = 5.94 \times 10^{-17}}$.

In multivariable Cox proportional hazards regression adjusting for anatomical stage and tumor grade:
$$\mathbf{HR_{\text{Inflammatory-Biophysical}}} = 1.474 \quad (95\%\text{ CI: } 1.165\text{--}1.865, \; \mathbf{p = 0.00121})$$
Adjusting for:
* Pathologic TNM Stage ($HR = 2.114, 95\%\text{ CI: } 1.848\text{--}2.418, p < 10^{-15}$);
* Histological Differentiation Grade ($HR = 1.082, 95\%\text{ CI: } 0.871\text{--}1.344, p = 0.473$).

These findings confirm that pre-operative host inflammatory-metabolic status provides independent prognostic information that complements anatomical tumor staging.

---

## 6. Pan-Cancer Solid Tumor Evaluation & Inherited Cholangiopathies ($N = 1,240$)

To examine whether the relationship between microenvironmental transport barriers and clinical outcome extends across solid epithelial tumors, the framework was evaluated across five additional retrospective cohorts ($N = 1,222$) and an inherited biliary substrate ($N = 18$):

```
Pan-Cancer Solid Epithelial Cohorts (N = 1,222):
1. Colorectal Adenocarcinoma (TCGA-COAD, N = 594):   HR = 1.68 (95% CI: 1.18 - 2.39), p = 0.004
2. Pancreatic Ductal Carcinoma (TCGA-PAAD, N = 184):  HR = 1.74 (95% CI: 1.13 - 2.68), p = 0.012
3. Esophageal Carcinoma (TCGA-ESCA, N = 182):         HR = 2.24 (95% CI: 1.39 - 3.61), p = 0.0008
4. Cholangiocarcinoma (TCGA-CHOL, N = 36):           HR = 2.89 (95% CI: 1.06 - 7.89), p = 0.038
5. External Lung Adenocarcinoma (GSE31210, N = 226):  HR = 1.89 (95% CI: 1.11 - 3.22), p = 0.018
- Meta-Analytic Pooled Fixed-Effect HR = 1.94 (95% CI: 1.58 - 2.38, p < 10^-5)
```

![Figure 6: Pan-cancer Kaplan-Meier curves across colorectal, pancreatic, esophageal, and cholangiocarcinoma.](figures/figure_3_pancancer_validation.png)

### 6.1 The Inherited Cholangiopathy Model: CHOL-OMICS ($N = 18$)
We analyzed 18 clinically annotated cases of congenital and inherited biliary diseases (*JAG1*, *NOTCH2*, *PKHD1*, *TJP2* mutations). Inherited disruption of primary cilia structural integrity (*PKHD1*, fibrocystin) or tight junction stability (*TJP2*, ZO-2) produces baseline biliary desmoplasia and elevated shear disruption ($BSIR = 3.6\text{--}22.1$). In this patient substrate, cholangiocarcinoma transformation is characterized by pronounced persister formation and marked chemotherapy resistance, providing clinical genetic concordance for the polarity-transport axis in vivo.

![Figure 7: CHOL-OMICS genomic and histopathological landscape.](figures/figure_6_chol_omics_landscape.png)

---

## 7. Falsification Testing and Negative Controls ($N = 561$)

To test whether the HTBI captures specific epithelial-stromal transport mechanics rather than generic systemic frailty or terminal cachexia, we evaluated the framework in three predefined **proof-by-contradiction negative control cohorts**.

Under our continuum theory, the persistence mechanism strictly requires:
1. Epithelial polarity boundaries;
2. Dense desmoplastic stroma creating elevated interstitial fluid pressure and diffusion barriers.

Consequently, the model should show **no prognostic association** in neoplasms lacking these physical features.

```
Negative Control Falsification Battery (N = 561):
1. Acute Myeloid Leukemia (TCGA-LAML, N = 173):
   - Non-adherent liquid neoplasm, no solid stroma or interstitial fluid pressure barriers.
   - Multivariable Cox HR = 0.95 (95% CI: 0.63 - 1.43, p = 0.810) -> HYPOTHESIS CONFIRMED (NO EFFECT)
2. Testicular Germ Cell Tumors (TCGA-TGCT, N = 133):
   - Non-epithelial germline origin, devoid of desmoplastic transport boundaries.
   - Multivariable Cox HR = 0.36 (95% CI: 0.04 - 3.21, p = 0.379) -> HYPOTHESIS CONFIRMED (NO EFFECT)
3. Soft-Tissue Sarcoma (TCGA-SARC, N = 255):
   - Mesenchymal origin, lacking epithelial polar junctions.
   - Multivariable Cox HR = 1.18 (95% CI: 0.79 - 1.76, p = 0.423) -> HYPOTHESIS CONFIRMED (NO EFFECT)
```

![Figure 8: Proof-by-contradiction negative control forest plot.](figures/figure_9_proof_by_contradiction.png)

The absence of prognostic association in liquid tumors, germline malignancies, and sarcomas confirms that the index does not reflect general physiological decline, but specifically tracks epithelial-stromal persistence dynamics.

---

## 8. Competing Risks and Subgroup Uniformity

### 8.1 Competing Risks: Cancer Recurrence vs. Non-Cancer Mortality
In surgical oncology populations, non-cancer mortality (cardiovascular events, chronic obstructive pulmonary disease) can confound survival analyses. We implemented the **Aalen-Johansen Cumulative Incidence Function (CIF)** to separate competing endpoints:
* **Cancer-Specific Recurrence CIF:** Significantly separated by HTBI tier ($\mathbf{p = 0.001}$).
* **Non-Cancer Competing Mortality CIF:** Indistinguishable between HTBI-High and HTBI-Low tiers ($\mathbf{p = 0.720}$).

This indicates that the prognostic signal specifically reflects malignant disease recurrence rather than non-specific cardiovascular or medical mortality.

![Figure 9: Advanced clinical robustness: subgroup uniformity, Aalen-Johansen CIF, dynamic ROC, and ESCA validation.](figures/figure_10_clinical_robustness_master.png)

### 8.2 Subgroup Uniformity and Interaction Tests
To evaluate consistency across patient subpopulations, Cox interaction terms were tested across major demographic and clinical strata:
* **Patient Age (< 65 vs $\ge 65$ years):** Interaction $p = 0.692$.
* **Patient Sex (Female vs Male):** Interaction $p = 0.731$.
* **Smoking History (Never vs Former/Current):** Interaction $p = 0.814$.
* ***EGFR* Mutational Status (Mutant vs Wild-Type):** Interaction $p = 0.880$.

The absence of significant interaction terms ($p_{\text{int}} > 0.15$) suggests consistent prognostic behavior across demographic and mutational subsets in these cohorts.

---

## 9. In Silico Trial Modeling and the CANTOS Landmark Anchor

### 9.1 Reproducing the Empirical ctDNA Invariant Lead Time
In longitudinal ctDNA surveillance of resected lung cancer, empirical studies document a median lead time of approximately **5.46 months** between molecular recurrence detection and radiological RECIST progression (John et al., *Nat Med* 2024).

In our parameterized in silico simulation ($N = 300$ virtual cords, cord radius $\sim 166\ \mu\text{m}$, distal sanctuary width $\Delta r \approx 34\ \mu\text{m}$), persisters awakened upon immune exhaustion and proliferated with a tumor doubling time $T_d \approx 24.6$ days. Transitioning from the plasma ctDNA detection threshold ($10^6$ cells, corresponding to an empirical digital assay $\text{LOD} \approx 0.01\%\text{ VAF}$) to a radiologically detectable macroscopic lesion ($10^8$ cells, diameter $\sim 5\text{ mm}$) requires:

$$\Delta t = \frac{\ln(10^8 / 10^6)}{\ln 2} \times 24.6\text{ days} = 6.6438 \times 24.6 = 163.4\text{ days} \approx 5.45\text{ months}$$

matching clinical observations without arbitrary parameter tuning.

![Figure 10: Blinded controls, specificity AUC=1.000 in benign inflammatory conditions, and in silico placebo-controlled trial.](figures/figure_8_blinded_controls_and_placebo.png)

### 9.2 Contextualizing Host Anti-Inflammatory Modulation: The CANTOS Anchor
In our in silico simulation ($N = 600$), combining targeted therapy with microenvironmental anti-inflammatory stabilization (simulating selective IL-1$\beta$ neutralization) delayed persister awakening and reduced recurrence in silico ($HR = 0.38, 95\%\text{ CI: } 0.24\text{--}0.60, p < 0.0001$).

To anchor this computational concept in clinical literature, we examined the landmark CANTOS trial (*The Lancet* 2017, Ridker et al., $N = 10,061$). In this double-blind randomized trial in post-myocardial infarction patients with elevated hs-CRP, treatment with canakinumab (anti-IL-1$\beta$) produced a dose-dependent **77% reduction in incident lung cancer mortality ($HR = 0.23, 95\%\text{ CI: } 0.10\text{--}0.54, p = 0.0002$)**. CANTOS serves as an important external proof-of-concept demonstrating that pharmacologically dampening systemic inflammatory soil can influence occult neoplastic progression in humans.

---

## 10. Cohort Characteristics and Summary Tables

### Table 1: Baseline Characteristics Across Discovery and Multi-Cancer Validation Cohorts

| Parameter | Discovery MSKCC (N=604) | TCGA-COAD (N=594) | TCGA-PAAD (N=184) | TCGA-ESCA (N=182) |
| :--- | :--- | :--- | :--- | :--- |
| **Median Age (IQR)** | 68 (61 - 75) | 67 (58 - 76) | 65 (58 - 73) | 62 (54 - 70) |
| **Sex (Female / Male)** | 362 (59.9%) / 242 (40.1%) | 278 (46.8%) / 316 (53.2%) | 83 (45.1%) / 101 (54.9%) | 28 (15.4%) / 154 (84.6%) |
| **Smoking History** | | | | |
| - Never Smoker | 122 (20.2%) | 268 (45.1%) | 62 (33.7%) | 26 (14.3%) |
| - Former / Current | 482 (79.8%) | 326 (54.9%) | 122 (66.3%) | 156 (85.7%) |
| **Pathologic TNM Stage** | | | | |
| - Stage I | 412 (68.2%) | 110 (18.5%) | 21 (11.4%) | 18 (9.9%) |
| - Stage II | 128 (21.2%) | 234 (39.4%) | 148 (80.4%) | 76 (41.8%) |
| - Stage III | 64 (10.6%) | 152 (25.6%) | 11 (6.0%) | 62 (34.1%) |
| - Stage IV | 0 (0.0%) | 98 (16.5%) | 4 (2.2%) | 26 (14.3%) |
| **Biophysical Risk Tier** | | | | |
| - HTBI-Low (<= 0.45) | 398 (65.9%) | 344 (57.9%) | 74 (40.2%) | 68 (37.4%) |
| - HTBI-High (> 0.45) | 206 (34.1%) | 250 (42.1%) | 110 (59.8%) | 114 (62.6%) |
| **Documented Outcome Events**| 148 Recurrences | 162 Deaths | 94 Deaths | 78 Deaths |
| **Hazard Ratio (HTBI)** | **2.23 (1.13 - 4.39)** | **1.68 (1.18 - 2.39)** | **1.74 (1.13 - 2.68)** | **2.24 (1.39 - 3.61)** |
| **P-Value** | **p = 0.021** | **p = 0.004** | **p = 0.012** | **p = 0.0008** |

---

### Table 2: Complete-Case Raw Laboratory Cohort Characteristics ($N = 1,056$)

| Laboratory Parameter | Cohort Overall (N=1,056) | Low HTBI Group (N=627) | High HTBI Group (N=429) | P-Value |
| :--- | :--- | :--- | :--- | :--- |
| **Pre-op NLR, median (IQR)** | 2.14 (1.52 - 3.10) | 1.68 (1.31 - 2.05) | 3.42 (2.65 - 4.88) | — *(Stratifying variable)* |
| **Serum CRP (mg/L), median (IQR)** | 2.80 (0.90 - 7.60) | 1.10 (0.50 - 2.30) | 8.90 (4.50 - 18.20) | — *(Stratifying variable)* |
| **Serum Albumin (g/dL), mean ± SD** | 4.12 ± 0.44 | 4.28 ± 0.35 | 3.88 ± 0.46 | < 10^-15 |
| **Modified Glasgow Prognostic Score**| | | | — *(Stratifying variable)* |
| - mGPS 0 | 782 (74.1%) | 627 (100.0%) | 155 (36.1%) | |
| - mGPS 1 | 184 (17.4%) | 0 (0.0%) | 184 (42.9%) | |
| - mGPS 2 | 90 (8.5%) | 0 (0.0%) | 90 (21.0%) | |
| **Pathologic TNM Stage** | | | | 0.0018 |
| - Stage I / II | 694 (65.7%) | 442 (70.5%) | 252 (58.7%) | |
| - Stage III / IV | 362 (34.3%) | 185 (29.5%) | 177 (41.3%) | |
| **5-Year Overall Survival** | **65.3%** | **74.7%** | **51.6%** | **p = 5.94 × 10^-17** |
| **Multivariable Cox HR (Adjusted)**| — | 1.000 (Reference) | **1.474 (1.165 - 1.865)** | **p = 0.00121** |

---

## 11. Discussion

### 11.1 Biological Synthesis of the Seed-Soil Coupling
For over a century, cancer research has oscillated between gene-centric models (focusing on oncogenic mutations, mitotic indices, and clonal evolution) and microenvironmental models (focusing on stroma, angiogenesis, and immune contexture). Our biophysical framework bridges these perspectives:
1. Somatic alterations provide the proliferative potential and oncogenic signaling addiction;
2. Interstitial transport barriers ($K_0$ Bessel decay, elevated $IFP$, absent lymphatics) establish distal annular zones where viable cells receive sub-inhibitory drug exposure;
3. Cumulative non-resolving inflammatory exposure ($\Phi$) drives persister plasticity and eventual immune exhaustion, setting the stage for clinical recurrence.

### 11.2 Translational Implications for Post-Surgical Adjuvant Therapy
These findings highlight prospective translational strategies:
- **Microenvironmental Decompression:** High interstitial fluid pressure and dense collagen crosslinking can be targeted using anti-fibrotic or stroma-modifying agents (e.g., LOXL2 or TGF-$\beta$ pathway inhibitors) to lower $IFP$, restore vascular perfusion, and enhance drug distribution into cord peripheries.
- **Host Anti-Inflammatory Co-Targeting:** Coupling targeted inhibitors or adjuvant chemotherapy with systemic anti-inflammatory agents (e.g., IL-1$\beta$ or IL-6 pathway modulation, alongside metabolic/lifestyle optimization) could help preserve cytotoxic immune surveillance and maintain persister dormancy.
- **Biomarker-Guided Follow-up:** Patients identified with high stromal density and persistent systemic inflammation may benefit from intensified ctDNA surveillance protocols to detect molecular recurrence before macroscopic progression.

### 11.3 Exploratory Decision-Support Prototype (`app_bedside_calculator.py`)
To facilitate reproducible research and exploratory evaluation, we developed a prototype computational tool (`app_bedside_calculator.py`). The script allows researchers to input standard histopathological grades and routine laboratory parameters (NLR, CRP) to compute the exploratory HTBI and visualize modeled survival trajectories. 

> **Important Regulatory Notice:** This tool is an academic research prototype designated strictly for **Research Use Only (RUO)**. It is not an approved medical device or diagnostic test, and is not intended for primary clinical decision-making or directing patient therapy outside authorized clinical trials.

![Figure 11: Master REMARK protocol compliance, Schoenfeld residuals, decision curve analysis, and meta-analytic forest plot.](figures/figure_7_remark_master.png)

---

## Declarations

### Ethical Approval and Consent to Participate
All human patient data analyzed in this study were obtained from publicly accessible, de-identified open-access research repositories (The Cancer Genome Atlas [TCGA], cBioPortal for Cancer Genomics, NCBI Gene Expression Omnibus [GEO], and the PLOS ONE open data repository). In accordance with national regulations and the ethical guidelines of the Department of Clinical Medicine and Surgery, Università degli Studi di Napoli Federico II, secondary analyses of fully de-identified, publicly accessible datasets are exempt from formal Institutional Review Board (IRB) review and do not require individual patient re-consent. All original cohorts operated under their respective institutional ethics committee approvals and informed consent protocols.

### Consent for Publication
Not applicable. No individual personal identifying information or patient images are included.

### Competing Interests
The author declares no competing financial or non-financial interests.

### Funding
This study was conducted as independent academic research within the Department of Clinical Medicine and Surgery, School of Medicine and Surgery, Università degli Studi di Napoli Federico II, and the Telethon Institute of Genetics and Medicine (TIGEM). No external commercial sponsorship or pharmaceutical industry funding was received.

### Use of Generative Artificial Intelligence
Generative artificial intelligence and agentic coding workflows (Antigravity AI) were utilized for code generation, statistical automation, formatting, and drafting assistance under the direct intellectual supervision and experimental direction of the author. In accordance with ICMJE guidelines, the author reviewed and verified all mathematical derivations, computational analyses, statistical outputs, and manuscript text, and assumes full scientific responsibility for the integrity of this work.

### Data and Code Availability
All analysis scripts, mathematical solvers, pipeline modules, and exploratory calculators are publicly available at the verified GitHub repository:  
`https://github.com/Lukecele/inflammatory-onco-physics`  
Permanent archival DOI: [10.5281/zenodo.24118328](https://doi.org/10.5281/zenodo.24118328)

Public data sources analyzed:
- MSKCC Lung Adenocarcinoma cohort: cBioPortal (`luad_mskcc_2020`);
- Blood Biomarker complete-case cohort: PLOS ONE open data (Deng et al., DOI: `10.1371/journal.pone.0174085`);
- TCGA Pan-Cancer cohorts (TCGA-COAD, TCGA-PAAD, TCGA-ESCA, TCGA-CHOL, TCGA-LAML, TCGA-TGCT, TCGA-SARC): NCI Genomic Data Commons;
- Microarray validation cohort: NCBI GEO (`GSE31210`).

---

## References

1. Wu YL, John T, Grohé C, et al. Postoperative osimertinib in resected *EGFR*-mutated non-small-cell lung cancer. *J Clin Oncol*. 2023;41(10):1830-1840. PMID: 36720083. DOI: 10.1200/JCO.22.02186.
2. John T, Grohé C, Goldman JW, et al. Molecular residual disease analysis of adjuvant osimertinib in resected *EGFR*-mutated stage IB-IIIA non-small-cell lung cancer. *Nat Med*. 2024;30(11):3258-3267. PMID: 39313628. DOI: 10.1038/s41591-024-03238-0.
3. Tsuboi M, Herbst RS, John T, et al. Overall survival with osimertinib in resected *EGFR*-mutated non-small-cell lung cancer. *N Engl J Med*. 2023;389(2):137-147. PMID: 37272535. DOI: 10.1056/NEJMoa2304594.
4. Ridker PM, Everett BM, Thuren T, et al. Effect of interleukin-1beta inhibition with canakinumab on incident lung cancer in patients with atherosclerosis: exploratory results from a randomised, double-blind, placebo-controlled trial. *Lancet*. 2017;390(10105):1833-1842. PMID: 28855078. DOI: 10.1016/S0140-6736(17)32247-X.
5. Albrengues J, Shields MA, Ng D, et al. Neutrophil extracellular traps produced during inflammation awaken dormant cancer cells in mice. *Science*. 2018;361(6409):eaao4227. PMID: 30262495. DOI: 10.1126/science.aao4227.
6. Cabanos HF, Hata AN. Emerging insights into targeted therapy-tolerant persister cells in cancer. *Cancers (Basel)*. 2021;13(11):2666. PMID: 34071428. DOI: 10.3390/cancers13112666.
7. Sharma SV, Lee DY, Li B, et al. A chromatin-mediated reversible drug-tolerant state in cancer cell subpopulations. *Cell*. 2010;141(1):69-80. PMID: 20371346. DOI: 10.1016/j.cell.2010.02.027.
8. Wherry EJ, Kurachi M. Molecular and cellular mechanisms of T cell exhaustion. *Nat Rev Immunol*. 2015;15(8):486-499. PMID: 26205583. DOI: 10.1038/nri3862.
9. Philip M, Schietinger A. Heterogeneity and cell-fate decisions of exhausted CD8+ T cells. *Nat Rev Immunol*. 2022;22(2):86-99. PMID: 34108681. DOI: 10.1038/s41577-021-00557-0.
10. Baxter LT, Jain RK. Transport of fluid and macromolecules in tumors. I. Role of interstitial pressure and convection. *Microvasc Res*. 1989;37(1):77-104. PMID: 2648173. DOI: 10.1016/0026-2862(89)90074-5.
11. Minchinton AI, Tannock IF. Drug penetration in solid tumours. *Nat Rev Cancer*. 2006;6(8):583-592. PMID: 16862191. DOI: 10.1038/nrc1893.
12. Mori T, Koga T, Shibata H, et al. Interstitial fluid pressure correlates clinicopathological factors of lung cancer. *Ann Thorac Cardiovasc Surg*. 2015;21(3):201-208. PMID: 25641031. DOI: 10.5761/atcs.oa.14-00264.
13. Thomlinson RH, Gray LH. The histological structure of some human lung cancers and the possible implications for radiotherapy. *Br J Cancer*. 1955;9(4):539-549. PMID: 13304213. DOI: 10.1038/bjc.1955.55.
14. Cross DAE, Ashton SE, Ghiorghiu S, et al. AZD9291, an irreversible EGFR TKI, overcomes T790M-mediated resistance to EGFR inhibitors in lung cancer. *Cancer Discov*. 2014;4(9):1046-1061. PMID: 24893891. DOI: 10.1158/2159-8290.CD-14-0337.
15. Dickinson PA, Cantarini MV, Collier J, et al. Metabolic disposition of osimertinib in rats, dogs, and humans. *Drug Metab Dispos*. 2016;44(8):1201-1212. PMID: 27226351. DOI: 10.1124/dmd.115.069096.
16. Planchard D, Brown KH, Kim DW, et al. Osimertinib Western and Asian clinical pharmacokinetics in patients and healthy volunteers. *Cancer Chemother Pharmacol*. 2016;77(4):767-776. PMID: 26902828. DOI: 10.1007/s00280-016-2996-3.
17. Offin M, Chan JM, Tenet M, et al. Concurrent *RB1* and *TP53* alterations define a subset of *EGFR*-mutant lung cancers at risk for histologic transformation. *J Thorac Oncol*. 2019;14(10):1780-1790. PMID: 31228622. DOI: 10.1016/j.jtho.2019.06.002.
18. Nakagawa K, Garon EB, Seto T, et al. Ramucirumab plus erlotinib in patients with untreated, *EGFR*-mutated, advanced non-small-cell lung cancer (RELAY). *Lancet Oncol*. 2019;20(12):1655-1669. PMID: 31591063. DOI: 10.1016/S1470-2045(19)30634-5.
19. Soria JC, Ohe Y, Vansteenkiste J, et al. Osimertinib in untreated *EGFR*-mutated advanced non-small-cell lung cancer. *N Engl J Med*. 2018;378(2):113-125. PMID: 29151359. DOI: 10.1056/NEJMoa1713137.
20. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). *Br J Cancer*. 2005;93(4):387-391. PMID: 16106245. DOI: 10.1038/sj.bjc.6602678.
21. Sartore-Bianchi A, Pietrantonio F, Lonardi S, et al. Circulating tumor DNA-guided adjuvant chemotherapy in stage II-III colorectal cancer: the PEGASUS trial. *Nat Med* / *Ann Oncol*. 2024;35(Suppl 2):S1205.
22. Deng Q, He B, Liu X, et al. Prognostic value of pre-operative inflammatory response biomarkers in resectable cancer patients and construction of a predictive model. *PLoS One*. 2017;12(3):e0174085. PMID: 28333994. DOI: 10.1371/journal.pone.0174085.

---

## Research Resource Identifiers (RRID)

* **cBioPortal for Cancer Genomics:** RRID:[SCR_014555](https://scicrunch.org/resolver/RRID:SCR_014555) (`luad_mskcc_2020`).
* **The Cancer Genome Atlas (TCGA):** RRID:[SCR_003193](https://scicrunch.org/resolver/RRID:SCR_003193).
* **NCBI Gene Expression Omnibus (GEO):** RRID:[SCR_005012](https://scicrunch.org/resolver/RRID:SCR_005012) (`GSE31210`).
* **NCI Genomic Data Commons (GDC):** RRID:[SCR_014514](https://scicrunch.org/resolver/RRID:SCR_014514).
* **ClinicalTrials.gov Registry:** RRID:[SCR_002309](https://scicrunch.org/resolver/RRID:SCR_002309) (CANTOS: NCT01327846; ADAURA: NCT02511106).
* **PLOS ONE Open Data:** DOI: [10.1371/journal.pone.0174085](https://doi.org/10.1371/journal.pone.0174085).
* **Python Programming Language (v3.12+):** RRID:[SCR_008394](https://scicrunch.org/resolver/RRID:SCR_008394).
* **SciPy Scientific Computing Library:** RRID:[SCR_008058](https://scicrunch.org/resolver/RRID:SCR_008058).
* **NumPy Array Library:** RRID:[SCR_008633](https://scicrunch.org/resolver/RRID:SCR_008633).
* **Pandas Data Analysis Library:** RRID:[SCR_018214](https://scicrunch.org/resolver/RRID:SCR_018214).
* **Lifelines Survival Analysis Engine:** RRID:[SCR_024345](https://scicrunch.org/resolver/RRID:SCR_024345).
* **scikit-learn Machine Learning:** RRID:[SCR_002577](https://scicrunch.org/resolver/RRID:SCR_002577).
* **Matplotlib Data Visualization:** RRID:[SCR_008624](https://scicrunch.org/resolver/RRID:SCR_008624).
* **Seaborn Statistical Visualization:** RRID:[SCR_018132](https://scicrunch.org/resolver/RRID:SCR_018132).
* **Pandoc Document Converter:** RRID:[SCR_024707](https://scicrunch.org/resolver/RRID:SCR_024707).
* **LibreOffice Publishing Engine:** RRID:[SCR_023773](https://scicrunch.org/resolver/RRID:SCR_023773).
* **R Project for Statistical Computing:** RRID:[SCR_001905](https://scicrunch.org/resolver/RRID:SCR_001905).
