# Universal Clinical & Biophysical Study Rubric v1.0
Adapted from the Evangelos Karatzas (vagkaratzas-review) scoring architecture.

## Scoring System
| Severity | Points | Impact |
| :--- | :---: | :--- |
| **gate** | — | A failing gate multiplies the entire score by 0.4 |
| **blocker** | 3 | Critical failure; invalidates clinical or mathematical claim |
| **major** | 2 | Significant methodological flaw or missing standard |
| **minor** | 1 | Minor presentation, formatting or documentation omission |
| **nit** | 0 | Stylistic comment only |

Formula: `Score = 10 * (Earned Points / Total Applicable Points) * (0.4 if any gate fails)`

---

## 1. Hard Execution Gates

| ID | Severity | Criterion | Verification Command |
| :---: | :---: | :--- | :--- |
| **G1** | gate | **Reproducible Code Execution:** All pipeline scripts execute with exit code 0 without unhandled exceptions. | `python3 simulation/simulate_relapse.py && python3 pipeline/01_multivariate_cox_engine.py` |
| **G2** | gate | **Statistical Assumption Integrity:** Schoenfeld residuals test confirms proportional hazards (p > 0.05). | `lifelines.statistics.proportional_hazard_test` |
| **G3** | gate | **Ethical & Legal Anonymization:** 100% GDPR-compliant pseudonymization of cohort data (zero PII/PHI). | Inspection of `chol_omics_anonymized_stratification.tsv` |

---

## 2. Mathematical & Biophysical Formulation (Theory)

| ID | Severity | Criterion |
| :---: | :---: | :--- |
| **T1** | blocker | **First-Principles Derivation:** Tumor persistence modeled as a coupled reaction-diffusion PDE system with modified Bessel $K_0$ drug penetration. |
| **T2** | blocker | **Analytical Lead Time Reproduction:** The continuum equations analytically reproduce the 5.46-month ctDNA lead time invariant. |
| **T3** | major | **History-Dependent Immune Exhaustion:** Exhaustion modeled as an exponential functional of cumulative inflammatory exposure ($\int \Phi d\tau$). |
| **T4** | minor | **Carrying Capacity Saturation:** Logistic saturation limit ($1 - (N_S + N_R)/K$) enforced on cell populations. |

---

## 3. Clinical & Statistical Validation (REMARK Standards)

| ID | Severity | Criterion |
| :---: | :---: | :--- |
| **C1** | blocker | **Independent Multivariable Significance:** HTBI achieves HR > 2.0 and p < 0.05 controlling for pathologic stage, sex, smoking, and mutations. |
| **C2** | major | **Multi-Cohort External Validation:** Evaluated across multiple independent patient cohorts ($N > 1,500$ patients). |
| **C3** | major | **Decision Curve Analysis (DCA):** Model demonstrates positive Net Clinical Benefit over "treat-all" and "treat-none" across clinical decision thresholds. |
| **C4** | major | **Empirical Calibration:** Calibration curve evaluated at 24 months shows agreement between predicted and observed recurrence along the 45° line. |
| **C5** | major | **Bootstrap Discrimination Gain:** 500-fold bootstrap confirms statistically significant incremental C-index ($\Delta C > 0, p < 0.05$). |
| **C6** | minor | **Bedside Parameter Feasibility:** Index computable from routine €5 complete blood count (NLR), standard biochemistry (CRP, albumin), and routine H&E histology. |
| **C7** | major | **Raw Data Integrity (Table 1):** Demographic and clinical baseline counts computed directly from raw certified patient files without data fabrication. |

---

## 4. Documentation & Translational Impact

| ID | Severity | Criterion |
| :---: | :---: | :--- |
| **D1** | major | **Complete Unified Manuscript Draft:** Full paper structured with Abstract, Intro, Methods, Results, Discussion, and Figures (`MANUSCRIPT_DRAFT.md`). |
| **D2** | minor | **REMARK 20-Point Audit Report:** Complete item-by-item compliance audit documented (`REMARK_COMPLIANCE_REPORT.md`). |
| **D3** | minor | **Interactive Decision Tool:** Standalone bedside clinical calculator implemented and functional (`app_bedside_calculator.py`). |
| **D4** | minor | **High-Resolution Figures:** Publication-grade figures (Figures 1 through 7) generated and verified. |
