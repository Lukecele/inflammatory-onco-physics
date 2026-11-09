# Table 1: Baseline Demographic and Clinical Characteristics of the Discovery Cohort (MSKCC 2020, N=604)
## Stratified by Host-Tumor Biophysical Index (HTBI)

**Data Source:** Memorial Sloan Kettering Cancer Center (MSKCC) Real-World Cohort  
**File:** `/home/luca/mskcc-cohort-hf/mskcc_2020_real_cohort.csv`  
**Statistical Testing:** Pearson's $\chi^2$ test for categorical variables; Mann-Whitney $U$ test for continuous variables.  

| Clinical Characteristic | Overall Cohort (N = 604) | Low HTBI (N = 172) | High HTBI (N = 432) | Statistical p-value |
| :--- | :---: | :---: | :---: | :---: |
| **Sex, n (%)** | | | | **p = 0.0074** ($\chi^2$) |
| • Female | 402 (66.6%) | 129 (75.0%) | 273 (63.2%) | |
| • Male | 202 (33.4%) | 43 (25.0%) | 159 (36.8%) | |
| **Smoking History, n (%)** | | | | **p < 0.0001** ($1.4 \times 10^{-48}$) |
| • Never Smoker | 138 (22.8%) | 108 (62.8%) | 30 (6.9%) | |
| • Ever Smoker | 466 (77.2%) | 64 (37.2%) | 402 (93.1%) | |
| **Pathologic Stage (TNM), n (%)** | | | | **p < 0.0001** ($6.2 \times 10^{-7}$) |
| • Stage I | 447 (74.0%) | 146 (84.9%) | 301 (69.7%) | |
| • Stage II | 95 (15.7%) | 26 (15.1%) | 69 (16.0%) | |
| • Stage III | 62 (10.3%) | 0 (0.0%) | 62 (14.4%) | |
| **Predominant Histology, n (%)** | | | | **p < 0.0001** |
| • Lepidic (Well differentiated / G1) | 88 (14.6%) | 88 (51.2%) | 0 (0.0%) | |
| • Acinar / Papillary (Intermediate / G2) | 411 (68.0%) | 84 (48.8%) | 327 (75.7%) | |
| • Solid / Micropapillary (G3 / Desmoplastic) | 105 (17.4%) | 0 (0.0%) | 105 (24.3%) | |
| **EGFR Mutation Status, n (%)** | | | | **p = 0.573** |
| • Wild Type | 422 (69.9%) | 123 (71.5%) | 299 (69.2%) | |
| • Mutant | 182 (30.1%) | 49 (28.5%) | 133 (30.8%) | |
| **Clinical Recurrence Event, n (%)** | | | | **p = 0.0010** ($\chi^2 = 10.74$) |
| • Recurrence-Free (Censored) | 506 (83.8%) | 158 (91.9%) | 348 (80.6%) | |
| • Documented Relapse (Event) | 98 (16.2%) | **14 (8.1%)** | **84 (19.4%)** | **2.40x Fold Increase** |
| **Relapse-Free Survival (RFS)** | | | | **p = 0.0051** (Mann-Whitney) |
| • Median Follow-Up, Months (IQR) | 21.0 (11.2 - 29.2) | 23.8 (14.5 - 30.2) | 19.9 (10.3 - 28.7) | |

---

*Note: All counts, percentages, medians, and IQR values were calculated directly from the raw MSKCC 2020 cohort without imputation or synthetic data generation.*
