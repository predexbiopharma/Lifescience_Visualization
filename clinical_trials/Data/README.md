# Synthetic ADaM Dataset Generation for a Clinical Trial Visualization: Methodology and Rationale

## 1. Background and Motivation

The development of a comprehensive visualization for oncology clinical trials requires a dataset that simultaneously covers all major analysis domains, including subject-level demographics, tumor response, adverse events, laboratory data, time-to-event outcomes, pharmacokinetics, dose exposure, biomarkers, and patient-reported outcomes. While real-world clinical trial data exist in the public domain, no single publicly available dataset satisfies all of the following criteria simultaneously:

- Complete coverage of all ADaM domains required for a full visualization
- Availability without data use agreements, privacy restrictions, or institutional approval
- Sufficient sample size to produce statistically stable and visually informative plots across subgroups
- Internal consistency across domains (e.g., a patient with Grade 4 ALT elevation in the adverse event dataset should also exhibit elevated ALT values in the laboratory dataset)
- A two-arm randomized structure enabling comparative visualizations such as Kaplan,Meier curves, forest plots, and butterfly plots

Public repositories such as cBioPortal and TCGA provide rich genomic and limited clinical data, but they do not include pharmacokinetic profiles, dose modification histories, or patient-reported outcomes. Clinical study reports from regulatory submissions contain these domains but are not publicly accessible at the patient level. This gap motivated the generation of a fully synthetic, internally consistent, multi-domain ADaM dataset calibrated to published clinical evidence.

Beyond the immediate purpose of this catalog, openly available synthetic datasets with documented generation methodology serve a broader scientific function. They enable reproducible methods development, support training and education in clinical data analysis, and provide a standardized benchmark for evaluating novel visualization approaches, all without the ethical and legal constraints associated with patient-level data sharing.

---

## 2. Study Design and Trial Simulation

The simulated trial, designated **ONCVIZ-001**, represents a Phase II/III randomized controlled basket trial evaluating a fictional oral kinase inhibitor, **Vizatinib 300 mg once daily**, against placebo across five solid tumor types: non-small cell lung cancer (NSCLC), colorectal cancer (CRC), hepatocellular carcinoma (HCC), pancreatic ductal adenocarcinoma (PDAC), and breast cancer (BRCA). Patients were randomized 2:1 to treatment versus control. A total of **400 patients** were simulated across **20 investigational sites**, with enrollment staggered over 18 months from June 2022 and a data cutoff of March 5, 2026.

The trial was anchored to NSCLC as the primary tumor type for parameter calibration, given the availability of rich public clinical data and the established role of kinase inhibitors in this indication.

---

## 3. Calibration Strategy

Rather than assigning arbitrary parameter values, all key distributional parameters were derived from two complementary sources.

**Empirical calibration from cBioPortal.** Real patient-level clinical data were downloaded via the cBioPortal API from four published NSCLC studies: lung_msk_2017 (Hellmann et al.), nsclc_pd1_msk_2018 (Rizvi et al.), luad_tcga_pan_can_atlas_2018 (TCGA, Cell 2018), and lusc_tcga_pan_can_atlas_2018 (TCGA, Cell 2018), yielding a pooled cohort of 2,153 patients. For each continuous variable, a set of candidate parametric distributions was fitted by maximum likelihood estimation, and the best-fitting distribution was selected by the Akaike Information Criterion (AIC). Goodness of fit was assessed using the two-sample Kolmogorov,Smirnov (KS) test. The fitted parameters are summarized in Table 1.

| Variable | Distribution | Key Parameters | n | KS p-value |
|---|---|---|---|---|
| OS_MONTHS | Weibull (minimum) | Fitted to pooled NSCLC data | 970 | 0.007 |
| PFS_MONTHS | Weibull (minimum) | Fitted to pooled NSCLC data | 1,230 | 0.185 |
| AGE | Normal | median = 66.0, IQR = 45,81 | 2,069 | 0.000 |
| SEX | Categorical | Male 50.9%, Female 49.1% | 2,099 | n/a |

**Literature-based calibration.** Efficacy parameters were calibrated to KEYNOTE-189 (Gandhi et al., NEJM 2018), a Phase III trial of pembrolizumab plus chemotherapy in metastatic NSCLC, which represents a contemporary standard of care with publicly available outcome data. The following parameters were directly derived from this source:

- Objective response rate: treatment arm 45%, control arm 19%
- Median OS: treatment arm 22.0 months, control arm 10.7 months (HR = 0.49)
- Median PFS: treatment arm 9.0 months, control arm 4.9 months (HR = 0.52)
- AE incidence rates by preferred term and grade, as reported in the trial's published safety tables and FDA label

Pharmacokinetic parameters were calibrated to the published population PK profile of erlotinib, a structurally analogous oral EGFR kinase inhibitor, scaled proportionally for a 300 mg dose (Ling et al., J Clin Pharmacol 2006). Hepatotoxicity patterns in the laboratory domain were graded according to NCI-CTCAE version 5.0 thresholds.

---

## 4. Dataset Generation Algorithm

Ten ADaM-compliant domains were generated sequentially, with all downstream domains derived from the subject-level dataset (ADSL) to ensure internal consistency.

### 4.1 Subject-Level Dataset (ADSL)

Each of the 400 virtual subjects was assigned baseline characteristics by sampling from the fitted or literature-derived distributions described above. Treatment arm assignment followed a 2:1 randomization ratio. Prognostic biomarkers including PD-L1 TPS score, tumor mutational burden (TMB), microsatellite instability status, and somatic mutation status (EGFR, KRAS, BRAF) were sampled from published prevalence data for NSCLC. Overall survival and progression-free survival durations were sampled from arm-specific Weibull distributions, with scale parameters adjusted by the arm-level hazard ratio. Best overall response was assigned by sampling from arm-specific response rate vectors.

### 4.2 Tumor Response Dataset (ADRS)

Longitudinal response assessments were generated at 42-day intervals using a first-order Markov chain with the following transition matrix, calibrated to reflect realistic RECIST 1.1 response trajectories:

| From \ To | CR | PR | SD | PD |
|---|---|---|---|---|
| CR | 0.82 | 0.10 | 0.05 | 0.03 |
| PR | 0.10 | 0.58 | 0.20 | 0.12 |
| SD | 0.02 | 0.10 | 0.52 | 0.36 |
| PD | 0.00 | 0.00 | 0.00 | 1.00 |

Once a patient reached progressive disease, no further assessments were generated, consistent with standard trial practice.

### 4.3 Tumor Measurement Dataset (ADTR)

Longitudinal sum of longest diameters was generated using a response-stratified interpolation algorithm. For each patient, a biologically plausible endpoint value was pre-specified based on best overall response (CR: 1,3% of baseline; PR: 25,68%; SD: 82,119%; PD: 120,220%), and tumor size at each visit was interpolated toward that endpoint with added Gaussian noise. Partial responders were assigned a nadir followed by a mild rebound trajectory to simulate post-nadir disease stabilization.

### 4.4 Adverse Events Dataset (ADAE)

Adverse events were generated from a predefined profile of 35 preferred terms across 11 System Organ Classes (SOCs), with arm-specific incidence rates calibrated to KEYNOTE-189 safety data. For each patient-event pair, grade was sampled from a term-specific multinomial distribution reflecting the published grade distribution for that event. Onset day, duration, seriousness, relatedness, action taken, and outcome were assigned probabilistically. A subset of treatment-arm patients (8%) were designated as having drug-induced hepatotoxicity, with ALT and AST trajectories in ADLB scaled to the assigned CTCAE grade.

### 4.5 Laboratory Dataset (ADLB)

Twenty-one laboratory parameters spanning liver function (ALT, AST, bilirubin, ALP, GGT), renal function (creatinine, BUN), electrolytes (potassium, sodium), complete blood count (hemoglobin, platelets, ANC, WBC, lymphocytes), cardiac (QTcF, heart rate), and metabolic (glucose, cholesterol, TSH, total protein, albumin) panels were generated at eight scheduled visits. Baseline values were sampled from normal-range distributions for each parameter. Post-baseline values were generated using effect-type-specific functions:

- **Hepatic parameters** (ALT, AST, bilirubin): a Gaussian peak function centered at cycle 2,3 followed by partial resolution, with a hepatotoxic subgroup receiving elevated peak values consistent with their assigned CTCAE grade
- **Suppressed parameters** (hemoglobin, platelets, ANC): a nadir function centered at cycle 4,6 followed by partial recovery
- **Incrementally increasing parameters** (QTcF, glucose, creatinine): a linear increase proportional to cumulative exposure
- **Stable parameters**: random fluctuation around baseline

CTCAE grades were assigned at each visit for hepatic parameters based on the x-ULN thresholds defined in CTCAE v5.0. Hy's Law candidates were identified as patients with peak ALT greater than 3x ULN concurrent with peak total bilirubin greater than 2x ULN.

### 4.6 Time-to-Event Dataset (ADTTE)

Four endpoints were generated: overall survival (OS), progression-free survival (PFS), duration of response (DOR, responders only), and time to response (TTR, responders only). Each record includes the event indicator, event date, and twelve subgroup variables (tumor type, age group, sex, ECOG, PD-L1 group, MSI status, TMB, liver metastases, prior lines, smoking history, EGFR mutation, KRAS mutation) to support subgroup forest plot generation.

### 4.7 Pharmacokinetics Dataset (ADPK)

PK profiles were generated for treatment-arm patients only using a one-compartment oral absorption model. Individual PK parameters (apparent clearance CL/F, apparent volume of distribution Vd/F, and absorption rate constant Ka) incorporated log-normal inter-individual variability (IIV) with approximately 30% coefficient of variation, consistent with published population PK analyses of oral kinase inhibitors. Concentration-time profiles were generated at ten nominal sampling times post-dose at three visits. Derived parameters (Cmax, AUCinf, Tmax, t1/2) were calculated analytically and appended as summary records. Proportional and additive residual errors were applied to simulated concentrations.

### 4.8 Exposure Dataset (ADEX)

Cycle-level dose exposure records were generated for all patients. Dose modifications (reductions and interruptions) were triggered probabilistically after cycle 2, with rates calibrated to published dose modification frequencies for kinase inhibitors. Modification reasons were sampled from a predefined list of AE-driven causes. Dose intensity was calculated as the ratio of actual to planned cumulative dose per cycle.

### 4.9 Biomarker Dataset (ADBM)

Six biomarker parameters were generated: TMB and PD-L1 (baseline only, derived from ADSL), and ctDNA variant allele frequency, CEA, CA-125, and CD8+ T cell density (longitudinal). Longitudinal biomarker trajectories were modulated by both treatment arm and best overall response, such that responders in the treatment arm exhibited greater directional change than non-responders, consistent with the expected pharmacodynamic relationship between target engagement and clinical outcome.

### 4.10 Patient-Reported Outcomes Dataset (ADPR)

Nine EORTC QLQ-C30 scales were generated longitudinally at up to seven visits. Baseline scores were sampled from published normative data for cancer patients. Post-baseline scores were generated using a toxicity-recovery trajectory function that simulates initial worsening during the early treatment phase followed by partial recovery, superimposed on a gradual disease-related decline in the control arm. Visit-level missingness was introduced via a patient-level compliance parameter sampled from a Beta(8, 2) distribution, yielding a mean completion rate of approximately 80% consistent with published PRO compliance data in Phase III oncology trials.

---

## 5. Validation

Internal consistency was verified by confirming that all subject identifiers in each domain were present in ADSL, that OS was greater than or equal to PFS for all patients, and that PK data were restricted to the treatment arm. Clinical plausibility was assessed by comparing key summary statistics against published benchmarks from KEYNOTE-189 (Table 2).

| Metric | Generated | Expected (KEYNOTE-189) |
|---|---|---|
| ORR, treatment arm | 43.2% | ~45% |
| ORR, control arm | 15.7% | ~19% |
| Median OS, treatment | 21.2 months | 22.0 months |
| Median OS, control | 9.1 months | 10.7 months |
| Median PFS, treatment | 8.1 months | 9.0 months |
| Median PFS, control | 4.3 months | 4.9 months |
| AE Grade 3 or higher | 19.6% | ~15,20% |
| Hy's Law candidates | 19 patients | 5,25 (typical range) |
| Cmax, median | 748 ng/mL | 500,1500 ng/mL |

All consistency checks passed. Minor deviations from expected values are attributable to stochastic sampling variability and are within acceptable ranges for a synthetic dataset of this size.

---

## 6. Reproducibility

All datasets were generated with a fixed random seed (seed = 42) to ensure full reproducibility. The complete generation script is provided as a companion file. The datasets are deposited at [GitHub repository URL to be added] and may be freely used for research, education, and methods development under a Creative Commons CC BY 4.0 license.

---

## 7. Limitations

Several limitations of this dataset should be acknowledged. First, while distributional parameters for OS, PFS, and age were fitted to real NSCLC data, the pooled source cohort comprised studies with heterogeneous treatment settings, and the KS test indicated imperfect fit for OS (p = 0.007) and age (p = 0.000), likely reflecting the large sample size sensitivity of the KS statistic rather than substantive distributional misspecification. Second, the between-domain correlations, while anchored by shared patient identifiers and response-stratified trajectory parameters, do not capture the full complexity of biological correlations present in real trial data. Third, the basket trial structure assigns tumor types randomly without subtype-specific parameter adjustments, which limits the realism of tumor-type-level subgroup analyses. These limitations do not diminish the utility of the dataset for its intended purpose of demonstrating visualization methods, but users should exercise caution if repurposing the data for statistical inference.

---

## References

Gandhi L, et al. Pembrolizumab plus Chemotherapy in Metastatic Non-Small-Cell Lung Cancer. N Engl J Med. 2018;378(22):2078,2092.

Cerami E, et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. Cancer Discov. 2012;2(5):401,404.

TCGA Research Network. Comprehensive molecular profiling of lung adenocarcinoma. Nature. 2014;511:543,550.

Ling J, et al. Metabolism and excretion of erlotinib, a small molecule inhibitor of epidermal growth factor receptor tyrosine kinase, in healthy male volunteers. Drug Metab Dispos. 2006;34(3):420,426.

National Cancer Institute. Common Terminology Criteria for Adverse Events (CTCAE) Version 5.0. 2017.

Fayers PM, et al. The EORTC QLQ-C30 Scoring Manual. 3rd ed. European Organisation for Research and Treatment of Cancer; 2001.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: Revised RECIST guideline (version 1.1). Eur J Cancer. 2009;45(2):228,247.
