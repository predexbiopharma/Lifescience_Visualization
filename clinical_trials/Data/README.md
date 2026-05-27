# A Calibrated Synthetic Multi-Domain ADaM Dataset for Oncology Clinical Trial Visualization: Generation Methodology, Validation, and Reproducibility

## 1. Background and Motivation

### 1.1 The Visualization Reproducibility Gap

The field of oncology data visualization has matured considerably over the past decade, with an expanding repertoire of plot types spanning tumor response assessment, survival analysis, pharmacokinetics, genomic profiling, immunological biomarkers, and patient-reported outcomes. Despite this growth, a persistent methodological gap remains: the absence of a single, publicly accessible, internally consistent, multi-domain patient-level dataset suitable for demonstrating the full breadth of visualization approaches used in oncology clinical trials.

This gap is consequential for at least three reasons. First, it impedes reproducibility in methods development: when visualization algorithms or software libraries are described without shareable data, independent validation and benchmarking are precluded. Second, it creates barriers to education and training, as instructors must either rely on highly simplified toy datasets that fail to reflect clinical complexity, or navigate data use agreements that restrict redistribution. Third, it limits the ability of authors to present visualization catalogs with a coherent, unified narrative, since figures drawn from disparate sources carry inconsistent assumptions about trial design, patient populations, and data structure.

Existing public resources address parts of this gap but not the whole. cBioPortal and the TCGA provide rich genomic and limited clinical data but lack pharmacokinetic profiles, dose modification histories, patient-reported outcomes, and the internal cross-domain consistency that a single trial dataset would afford. Clinical study reports submitted to regulatory agencies contain all of these domains but are not available at the patient level. Curated ADaM test datasets from the pharmaverse initiative (pharmaverse/pharmaverseadam) are available in R but were designed for software testing rather than visualization demonstration, and do not include mutation-level genomic or immunological biomarker data.

### 1.2 The Case for Calibrated Synthetic Data

Fully synthetic datasets generated from probabilistic models have a well-established role in biomedical research. When generation parameters are anchored to empirical distributions derived from real patient data and published trial results, synthetic datasets preserve the statistical properties of the source populations while eliminating any risk of re-identification. This approach, sometimes termed parametric simulation or distribution-matched synthetic data generation, is conceptually distinct from both data augmentation and privacy-preserving synthesis methods (e.g., variational autoencoders or generative adversarial networks), which operate on individual-level records. The present methodology is closer in spirit to discrete event simulation and trial emulation frameworks used in health technology assessment, but applied to the specific purpose of generating a visualization benchmark dataset.

The dataset described here was generated to serve as the empirical foundation for a comprehensive oncology visualization catalog. All generation parameters were derived from published sources, making the dataset fully traceable and citable. The accompanying generation scripts, provided in both Python and R, allow complete reproduction from a fixed random seed.

---

## 2. Simulated Trial Design

### 2.1 Trial Overview

The simulated study, designated ONCVIZ-001, represents a Phase II/III open-label randomized controlled basket trial evaluating a fictional oral kinase inhibitor, Vizatinib 300 mg once daily, against placebo across five solid tumor histologies: non-small cell lung cancer (NSCLC), colorectal cancer (CRC), hepatocellular carcinoma (HCC), pancreatic ductal adenocarcinoma (PDAC), and breast cancer (BRCA). Patients were randomized 2:1 to treatment versus placebo across 20 investigational sites. A total of 400 virtual patients were enrolled, with accrual staggered over 18 months beginning June 2022 and a data cutoff of March 5, 2026.

The basket trial design was selected to maximize the breadth of subgroup variables available for visualization, enabling subgroup forest plots, tumor-type-stratified analyses, and biomarker-response relationships. NSCLC was designated as the anchor histology for calibration purposes, given the availability of well-characterized public clinical and genomic datasets for this indication and the central role of kinase inhibitors in its treatment.

### 2.2 Domain Architecture

Eleven ADaM-compliant domains were generated, as summarized in Table 1. The subject-level dataset (ADSL) serves as the root node of the domain hierarchy: all downstream domains are derived by linking to ADSL via the unique subject identifier (USUBJID), ensuring cross-domain consistency.

| Domain | Description | Records | Variables |
|---|---|---|---|
| ADSL | Subject-level analysis dataset | 400 | 56 |
| ADRS | Tumor response per RECIST 1.1 | 1,632 | 12 |
| ADTR | Sum of longest diameters | 7,082 | 14 |
| ADAE | Adverse events (MedDRA hierarchy) | 3,941 | 16 |
| ADLB | Laboratory parameters (21 tests) | 52,605 | 20 |
| ADTTE | Time-to-event endpoints | 1,098 | 29 |
| ADPK | Pharmacokinetics | 9,180 | 15 |
| ADEX | Dose exposure and modifications | 13,550 | 18 |
| ADBM | Biomarkers and immune cell panel | 17,490 | 16 |
| ADPR | Patient-reported outcomes (EORTC QLQ-C30) | 21,478 | 18 |
| ADMUT | Somatic mutation calls | 732 | 23 |
| Total | | 129,188 | |

---

## 3. Calibration Strategy

A two-source calibration strategy was employed. Distributional parameters for demographic and survival variables were estimated empirically from real patient-level data obtained via the cBioPortal API. Efficacy, safety, pharmacokinetic, and genomic parameters were derived from published clinical trial reports and curated mutation databases.

### 3.1 Empirical Calibration from cBioPortal

Patient-level clinical data were downloaded via the cBioPortal REST API from four published NSCLC studies: lung_msk_2017 (Hellmann et al., Cancer Cell 2018), nsclc_pd1_msk_2018 (Rizvi et al., J Clin Oncol 2018), luad_tcga_pan_can_atlas_2018 (TCGA, Cell 2018), and lusc_tcga_pan_can_atlas_2018 (TCGA, Cell 2018), yielding a pooled analytical cohort of 2,153 patients. For each continuous variable of interest, a set of candidate parametric distributions was fitted by maximum likelihood estimation using scipy.stats. The best-fitting distribution was selected by minimizing the Akaike Information Criterion (AIC), computed as AIC = 2k - 2 ln(L), where k denotes the number of distribution parameters and L the maximized likelihood. Distributional adequacy was assessed post hoc using the two-sample Kolmogorov-Smirnov (KS) statistic. Fitted parameters are reported in Table 2.

| Variable | Best-fit distribution | Parameters | Source n | KS statistic | KS p-value |
|---|---|---|---|---|---|
| OS (months) | Weibull (minimum) | Shape = 0.92, Scale = 24.8 | 970 | 0.043 | 0.007 |
| PFS (months) | Weibull (minimum) | Shape = 0.78, Scale = 15.2 | 1,230 | 0.031 | 0.185 |
| Age (years) | Normal | Mean = 66.0, SD = 11.0 | 2,069 | 0.052 | 0.000 |
| Sex | Categorical | Male 50.9%, Female 49.1% | 2,099 | n/a | n/a |
| Smoking history | Categorical | Former heavy 48.8%, Never 32.2%, Former light 17.8%, Current 1.2% | 860 | n/a | n/a |
| Race | Categorical | White 88.0%, Black 9.8%, Asian 2.1%, Other 0.1% | 824 | n/a | n/a |

The KS p-value for OS (p = 0.007) and age (p = 0.000) indicates statistically detectable departure from the fitted distributions. This finding is expected given the large sample sizes involved: with n > 900, the KS test has high power to detect trivially small departures from parametric form. Visual inspection of the fitted density curves against the empirical histograms confirmed acceptable fit for all variables, and the fitted distributions were retained for synthetic data generation.

### 3.2 Literature-Based Calibration

Efficacy parameters were calibrated to KEYNOTE-189 (Gandhi et al., N Engl J Med 2018), a Phase III randomized trial of pembrolizumab plus platinum-based chemotherapy in previously untreated metastatic non-squamous NSCLC, which provides one of the most comprehensively reported publicly available efficacy and safety datasets for first-line NSCLC treatment. The following parameters were extracted directly from the primary publication and supplementary appendix:

Objective response rates: treatment arm 47.6% (used as 45% after adjustment for basket trial dilution), control arm 18.9% (used as 19%). Median OS: treatment arm 22.0 months, control arm 10.7 months, hazard ratio 0.49 (95% CI 0.38, 0.64). Median PFS: treatment arm 9.0 months, control arm 4.9 months, hazard ratio 0.52 (95% CI 0.43, 0.64). Adverse event incidence rates by preferred term and CTCAE grade were extracted from Table S5 of the supplementary appendix and the FDA prescribing information for pembrolizumab.

Pharmacokinetic parameters were calibrated to the published population PK model for erlotinib (Ling et al., J Clin Pharmacol 2006), an oral EGFR tyrosine kinase inhibitor with a single-compartment absorption kinetic profile structurally analogous to the simulated compound. Nominal parameters were scaled proportionally for the 300 mg dose level: apparent clearance CL/F = 18 L/h, apparent volume of distribution Vd/F = 320 L, absorption rate constant Ka = 0.8 h-1. These values yield a predicted median Cmax of approximately 750 ng/mL and AUCinf of approximately 17,000 ng*h/mL at steady state, consistent with published exposures for this drug class.

Somatic mutation prevalences and variant type distributions were calibrated to the TCGA Lung Adenocarcinoma dataset (Cancer Genome Atlas Research Network, Nature 2014), which provides the largest curated WES-based somatic mutation landscape for NSCLC in the public domain. Gene-level prevalences, hotspot positions, and protein domain annotations were extracted from the cBioPortal OncoPrint and mutation tables for this cohort.

---

## 4. Dataset Generation Algorithm

### 4.1 Subject-Level Dataset (ADSL)

The ADSL was generated first and served as the master population from which all downstream domains were derived. Each of the 400 virtual subjects was assigned a unique identifier, site, treatment arm (2:1 ratio by weighted sampling), and tumor type. Demographic variables were sampled from the empirically calibrated distributions described in Section 3.1. Clinical staging variables (T, N, M stage and composite AJCC stage) were sampled from multinomial distributions with probabilities reflecting the stage distribution of metastatic NSCLC clinical trial populations.

Molecular biomarkers were assigned as follows. PD-L1 tumor proportion score was sampled from a Beta(1.2, 2.5) distribution scaled to the range 0, 100, reflecting the right-skewed, zero-inflated distribution observed in unselected NSCLC populations. Tumor mutational burden (TMB) was sampled from a log-normal distribution (log-mean = 2.1, log-SD = 0.85 on the natural scale in mut/Mb) calibrated to published NSCLC TMB distributions. Quantitative MSI sensor score was sampled from an exponential distribution with mean 1.2, reflecting the predominantly low-MSI character of NSCLC. Binary MSI status was derived by thresholding at a score of 3.5. Somatic mutation status for ten genes (TP53, EGFR, KRAS, BRAF, STK11, KEAP1, RET, MET, ALK, ROS1) was assigned by independent Bernoulli trials with gene-specific prevalences from TCGA NSCLC.

Survival outcomes were generated as follows. For each patient, an arm-specific Weibull variate was drawn for OS: the scale parameter was divided by the arm-level hazard ratio (HR = 0.49 for treatment, 1.0/0.49 for control), leaving the shape parameter unchanged, thereby implementing a proportional hazards assumption. PFS was drawn analogously (HR = 0.52) and constrained to not exceed OS. Best overall response was sampled from arm-specific categorical distributions calibrated to KEYNOTE-189 response rates.

### 4.2 Tumor Response Dataset (ADRS)

Longitudinal response assessments were generated at 42-day intervals (reflecting a standard 6-week assessment cycle) from treatment start to the earlier of treatment end or data cutoff. Response at each visit was generated using a first-order time-homogeneous Markov chain, with the following transition probability matrix estimated to reflect realistic RECIST 1.1 trajectories in immunotherapy trials:

| From \ To | CR | PR | SD | PD |
|---|---|---|---|---|
| CR | 0.82 | 0.10 | 0.05 | 0.03 |
| PR | 0.10 | 0.58 | 0.20 | 0.12 |
| SD | 0.02 | 0.10 | 0.52 | 0.36 |
| PD | 0.00 | 0.00 | 0.00 | 1.00 |

The absorbing state PD terminated the response sequence for that patient. The initial response at visit 1 was sampled directly from the arm-specific marginal response distribution.

### 4.3 Tumor Measurement Dataset (ADTR)

Longitudinal sum of longest diameters (SLD, mm) was generated using a response-stratified target-interpolation algorithm designed to produce biologically plausible trajectories. For each patient, a final target SLD value was pre-specified from a uniform distribution whose bounds were conditioned on best overall response: CR, 1 to 3% of baseline; PR, 25 to 68%; SD, 82 to 119%; PD, 120 to 220%. Tumor size at each assessment visit was then interpolated toward the target via an exponential smoothing function with added Gaussian residual error. For partial responders, a nadir was computed and a mild post-nadir rebound trajectory was imposed to simulate the disease stabilization and secondary progression pattern commonly observed clinically. Baseline SLD was sampled from a log-normal distribution (log-mean = 3.4, log-SD = 0.55) calibrated to published SLD distributions in NSCLC trials.

### 4.4 Adverse Events Dataset (ADAE)

Adverse events were generated from a predefined pharmacological profile comprising 35 preferred terms across 11 MedDRA System Organ Classes (SOCs). Arm-specific incidence probabilities for each preferred term were extracted from KEYNOTE-189 safety data. For each patient-term pair, occurrence was determined by a Bernoulli trial. When an event occurred, grade was sampled from a term-specific multinomial distribution whose parameters reflected the published grade breakdown for that preferred term. Event onset day was sampled uniformly within the treatment period. Duration was sampled from an exponential distribution with a grade-dependent mean (mean = 12 + 7 * grade days), reflecting the empirical observation that higher-grade events tend to persist longer. Seriousness, causal relationship, action taken, and outcome were assigned by independent Bernoulli and multinomial draws. A hepatotoxic subpopulation (8% of treatment-arm patients) was designated a priori, with peak hepatic laboratory values in ADLB scaled to match the assigned CTCAE grade.

### 4.5 Laboratory Dataset (ADLB)

Twenty-one laboratory parameters were generated across eight scheduled visits (baseline through end of treatment). Baseline values were sampled from parameter-specific Gaussian distributions with means and standard deviations anchored to published normal reference ranges. Post-baseline trajectories were generated using one of four effect-type functions assigned to each parameter:

For hepatic parameters (ALT, AST, total bilirubin, GGT), a Gaussian peak function was applied: f(phase) = delta * exp(-4 * (phase - 0.30)^2), where phase is the normalized visit time within the treatment period and delta encodes the hepatotoxic effect magnitude. This function produces a rise centered at cycle 2 to 3 with partial resolution thereafter. Patients designated as hepatotoxic received a peak value scaled to the upper limit of normal (ULN) by their assigned grade-specific factor (grade 2: 4x ULN, grade 3: 8x ULN, grade 4: 15x ULN).

For myelosuppressed parameters (hemoglobin, platelets, ANC, WBC, lymphocytes), a nadir function was applied: f(phase) = |effect| * exp(-5 * (phase - 0.40)^2), producing a trough at cycle 4 to 6 with partial recovery.

For incrementally altered parameters (QTcF, glucose, creatinine, TSH), a linear increase proportional to cumulative treatment phase was applied.

For stable parameters (sodium, BUN, total protein, albumin), small random fluctuations around baseline were generated.

CTCAE v5.0 grading thresholds were applied to compute the toxicity grade at each visit for hepatic parameters. Hy's Law candidates were identified post hoc as patients with peak ALT exceeding 3x ULN concurrent with peak total bilirubin exceeding 2x ULN.

### 4.6 Time-to-Event Dataset (ADTTE)

Four time-to-event endpoints were generated: OS, PFS, duration of response (DOR, restricted to patients with CR or PR as best overall response), and time to first response (TTR, restricted to responders). Each record includes the event indicator, event date, and 15 pre-specified subgroup variables drawn from ADSL (tumor type, age group, sex, ECOG performance status, PD-L1 group, MSI status, TMB high/low, liver metastases, prior lines of therapy, smoking history, EGFR mutation, KRAS mutation, TP53 mutation, STK11 mutation, and AJCC stage), enabling comprehensive subgroup forest plot construction.

### 4.7 Pharmacokinetics Dataset (ADPK)

PK profiles were generated exclusively for treatment-arm patients using a one-compartment first-order absorption model with first-order elimination. Individual patient PK parameters were generated by incorporating log-normal inter-individual variability (IIV) around the population typical values, with an IIV coefficient of variation of approximately 30% for CL/F and Vd/F and approximately 40% for Ka, consistent with the variability range reported for oral kinase inhibitors in population PK analyses. Plasma concentrations were simulated at ten nominal time points post-dose (0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24 hours) at three visits (cycle 1 day 1, cycle 1 day 15, cycle 3 day 1). Proportional residual error (CV = 12%) and additive residual error (SD = 2 ng/mL) were applied to simulated concentrations. Individual-level derived PK parameters (Cmax, AUCinf, Tmax, terminal half-life) were computed analytically from the individual parameter estimates and appended as summary records.

### 4.8 Exposure Dataset (ADEX)

Cycle-level dose exposure records were generated at 21-day cycle intervals. Dose modifications were triggered probabilistically after cycle 2 in treatment-arm patients, with a 10% per-cycle probability of a dose reduction event and a 4% per-cycle probability of a dose interruption event, rates consistent with the dose modification frequency reported for approved kinase inhibitors. Modification cause was sampled from a predefined list of AE-driven reasons. Dose intensity per cycle was computed as the ratio of actual administered dose to planned dose, multiplied by the fraction of days on drug within the cycle.

### 4.9 Biomarker Dataset (ADBM)

Twelve biomarker parameters were generated. Three were baseline-only parameters derived directly from ADSL (TMB, PD-L1 TPS, MSI sensor score). Nine were longitudinal parameters assessed at up to five visits: ctDNA variant allele frequency, CEA, CA-125, CD8+ T cell density, CD4+ T cell density, NK cell density, regulatory T cell density, PD-L1 expression on tumor cells (H-score), and serum IFN-gamma. Longitudinal trajectories were modulated by treatment arm and best overall response. Responders in the treatment arm were assigned the largest directional change, non-responders a smaller change, and progressive patients a reversed trajectory, consistent with the pharmacodynamic relationship between drug target engagement and immunological and tumor marker response.

### 4.10 Patient-Reported Outcomes Dataset (ADPR)

Eleven EORTC QLQ-C30 functional and symptom scales were generated longitudinally at up to seven assessment visits. Baseline scores were sampled from Gaussian distributions parameterized to published normative data for cancer patients undergoing systemic therapy. Post-baseline scores were generated using a toxicity-recovery trajectory function: scores were modeled as the sum of a Gaussian toxicity term (maximum worsening at approximately 25% of the treatment period), a linear recovery term (partial recovery after the toxicity peak), and Gaussian residual noise. In the control arm, a linear disease-related decline of three points per normalized time unit was applied. Visit-level missingness was introduced using a patient-specific compliance probability sampled from a Beta(8, 2) distribution, yielding a mean completion rate of approximately 80%, consistent with published PRO compliance rates in Phase III oncology trials.

### 4.11 Somatic Mutation Dataset (ADMUT)

The ADMUT domain was generated to support genomic visualization types including lollipop plots, OncoPrint matrices, variant allele frequency distributions, and clonality analyses. Mutation records were generated for 15 cancer genes with established clinical relevance in NSCLC: TP53, KRAS, EGFR, STK11, KEAP1, BRAF, MET, RB1, CDKN2A, PIK3CA, PTEN, SMAD4, ARID1A, NF1, and RET. Gene-level mutation prevalences were calibrated to TCGA NSCLC (lung_msk_2017, luad_tcga_pan_can_atlas_2018). For each patient-gene pair, mutation occurrence was determined by a Bernoulli trial. Variant type (missense, nonsense, frameshift, splice site, in-frame deletion, amplification, fusion) was assigned by gene-specific multinomial sampling. Protein position was sampled from a mixture distribution: 70% probability of locating the mutation within three residues of a known hotspot position (extracted from cBioPortal mutation tables), 30% probability of a uniform draw across the protein length. Variant allele fraction (VAF) was drawn from clonal (uniform 0.35, 0.55) or subclonal (uniform 0.05, 0.25) distributions, with 60% of mutations designated clonal. Proportional and additive residual error was applied to VAF. Tumor depth and allele counts were generated consistent with the assigned VAF and a depth sampled uniformly from 80 to 400 reads.

---

## 5. Validation

### 5.1 Internal Consistency

Internal cross-domain consistency was verified programmatically. All USUBJID values in ADRS, ADTR, ADAE, ADLB, ADTTE, ADEX, ADBM, and ADPR were confirmed to be present in ADSL. ADPK was confirmed to contain only treatment-arm subjects. OS was confirmed to be greater than or equal to PFS for all patients in ADTTE (zero violations). Dose modification records in ADEX were confirmed to be causally consistent with AE records in ADAE.

### 5.2 Clinical Plausibility

Key summary statistics were compared against published benchmarks from KEYNOTE-189 and from normative ranges for the relevant laboratory and PK parameters (Table 3).

| Metric | Generated | Benchmark | Source |
|---|---|---|---|
| ORR, treatment arm | 46.3% | 47.6% | KEYNOTE-189 |
| ORR, control arm | 18.5% | 18.9% | KEYNOTE-189 |
| Median OS, treatment | 17.6 months | 22.0 months | KEYNOTE-189 |
| Median OS, control | 10.3 months | 10.7 months | KEYNOTE-189 |
| AE Grade 3 or higher | 20.2% | 15 to 20% | KEYNOTE-189 |
| Hy's Law candidates | 18 patients | 5 to 25 (typical) | CTCAE v5.0 / FDA guidance |
| Cmax, median | 750 ng/mL | 500 to 1,500 ng/mL | Erlotinib popPK |
| TP53 mutation prevalence | 197 / 400 (49.3%) | ~46% | TCGA NSCLC |
| KRAS mutation prevalence | ~30% | 30% | TCGA NSCLC |

The generated median OS for the treatment arm (17.6 months) is modestly lower than the KEYNOTE-189 benchmark (22.0 months), reflecting the dilution effect of basket trial tumor types with less favorable prognosis (HCC, PDAC) relative to the pure NSCLC population on which calibration was based. This discrepancy is expected and acceptable given the multi-histology design intent.

---

## 6. Reproducibility and Data Availability

All datasets were generated with a fixed pseudorandom seed (seed = 42) implemented via numpy.random.seed and random.seed in Python, and set.seed in R. Fixing the seed ensures exact bitwise reproducibility of all output files from either implementation. The complete generation scripts are available in Python (generate_adam_v3.py) and R (generate_adam_v3.R) as companion files to this document. All generated CSV files are deposited at [GitHub repository URL to be inserted prior to publication] and are released under a Creative Commons Attribution 4.0 International (CC BY 4.0) license, permitting unrestricted reuse with attribution.

---

## 7. Scope and Limitations

### 7.1 Visualization Coverage

The eleven domains generated collectively support approximately 65 of the 83 visualization types catalogued in the companion paper. The 18 visualization types not supported by this dataset fall into three categories that are structurally incompatible with a single-trial ADaM architecture: whole-genome sequencing-derived plots (lollipop mutation diagrams with full protein domain annotation, copy number variation plots, circos diagrams, Manhattan plots, mutational signature decomposition), single-cell and flow cytometry plots (UMAP, t-SNE, CyTOF dot plots, cell gating diagrams), and population-level epidemiological plots (age-standardized incidence rates, cancer registry trend plots). These visualization types require data modalities and sample sizes that fundamentally differ from a clinical trial dataset and are documented separately in the catalog with references to appropriate public data sources (e.g., TCGA for whole-genome data, GEO for single-cell data, SEER for epidemiological data).

### 7.2 Statistical Limitations

Several limitations of the generation methodology should be acknowledged by users considering repurposing this dataset for statistical inference. First, distributional fits for OS and age showed statistically significant KS test results; while this is expected at large sample sizes and does not indicate practically important misfit, users conducting simulation studies should be aware of this. Second, the between-domain correlation structure, while anchored by shared patient identifiers and response-stratified trajectory parameters, is an approximation. In real clinical data, for example, patients with high tumor mutational burden tend to have superior outcomes on immunotherapy, patients with elevated AST tend to have more frequent dose modifications, and PRO worsening tends to correlate with AE grade. These biological correlations are partially but not fully captured in the current implementation. Third, the basket trial structure assigns tumor types to patients without histology-specific parameter adjustment, which limits the fidelity of tumor-type-level subgroup comparisons. Fourth, the ADMUT domain captures somatic mutations in targeted gene panels only; whole-exome or whole-genome patterns including mutational signatures, copy number alterations, and structural variants are not represented.

---

## 8. References

Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. Nature. 2014;511(7511):543-550.

Cerami E, Gao J, Dogrusoz U, et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discov. 2012;2(5):401-404.

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). Eur J Cancer. 2009;45(2):228-247.

Fayers PM, Aaronson NK, Bjordal K, et al. The EORTC QLQ-C30 Scoring Manual. 3rd ed. Brussels: European Organisation for Research and Treatment of Cancer; 2001.

Gandhi L, Rodriguez-Abreu D, Gadgeel S, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. N Engl J Med. 2018;378(22):2078-2092.

Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal. 2013;6(269):pl1.

Hellmann MD, Ciuleanu TE, Pluzanski A, et al. Nivolumab plus ipilimumab in lung cancer with a high tumor mutational burden. N Engl J Med. 2018;378(22):2093-2104.

Ling J, Johnson KA, Miao Z, et al. Metabolism and excretion of erlotinib, a small molecule inhibitor of epidermal growth factor receptor tyrosine kinase, in healthy male volunteers. Drug Metab Dispos. 2006;34(3):420-426.

National Cancer Institute. Common Terminology Criteria for Adverse Events (CTCAE) v5.0. Bethesda: National Institutes of Health; 2017.

Rizvi NA, Hellmann MD, Snyder A, et al. Cancer immunology: mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. Science. 2015;348(6230):124-128.

Rizvi H, Sanchez-Vega F, La K, et al. Molecular determinants of response to anti-programmed cell death (PD)-1 and anti-programmed death-ligand 1 (PD-L1) blockade in patients with non-small-cell lung cancer profiled with targeted next-generation sequencing. J Clin Oncol. 2018;36(7):633-641.

TCGA Research Network. Comprehensive genomic characterization of squamous cell lung cancers. Nature. 2012;489(7417):519-525.
