# Lifescience Visualization

**A systematic, production-ready catalog of clinical trial visualizations for oncology research, grounded in calibrated synthetic ADaM data, implemented in R, and structured for direct deployment in regulatory and scientific contexts.**

---

## Overview

Clinical trials generate some of the most consequential data in medicine. The analytical standards governing how that data is collected, modeled, and reported have advanced considerably over the past two decades, yet the standards governing how it is *visualized* have lagged behind. Regulatory submissions, peer-reviewed publications, and clinical team reviews continue to rely on dense summary tables that demand substantial interpretive effort and offer little intuition about what is actually happening at the level of individual patients or biological mechanisms.

This repository addresses that gap systematically. It is a comprehensive, domain-stratified implementation of the visualization methods used across the full spectrum of oncology clinical trial research, from individual patient response timelines to population-level pharmacokinetics, from adverse event profiling to genomic landscape characterization, from single-cell immunophenotyping to competing-risks survival analysis. Most implementations are built on a shared synthetic dataset, ONCVIZ-001, that follows ADaM regulatory standards and is fully traceable to published clinical trial data; a subset of the genomic and cellular visualizations (domains 03–04) instead use plot-local synthetic data calibrated to realistic biological ranges, since their source modalities (TCGA-style variant calls, single-cell panels) are not part of the ONCVIZ-001 ADaM schema. Data provenance is stated per plot in the tables below.

The goal is not a visualization gallery. It is a rigorous reference implementation: grounded in regulatory methodology, calibrated to published efficacy and safety benchmarks, implemented cleanly in R, and structured so that each visualization can be understood, validated, and deployed independently.

---

## Repository Status

> **Current version:** Active development, approaching v1.0 (Oncology)
>
> The shared ADaM synthetic dataset (`Data/V1/`) is complete and fully validated across 13 domains and 26,723 records. **Version 1 scope has been narrowed to the most commonly used visualizations** across domains 01–09, so v1 can ship as a tight, high-value reference implementation rather than an exhaustive gallery. Niche, low-frequency, or highly specialized plot types (including the entirety of domains 10–12) have been moved to the [Deferred to Future Versions](#deferred-to-future-versions) list below and are out of scope for v1. See the table below for current implementation status.

---

## Visualization Catalog

The catalog is organized into twelve analytical domains reflecting the distinct scientific questions that each class of visualization is designed to answer.

### 01 · Response Assessment

The central question of any oncology efficacy trial is whether the treatment works — and for whom. Response visualizations operate at the level of individual patients, making inter-patient heterogeneity legible and enabling the nuanced interpretation of efficacy that aggregate statistics systematically obscure.

| Visualization | Primary Data | Status |
|---|---|---|
| Waterfall Plot | ADTR, ADRS | ✅ Implemented |
| Spider Plot | ADTR | ✅ Implemented |
| Swimmer Plot | ADSL, ADRS, ADTR | ✅ Implemented |
| Best Overall Response (BOR) Plot | ADRS | ✅ Implemented |
| Tumor Burden Plot | ADTR | ✅ Implemented |
| Response Heatmap (ORR by Dose-Cohort × Tumor-Type) | ADRS | ✅ Implemented |

---

### 02 · Survival & Time-to-Event

Most Phase III oncology trials are powered to detect a survival endpoint. The Kaplan-Meier curve has been the standard representation of time-to-event data since Kaplan and Meier (1958) and remains indispensable because it conveys both the magnitude and the temporal uncertainty of a treatment effect simultaneously. This domain also includes model diagnostics, newer estimands increasingly required under the proportional hazards assumption (RMST), and composite endpoint frameworks gaining traction in regulatory submissions.

| Visualization | Primary Data | Status |
|---|---|---|
| Kaplan-Meier Curve | ADTTE | ✅ Implemented |
| Overall Survival (OS) Curve | ADTTE | ✅ Implemented |
| Progression-Free Survival (PFS) Curve | ADTTE | ✅ Implemented |
| Event-Free Survival (EFS) Curve | ADTTE | ✅ Implemented |
| Disease-Free Survival (DFS) Curve | ADTTE | ✅ Implemented |
| Time to Response (TTR) Plot | ADTTE | ✅ Implemented |
| Time to Progression (TTP) Plot | ADTTE | ✅ Implemented |
| Duration of Response (DOR) Plot | ADTTE | ✅ Implemented |
| Landmark Analysis Plot | ADTTE | ✅ Implemented |
| Competing Risks Curve (Cumulative Incidence) | ADTTE | ✅ Implemented |
| Restricted Mean Survival Time (RMST) Plot | ADTTE | ✅ Implemented |

---

### 03 · Biomarker & Genomics

Genomic and molecular biomarker visualizations span a wide methodological range, from individual variant annotation (lollipop plots) to population-level mutational pattern decomposition (signature plots), and from pairwise statistical testing (volcano plots) to whole-genome structural visualization (circos diagrams). ONCVIZ-001's ADaM schema does not include raw genomic call sets at the resolution these visualization types require, so — with the exception of the TMB plot, which reads live mutation calls from ADMUT — this domain's implementations use plot-local synthetic data, generated with a fixed `set.seed()` per script and calibrated to realistic variant frequencies, chromosomal distributions, and effect-size ranges. This keeps every figure fully reproducible while making explicit that it is illustrative rather than drawn from the shared trial dataset.

| Visualization | Primary Data | Status |
|---|---|---|
| Lollipop Plot (Mutation) | Synthetic (self-contained) | ✅ Implemented |
| OncoPrint / Oncoprint Heatmap | Synthetic (self-contained) | ✅ Implemented |
| Volcano Plot | Synthetic (self-contained) | ✅ Implemented |
| Forest Plot (Subgroup Analysis) | ADTTE | ✅ Implemented (see 08 · Meta-Analysis & Comparison) |
| Mutation Landscape Plot | Synthetic (self-contained) | ✅ Implemented |
| Copy Number Variation (CNV) Plot | Synthetic (self-contained) | ✅ Implemented |
| Circos Plot | Synthetic (self-contained) | ✅ Implemented |
| Manhattan Plot | Synthetic (self-contained) | ✅ Implemented |
| Miami Plot | Synthetic (self-contained) | ✅ Implemented |
| Rainfall Plot | Synthetic (self-contained) | ✅ Implemented |
| Mutational Signature Plot | Synthetic (self-contained) | ✅ Implemented |
| VAF (Variant Allele Frequency) Plot | Synthetic (self-contained) | ✅ Implemented |
| ctDNA Dynamics Plot | Synthetic (self-contained) | ✅ Implemented |
| TMB (Tumor Mutational Burden) Plot | ADMUT | ✅ Implemented |
| MSI (Microsatellite Instability) Plot | Synthetic (self-contained) | ✅ Implemented |

---

### 04 · Immunology & Cellular

Single-cell and flow cytometry-based visualizations represent a methodologically distinct class that typically requires data modalities differing structurally from standard ADaM clinical trial datasets. In this repository, five of the six implementations run on the longitudinal immune cell panel in ADBM (joined to ADSL for patient-level covariates) rather than external single-cell files, since ONCVIZ-001's ADBM panel is sufficient to demonstrate the visual grammar of each plot type; the CyTOF Dot Plot uses plot-local synthetic data.

| Visualization | Primary Data | Status |
|---|---|---|
| Flow Cytometry Plot (Scatter / Gating) | ADSL, ADBM | ✅ Implemented |
| UMAP Plot | ADSL, ADBM | ✅ Implemented |
| t-SNE Plot | ADSL, ADBM | ✅ Implemented |
| Cell Composition Bar Plot | ADSL, ADBM | ✅ Implemented |
| Immune Cell Infiltration Heatmap | ADSL, ADBM | ✅ Implemented |
| CyTOF Dot Plot | Synthetic (self-contained) | ✅ Implemented |

---

### 05 · Safety & Toxicity

Safety is not a secondary concern in clinical trials; in Phase I, it is the primary endpoint. Adverse event data is high-dimensional, hierarchically coded via MedDRA, and graded by CTCAE severity, making it one of the most analytically demanding domains to visualize rigorously. This domain is fully supported by ADAE, joined to ADSL for patient-level covariates.

| Visualization | Primary Data | Status |
|---|---|---|
| Adverse Event (AE) Bar Chart | ADAE, ADSL | ✅ Implemented |
| Toxicity Heatmap | ADAE, ADSL | ✅ Implemented |

*Dose-Limiting Toxicity (DLT) Plot, Time-to-Toxicity Plot, Exposure-Response Plot, and Dose Escalation Plot (3+3/BOIN) are Phase I dose-finding niche visualizations — deferred to v2 (see below).*

---

### 06 · Pharmacokinetics / Pharmacodynamics (PK/PD)

PK/PD visualizations underpin dose selection, exposure-response characterization, and the regulatory justification of the proposed dose. Population PK modeling and visual predictive checks are standard components of NDA and BLA submissions. This domain is fully supported by ADPK, ADBM, ADTR, and ADSL.

| Visualization | Primary Data | Status |
|---|---|---|
| PK Concentration-Time Curve | ADPK | ✅ Implemented |
| Trough Level Plot | ADPK | ✅ Implemented |
| PD Biomarker Plot | ADBM | ✅ Implemented |
| Exposure-Efficacy Plot | ADPK, ADTR | ✅ Implemented |
| Waterfall + PK Overlay | ADSL, ADTR, ADPK | ✅ Implemented |

---

### 07 · Imaging & Tumor Measurement

Target lesion measurement trajectories and scan-level timeline visualizations bridge the gap between radiology data and clinical interpretation, supporting the assessment of response kinetics and measurement variability across sites and time points. Supported by ADTR, joined to ADSL for best-response labeling.

| Visualization | Primary Data | Status |
|---|---|---|
| Sum of Longest Diameters (SLD) Over Time | ADTR, ADSL | ✅ Implemented |
| Target Lesion Change Plot | ADTR, ADSL | ✅ Implemented |

*Scan Timeline Plot is a lower-frequency visualization — deferred to v2 (see below).*

---

### 08 · Meta-Analysis & Comparison

Meta-analytic visualizations synthesize evidence across trials, treatment arms, or subgroups. The forest plot is among the most common figure types in published randomized controlled trials and is central to regulatory label negotiations. Network meta-analysis plots, benefit-risk visualizations, and funnel plots address the increasingly structured quantitative frameworks used in health technology assessment.

| Visualization | Primary Data | Status |
|---|---|---|
| Forest Plot | ADTTE | ✅ Implemented |
| Benefit-Risk Plot | ADSL, ADAE | ✅ Implemented |

*Funnel Plot, Network Meta-Analysis (NMA) Plot, and Tornado Plot require external multi-trial data and are lower-frequency in single-trial reporting — deferred to v2 (see below).*

---

### 09 · Trial Design & Patient Flow

Trial flow and exposure visualizations document the operational execution of a study: enrollment trajectories, treatment compliance, dose intensity, and patient disposition. CONSORT diagrams are required components of randomized trial publications under ICMJE reporting standards. Supported by ADRAND, ADSL, ADRS, and ADEX.

| Visualization | Primary Data | Status |
|---|---|---|
| CONSORT Diagram | ADRAND, ADSL | ✅ Implemented |
| Enrollment Over Time Plot | ADRAND, ADSL | ✅ Implemented |
| Treatment Exposure Plot | ADSL, ADRS | ✅ Implemented |
| Dose Intensity Plot | ADEX, ADSL | ✅ Implemented |

*Relative Dose Intensity (RDI) Plot is a more granular variant of Dose Intensity Plot — deferred to v2 (see below).*

---

### 10 · Cell Therapy / CAR-T Specific

Cellular immunotherapy trials produce visualization challenges that have no analog in conventional pharmacological trials: CAR-T cell expansion kinetics spanning orders of magnitude, cytokine release syndrome timelines requiring hour-level resolution, and bone marrow response assessments reflecting a distinct efficacy paradigm from RECIST-based solid tumor measurement. Partially supported by ADBM and ADRS.

> **Out of scope for v1** — modality-specific to cell therapy trials, not general-purpose. All three plot types in this domain are deferred to v2 (see [Deferred to Future Versions](#deferred-to-future-versions)).

---

### 11 · Radiomics & Imaging Analytics

Radiomics-based visualizations extract quantitative imaging features from radiology scans and connect them to clinical outcomes. These approaches require imaging-derived feature matrices that extend beyond standard ADaM data structures; where ADTR cannot serve as a proxy, appropriate public imaging datasets (TCIA) are referenced.

> **Out of scope for v1** — requires external imaging feature data and is not part of standard clinical trial reporting. All plot types in this domain are deferred to v2 (see [Deferred to Future Versions](#deferred-to-future-versions)).

---

### 12 · Epidemiology & Incidence

Population-level cancer epidemiology visualizations require registry-level datasets that are structurally incompatible with the single-trial ADaM architecture. This domain references SEER and GLOBOCAN as authoritative public data sources.

> **Out of scope for v1** — population-registry visualizations, not single-trial reporting outputs. All plot types in this domain are deferred to v2 (see [Deferred to Future Versions](#deferred-to-future-versions)).

---

## Deferred to Future Versions

The following plot types are lower-frequency, modality-specific, or dependent on external data outside the core ADaM/ONCVIZ-001 scope. They remain part of the long-term catalog vision but are **not** part of v1 and will be added in v2 or later, prioritized roughly in the order listed.

**05 · Safety & Toxicity**
- Dose-Limiting Toxicity (DLT) Plot
- Time-to-Toxicity Plot
- Exposure-Response Plot
- Dose Escalation Plot (3+3 / BOIN)

**07 · Imaging & Tumor Measurement**
- Scan Timeline Plot

**08 · Meta-Analysis & Comparison**
- Funnel Plot
- Network Meta-Analysis (NMA) Plot
- Tornado Plot

**09 · Trial Design & Patient Flow**
- Relative Dose Intensity (RDI) Plot

**10 · Cell Therapy / CAR-T Specific**
- CAR-T Cell Expansion Curve
- Cytokine Release Syndrome (CRS) Timeline
- Bone Marrow Response Plot

**11 · Radiomics & Imaging Analytics**
- Heatmap Overlay (Tumor Heterogeneity)
- Radiomics Feature Importance Plot
- Lesion Size Over Time Plot

**12 · Epidemiology & Incidence**
- Age-Standardized Incidence Rate (ASR) Plot
- Cancer Incidence Trend Line
- Mortality-to-Incidence Ratio Plot
- Stage Distribution Bar Chart
- Prevalence Pie / Donut Chart

---

## Phase 1 Benchmark Artifacts

`Phase1_Benchmark/` contains the supporting artifacts for the calibration study "Human-Calibrated Validation of Large Language Model-Generated Oncology Clinical Trial Visualizations" (5 plot types, 20 injected errors, 228 scored responses from Claude Opus 4.8 and GPT-5.6).

> **Provenance notice:** `Phase1_Benchmark/raw_outputs/combined_final_strict_prompt.csv` (the 228 raw model responses) is the exact, byte-identical data underlying every reported metric. The injector, checker, and response-parsing scripts in that folder are **faithful reconstructions of the original implementation's documented specification, not the byte-identical original files**, which are no longer available; three checker functions are verbatim from the original session and are marked as such in the code. See [`Phase1_Benchmark/README.md`](Phase1_Benchmark/README.md) for the full breakdown of what is exact versus reconstructed.

---

## Synthetic Dataset: ONCVIZ-001

All visualizations in domains 01–02 and 05–09 are demonstrated on a single shared synthetic ADaM dataset — **ONCVIZ-001** — simulating a Phase I/II open-label dose-escalation and randomized basket trial of a fictional oral kinase inhibitor (Vizatinib 300 mg QD) across five solid tumor histologies (NSCLC, CRC, HCC, PDAC, BRCA). Domains 03–04 (Biomarker/Genomics and Immunology/Cellular) draw partly on this dataset (ADBM, ADSL) and partly on plot-local synthetic data calibrated separately, as noted per plot above. The dataset comprises 80 virtual patients (Phase I: 20 treatment-only; Phase II: 42 treatment + 18 control) with a data cutoff of March 5, 2026. Across both phases, 62 patients are in the TREATMENT arm and 18 are in the CONTROL arm (verified against `ADSL.csv`, `ARM` variable); this matches the n = 62 cited in individual plot READMEs (e.g., Treatment Exposure Plot).

Each histology carries independently calibrated response rates, survival parameters, mutation prevalences, toxicity profiles, and mutational signatures, enabling biologically valid per-tumor subgroup analyses. Generation parameters were anchored to empirical distributions from cBioPortal and to published trial benchmarks from KEYNOTE-189 (Gandhi et al., *N Engl J Med* 2018), IMbrave150 (Finn et al., *N Engl J Med* 2020), OlympiAD (Robson et al., *N Engl J Med* 2017), the TCGA PanCancer Atlas, COSMIC v3.3, and the erlotinib population PK model (Ling et al., *J Clin Pharmacol* 2006). All outputs are exactly reproducible from `seed = 42`.

| Domain | File | Records | Description |
|---|---|---|---|
| ADSL | `Data/V1/ADSL.csv` | 80 | Subject-level: demographics, treatment arm, survival outcomes, biomarker and mutation status |
| ADRS | `Data/V1/ADRS.csv` | 769 | RECIST 1.1 response assessments per visit, Markov-chain trajectories |
| ADTR | `Data/V1/ADTR.csv` | 769 | Sum of longest diameters (SLD) over time |
| ADAE | `Data/V1/ADAE.csv` | 752 | Adverse events with MedDRA coding, CTCAE v5 grading, histology-specific incidence |
| ADLB | `Data/V1/ADLB.csv` | 11,880 | Laboratory parameters (20 tests incl. immune cell panel, CTCAE grading for all applicable parameters) |
| ADTTE | `Data/V1/ADTTE.csv` | 377 | Time-to-event: OS, PFS, EFS, TTP, DOR, TTR, DFS with 15 subgroup variables and landmark flags |
| ADPK | `Data/V1/ADPK.csv` | 1,984 | Plasma PK profiles (1-compartment model, treatment arm only) |
| ADEX | `Data/V1/ADEX.csv` | 2,361 | Dose exposure and cycle-level modifications including re-escalation records |
| ADBM | `Data/V1/ADBM.csv` | 3,321 | Longitudinal biomarkers and immune cell panel |
| ADPR | `Data/V1/ADPR.csv` | 3,888 | Patient-reported outcomes (EORTC QLQ-C30, 8 scales) |
| ADMUT | `Data/V1/ADMUT.csv` | 187 | Somatic mutation calls (15 cancer genes, tumor-stratified prevalence) |
| ADRAND | `Data/V1/ADRAND.csv` | 92 | Screening and randomization log (92 screened, 12 screen failures, 80 randomized) |
| ADSIG | `Data/V1/ADSIG.csv` | 263 | Mutational signatures (SBS, COSMIC v3.3, 12 unique SBS signatures) |
| **Total** | | **26,723** | |

The CSVs above ship compressed as `Data/V1/oncoviz_data_v1.zip` and must be extracted before running any script (`unzip Data/V1/oncoviz_data_v1.zip -d Data/V1/`); they are regenerated deterministically by `Data/V1/generate_adam_synthetic.R`. For complete dataset documentation including calibration methodology, domain architecture, internal consistency validation, and full reference list, see [`Data/README.md`](Data/README.md).

---

## Repository Structure

```
Lifescience_Visualization/
│
├── README.md
├── Data/
│   └── V1/                            ← Shared synthetic ADaM dataset (ONCVIZ-001)
│       ├── oncoviz_data_v1.zip        ← ADSL, ADRS, ADTR, ADAE, ADLB, ADTTE, ADPK,
│       │                                 ADEX, ADBM, ADPR, ADMUT, ADRAND, ADSIG (.csv)
│       └── generate_adam_synthetic.R  ← Deterministic generator (seed = 42)
│
└── clinical_trials/
    └── Oncology/
        ├── Response_Assessment/       ← 6 plot types
        │   ├── Bor_plot/
        │   ├── Response_Heatmap_plot/
        │   ├── Spider_plot/
        │   ├── swimmer_plot/
        │   ├── Tumor_Burden_plot/
        │   └── Waterfall_plot/
        ├── Survival_TimeToEvent/      ← 11 plot types (2 scripts)
        │   ├── SurvivalCurves_plot/
        │   └── TimetoEvent_plot/
        ├── Biomarker_Genomics/        ← 14 plot types
        ├── Immunology_Cellular/       ← 6 plot types
        ├── Meta_Analysis/             ← Forest Plot, Benefit-Risk Plot (v1 scope only)
        │   ├── forest_plot/
        │   └── Benefit_Risk_plot/
        ├── PK_PD/                     ← 5 plot types
        ├── Safety_Toxicity/           ← 2 plot types
        ├── Imaging_Tumor_Measurement/ ← 2 plot types
        └── Trial_Design_Patient_Flow/ ← 4 plot types
# Cell_Therapy_CART, Radiomics_Imaging, and Epidemiology_Incidence are deferred to
# v2+ and intentionally do not exist in the repo yet — see "Deferred to Future Versions" above.
```

Each visualization folder follows a uniform structure:

```
Plot_Name/
├── README.md          ← analytical rationale, design decisions, regulatory context, references
├── plot_name.R        ← standalone R script (reads from ../../Data/V1/)
└── output/             ← example outputs at publication resolution (300 dpi)
```

---

## Standards & Methodological Grounding

The visualizations in this repository are implemented in accordance with published regulatory and analytical standards:

- **RECIST 1.1** — Eisenhauer EA et al. *Eur J Cancer* 2009;45(2):228–247.
- **iRECIST** — Seymour L et al. *Lancet Oncol* 2017;18(3):e143–e152.
- **ICH E9 (R1)** — Addendum on Estimands and Sensitivity Analysis in Clinical Trials. 2019.
- **ICH E3** — Structure and Content of Clinical Study Reports. 1995.
- **FDA Exposure-Response Guidance** — Study Design, Data Analysis, and Regulatory Applications. 2019.
- **CTCAE v5.0** — National Cancer Institute. 2017.
- **ADaM Implementation Guide** — CDISC ADaM Team. Version 1.3. 2021.
- **CONSORT 2010** — Schulz KF et al. *JAMA* 2010;303(7):681–683.
- Chia PL et al. Current and evolving methods to visualize biological data in cancer research. *JNCI* 2016;108(8).

---

## Requirements

```r
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2",
                   "survival", "survminer", "scales", "ggrepel"))
```

---

## License

Dataset released under **Creative Commons Attribution 4.0 International (CC BY 4.0)** — see [`Data/LICENSE`](Data/LICENSE).
Code released under **MIT License** — see [`LICENSE`](LICENSE).
