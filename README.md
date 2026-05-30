# Lifescience Visualization

**A systematic, production-ready catalog of clinical trial visualizations for oncology research, grounded in calibrated synthetic ADaM data, implemented in Python, and structured for direct deployment in regulatory and scientific contexts.**

---

## Overview

Clinical trials generate some of the most consequential data in medicine. The analytical standards governing how that data is collected, modeled, and reported have advanced considerably over the past two decades, yet the standards governing how it is *visualized* have lagged behind. Regulatory submissions, peer-reviewed publications, and clinical team reviews continue to rely on dense summary tables that demand substantial interpretive effort and offer little intuition about what is actually happening at the level of individual patients or biological mechanisms.

This repository addresses that gap systematically. It is a comprehensive, domain-stratified implementation of the visualization methods used across the full spectrum of oncology clinical trial research, from individual patient response timelines to population-level pharmacokinetics, from adverse event profiling to genomic landscape characterization, from single-cell immunophenotyping to competing-risks survival analysis. Each implementation is built on a shared synthetic dataset, ONCVIZ-001, that follows ADaM regulatory standards and is fully traceable to published clinical trial data, making every figure in this catalog reproducible, citable, and adaptable to real study data with minimal modification.

The goal is not a visualization gallery. It is a rigorous reference implementation: grounded in regulatory methodology, calibrated to published efficacy and safety benchmarks, implemented cleanly in Python, and structured so that each visualization can be understood, validated, and deployed independently.

---

## Repository Status

> **Current version:** Active development, Phase I (Oncology)
>
> The shared ADaM synthetic dataset (`Data/`) is complete and fully validated across 11 domains and 129,188 records. Visualization modules are being re-implemented progressively against this unified dataset, replacing earlier per-folder data copies. See the table below for current implementation status.

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
| Tumor Burden Plot | ADTR | 🔄 In progress |

---

### 02 · Survival & Time-to-Event

Most Phase III oncology trials are powered to detect a survival endpoint. The Kaplan-Meier curve has been the standard representation of time-to-event data since Kaplan and Meier (1958) and remains indispensable because it conveys both the magnitude and the temporal uncertainty of a treatment effect simultaneously. This domain also includes model diagnostics (log-log plots, Schoenfeld residuals), newer estimands increasingly required under the proportional hazards assumption (RMST), and composite endpoint frameworks gaining traction in regulatory submissions.

| Visualization | Primary Data | Status |
|---|---|---|
| Kaplan-Meier Curve | ADTTE | 🔜 Planned |
| Overall Survival (OS) Curve | ADTTE | 🔜 Planned |
| Progression-Free Survival (PFS) Curve | ADTTE | 🔜 Planned |
| Event-Free Survival (EFS) Curve | ADTTE | 🔜 Planned |
| Disease-Free Survival (DFS) Curve | ADTTE | 🔜 Planned |
| Time to Response (TTR) Plot | ADTTE | 🔜 Planned |
| Time to Progression (TTP) Plot | ADTTE | 🔜 Planned |
| Duration of Response (DOR) Plot | ADTTE | 🔜 Planned |
| Landmark Analysis Plot | ADTTE | 🔜 Planned |
| Competing Risks Curve (Cumulative Incidence) | ADTTE | 🔜 Planned |
| Restricted Mean Survival Time (RMST) Plot | ADTTE | 🔜 Planned |

---

### 03 · Biomarker & Genomics

Genomic and molecular biomarker visualizations span a wide methodological range, from individual variant annotation (lollipop plots) to population-level mutational pattern decomposition (signature plots), and from pairwise statistical testing (volcano plots) to whole-genome structural visualization (circos diagrams). This domain draws primarily on ADMUT and ADBM, supplemented by ADSL for co-variates.

| Visualization | Primary Data | Status |
|---|---|---|
| Lollipop Plot (Mutation) | ADMUT | 🔜 Planned |
| OncoPrint / Oncoprint Heatmap | ADMUT, ADSL | 🔜 Planned |
| Volcano Plot | ADBM | 🔜 Planned |
| Forest Plot (Subgroup Analysis) | ADTTE, ADSL | 🔜 Planned |
| Mutation Landscape Plot | ADMUT | 🔜 Planned |
| Copy Number Variation (CNV) Plot | External (TCGA) | 🔜 Planned |
| Circos Plot | External (TCGA) | 🔜 Planned |
| Manhattan Plot | External (TCGA) | 🔜 Planned |
| Miami Plot | External (TCGA) | 🔜 Planned |
| Rainfall Plot | ADMUT | 🔜 Planned |
| Mutational Signature Plot | External (TCGA / COSMIC) | 🔜 Planned |
| VAF (Variant Allele Frequency) Plot | ADMUT | 🔜 Planned |
| ctDNA Dynamics Plot | ADBM | 🔜 Planned |
| TMB (Tumor Mutational Burden) Plot | ADSL, ADMUT | 🔜 Planned |
| MSI (Microsatellite Instability) Plot | ADSL | 🔜 Planned |

---

### 04 · Immunology & Cellular

Single-cell and flow cytometry-based visualizations represent a methodologically distinct class that requires data modalities, cell-level event tables, spectral unmixing outputs, single-cell transcriptomic count matrices, that differ structurally from standard ADaM clinical trial datasets. This domain is partially supported by longitudinal immune cell panel data in ADBM; plot types requiring single-cell resolution reference appropriate public data sources (GEO, Human Cell Atlas).

| Visualization | Primary Data | Status |
|---|---|---|
| Flow Cytometry Plot (Scatter / Gating) | External (GEO) | 🔜 Planned |
| UMAP Plot | External (GEO / HCA) | 🔜 Planned |
| t-SNE Plot | External (GEO / HCA) | 🔜 Planned |
| Cell Composition Bar Plot | ADBM | 🔜 Planned |
| Immune Cell Infiltration Heatmap | ADBM | 🔜 Planned |
| CyTOF Dot Plot | External (GEO) | 🔜 Planned |

---

### 05 · Safety & Toxicity

Safety is not a secondary concern in clinical trials, in Phase I, it is the primary endpoint. Adverse event data is high-dimensional, hierarchically coded via MedDRA, and graded by CTCAE severity, making it one of the most analytically demanding domains to visualize rigorously. This domain is fully supported by ADAE and ADLB.

| Visualization | Primary Data | Status |
|---|---|---|
| Adverse Event (AE) Bar Chart | ADAE | 🔜 Planned |
| Dose-Limiting Toxicity (DLT) Plot | ADAE, ADEX | 🔜 Planned |
| Toxicity Heatmap | ADAE, ADLB | ✅ Implemented |
| Time-to-Toxicity Plot | ADAE | 🔜 Planned |
| Exposure-Response Plot | ADPK, ADAE | 🔜 Planned |
| Dose Escalation Plot (3+3 / BOIN) | ADEX, ADAE | 🔜 Planned |

---

### 06 · Pharmacokinetics / Pharmacodynamics (PK/PD)

PK/PD visualizations underpin dose selection, exposure-response characterization, and the regulatory justification of the proposed dose. Population PK modeling and visual predictive checks are standard components of NDA and BLA submissions. This domain is fully supported by ADPK and ADBM.

| Visualization | Primary Data | Status |
|---|---|---|
| PK Concentration-Time Curve | ADPK | 🔜 Planned |
| Trough Level Plot | ADPK | 🔜 Planned |
| PD Biomarker Plot | ADBM, ADPK | 🔜 Planned |
| Exposure-Efficacy Plot | ADPK, ADRS | 🔜 Planned |
| Waterfall + PK Overlay | ADTR, ADPK | 🔜 Planned |

---

### 07 · Imaging & Tumor Measurement

Target lesion measurement trajectories and scan-level timeline visualizations bridge the gap between radiology data and clinical interpretation, supporting the assessment of response kinetics and measurement variability across sites and time points. Supported by ADTR and ADRS.

| Visualization | Primary Data | Status |
|---|---|---|
| Sum of Longest Diameters (SLD) Over Time | ADTR | 🔜 Planned |
| Target Lesion Change Plot | ADTR | 🔜 Planned |
| Scan Timeline Plot | ADRS, ADTR | 🔜 Planned |

---

### 08 · Meta-Analysis & Comparison

Meta-analytic visualizations synthesize evidence across trials, treatment arms, or subgroups. The forest plot is among the four most common figure types in published randomized controlled trials and is central to regulatory label negotiations. Network meta-analysis plots, benefit-risk visualizations, and funnel plots address the increasingly structured quantitative frameworks used in health technology assessment.

| Visualization | Primary Data | Status |
|---|---|---|
| Forest Plot | ADTTE, ADSL | 🔜 Planned |
| Funnel Plot | External (published studies) | 🔜 Planned |
| Network Meta-Analysis (NMA) Plot | External (published studies) | 🔜 Planned |
| Benefit-Risk Plot | ADTTE, ADAE | 🔜 Planned |
| Tornado Plot | ADTTE, ADAE | 🔜 Planned |

---

### 09 · Trial Design & Patient Flow

Trial flow and exposure visualizations document the operational execution of a study, enrollment trajectories, treatment compliance, dose intensity, and patient disposition. CONSORT diagrams are required components of randomized trial publications under ICMJE reporting standards. Supported by ADSL and ADEX.

| Visualization | Primary Data | Status |
|---|---|---|
| CONSORT Diagram | ADSL | 🔜 Planned |
| Enrollment Over Time Plot | ADSL | 🔜 Planned |
| Treatment Exposure Plot | ADEX | 🔜 Planned |
| Dose Intensity Plot | ADEX | 🔜 Planned |
| Relative Dose Intensity (RDI) Plot | ADEX | 🔜 Planned |

---

### 10 · Cell Therapy / CAR-T Specific

Cellular immunotherapy trials, particularly CAR-T, produce visualization challenges that have no analog in conventional pharmacological trials: CAR-T cell expansion kinetics spanning orders of magnitude, cytokine release syndrome timelines requiring hour-level resolution, and bone marrow response assessments reflecting a distinct efficacy paradigm from RECIST-based solid tumor measurement. Partially supported by ADBM and ADRS.

| Visualization | Primary Data | Status |
|---|---|---|
| CAR-T Cell Expansion Curve | ADBM | 🔜 Planned |
| Cytokine Release Syndrome (CRS) Timeline | ADAE, ADBM | 🔜 Planned |
| Bone Marrow Response Plot | ADBM | 🔜 Planned |

---

### 11 · Radiomics & Imaging Analytics

Radiomics-based visualizations extract quantitative imaging features from radiology scans and connect them to clinical outcomes. These approaches require imaging-derived feature matrices that extend beyond standard ADaM data structures; where ADTR cannot serve as a proxy, appropriate public imaging datasets (TCIA) are referenced.

| Visualization | Primary Data | Status |
|---|---|---|
| Heatmap Overlay (Tumor Heterogeneity) | External (TCIA) | 🔜 Planned |
| Radiomics Feature Importance Plot | External (TCIA) | 🔜 Planned |
| Lesion Size Over Time Plot | ADTR | 🔜 Planned |

---

### 12 · Epidemiology & Incidence

Population-level cancer epidemiology visualizations, incidence trends, mortality-to-incidence ratios, stage distribution patterns — require registry-level datasets that are structurally incompatible with the single-trial ADaM architecture. This domain references SEER (Surveillance, Epidemiology, and End Results) and GLOBOCAN as authoritative public data sources.

| Visualization | Primary Data | Status |
|---|---|---|
| Age-Standardized Incidence Rate (ASR) Plot | External (SEER / GLOBOCAN) | 🔜 Planned |
| Cancer Incidence Trend Line | External (SEER / GLOBOCAN) | 🔜 Planned |
| Mortality-to-Incidence Ratio Plot | External (SEER / GLOBOCAN) | 🔜 Planned |
| Stage Distribution Bar Chart | ADSL / External | 🔜 Planned |
| Prevalence Pie / Donut Chart | External (SEER / GLOBOCAN) | 🔜 Planned |

---

## Synthetic Dataset: ONCVIZ-001

All visualizations in domains 01–11 are demonstrated on a single shared synthetic ADaM dataset — **ONCVIZ-001** — simulating a Phase II/III open-label randomized basket trial of a fictional oral kinase inhibitor (Vizatinib 300 mg QD) across five solid tumor histologies (NSCLC, CRC, HCC, PDAC, BRCA) at 20 investigational sites. The dataset comprises 400 virtual patients enrolled over 18 months with a data cutoff of March 5, 2026.

Generation parameters were anchored to empirical distributions from cBioPortal (n = 2,153 pooled NSCLC patients) and to published trial benchmarks from KEYNOTE-189 (Gandhi et al., *N Engl J Med* 2018), the TCGA Lung Adenocarcinoma dataset, and the erlotinib population PK model (Ling et al., *J Clin Pharmacol* 2006). All outputs are exactly reproducible from `seed = 42`.

| Domain | File | Records | Description |
|---|---|---|---|
| ADSL | `Data/ADSL.csv` | 400 | Subject-level: demographics, treatment arm, survival outcomes, biomarker status |
| ADRS | `Data/ADRS.csv` | 1,632 | RECIST 1.1 response assessments per visit |
| ADTR | `Data/ADTR.csv` | 7,082 | Sum of longest diameters (SLD) over time |
| ADAE | `Data/ADAE.csv` | 3,941 | Adverse events with MedDRA coding and CTCAE grading |
| ADLB | `Data/ADLB.csv` | 52,605 | Laboratory parameters (21 tests, 8 visits) |
| ADTTE | `Data/ADTTE.csv` | 1,098 | Time-to-event: OS, PFS, DOR, TTR with 15 subgroup variables |
| ADPK | `Data/ADPK.csv` | 9,180 | Plasma PK profiles (1-compartment model, treatment arm only) |
| ADEX | `Data/ADEX.csv` | 13,550 | Dose exposure and cycle-level modifications |
| ADBM | `Data/ADBM.csv` | 17,490 | Longitudinal biomarkers and immune cell panel |
| ADPR | `Data/ADPR.csv` | 21,478 | Patient-reported outcomes (EORTC QLQ-C30) |
| ADMUT | `Data/ADMUT.csv` | 732 | Somatic mutation calls (15 cancer genes) |
| **Total** | | **129,188** | |

For complete dataset documentation including calibration methodology, domain architecture, internal consistency validation, and full reference list, see [`Data/README.md`](Data/README.md).

---

## Repository Structure

```
Lifescience_Visualization/
│
├── README.md
├── Data/                              ← Shared synthetic ADaM dataset (ONCVIZ-001)
│   ├── ADSL.csv  ADRS.csv  ADTR.csv
│   ├── ADAE.csv  ADLB.csv  ADTTE.csv
│   ├── ADPK.csv  ADEX.csv  ADBM.csv
│   ├── ADPR.csv  ADMUT.csv
│   ├── generate_adam_synthetic.R
│   └── README.md                     ← Full dataset methodology paper
│
└── oncology/
    ├── 01_Response_Assessment/
    │   ├── Waterfall_Plot/
    │   ├── Spider_Plot/
    │   ├── Swimmer_Plot/
    │   ├── BOR_Plot/
    │   └── Tumor_Burden_Plot/
    ├── 02_Survival_Time_to_Event/     ← 11 plot types
    ├── 03_Biomarker_Genomics/         ← 15 plot types
    ├── 04_Immunology_Cellular/        ← 6 plot types
    ├── 05_Safety_Toxicity/            ← 6 plot types
    ├── 06_PK_PD/                      ← 5 plot types
    ├── 07_Imaging_Tumor_Measurement/  ← 3 plot types
    ├── 08_Meta_Analysis_Comparison/   ← 5 plot types
    ├── 09_Trial_Design_Patient_Flow/  ← 5 plot types
    ├── 10_Cell_Therapy_CART/          ← 3 plot types
    ├── 11_Radiomics_Imaging/          ← 3 plot types
    └── 12_Epidemiology_Incidence/     ← 5 plot types
```

Each visualization folder follows a uniform structure:

```
Plot_Name/
├── README.md          ← analytical rationale, design decisions, regulatory context, references
├── plot_name.py       ← standalone Python script (reads from ../../Data/)
└── output/            ← example outputs at publication resolution (300 dpi)
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

```
pandas
numpy
matplotlib
scipy
lifelines       # survival analysis
```

---

## License

Dataset released under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.  
Code released under **MIT License**.  
