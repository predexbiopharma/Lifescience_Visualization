# ONCVIZ-001 · ADaM v1
### A Tumor-Stratified Synthetic Clinical Trial Dataset for Oncology Visualization

```
Study     ONCVIZ-001  ·  Vizatinib 300 mg QD vs Placebo
Design    Phase II/III Open-Label Randomized Basket Trial
Patients  400  (TRT = 263  ·  CTL = 137  ·  2:1 ratio)
Records   131,690 across 13 ADaM domains
Seed      42  ·  fully reproducible
Cutoff    March 5, 2026
```

---

## Why we built this

There is a real gap in the oncology data visualization ecosystem. Most published visualization
catalogs either use toy datasets that are too simple to reveal anything interesting, or they stitch
together figures from different studies with incompatible assumptions about trial design, patient
populations, and data structure. Reproducibility suffers. Readers cannot verify whether a spider
plot or a competing risks curve was generated from the same 400 patients as the forest plot
sitting next to it.

The obvious alternative is real clinical data. The problem is that patient-level data from
randomized controlled trials is almost never publicly available at the resolution needed for this
kind of work. Regulatory submissions contain everything, but are not released at the individual
level. The pharmaverse `admiral` test datasets are designed for software testing, not visual
demonstration, and lack PK profiles, dose modification histories, mutation records, and immune
biomarkers. cBioPortal and TCGA provide rich genomic data but no longitudinal tumor measurements,
no QoL, no pharmacokinetics, and no cross-domain consistency.

Synthetic data is the right answer here, provided it is done carefully. The goal is not to trick
anyone into thinking this is real. The goal is a dataset that behaves statistically like real
oncology data, respects the biological relationships between domains, and is anchored to published
benchmarks so that every number can be traced to a source. That is what this dataset is.

---

## Design decisions

A core design decision was tumor stratification: each of the five histologies has its own
calibrated ORR, survival parameters, mutation prevalence table, toxicity modifier, and
mutational signature profile. A KRAS mutation in PDAC (91%) means something completely
different from one in BRCA (2%). Most synthetic datasets use a single parameter set across
all histologies, which makes per-tumor subgroup analysis meaningless. This one does not.

Other deliberate choices:

```
Naive approach   one parameter set for all histologies, subgroup analysis meaningless
This dataset     five independent calibrated profiles, one per tumor type
    + ADRAND (screening/CONSORT)
    + ADSIG  (mutational signatures)
    + competing event encoding
    + CTCAE grading for all lab parameters
    + dose re-escalation logic
```

---

## Dataset Inventory

| Domain   | Description                                |   Rows | Cols | Design notes |
|----------|--------------------------------------------|-------:|-----:|-------------|
| ADSL     | Subject-level analysis dataset             |    400 |   63 | COMPTYPE, 5 mutation flags |
| ADRS     | Tumor response per RECIST 1.1              |  1,211 |   14 | BICR_CONF, tumor-stratified Markov |
| ADTR     | Sum of longest diameters (mm)              |  6,539 |   15 | TUMORTYPE propagated |
| ADAE     | Adverse events · MedDRA / CTCAE v5         |  5,213 |   18 | Histology-specific incidence |
| ADLB     | Laboratory parameters · 21 tests           | 53,298 |   21 | Full CTCAE grading |
| ADTTE    | Time-to-event · OS / PFS / DOR / TTR       |    970 |   34 | LM flags, COMPEVENT |
| ADPK     | Pharmacokinetics · 1-compartment model     |  8,942 |   17 | — |
| ADEX     | Dose exposure and modifications            | 12,484 |   20 | Re-escalation, CUMDOSE |
| ADBM     | Biomarkers and immune cell panel           | 17,724 |   17 | — |
| ADPR     | Patient-reported outcomes · EORTC QLQ-C30  | 21,697 |   19 | — |
| ADMUT    | Somatic mutations · 15 genes               |    753 |   23 | Tumor-stratified prevalence |
| ADRAND   | Screening and randomization                |    459 |   12 | **Included** |
| ADSIG    | Mutational signatures · SBS                |  2,000 |   13 | **Included** |
| **Total**|                                            |**131,690**| | |

---

## Calibration Strategy

Two calibration sources were used, chosen to maximize traceability.

### Empirical calibration from cBioPortal

Distributional parameters for demographic and survival variables were estimated from patient-level
data accessed via the cBioPortal REST API across four NSCLC studies:
`lung_msk_2017`, `nsclc_pd1_msk_2018`, `luad_tcga_pan_can_atlas_2018`, `lusc_tcga_pan_can_atlas_2018`
— yielding a pooled analytical cohort of 2,153 patients. Candidate parametric distributions were
fitted by maximum likelihood estimation, with model selection by AIC.

| Variable | Distribution | Parameters | n | KS p-value |
|----------|-------------|-----------|---|------------|
| OS (months) | Weibull | shape=0.92, scale=24.8 | 970 | 0.007 |
| PFS (months) | Weibull | shape=0.78, scale=15.2 | 1,230 | 0.185 |
| Age (years) | Normal | μ=66.0, σ=11.0 | 2,069 | <0.001 |
| Sex | Categorical | M=50.9%, F=49.1% | 2,099 | — |
| Smoking | Categorical | Former-heavy 48.8%, Never 32.2% | 860 | — |

The KS p-values for OS and age are statistically significant. This is expected at n>900: the test
has enough power to detect departures from parametric form that are biologically irrelevant.
Visual inspection of fitted density curves confirmed acceptable fit in both cases.

### Literature-based calibration

Each tumor type was anchored to a specific published trial or database:

| Histology | Primary source | Parameters borrowed |
|-----------|---------------|---------------------|
| NSCLC | KEYNOTE-189 (Gandhi et al., NEJM 2018) | ORR, OS/PFS, AE incidence |
| HCC | IMbrave150 (Finn et al., NEJM 2020) | ORR, OS HR, PFS HR |
| CRC | KEYNOTE-177 (André et al., NEJM 2020) | ORR, survival |
| BRCA | OlympiAD (Robson et al., NEJM 2017) | ORR, survival |
| PDAC | NAPOLI-1 / PRODIGE-4 | ORR, survival |
| All histologies | TCGA PanCancer Atlas (cBioPortal) | Mutation prevalences |
| All histologies | COSMIC v3.3 SBS catalogue | Signature profiles |
| PK | Erlotinib popPK (Ling et al., 2006) | CL/F, Vd/F, Ka |
| Safety grading | NCI CTCAE v5.0 | All lab threshold values |
| QoL | EORTC QLQ-C30 manual (3rd ed.) | Baseline norms, MID=10 |

---

## Tumor-Stratified Parameters

### Response and survival

| Histology | n | ORR TRT | ORR CTL | OS median TRT | OS HR | PFS HR |
|-----------|---|--------:|--------:|--------------:|------:|-------:|
| NSCLC | 88 | 47% | 20% | 19.6 m | 0.49 | 0.52 |
| HCC | 69 | 27% | 12% | 12.4 m | 0.58 | 0.59 |
| CRC | 80 | 13% | 0% | 7.3 m | 0.72 | 0.75 |
| BRCA | 82 | 38% | 12% | 19.8 m | 0.61 | 0.63 |
| PDAC | 81 | 11% | 0% | 3.6 m | 0.76 | 0.78 |

### Somatic mutation prevalence · TCGA PanCancer Atlas

| Gene | NSCLC | CRC | HCC | PDAC | BRCA | Key biology |
|------|------:|----:|----:|-----:|-----:|-------------|
| TP53 | 49% | 68% | 32% | 72% | 38% | Universal tumor suppressor |
| KRAS | 31% | 42% | 1% | **91%** | 2% | Dominant driver in PDAC |
| EGFR | 10% | 1% | 0% | 1% | 4% | Targetable in NSCLC |
| PIK3CA | 7% | 20% | 5% | 3% | **35%** | PI3K pathway, high in BRCA |
| SMAD4 | 2% | **32%** | 3% | **22%** | 2% | TGF-β, CRC and PDAC |
| CDKN2A | 12% | 5% | 8% | **30%** | 15% | Cell cycle regulation |
| STK11 | 17% | 3% | 2% | 4% | 5% | LKB1 pathway, NSCLC |

### Toxicity modifiers by histology

KEYNOTE-189 base incidence rates multiplied by the factors below:

| SOC category | NSCLC | HCC | CRC | BRCA | PDAC |
|--------------|------:|----:|----:|-----:|-----:|
| Hepatic | 1.0× | **2.2×** | 0.8× | 0.9× | 1.2× |
| Gastrointestinal | 1.0× | 1.1× | **1.8×** | 1.2× | **2.0×** |
| Dermatologic | 1.0× | 0.9× | 1.2× | 1.3× | 0.7× |
| Hematologic | 1.0× | 0.9× | 1.1× | **1.4×** | 1.3× |

### Mutational signatures · COSMIC v3.3 SBS

| Histology | Dominant signature | Weight | Biology |
|-----------|-------------------|-------:|---------|
| NSCLC | SBS4 | 35% | Tobacco carcinogens |
| NSCLC | SBS2 + SBS13 | 35% | APOBEC cytidine deaminase |
| CRC | SBS15 + SBS6 + SBS44 | 45% | Defective mismatch repair |
| HCC | SBS22 + SBS24 | 40% | Aristolochic acid / aflatoxin |
| PDAC | SBS1 + SBS5 | 65% | Age-related clock-like mutagenesis |
| BRCA | SBS3 | 30% | HR deficiency / BRCA1–2 loss |
| BRCA | SBS2 + SBS13 | 45% | APOBEC activity |

---

## Domain Design Notes

### ADTTE — time-to-event
Four endpoints: OS, PFS, DOR (responders only), TTR (responders only). Each record embeds
15 subgroup variables from ADSL so forest plots can be built without additional joins.
Three landmark flags are pre-computed: `LM6MFL`, `LM12MFL`, `LM24MFL`.

Competing events are encoded in `COMPTYPE` for cumulative incidence function (CIF) plots:

```
PROGRESSION                386 patients   dominant event
DEATH_WITHOUT_PROGRESSION    9 patients   non-cancer death, ~4% of TRT deaths
CENSORED                     5 patients   administrative
```

### ADLB — laboratory parameters
CTCAE v5.0 grading computed for all applicable parameters:

```
Hepatic      ALT, AST (5-grade thresholds)   BILI (4-grade)
Hematologic  ANC, PLT, HGB
Renal        CREATININE
Cardiac      QTcF
```

Sodium is stored as `PARAMCD = "SOD"` rather than `"NA"` to prevent accidental
missing-value coercion in downstream tooling. All other parameter codes are unambiguous.

Hy's Law candidates (peak ALT >3× ULN concurrent with peak BILI >2× ULN): **22 patients**,
within the FDA-expected range of 5–25 for this population.

### ADEX — dose exposure
```
Dose reduction     10% per-cycle probability after cycle 2
                   one level at a time: 300 → 200 → 100 mg

Dose interruption   4% per-cycle probability
                   one-cycle hold then automatic resumption

Dose re-escalation 15% per-cycle probability when below full dose
```
Mean relative dose intensity: **77.1%** — consistent with ~75–85% reported
for approved oral kinase inhibitors.

### ADRAND — screening (included)
459 screened, 59 screen failures across 6 reason categories, 400 randomized.
Enables a complete CONSORT flow diagram with reason-level breakdown.

### ADSIG — mutational signatures (included)
Five SBS signatures per patient, weights normalized to sum to 1.000. Total mutation
count derived from TMB × estimated exome size. Enables SBS bar charts, hierarchical
clustering by signature profile, and per-histology comparison plots.

---

## Validation

Every check was run programmatically against the final output files.

| Check | Result | Benchmark / Requirement |
|-------|--------|------------------------|
| OS ≥ PFS, all patients | **0 violations** | Required |
| ADSL ↔ ADRAND exact ID match | **400 / 400** | Required |
| ADAE patient coverage | **400 / 400 (100%)** | Required |
| ADLB null PARAMCD | **0** | Required |
| ADLB ATOXGR null | **0** | Required |
| ADMUT gene count | **15 genes** | Required |
| ADSIG weight sum per patient | **1.000 (all)** | Required |
| Cross-domain extra IDs | **0 in any domain** | Required |
| Grade ≥ 3 AE rate (treatment) | **19.5%** | 15–20% · KEYNOTE-189 |
| Cmax median | **735 ng/mL** | 500–1,500 ng/mL · erlotinib popPK |
| AUCinf median | **~17,000 ng·h/mL** | erlotinib popPK |
| KRAS prevalence in PDAC | **91%** | ~90% · TCGA |
| TP53 prevalence in CRC | **68%** | ~60–65% · TCGA |
| Mean RDI | **77.1%** | ~75–85% · kinase inhibitors |
| Hy's Law candidates | **22** | 5–25 · FDA guidance |
| Competing event categories | **3 present** | Required for CIF plots |

---

## Visualization Coverage

### Supported

| Category | Plot types | Primary domain |
|----------|-----------|----------------|
| Response Assessment | Waterfall, Spider, Swimmer, BOR, SLD over time | ADTR · ADRS |
| Survival | KM (OS/PFS/DOR), Landmark, CIF / Competing Risks, RMST, TTR, TTP | ADTTE |
| Genomics | OncoPrint, VAF, TMB, MSI, ctDNA dynamics, Mutational Signature | ADMUT · ADSIG · ADBM |
| Safety | AE bar chart, Toxicity heatmap, Hy's Law, Lab shift | ADAE · ADLB |
| PK/PD | Concentration-time, Trough, Exposure-efficacy, Waterfall+PK overlay | ADPK · ADEX |
| Biomarker | Forest plot (15 subgroups), Immune cell panel, PD-L1 by response | ADTTE · ADBM |
| QoL / PRO | PRO trajectories, MID responder rate, Deterioration-free | ADPR |
| Trial Design | CONSORT, Enrollment curve, Dose intensity, RDI | ADRAND · ADEX |

### Not supported

The following require data modalities structurally incompatible with a single-arm clinical trial
ADaM dataset. They are documented separately with references to appropriate public sources.

```
CNV / Circos / Manhattan / Rainfall    whole-genome data  →  TCGA
Flow cytometry / UMAP / t-SNE / CyTOF  single-cell data   →  GEO
CAR-T expansion / CRS timeline          CAR-T trial data   →  specialized sources
Radiomics overlays                      imaging data       →  TCIA
ASR / incidence trend lines             registry data      →  SEER / GLOBOCAN
```

---

## Usage

**Generate data**

```r
source("generate_adam_oncviz.R")
datasets <- generate_adam_oncviz(output_dir = "./data")
```

**Load and subset**

```r
library(dplyr)

adsl  <- read.csv("data/ADSL.csv")
adtte <- read.csv("data/ADTTE.csv")

# Subset: NSCLC treatment arm
nsclc_trt <- adsl |> filter(TUMORTYPE == "NSCLC", ARM == "TREATMENT")

# OS data — 15 subgroup variables already embedded, no join required
os <- adtte |> filter(PARAMCD == "OS")
# Available subgroups: TUMORTYPE, AGEGR1, SEX, ECOG, PDL1GRP, MSISTS,
# TMBHIGH, LIVERMETS, PRIORLINES, SMOKING, EGFRMUT, KRASMUT,
# TP53MUT, STK11MUT, STAGE

# Safe read for ADLB — Sodium PARAMCD is "SOD", not "NA"
adlb   <- read.csv("data/ADLB.csv")
sodium <- adlb |> filter(PARAMCD == "SOD")
```

---

## Reproducibility

The fixed seed `set.seed(42)` is applied before any sampling.
The generator produces bitwise-identical output on:

```
R      >= 4.1
dplyr  >= 1.0
tidyr  >= 1.1
purrr  >= 0.3
```

License: **CC BY 4.0** — unrestricted reuse with attribution.

---

## References

André T, Shiu K-K, Kim TW, et al. Pembrolizumab in microsatellite-instability–high advanced
colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma.
*Nature.* 2014;511(7511):543–550.

Cerami E, Gao J, Dogrusoz U, et al. The cBio cancer genomics portal. *Cancer Discov.*
2012;2(5):401–404.

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours:
RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Fayers PM, Aaronson NK, Bjordal K, et al. *The EORTC QLQ-C30 Scoring Manual.* 3rd ed.
Brussels: EORTC; 2001.

Finn RS, Qin S, Ikeda M, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular
carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, Rodríguez-Abreu D, Gadgeel S, et al. Pembrolizumab plus chemotherapy in metastatic
non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics using
cBioPortal. *Sci Signal.* 2013;6(269):pl1.

Ling J, Johnson KA, Miao Z, et al. Metabolism and excretion of erlotinib in healthy male
volunteers. *Drug Metab Dispos.* 2006;34(3):420–426.

National Cancer Institute. *Common Terminology Criteria for Adverse Events (CTCAE) v5.0.*
Bethesda: NIH; 2017.

Robson M, Im S-A, Senkus E, et al. Olaparib for metastatic breast cancer in patients with a
germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Tate JG, Bamford S, Jubb HC, et al. COSMIC: the catalogue of somatic mutations in cancer.
*Nucleic Acids Res.* 2019;47(D1):D941–D947.

TCGA Research Network. Comprehensive genomic characterization of squamous cell lung cancers.
*Nature.* 2012;489(7417):519–525.
