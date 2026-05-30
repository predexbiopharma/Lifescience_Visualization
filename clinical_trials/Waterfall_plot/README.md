# Waterfall Plot

A patient-level visualization displaying best percentage change from baseline in sum of longest diameters (SLD) per RECIST 1.1 criteria, stratified by histology and best overall response, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 400 | **Cutoff:** 05 March 2026 | **Language:** R · base graphics | **License:** CC BY 4.0

---

## The gap this fills

Most published oncology visualization catalogs either rely on toy datasets too simple to reveal anything interesting, or stitch together figures from different studies with incompatible assumptions. Readers cannot verify whether a waterfall plot was generated from the same patients as the forest plot next to it.

| | Prior art | This work |
|---|---|---|
| Histology parameters | Single set — subgroup analysis meaningless | Five independently calibrated profiles |
| Cross-domain consistency | PK, QoL, mutations from different studies | 13 ADaM domains, 131,690 records, same 400 patients |
| Competing events / landmarks | Not available | Encoded in ADTTE; LM6/12/24 flags pre-computed |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, TCGA PanCancer Atlas |

---

## Visual anatomy

```
 60% |
     |
  0% |---|---|---|---|---|---|---|---|---|---  ← zero line
     |
-30% | - - - - - - - - - - - - - - - - - -   ← PR threshold (RECIST 1.1)
     |
     Bars sorted: greatest reduction (left) → greatest increase (right)
```

| Element | Description |
|---|---|
| Bar below 0% | Tumor reduction from baseline in SLD |
| Bar above 0% | Tumor increase from baseline in SLD |
| Dashed line at −30% | PR threshold: ≥30% reduction qualifies as Partial Response |
| Dashed line at +20% | PD threshold: ≥20% increase indicates Progressive Disease |
| Dark green bar | Complete Response (CR) |
| Blue bar | Partial Response (PR) |
| Yellow bar | Stable Disease (SD) |
| Red bar | Progressive Disease (PD) |
| Grey bar at 0% | Not Evaluable (NE) — no valid baseline or post-baseline measurement |
| Tumor type strip | Color-coded strip below bars — NSCLC / CRC / HCC / PDAC / BRCA |

---

## Evaluability criteria

A patient is **Evaluable** if they have:
- A valid ADTR baseline record with `AVISIT = "BASELINE"` and `AVAL > 0`
- At least one post-baseline ADTR record with a non-missing `PCHG` value

Patients failing either criterion are classified as **Not Evaluable** and displayed as a 0% grey bar. They are excluded from the ORR denominator in evaluable-population analyses.

```r
# Evaluability logic — waterfall.R
ev <- nrow(bl) > 0 &&
      !is.na(bl$AVAL[1]) &&
      suppressWarnings(as.numeric(bl$AVAL[1])) > 0 &&
      nrow(post) > 0

# Best % change = minimum PCHG across all post-baseline visits
pct <- if (length(pchg)) min(pchg) else 0
```

---

## Plot variants

### `waterfall_5panel.png` — Histology-stratified 5-panel

Five side-by-side panels (NSCLC, CRC, HCC, PDAC, BRCA), bars colored by best overall response. Because each tumor type was calibrated to a separate published trial, the panels behave differently: NSCLC and BRCA show deep responders, CRC and PDAC show predominantly shallow or absent response.

- Panel header: tumor type · N · CR/PR/SD/PD counts · ORR
- Response legend in right margin
- Tumor type color strip below each panel
- Output: 26 × 13 in @ 300 DPI

### `waterfall_treatment_arm.png` — Treatment arm, all histologies

All 263 treatment-arm patients sorted by best % change, bars colored by response. Dual legend shows response distribution (n per category) and tumor type breakdown.

- ORR annotation in top-left legend box
- Output: 30 × 13 in @ 300 DPI

### `waterfall_all_patients.png` — Control vs Treatment

All 400 patients with arms separated by a dashed vertical line. Each arm annotated with its ORR. Directly visualizes the treatment effect: the treatment arm shows a clear leftward shift in bar distribution.

- Per-arm ORR annotation boxes
- Output: 38 × 16 in @ 300 DPI

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Placebo
Design   Phase II/III Open-Label Randomized Basket Trial
N        400  (TRT = 263 · CTL = 137 · 2:1 ratio)
Records  131,690 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   March 5, 2026
```

### Tumor-stratified response parameters (treatment arm)

| Histology | n | ORR TRT | ORR CTL | OS median | Calibration source |
|-----------|---|--------:|--------:|----------:|---|
| NSCLC | 88 | 47% | 20% | 19.6 m | KEYNOTE-189 (Gandhi et al., NEJM 2018) |
| HCC | 69 | 27% | 12% | 12.4 m | IMbrave150 (Finn et al., NEJM 2020) |
| CRC | 80 | 13% | 0% | 7.3 m | KEYNOTE-177 (André et al., NEJM 2020) |
| BRCA | 82 | 38% | 12% | 19.8 m | OlympiAD (Robson et al., NEJM 2017) |
| PDAC | 81 | 11% | 0% | 3.6 m | NAPOLI-1 / PRODIGE-4 |

### Domains used by this script

| Domain | Description | Rows | Key variables |
|--------|-------------|-----:|---|
| ADSL | Subject-level | 400 | `USUBJID`, `ARM`, `TUMORTYPE` |
| ADRS | Tumor response per RECIST 1.1 | 1,211 | `PARAMCD="OVRLRESP"`, `AVALC` (CR/PR/SD/PD/NE) |
| ADTR | Sum of longest diameters (mm) | 6,539 | `AVISIT`, `AVAL` (SLD mm), `PCHG` (% change) |

```r
# Load the three required ADaM domains
adsl <- read.csv("ADSL.csv", stringsAsFactors=FALSE)
adrs <- read.csv("ADRS.csv", stringsAsFactors=FALSE)
adtr <- read.csv("ADTR.csv", stringsAsFactors=FALSE)

# Subset to treatment arm
df_trt <- df_all[df_all$arm == "TREATMENT", ]
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `waterfall.R` | Main R script — all three plot variants | |
| `ADSL.csv` | Subject-level dataset | 400 rows |
| `ADRS.csv` | Overall response assessments | 1,211 rows |
| `ADTR.csv` | Tumor measurements (SLD) | 6,539 rows |
| `waterfall_5panel.png` | Histology-stratified 5-panel | 26×13 in · 300 DPI |
| `waterfall_treatment_arm.png` | Treatment arm · all histologies | 30×13 in · 300 DPI |
| `waterfall_all_patients.png` | All patients · control vs treatment | 38×16 in · 300 DPI |

---

## When to use

**Appropriate:**
- Phase I/II efficacy signal reporting
- Communicating RECIST response distribution across a cohort
- Comparing histology-specific response profiles in basket trials
- Supplementary figures in clinical publications

**Limitations:**
- Does not show duration of response — use swimmer plots
- Does not show time to response — use TTR Kaplan-Meier curves
- Sorting by % change makes individual patient trajectories untrackable
- Not suitable for non-solid tumor endpoints (e.g. hematologic malignancies)

---

## Requirements

```
R >= 4.1
```

Base R graphics only — no external packages required.

---

## References

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Cerami E, et al. The cBio cancer genomics portal. *Cancer Discov.* 2012;2(5):401–404.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Tate JG, et al. COSMIC: the catalogue of somatic mutations in cancer. *Nucleic Acids Res.* 2019;47(D1):D941–D947.
