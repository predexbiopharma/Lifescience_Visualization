# Swimmer Plot

A patient-level visualization displaying individual treatment duration, longitudinal response assessments, post-treatment follow-up, and outcome status (ongoing, death, discontinuation) per RECIST 1.1 criteria, stratified by histology and best overall response, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 263 (Treatment Arm) | **Cutoff:** 31 May 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Swimmer plots are well established in oncology publications, yet existing open-source implementations share a consistent set of structural limitations. Most public examples use toy datasets of 15–36 patients with a single tumor type and no calibration to real trial benchmarks. The few multi-histology examples available draw patients from incompatible studies, making cross-domain comparison meaningless. Purpose-built packages such as `ggswim` (CHOP, 2025) and `swimplot` (PMH, 2022) solve the grammar problem elegantly, but ship with de-identified illustrative data that behaves nothing like a real oncology trial.

| | Prior art | This work |
|---|---|---|
| Dataset size | 15–36 synthetic patients (swimplot, ggswim vignettes) | 263 treatment-arm patients from ONCVIZ-001 |
| Tumor histologies | Single histology or uncalibrated mix | Five independently calibrated profiles |
| Cross-domain consistency | Swimmer and waterfall from different studies | Same 400 patients, 13 ADaM domains, 131,690 records |
| Follow-up encoding | Treatment bar only | Separate follow-up bar extending to OS date |
| Outcome symbols | Rarely distinguished | Arrow (ongoing) · ✕ (death) · \| (discontinuation) |
| Response markers | Color only, or not time-located | Shape + color markers at exact assessment month |
| Sorting variants | Duration only | Three purpose-built variants (BOR · duration · tumor type) |
| Reproducibility | No fixed seed, or seed not stated | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why the swimmer plot complements the waterfall plot

The waterfall plot answers *how much*, depth of tumor response at best assessment. The swimmer plot answers *how long* and *what happened next*, the temporal arc of each patient's journey from treatment start through follow-up or death. Used together from the same dataset they provide the complete RECIST efficacy picture that neither plot can deliver alone.

| Question | Waterfall | Swimmer |
|---|---|---|
| Depth of response (% SLD change) | ✓ | — |
| Duration of response | — | ✓ |
| Time to first response | — | ✓ |
| Post-treatment follow-up | — | ✓ |
| Outcome (ongoing / death) | — | ✓ |
| Responders vs. non-responders at a glance | ✓ | ✓ |

A key limitation noted in the oncology visualization literature is that classic swimmer plots display only patients with favorable responses (CR or PR), leaving SD and PD patients invisible. The three variants here address this explicitly — Variants A and B are responder-focused (suitable for publication), while Variant C adds histology stratification to reveal which tumor types drove the response signal.

---

## Visual anatomy

```
Patient ─┤
  ID     │
         │  [████████████████░░░░░]──────────────────────── ▶  ongoing
         │  [████████████░░░░░░░░░░░░░░] ✕                      death
         │  [████████████████████] |                             discontinuation
         │
         └──────────────────────────────────────────────────────
              0m     6m     12m    18m    24m    36m
                      Months from Treatment Start
```

| Element | Description |
|---|---|
| Colored solid bar | Time on treatment - color encodes tumor type |
| Grey translucent bar | Post-treatment follow-up to OS date (off-treatment, alive) |
| ▶ Arrow | Patient ongoing at data cutoff |
| ✕ Symbol | Death |
| \| Symbol | Discontinuation (not death) |
| Filled circle (●) | Complete Response (CR) marker at assessment month |
| Filled circle (●) | Partial Response (PR) marker at assessment month |
| Filled triangle (▲) | Stable Disease (SD) marker |
| Filled diamond (◆) | Progressive Disease (PD) marker |
| × | Not Evaluable (NE) |
| Left color strip | Tumor type - NSCLC / CRC / HCC / PDAC / BRCA |
| BOR column | Best Overall Response label at right margin |
| Milestone lines | Dashed verticals at 6, 12, 18, 24, 36 months |
| Alternating row shading | Guides the eye across long patient IDs |

---

## Plot variants

### `swimmer_A_responders.png` - Responders (CR + PR), sorted by BOR then duration

Shows all treatment-arm patients with a best overall response of CR or PR. Bars sorted first by BOR group (CR above PR), then by descending treatment duration within each group. A dashed horizontal separator and group labels distinguish the two BOR tiers. Best for Phase I/II efficacy signal figures and supplementary materials in publications.

- BOR group labels in right margin
- CR/PR separator line
- Output: 16 × dynamic height @ 300 DPI

### `swimmer_B_12mo_responders.png` - Responders with ≥ 12 months on treatment

Subset of Variant A filtered to patients who remained on treatment for at least 12 months. Highlights durable responders — a clinically meaningful cut point for regulatory submissions and dossiers. Directly addresses the question raised by reviewers: *"How many patients maintained response beyond one year?"*

- Same layout as Variant A
- Output: 16 × dynamic height @ 300 DPI

### `swimmer_C_responders_by_histology.png` - Responders stratified by tumor type

All CR + PR patients sorted first by tumor type (NSCLC → BRCA → HCC → CRC → PDAC), then by BOR and duration within each histology block. Left-margin colored rectangles and rotated tumor-type labels act as panel headers without the whitespace cost of true facets. Horizontal separator lines delineate histology groups. Best for basket trial figures comparing response depth and duration across tumor types.

- Left-margin tumor-type strip (colored rect + rotated label)
- Histology separator lines
- Output: 20 × dynamic height @ 300 DPI

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Placebo
Design   Phase II/III Open-Label Randomized Basket Trial
N        400  (TRT = 263 · CTL = 137 · 2:1 ratio)
Records  131,690 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   31 May 2026
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

| Domain | Description | Rows | Key variables used |
|--------|-------------|-----:|---|
| ADSL | Subject-level | 400 | `USUBJID`, `ARM`, `TUMORTYPE`, `TRTSDT`, `TRTEDT`, `TRTDURD`, `EOSSTT`, `DCSREAS`, `BESTRSPC`, `OSCR`, `OSDTC` |
| ADRS | Tumor response per RECIST 1.1 | 1,211 | `PARAMCD="OVRLRESP"`, `ADT`, `AVALC` (CR/PR/SD/PD/NE) |

```r
# Load required ADaM domains
ADSL <- read.csv("ADSL.csv") |>
  mutate(across(c(TRTSDT, TRTEDT, CUTDTC, OSDTC, PFSDTC), as.Date))

ADRS <- read.csv("ADRS.csv") |>
  mutate(ADT = as.Date(ADT))

# Subset to treatment arm and compute treatment duration in months
trt <- ADSL |>
  filter(ARM == "TREATMENT") |>
  mutate(TRTDUR_MO = TRTDURD / DAY2MO)   # DAY2MO = 30.4375
```

---

## Key implementation details

### Treatment bar vs. follow-up bar

Two separate bar layers encode distinct clinical phases. The treatment bar (colored, tumor-type fill) runs from day 0 to the earlier of treatment end or data cutoff. The follow-up bar (grey, narrower) runs from treatment end to the OS date when a patient discontinued but survived beyond treatment end. This separation makes it visually immediate whether a patient died on treatment, died after treatment, or is still alive.

```r
# Treatment duration (capped at cutoff for ongoing patients)
trt_end_dt <- if (ongoing) cutoff else trtedt
trt_mo     <- max(as.numeric(trt_end_dt - trtsdt) / DAY2MO, 0.2)

# Follow-up bar (only when patient died after treatment end)
fu_end_mo <- NA_real_
if (!ongoing && !is.na(s$OSDTC)) {
  osdtc <- as.Date(s$OSDTC)
  if (!is.na(osdtc) && osdtc > trtedt)
    fu_end_mo <- as.numeric(osdtc - trtsdt) / DAY2MO
}
```

### Response markers

Each ADRS record for a patient is placed as a shaped, filled marker at its exact assessment month relative to treatment start. Marker shape encodes response category (circle = CR/PR, triangle = SD, diamond = PD, × = NE), fill color encodes the same, and stroke weight is held constant for legibility at small sizes.

```r
resps <- pt_rs |>
  mutate(mo = as.numeric(ADT - trtsdt) / DAY2MO) |>
  filter(mo >= 0, mo <= trt_mo + 2) |>
  select(mo, resp = AVALC)
```

### Outcome symbol placement

End-of-bar symbols are positioned 0.5 months past the rightmost bar edge (treatment bar or follow-up bar, whichever extends further), ensuring they never overlap with bar content regardless of scale.

```r
sym_df <- df |>
  mutate(sx = ifelse(!is.na(fu_end_mo) & fu_end_mo > trt_mo,
                     fu_end_mo + 0.5, trt_mo + 0.5))
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `swimmer_plot.R` | Main R script — all three plot variants | |
| `ADSL.csv` | Subject-level dataset | 400 rows |
| `ADRS.csv` | Overall response assessments | 1,211 rows |
| `swimmer_A_responders.png` | CR + PR · sorted by BOR | 16 × dynamic · 300 DPI |
| `swimmer_B_12mo_responders.png` | CR + PR · ≥12 months on treatment | 16 × dynamic · 300 DPI |
| `swimmer_C_responders_by_histology.png` | CR + PR · stratified by tumor type | 20 × dynamic · 300 DPI |

Height is computed dynamically as `max(n_patients × 0.42 + 3, 10)` inches, capped at 28–30 inches, so the plot scales gracefully from small subsets to the full cohort.

---

## When to use

**Appropriate:**
- Displaying duration of response alongside BOR in Phase I/II efficacy readouts
- Communicating the proportion of patients who remain on treatment at landmark timepoints
- Supplementary figures in clinical publications (commonly paired with waterfall plots)
- Regulatory dossiers requiring patient-level response timelines
- Basket trial figures where histology-specific duration patterns need to be visible

**Limitations:**
- Does not show magnitude of tumor shrinkage, use waterfall plots
- Does not show continuous tumor size trajectories, use spider plots
- Becomes hard to read beyond ~100 patients per panel; consider subgroup filtering
- Sorting by duration means individual patients cannot be tracked across multiple figures

---

## Requirements

```
R >= 4.1
ggplot2 >= 3.4
dplyr >= 1.1
lubridate >= 1.9
patchwork >= 1.2   # used for multi-panel assembly in extended variants
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Hanna R. ggswim: Create swimmer plots with ggplot2. Presented at R/Medicine 2025; Children's Hospital of Philadelphia. https://github.com/CHOP-CGTInformatics/ggswim

Kaplan EL, Meier P. Nonparametric estimation from incomplete observations. *J Am Stat Assoc.* 1958;53(282):457–481.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Shalhout SZ, Miller DM. Graphical representation of survival: swimmer plots for clinical trials in oncology. *The Miller Lab.* 2020. https://www.themillerlab.io/posts/swimmer_plots/

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
