# Response Heatmap

A response heatmap visualizes Objective Response Rate (ORR) across two categorical dimensions simultaneously, typically dose cohort and tumor type, allowing rapid identification of which patient subgroups respond best to treatment.

## When to use

- Phase I/II oncology trials with multiple cohorts or dose levels
- Cross-tabulation of response by tumor type, biomarker status, or subgroup
- Communicating dose-response patterns to clinical and regulatory audiences
- Identifying subgroups with high or low ORR at a glance

## What this plot shows

| Element | Description |
|---|---|
| Tile color | ORR (%) — blue = low, red = high |
| Large number | ORR % for that cell |
| CR/PR/SD/PD counts | Response breakdown per cell |
| `n=` | Number of patients in that cell |
| Grey tile | No patients in that subgroup |

## Evaluability

ORR is calculated as (CR + PR) / n × 100 for all patients in each cohort × tumor type cell. Cells with zero patients are shown as grey.

## Files

| File | Description |
|---|---|
| `response_heatmap.py` | Study-agnostic Python script |
| `data/ADSL.csv` | Synthetic ADaM subject-level dataset (N=50) |
| `output/response_heatmap.png` | Example output |

## Requirements

```
pandas
numpy
matplotlib
```

## References

- Jänne PA, et al. (2015). AZD9291 in EGFR inhibitor–resistant non–small-cell lung cancer. *NEJM*.
- Seymour L, et al. (2017). iRECIST. *Lancet Oncology*.
