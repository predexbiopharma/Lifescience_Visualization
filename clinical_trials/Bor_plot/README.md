# Best Overall Response (BOR) Bar Chart

A Best Overall Response bar chart summarizes the distribution of antitumor response categories across a trial population, expressed as the percentage of patients achieving each RECIST-defined response: Complete Response (CR), Partial Response (PR), Stable Disease (SD), or Progressive Disease (PD). It is among the most commonly reported efficacy figures in oncology clinical trial publications and regulatory submissions.

## When to use

- Phase I/II oncology trials with RECIST-based response assessment
- Communicating the efficacy profile of a new agent across dose cohorts
- Comparing response distributions between treatment arms or patient subgroups
- Complementing the waterfall plot with aggregate response statistics

## What this plot shows

| Element | Description |
|---|---|
| Stacked bar (left panel) | Proportion of patients in each response category per cohort |
| Grouped bar (right panel) | Cross-cohort comparison of each response category |
| Summary box | Overall ORR and per-category counts for the full population |
| `n=` label | Number of patients per cohort |

## Objective Response Rate (ORR)

ORR is defined as the proportion of patients achieving CR or PR:

```
ORR = (CR + PR) / N × 100
```

This is the primary efficacy endpoint in most Phase II oncology trials and a key secondary endpoint in Phase I dose-escalation studies.

## Files

| File | Description |
|---|---|
| `bor_chart.py` | Study-agnostic Python script |
| `data/ADSL.csv` | Synthetic ADaM subject-level dataset (N=50) |
| `output/bor_chart.png` | Example output |

## Requirements

```
pandas
numpy
matplotlib
```

## References

- Eisenhauer EA, et al. (2009). New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). *European Journal of Cancer*, 45(2), 228–247.
- Seymour L, et al. (2017). iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncology*, 18(3), e143–e152.
