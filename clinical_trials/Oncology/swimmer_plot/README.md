# Swimmer Plot

A swimmer plot is a patient-level visualization used in oncology clinical trials to display treatment duration, response assessments, and study endpoints for each individual patient on a single chart.

## When to use

- Phase I/II oncology trials
- Displaying individual patient timelines (treatment start → end → follow-up)
- Showing response assessments over time (CR, PR, SD, PD)
- Communicating study endpoints (death, withdrawal, still on treatment)

## What this plot shows

| Element | Description |
|---|---|
| Blue bar | Duration on treatment |
| Grey bar | Off treatment, still on study follow-up |
| `▶` arrow | Patient still on treatment at data cutoff |
| `✕` black | Death or early termination |
| `✕` red | Withdrawal of consent / lost to follow-up |
| `●` green | Complete Response (CR) |
| `■` blue | Partial Response (PR) |
| `▲` yellow | Stable Disease (SD) |
| `◆` red | Progressive Disease (PD) |

## Files

| File | Description |
|---|---|
| `swimmer_plot.py` | Study-agnostic Python script |
| `data/ADSL.csv` | Synthetic ADaM subject-level dataset (N=50) |
| `data/ADRS.csv` | Synthetic ADaM response dataset |
| `data/ADTR.csv` | Synthetic ADaM tumor measurement dataset |
| `output/swimmer_plot_all.png` | All cohorts output |
| `output/swimmer_plot_cohort1.png` | Cohort 1 output |
| `output/swimmer_plot_cohort2.png` | Cohort 2 output |
| `output/swimmer_plot_cohort3.png` | Cohort 3 output |
| `output/swimmer_plot_cohort4.png` | Cohort 4 output |

## Requirements

```
pandas
numpy
matplotlib
```

## References

- Matange S. (2014). Swimmer Plot. SAS Graphically Speaking.
- JNCI (2016). Current and Evolving Methods to Visualize Biological Data in Cancer Research.
