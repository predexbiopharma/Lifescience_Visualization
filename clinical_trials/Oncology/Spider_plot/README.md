# Spider Plot

A spider plot (also known as a spaghetti plot) is a patient-level visualization used in oncology clinical trials to display the change in target lesion size over time for each individual patient, allowing assessment of tumor response trajectories throughout treatment.

## When to use

- Phase I/II oncology trials
- Displaying individual tumor burden trajectories over time
- Visualizing durability of response alongside waterfall plot
- Identifying patients with delayed or sustained responses

## What this plot shows

| Element | Description |
|---|---|
| Each line | One patient's tumor burden change over time |
| Y = 0% | Baseline (no change) |
| `--` at -30% | PR threshold (≥30% reduction) |
| `--` at +20% | PD threshold (≥20% increase) |
| Values clipped at +200% | Extreme progressors capped for display |

## Evaluability criteria

A patient is included if they have:
- A valid baseline sum of target lesion diameters (`> 0`)
- At least one post-baseline tumor measurement with a valid assessment date

## Files

| File | Description |
|---|---|
| `spider_plot.py` | Study-agnostic Python script |
| `data/ADSL.csv` | Synthetic ADaM subject-level dataset (N=50) |
| `data/ADTR.csv` | Synthetic ADaM tumor measurement dataset |
| `output/spider_all.png` | All cohorts output |
| `output/spider_cohort1.png` | Cohort 1 output |
| `output/spider_cohort2.png` | Cohort 2 output |
| `output/spider_cohort3.png` | Cohort 3 output |
| `output/spider_cohort4.png` | Cohort 4 output |

## Requirements

```
pandas
numpy
matplotlib
```

## References

- Seymour L, et al. (2017). iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncology*.
- JNCI (2016). Current and Evolving Methods to Visualize Biological Data in Cancer Research.
