"""
ONCVIZ-001 Phase 1 -- unique-sample-level scoring (majority vote).

Built in response to a specific question: the primary metrics in score.py
treat each of the 3 repeated runs per sample as an independent
observation. This script instead collapses the 3 repeats per (model,
sample_id, input_condition) into one majority-vote verdict, then computes
sensitivity/specificity/critical-FNR at that level.

Majority-vote rule (chosen and justified explicitly, not arbitrarily):
simple majority on VERDICT (CORRECT/FLAWED) only -- ignore
predicted_error_id for this vote. Rationale:
  - VERDICT is binary, so with 3 valid observations a true tie is
    mathematically impossible (only 3-0 or 2-1 can occur). This makes the
    rule well-defined without needing a tie-breaker, except in the edge
    case of missing/unparsed repeats (handled below).
  - Requiring predicted_error_id agreement too (a stricter alternative)
    would conflate two different capabilities (detection vs. category
    naming) that this study deliberately treats as separate axes -- see
    semantic_scoring/semantic_score.py and the manuscript's three-axis
    framing (detection, semantic recognition, format compliance).
  - A rule that counts any single FLAWED vote among 3 as FLAWED (a
    maximally sensitive alternative) is not a "majority" rule at all --
    it introduces a one-directional bias (inflates flagged-as-FLAWED
    rate) rather than resolving disagreement neutrally.

Tie-breaking (only relevant if a repeat is missing/unparsed, reducing
effective n below 3): unresolvable ties are reported as their own
category and EXCLUDED from the metric denominators, rather than assigned
a default value -- exclusion is more transparent than an implicit
assumption.

Reconstructed from the documented build session for this benchmark
(originally built live, in direct response to a user question, and
verified against the real raw-output data before being written up here).

Usage:
    python score_unique_sample.py --results your_results.csv [--condition image_only]
"""

import argparse
import pandas as pd


def majority(x: pd.Series) -> str:
    vc = x.value_counts()
    if vc.empty:
        return None
    top = vc.iloc[0]
    n_tied_at_top = (vc == top).sum()
    return "TIE" if n_tied_at_top > 1 else vc.index[0]


def compute(df: pd.DataFrame, condition: str | None) -> None:
    df = df.copy()
    df["ground_truth_severity"] = df["ground_truth_severity"].fillna("none")
    df["gt_flawed"] = df.ground_truth_error_id != "NONE"
    df["pred_flawed"] = df.predicted_verdict.str.upper() == "FLAWED"

    if condition:
        df = df[df.input_condition == condition]

    # dropna=False is required: reference (NONE) samples have a NaN/blank
    # ground_truth_severity, and pandas' default groupby silently drops
    # rows with NaN in any grouping key -- this was a real bug caught
    # during the original analysis (it made specificity come out as NaN
    # for every model because all reference-sample rows vanished).
    grp = df.groupby(
        ["model_name", "sample_id", "input_condition", "gt_flawed", "ground_truth_severity"],
        dropna=False,
    )
    maj = grp["pred_flawed"].apply(
        lambda x: majority(x.map({True: "FLAWED", False: "CORRECT"}))
    ).reset_index()
    maj.columns = list(maj.columns[:-1]) + ["majority_verdict"]

    n_ties = (maj.majority_verdict == "TIE").sum()
    print(f"Unresolvable ties (excluded from metrics): {n_ties}")
    maj = maj[maj.majority_verdict != "TIE"].copy()
    maj["maj_flawed"] = maj.majority_verdict == "FLAWED"

    for model, g in maj.groupby("model_name"):
        tp = ((g.gt_flawed) & (g.maj_flawed)).sum()
        fn = ((g.gt_flawed) & (~g.maj_flawed)).sum()
        tn = ((~g.gt_flawed) & (~g.maj_flawed)).sum()
        fp = ((~g.gt_flawed) & (g.maj_flawed)).sum()
        sens = tp / (tp + fn) if (tp + fn) else float("nan")
        spec = tn / (tn + fp) if (tn + fp) else float("nan")
        crit = g[g.ground_truth_severity == "critical"]
        crit_fn = ((crit.gt_flawed) & (~crit.maj_flawed)).sum()
        crit_tot = (crit.gt_flawed).sum()
        crit_fnr = crit_fn / crit_tot if crit_tot else float("nan")
        print(f"{model}: sensitivity={sens:.3f} (n={tp + fn})  "
              f"specificity={spec:.3f} (n={tn + fp})  "
              f"critical_FNR={crit_fnr:.3f} (n={crit_tot})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--condition", default=None,
                     help="Restrict to one input_condition (e.g. image_only). "
                          "Default: pool across all conditions present.")
    args = ap.parse_args()
    df = pd.read_csv(args.results)
    compute(df, args.condition)


if __name__ == "__main__":
    main()
