"""
ONCVIZ-001 Phase 1 -- primary scoring script.

Expected input: a CSV with columns
    sample_id, model_name, input_condition, ground_truth_error_id,
    ground_truth_severity, predicted_verdict, predicted_error_id, confidence

ground_truth_error_id is "NONE" for the reference (correct) samples.
predicted_verdict is CORRECT/FLAWED as returned by the model.
predicted_error_id is the closed-list code the model returned (or OTHER/NONE).

BUG HISTORY (kept for transparency): an earlier version of this script
silently counted unparseable/NaN predicted_verdict rows as "not flawed",
artificially inflating specificity and lowering the false-negative rate.
Fixed to explicitly exclude and report unparsed rows instead. This is the
corrected version.

Note on unit of analysis: this script's default behavior treats each of
the 3 repeated runs per sample as an independent observation (matching
what the original manuscript reported as its primary metric). For
unique-sample-level metrics (majority vote across the 3 repeats), see
score_unique_sample.py in this directory, which was added after a
reviewer asked for it specifically -- the two give different numbers
(most notably for GPT-5.6, whose repeat-level inconsistency was partly
masked by pseudo-replication in the repeat-level numbers).

Reconstructed from the documented build session for this benchmark.

Usage:
    python score.py --results your_results.csv
"""

import argparse
import pandas as pd


def binary_labels(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["unparsed"] = df.predicted_verdict.isna()
    df["gt_flawed"] = df.ground_truth_error_id != "NONE"
    df["pred_flawed"] = df.predicted_verdict.str.upper() == "FLAWED"
    return df


def parsing_report(df: pd.DataFrame) -> None:
    n_unparsed = df.unparsed.sum()
    if n_unparsed:
        print(f"WARNING: {n_unparsed}/{len(df)} rows have an unparseable "
              f"predicted_verdict (blank/NaN). These are EXCLUDED from "
              f"sensitivity/specificity below, not silently counted as "
              f"correct. Inspect raw_description for these rows -- likely "
              f"the model didn't follow the fixed output format.")
        cols = [c for c in ["sample_id", "model_name", "input_condition", "repeat"] if c in df.columns]
        print(df[df.unparsed][cols].to_string(index=False))
        print()


def detection_metrics(df: pd.DataFrame) -> dict:
    tp = ((df.gt_flawed) & (df.pred_flawed)).sum()
    fn = ((df.gt_flawed) & (~df.pred_flawed)).sum()
    tn = ((~df.gt_flawed) & (~df.pred_flawed)).sum()
    fp = ((~df.gt_flawed) & (df.pred_flawed)).sum()
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    fnr = fn / (tp + fn) if (tp + fn) else float("nan")
    return {"sensitivity": sens, "specificity": spec,
            "false_negative_rate": fnr, "n": len(df)}


def category_match_accuracy(df: pd.DataFrame) -> float:
    flawed = df[df.gt_flawed]
    return (flawed.predicted_error_id == flawed.ground_truth_error_id).mean()


def critical_fnr(df: pd.DataFrame) -> float:
    crit = df[(df.ground_truth_severity == "critical")]
    if crit.empty:
        return float("nan")
    missed = crit[crit.pred_flawed == False]  # noqa: E712
    return len(missed) / len(crit)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    args = ap.parse_args()

    df = binary_labels(pd.read_csv(args.results))
    parsing_report(df)
    df = df[~df.unparsed].copy()  # exclude unparseable rows from all metrics below

    print("=== Overall (pooled across models) ===")
    print(detection_metrics(df))
    print("Category-match accuracy:", round(category_match_accuracy(df), 3))
    print("Critical-error false-negative rate:", round(critical_fnr(df), 3))

    print("\n=== Per model ===")
    for model, g in df.groupby("model_name"):
        print(f"--- {model} ---")
        print(detection_metrics(g))
        print("Category-match accuracy:", round(category_match_accuracy(g), 3))
        print("Critical FNR:", round(critical_fnr(g), 3))

    print("\n=== Per input condition ===")
    for cond, g in df.groupby("input_condition"):
        print(f"--- {cond} ---")
        print(detection_metrics(g))

    print("\n=== Pass/fail vs pre-specified calibration threshold ===")
    overall = detection_metrics(df)
    passed = (overall["sensitivity"] >= 0.90) and (critical_fnr(df) <= 0.05)
    print("PASS" if passed else "FAIL",
          "(threshold: sensitivity >= 0.90 AND critical FNR <= 0.05)")


if __name__ == "__main__":
    main()
