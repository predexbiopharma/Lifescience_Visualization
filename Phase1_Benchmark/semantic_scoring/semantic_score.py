"""
ONCVIZ-001 Phase 1 -- semantic (soft) scoring: does the model's free-text
description semantically identify the correct error, even if it didn't
output the exact closed-list code?

Deliberately rule-based (keyword matching), NOT another LLM call -- using
an LLM to judge "did the model understand the error" would just
reintroduce the same bias question this script exists to address (whose
semantic judgment do we trust?). A fixed, published keyword list is fully
auditable by a human or any third party instead.

BUG HISTORY (important, kept for transparency): the first version of this
scorer matched keywords regardless of polarity/negation -- a description
like "the curve IS monotonic, no issue found" (i.e. explicitly denying
the KM01 error) still counted as a match because it contained the word
"monotonic". This inflated GPT-5.6's apparent semantic accuracy from
76.8% to a wrongly-high 86.9% in one intermediate result. It was caught
by a user challenge ("are you sure?"), not by the original author, and
fixed by requiring predicted_verdict == "FLAWED" as a precondition for
any semantic match. This is the corrected version.

Reconstructed from the documented build session for this benchmark.
"""

import argparse
import re
import pandas as pd

# One entry per taxonomy error: a list of keyword patterns (regex, case-
# insensitive) where at least MIN_HITS must appear in raw_description for
# the model's free-text answer to count as a semantic match. Keywords are
# kept close to the plain-English error name, not to the closed-list code
# itself, so this doesn't just re-test string matching in disguise.
KEYWORDS = {
    "KM01": (["monoton", "increas", "goes up", "rises"], 1),
    "KM02": (["censor"], 1),
    "KM03": (["start", "100%", "begin"], 1),
    "KM04": (["at.?risk", "risk table", "number at risk"], 1),
    "FP01": (["confidence interval", r"\bci\b", "does not (contain|span|include)", "point estimate"], 1),
    "FP02": (["null", "reference line", "vertical line"], 1),
    "FP03": (["log", "linear", "scale", "axis"], 1),
    "FP04": ([r"\bn\s*=", "sample size", "subgroup"], 1),
    "WF01": (["color", "recist", "response categor"], 1),
    "WF02": (["pd threshold", r"\+?20%", "reference line", "threshold"], 1),
    "WF03": (["sort", "order", "magnitude", "not sorted", "not ordered"], 1),
    "WF04": (["best overall response", "bor", "mismatch", "inconsist"], 1),
    "SW01": (["marker", "symbol", "misalign", "timepoint", "assessment date"], 1),
    "SW02": (["bar length", "duration", "inconsist"], 1),
    "SW03": (["arrow", "ongoing"], 1),
    "SW04": (["sort", "order", "unsorted", "not ordered"], 1),
    "CN01": (["sum", "add up", "arithmetic", "don'?t match", "doesn'?t match", "inconsist"], 1),
    "CN02": (["exclusion", "reason"], 1),
    "CN03": (["arm", "randomiz", "total", "match"], 1),
    "CN04": ([r"follow.?up", "loss", "discontinu"], 1),
}


def semantic_match(row) -> bool:
    gt = row.ground_truth_error_id
    if gt == "NONE" or gt not in KEYWORDS:
        return None  # not scoreable (reference sample, or unknown code)
    # CRITICAL (the fixed bug): only count a semantic match if the model
    # actually said the plot was FLAWED. Keyword presence alone is not
    # sufficient; it must co-occur with a positive flaw verdict, or a
    # denial containing the same keyword gets wrongly credited.
    if str(row.predicted_verdict).strip().upper() != "FLAWED":
        return False
    text = str(row.raw_description).lower()
    patterns, min_hits = KEYWORDS[gt]
    hits = sum(1 for p in patterns if re.search(p, text))
    return hits >= min_hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.results)
    df = df[df.ground_truth_error_id != "NONE"].copy()
    df["semantic_match"] = df.apply(semantic_match, axis=1)
    scoreable = df[df.semantic_match.notna()]

    print("=== Semantic (keyword-based) match rate, per model ===")
    for model, g in scoreable.groupby("model_name"):
        rate = g.semantic_match.mean()
        print(f"{model}: {rate:.3f}  (n={len(g)})")

    print("\n=== Compare: exact code match vs semantic match, per model ===")
    for model, g in scoreable.groupby("model_name"):
        exact = (g.predicted_error_id == g.ground_truth_error_id).mean()
        semantic = g.semantic_match.mean()
        print(f"{model}: exact={exact:.3f}  semantic={semantic:.3f}  gap={semantic - exact:+.3f}")

    print("\n=== Per error type (pooled across models) ===")
    for eid, g in scoreable.groupby("ground_truth_error_id"):
        print(f"{eid}: semantic_match={g.semantic_match.mean():.2f}  (n={len(g)})")


if __name__ == "__main__":
    main()
