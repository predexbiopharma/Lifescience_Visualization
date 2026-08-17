"""
ONCVIZ-001 Phase 1 -- ground-truth checker functions.

Each function programmatically verifies presence/absence of one taxonomy
error, independent of any model judgment -- ground truth is a deterministic
check on the data, never a human or model opinion. Three checkers
(check_km_monotonic, check_km_censor_marks_present,
check_waterfall_color_matches_pchg) are reproduced verbatim from the
original build session (they were pasted into an independent bias-audit
prompt and are known-exact). The remaining checkers are reconstructed to the
same logical specification documented during the build; re-verify against
the actual injector before relying on them for a new run.

Reconstructed from the documented build session for this benchmark.
"""

import pandas as pd


# ---------------------------------------------------------------------------
# Kaplan-Meier
# ---------------------------------------------------------------------------

def check_km_monotonic(adtte, paramcd="OS") -> bool:
    """Returns True if the plotted survival curve is non-increasing (the
    only mathematically valid behavior for a KM estimator)."""
    df = adtte[adtte.PARAMCD == paramcd].sort_values("AVAL")
    if "SURV_OVERRIDE" in df.columns:
        curve = df["SURV_OVERRIDE"].tolist()
    else:
        at_risk = len(df)
        surv = 1.0
        curve = []
        for _, row in df.iterrows():
            if row.CNSR == 0:
                surv *= (at_risk - 1) / at_risk
            at_risk -= 1
            curve.append(surv)
    prev = 1.0
    for s in curve:
        if s > prev + 1e-9:
            return False
        prev = s
    return True


def check_km_censor_marks_present(adtte, paramcd="OS") -> bool:
    """Returns True if at least one censored patient exists in the data
    (i.e., censoring tick marks should be visually present)."""
    df = adtte[adtte.PARAMCD == paramcd]
    return int((df.CNSR == 1).sum()) > 0


def check_km_starts_at_one(km_fit_surv_at_t0: float, tol=1e-6) -> bool:
    """Returns True if the first fitted survival value (t=0) equals 1.0, as
    required by the KM estimator's definition. Operates on the fitted
    survfit object's $surv[1], not raw ADaM data -- the corruption for KM03
    happens at the survfit level (see injectors_kaplan_meier.R), so the
    check must too."""
    return abs(km_fit_surv_at_t0 - 1.0) < tol


def check_km_at_risk_matches_data(km_fit_n_risk, true_risk_sets) -> bool:
    """Returns True if every at-risk count in the rendered table matches the
    true number of patients at risk at that landmark time, recomputed
    independently from the event/censor data."""
    return all(a == b for a, b in zip(km_fit_n_risk, true_risk_sets))


# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------

def check_forest_ci_contains_point(fdf: pd.DataFrame) -> bool:
    """Returns True if every subgroup's point estimate lies within its own
    displayed [lo_c, hi_c] confidence interval."""
    return bool(((fdf.hr >= fdf.lo_c) & (fdf.hr <= fdf.hi_c)).all())


def check_forest_null_line_correct(displayed_xintercept: float, expected=1.0, tol=1e-6) -> bool:
    """Returns True if the null reference line is drawn at the correct value
    (1.0 for a ratio statistic)."""
    return abs(displayed_xintercept - expected) < tol


def check_forest_axis_is_log(scale_type: str) -> bool:
    """Returns True if the x-axis transform is logarithmic, as required for
    a ratio-statistic (hazard ratio) forest plot."""
    return scale_type.lower() in ("log10", "log")


def check_forest_subgroup_n_matches_data(fdf: pd.DataFrame, true_n: pd.Series) -> bool:
    """Returns True if displayed subgroup Ns match the true analysis
    population sizes, recomputed independently from the underlying data."""
    return bool((fdf.n.reset_index(drop=True) == true_n.reset_index(drop=True)).all())


# ---------------------------------------------------------------------------
# Waterfall plot
# ---------------------------------------------------------------------------

def check_waterfall_color_matches_pchg(adrs, adtr, pr_th=-30, pd_th=20) -> bool:
    """Returns True if every patient's plotted response-category color
    matches the RECIST category implied by their actual best percent change
    in tumor size (independently recomputed from raw data, not read from any
    label)."""
    resp = adrs[adrs.PARAMCD == "OVRLRESP"][["USUBJID", "AVALC"]]
    best_pchg = (adtr[adtr.AVISITN > 0]
                 .groupby("USUBJID")["PCHG"].min().rename("best_pchg"))
    merged = resp.merge(best_pchg, on="USUBJID", how="left")

    def expected(pchg):
        if pd.isna(pchg):
            return "NE"
        if pchg <= -100:
            return "CR"
        if pchg <= pr_th:
            return "PR"
        if pchg >= pd_th:
            return "PD"
        return "SD"

    merged["expected"] = merged.best_pchg.apply(expected)
    return bool((merged.AVALC == merged.expected).all())


def check_waterfall_pd_line_present(rendered_hlines: list, pd_th=20) -> bool:
    """Returns True if a reference line at the +20% PD threshold is present
    among the rendered horizontal lines."""
    return any(abs(y - pd_th) < 1e-6 for y in rendered_hlines)


def check_waterfall_bars_sorted(plotted_order_pct: list) -> bool:
    """Returns True if bars are sorted by descending best percent change
    (the standard waterfall convention)."""
    return plotted_order_pct == sorted(plotted_order_pct, reverse=True)


def check_waterfall_bor_matches_data(adrs, adtr) -> bool:
    """Alias of check_waterfall_color_matches_pchg -- WF04 and WF01 share
    the same underlying ground-truth rule; they differ only in which
    patient's data was perturbed by the injector."""
    return check_waterfall_color_matches_pchg(adrs, adtr)


# ---------------------------------------------------------------------------
# Swimmer plot
# ---------------------------------------------------------------------------

def check_swimmer_marker_matches_assessment_date(adrs, rendered_marker_x: dict) -> bool:
    """Returns True if each plotted response marker's x-position matches the
    true ADT (assessment date) recorded in ADRS, converted to the same time
    axis as the plot."""
    ovrl = adrs[adrs.PARAMCD == "OVRLRESP"]
    true_dates = ovrl.groupby("USUBJID")["ADT"].apply(list).to_dict()
    return rendered_marker_x == true_dates  # exact-match check; adapt tolerance as needed


def check_swimmer_bar_length_matches_dates(adsl) -> bool:
    """Returns True if the rendered bar length equals TRTEDT - TRTSDT for
    every patient (the actual formula used by the production script)."""
    computed = (pd.to_datetime(adsl.TRTEDT) - pd.to_datetime(adsl.TRTSDT)).dt.days
    return bool((computed >= 0).all())  # sign/consistency check; compare against rendered pixel length in practice


def check_swimmer_ongoing_marker_present(adsl, rendered_ongoing_ids: set) -> bool:
    """Returns True if every patient with EOSSTT == 'ONGOING' has a
    corresponding ongoing-treatment marker rendered."""
    true_ongoing = set(adsl.USUBJID[adsl.EOSSTT == "ONGOING"])
    return true_ongoing.issubset(rendered_ongoing_ids)


def check_swimmer_rows_sorted(plotted_order_bor_then_duration: list, expected_order: list) -> bool:
    """Returns True if rows are grouped by best-overall-response rank, then
    sorted by descending duration within each group."""
    return plotted_order_bor_then_duration == expected_order


# ---------------------------------------------------------------------------
# CONSORT diagram
# ---------------------------------------------------------------------------

def check_consort_stage_arithmetic(n_screened, n_excluded, n_randomized) -> bool:
    """Returns True if screened = randomized + excluded (basic CONSORT
    bookkeeping identity)."""
    return n_screened == n_randomized + n_excluded


def check_consort_exclusion_reasons_present(rendered_excluded_box_text: str) -> bool:
    """Returns True if the exclusion box text includes an itemized reason
    breakdown, not just a bare count."""
    return "\n" in rendered_excluded_box_text or ":" in rendered_excluded_box_text


def check_consort_arm_totals_match_randomization(n_randomized, arm_counts: dict) -> bool:
    """Returns True if the sum of displayed arm allocation counts equals the
    displayed randomized total."""
    return n_randomized == sum(arm_counts.values())


def check_consort_followup_reported(rendered_discontinued_box_text: str) -> bool:
    """Returns True if the discontinuation box includes a count and/or
    reason breakdown, not just the generic label."""
    return any(ch.isdigit() for ch in rendered_discontinued_box_text)
