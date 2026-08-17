# =============================================================================
# ONCVIZ-001 Phase 1 -- Waterfall plot error injectors
# Reconstructed from the documented build session for this benchmark.
# =============================================================================

# ---- WF01: incorrect RECIST response color-coding (data_injection) --------
# Relabel a TREATMENT-arm patient's OVRLRESP row (whose best response is PR,
# no CR) from PR to SD, so the plotted color no longer matches their actual
# best percent change. NOTE: must filter to the TREATMENT arm specifically --
# an earlier version of this injector picked a CONTROL-arm patient by
# mistake, which had zero visible effect since the waterfall script only
# plots the TREATMENT arm (caught via md5 hash comparison against the
# reference image, not by inspection).
wf01_miscode_response_color <- function(adrs, adsl) {
  trt_ids <- adsl$USUBJID[adsl$ARM == "TREATMENT"]
  ovrl <- adrs[adrs$PARAMCD == "OVRLRESP" & adrs$USUBJID %in% trt_ids, ]
  candidates <- unique(ovrl$USUBJID[ovrl$AVALC == "PR"])
  # exclude patients who also have a CR row (best response would stay PR-only)
  has_cr <- unique(ovrl$USUBJID[ovrl$AVALC == "CR"])
  candidates <- setdiff(candidates, has_cr)
  target <- candidates[1]
  mask <- adrs$USUBJID == target & adrs$PARAMCD == "OVRLRESP" & adrs$AVALC == "PR"
  adrs$AVALC[mask] <- "SD"
  adrs
}

# ---- WF02: missing PD threshold reference line (code_patch) ---------------
# Remove the +20% PD threshold geom_hline call from make_waterfall_single()
# in a patched copy of waterfall_plot.R, while leaving the -30% PR line intact.
# (Diff shown as pseudocode; apply as a patched script copy.)
wf02_remove_pd_line_patch <- '
# In make_waterfall_single(), delete this line:
#   geom_hline(yintercept = PD_TH, linewidth = 0.8, color = "#444444", linetype = "dashed") +
# while keeping the PR_TH line.
'

# ---- WF03: bars not sorted by magnitude (code_patch) -----------------------
wf03_shuffle_bar_order_patch <- '
# In make_waterfall_single(), replace:
#   df <- df_in |> arrange(desc(pct)) |> mutate(x = row_number())
# with:
#   df <- df_in |> { \\(d) d[sample(nrow(d)), ] }() |> mutate(x = row_number())
'

# ---- WF04: best overall response mismatch (data_injection) ----------------
# Relabel a different TREATMENT-arm PR patient's row to CR, contradicting
# their actual tumor-shrinkage data (best response should still be PR).
wf04_mismatch_bor <- function(adrs, adsl, exclude_usubjid = NULL) {
  trt_ids <- adsl$USUBJID[adsl$ARM == "TREATMENT"]
  ovrl <- adrs[adrs$PARAMCD == "OVRLRESP" & adrs$USUBJID %in% trt_ids, ]
  has_cr <- unique(ovrl$USUBJID[ovrl$AVALC == "CR"])
  candidates <- setdiff(unique(ovrl$USUBJID[ovrl$AVALC == "PR"]), has_cr)
  candidates <- setdiff(candidates, exclude_usubjid)
  target <- candidates[1]
  mask <- adrs$USUBJID == target & adrs$PARAMCD == "OVRLRESP" & adrs$AVALC == "PR"
  adrs$AVALC[mask] <- "CR"
  adrs
}
