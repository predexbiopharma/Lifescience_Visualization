# =============================================================================
# ONCVIZ-001 Phase 1 -- Swimmer plot error injectors
# =============================================================================

# ---- SW01: response symbol misaligned with timepoint (data_injection) -----
# Shift a TREATMENT-arm NSCLC responder's (CR/PR) assessment date ~200 days
# earlier than actually recorded, misaligning the plotted marker from the bar.
sw01_shift_response_date <- function(adrs, adsl, shift_days = -200) {
  nsclc_trt <- adsl$USUBJID[adsl$ARM == "TREATMENT" & adsl$TUMORTYPE == "NSCLC"]
  ovrl <- adrs[adrs$PARAMCD == "OVRLRESP" & adrs$USUBJID %in% nsclc_trt, ]
  responders <- unique(ovrl$USUBJID[ovrl$AVALC %in% c("CR", "PR")])
  target <- responders[1]
  mask <- adrs$USUBJID == target & adrs$PARAMCD == "OVRLRESP" & adrs$AVALC %in% c("CR", "PR")
  adrs$ADT[mask] <- as.character(as.Date(adrs$ADT[mask]) + shift_days)
  adrs
}

# ---- SW02: bar length inconsistent with duration data (data_injection) ----
# IMPORTANT: the swimmer script computes bar length from TRTEDT - TRTSDT
# directly, NOT from the ADSL$TRTDURD column. An earlier version of this
# injector modified TRTDURD and had zero visible effect (caught via md5 hash
# comparison against the reference image). Modify TRTEDT instead.
sw02_shift_trtedt <- function(adsl, shrink_factor = 0.4) {
  trt_nsclc <- which(adsl$ARM == "TREATMENT" & adsl$TUMORTYPE == "NSCLC")
  idx <- trt_nsclc[1]
  start <- as.Date(adsl$TRTSDT[idx])
  true_end <- as.Date(adsl$TRTEDT[idx])
  new_end <- start + (true_end - start) * shrink_factor
  adsl$TRTEDT[idx] <- as.character(new_end)
  adsl
}

# ---- SW03: missing ongoing-treatment arrow (code_patch) --------------------
sw03_drop_ongoing_markers_patch <- '
# Replace:
#   ong_df <- end_sym[end_sym$end_t == "ongoing", ]
# with:
#   ong_df <- end_sym[end_sym$end_t == "ongoing", ][0, ]
'

# ---- SW04: patient rows not consistently sorted (code_patch) --------------
sw04_randomize_row_order_patch <- '
# In sort_swimmer(), replace both:
#   ev <- df[df$eval, ]  |> arrange(bor_n, desc(trt_mo))
#   ne <- df[!df$eval, ] |> arrange(desc(trt_mo))
# with:
#   ev <- df[df$eval, ]  |> { \\(d) d[sample(nrow(d)), ] }()
#   ne <- df[!df$eval, ] |> { \\(d) d[sample(nrow(d)), ] }()
'
