# =============================================================================
# ONCVIZ-001 Phase 1 -- CONSORT diagram error injectors
# =============================================================================

# ---- CN01: numbers don't sum across CONSORT stages (data_injection) -------
# Add a phantom patient row that is neither randomized (RANDFL='N') nor
# screen-failed (SFFL='N'), which increases n_screened by 1 without changing
# n_randomized or n_failed -- breaking the screened = randomized + excluded
# identity the diagram displays.
cn01_add_unaccounted_patient <- function(adrand) {
  fake <- adrand[1, ]
  fake$USUBJID <- "SYNTH-CN01"
  fake$RANDFL <- "N"
  fake$SFFL <- "N"
  rbind(adrand, fake)
}

# ---- CN02: missing exclusion reasons (code_patch) ---------------------------
cn02_drop_exclusion_reasons_patch <- '
# Replace:
#   sprintf("Excluded (n = %d)\\n%s", n_failed, reasons_txt(sf_reasons, "SFREASN"))
# with:
#   sprintf("Excluded (n = %d)", n_failed)
'

# ---- CN03: randomization N does not match arm totals (data_injection) -----
# Relabel one CONTROL-arm patient as TREATMENT in ADSL, so the displayed arm
# counts (Treatment/Control boxes) no longer sum consistently against the
# actual randomization split recorded elsewhere.
cn03_relabel_arm <- function(adsl) {
  idx <- which(adsl$ARM == "CONTROL")[1]
  adsl$ARM[idx] <- "TREATMENT"
  adsl
}

# ---- CN04: missing follow-up/loss-to-follow-up reporting (code_patch) -----
cn04_drop_discontinuation_detail_patch <- '
# Replace both:
#   bDL <- box_df(5,  disc_top - h4, 40, h4,
#                 sprintf("Discontinued intervention (n = %d)\\n%s",
#                         nrow(disc_trt), reasons_txt(disc_trt_reasons, "DCSREAS")))
#   bDR <- box_df(55, disc_top - h4, 40, h4,
#                 sprintf("Discontinued intervention (n = %d)\\n%s",
#                         nrow(disc_ctrl), reasons_txt(disc_ctrl_reasons, "DCSREAS")))
# with:
#   bDL <- box_df(5,  disc_top - h4, 40, h4, "Discontinued intervention")
#   bDR <- box_df(55, disc_top - h4, 40, h4, "Discontinued intervention")
'
