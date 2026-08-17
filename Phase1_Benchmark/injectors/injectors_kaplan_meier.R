# =============================================================================
# ONCVIZ-001 Phase 1 -- Kaplan-Meier error injectors
# Reconstructed from the documented build session for this benchmark.
# Two injection methods are used, matching what the real survivalcurves_plot.R
# production script required:
#   - data_injection: perturb the input ADaM CSV only, script untouched
#   - code_patch: modify a copy of the real script, because survminer's
#     print.ggsurvplot() rebuilds the plot from the underlying survfit object
#     on every call -- post-hoc edits to $plot$data or $table$data are
#     silently ignored, so the corruption must happen on km_fit itself,
#     BEFORE ggsurvplot() is called. This was discovered empirically during
#     the original build and is the reason KM01/KM03/KM04 are code_patch
#     rather than data_injection despite superficially looking like data-level
#     changes.
# =============================================================================

# ---- KM02: missing censoring marks (data_injection) -----------------------
# Set CNSR=0 for all subjects in the ADaM survival dataset so no censoring
# ticks are drawn on the curve, even though the underlying event process
# should have produced some.
km02_zero_censor_flag <- function(adtte) {
  adtte$CNSR <- 0
  adtte
}

# ---- KM01: non-monotonic survival curve (code_patch) -----------------------
# Insert immediately after km_fit <- survfit(...) and BEFORE ggsurvplot()
# is called on it, in a patched copy of survivalcurves_plot.R.
km01_swap_adjacent_survival <- function(km_fit) {
  diffs <- which(diff(km_fit$surv) != 0)
  if (length(diffs) >= 1) {
    k <- diffs[ceiling(length(diffs) / 2)]
    tmp <- km_fit$surv[k]
    km_fit$surv[k] <- km_fit$surv[k + 1]
    km_fit$surv[k + 1] <- tmp
  }
  km_fit
}

# ---- KM03: curve does not start at 1.0 (code_patch) ------------------------
km03_shift_survival_down <- function(km_fit, shift = 0.08) {
  km_fit$surv <- pmax(0, km_fit$surv - shift)
  km_fit
}

# ---- KM04: at-risk table mismatch (code_patch) ------------------------------
km04_corrupt_n_risk <- function(km_fit, seed = 4) {
  set.seed(seed)
  km_fit$n.risk <- pmax(0, km_fit$n.risk - sample(3:9, length(km_fit$n.risk), replace = TRUE))
  km_fit
}
