# =============================================================================
# ONCVIZ-001 Phase 1 -- Forest plot error injectors
# All 4 forest-plot errors are code_patch: the real forest_plot.R has no
# configurable "flag" hooks, so each requires a small, commented change to a
# copy of the production script.
# =============================================================================

# ---- FP01: CI does not span point estimate ---------------------------------
fp01_push_point_outside_ci_patch <- '
# Insert immediately after fdf is constructed with lo_c/hi_c columns:
fp01_idx <- which(!fdf$overall & !is.na(fdf$hi_c))[1]
fdf$hr[fp01_idx] <- fdf$hi_c[fp01_idx] * 1.6
'

# ---- FP02: wrong null reference line ---------------------------------------
fp02_wrong_null_line_patch <- '
# Replace:
#   geom_vline(xintercept = 1, linetype = "dashed", color = "#888888") +
# with:
#   geom_vline(xintercept = 2, linetype = "dashed", color = "#888888") +
# NOTE: use a visibly-wrong-but-plausible value (2), not 0, which would be
# invisible/off the log-scale axis and undermine the "wrong null line" signal.
'

# ---- FP03: asymmetric CI on log scale plotted linearly ---------------------
fp03_linear_axis_for_ratio_patch <- '
# Replace:
#   scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4), limits = c(0.08, 6)) +
# with:
#   scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 1, 2, 4), limits = c(0.08, 6)) +
'

# ---- FP04: missing/incorrect subgroup N -------------------------------------
fp04_inflate_subgroup_n_patch <- '
# Insert immediately before the fdf$label <- factor(paste0(...)) line:
fdf$n_true <- fdf$n
fdf$n <- fdf$n + 25
# then let the existing label-construction code use the inflated fdf$n
'

# -----------------------------------------------------------------------------
# KNOWN SHARED DEFECT (not part of the injected taxonomy, found via free-format
# audit and fixed post hoc): lo_c/hi_c were clipped to (0.05, 10) while the
# plot's own axis limits are (0.08, 6). Because the clip range exceeded the
# axis range, ggplot2 silently dropped out-of-range CIs entirely (as NA)
# rather than truncating them at the axis edge, for small-n subgroups whose
# true CI exceeded the clip bounds (NSCLC, HCC, PDAC, PD-L1 High in the
# reference image). Fix: match clip bounds to the axis exactly.
fp_fix_ci_clipping_patch <- '
# Replace:
#   lo_c = pmax(lo, 0.05), hi_c = pmin(hi, 10)
# with:
#   lo_c = pmax(lo, 0.08), hi_c = pmin(hi, 6)
'
