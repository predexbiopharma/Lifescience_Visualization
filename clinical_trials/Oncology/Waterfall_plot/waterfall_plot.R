library(grDevices)
library(graphics)

ADSL_PATH <- "ADSL.csv"
ADRS_PATH <- "ADRS.csv"
ADTR_PATH <- "ADTR.csv"
OUTPUT_DIR <- ""
out <- function(f) if (OUTPUT_DIR == "") f else file.path(OUTPUT_DIR, f)

RESP_COLORS <- c(CR="#1a6b3c", PR="#2b83ba", SD="#fdae61", PD="#d7191c", NE="#7f7f7f")
TUMOR_COLORS <- c(NSCLC="#2b83ba", CRC="#1a9641", HCC="#e6007e", PDAC="#d95f02", BRCA="#7b2d8b")
TUMOR_ORDER <- c("NSCLC","CRC","HCC","PDAC","BRCA")
PR_TH <- -30
PD_TH <- 20
DPI <- 300
CUTOFF <- "Data cutoff: 05 Mar 2026"

adsl <- read.csv(ADSL_PATH, stringsAsFactors=FALSE)
adrs <- read.csv(ADRS_PATH, stringsAsFactors=FALSE)
adtr <- read.csv(ADTR_PATH, stringsAsFactors=FALSE)

build_df <- function(adsl_sub) {
  rows <- vector("list", nrow(adsl_sub))
  rank <- c(CR=1, PR=2, SD=3, PD=4)

  for (i in seq_len(nrow(adsl_sub))) {
    s <- adsl_sub[i, ]
    uid <- s$USUBJID
    tr <- adtr[adtr$USUBJID == uid, ]
    bl <- tr[tr$AVISIT == "BASELINE", ]
    post <- tr[tr$AVISIT != "BASELINE", ]

    ev <- nrow(bl) > 0 &&
      !is.na(bl$AVAL[1]) &&
      suppressWarnings(as.numeric(bl$AVAL[1])) > 0 &&
      nrow(post) > 0

    pct <- 0
    if (ev) {
      pchg <- suppressWarnings(as.numeric(post$PCHG))
      pchg <- pchg[!is.na(pchg)]
      pct <- if (length(pchg)) min(pchg) else 0
    }

    rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP", ]
    resp <- "NE"
    if (nrow(rs) > 0) {
      found <- rs$AVALC[rs$AVALC %in% names(rank)]
      if (length(found)) resp <- found[which.min(rank[found])]
    }

    rows[[i]] <- data.frame(
      uid=uid,
      arm=s$ARM,
      tumor=s$TUMORTYPE,
      pct=pct,
      resp=resp,
      ev=ev,
      stringsAsFactors=FALSE
    )
  }

  do.call(rbind, rows)
}

df_all <- build_df(adsl)
df_trt <- df_all[df_all$arm == "TREATMENT", ]

orr_vals <- function(sub) {
  n <- nrow(sub)
  norr <- sum(sub$resp %in% c("CR","PR"))
  pct <- if (n) norr / n * 100 else 0
  c(n=n, norr=norr, pct=pct)
}

style_ax <- function(n_bars, ymin, ylim_bottom, ylim_top, tick_max, ylabel=TRUE) {
  plot(
    NA, NA,
    xlim=c(-0.5, n_bars + 0.5),
    ylim=c(ylim_bottom, ylim_top),
    xaxt="n", yaxt="n", xlab="", ylab="", bty="n"
  )
  ytk <- seq(floor(ymin / 20) * 20, tick_max, by=20)
  axis(2, at=ytk, labels=paste0(ytk, "%"), las=1, cex.axis=0.9)
  abline(h=ytk, col="#eeeeee", lwd=0.5)
  axis(2, at=ytk, labels=FALSE, col="#444444", lwd=0.8)
  if (ylabel) {
    mtext("Best % Change from Baseline in SLD (%)", side=2, line=4.6, font=2, cex=1.0)
  }
}

reflines <- function(n_bars, thr_labels=TRUE, fs=0.75) {
  abline(h=0, col="#111111", lwd=1.5)
  abline(h=PD_TH, col="#777777", lwd=1.0, lty=2)
  abline(h=PR_TH, col="#777777", lwd=1.0, lty=2)
  if (thr_labels) {
    text(n_bars + 1.0, PD_TH + 1.8, "+20% PD", cex=fs, col="#666666", adj=c(0,0), font=3, xpd=NA)
    text(n_bars + 1.0, PR_TH - 1.8, "\u221230% PR", cex=fs, col="#666666", adj=c(0,1), font=3, xpd=NA)
  }
}

draw_bars <- function(sub, xp, bw, color_by="resp") {
  for (i in seq_len(nrow(sub))) {
    row <- sub[i, ]
    col <- if (color_by == "resp") RESP_COLORS[row$resp] else TUMOR_COLORS[row$tumor]
    if (is.na(col)) col <- "#888888"
    rect(xp[i] - bw / 2, min(0, row$pct), xp[i] + bw / 2, max(0, row$pct), col=col, border=NA)
  }
}

tumor_strip <- function(sub, xp, bw, label=TRUE) {
  n <- nrow(sub)
  plot(NA, NA, xlim=c(-0.5, n + 0.5), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  for (i in seq_len(n)) {
    col <- TUMOR_COLORS[sub$tumor[i]]
    if (is.na(col)) col <- "#888888"
    rect(xp[i] - bw / 2, 0.05, xp[i] + bw / 2, 0.95, col=col, border=NA)
  }
  if (label) {
    text(-0.8, 0.5, "Tumor\ntype", adj=c(1,0.5), cex=0.75, font=2, col="#444444", xpd=NA)
  }
}

panel_header <- function(tumor, stats_txt) {
  usr <- par("usr")
  x <- usr[1]
  y_rng <- diff(usr[3:4])
  text(x, usr[4] + 0.075 * y_rng, tumor, adj=c(0,0), cex=1.15, font=2, col=TUMOR_COLORS[tumor], xpd=NA)
  text(
    x, usr[4] + 0.015 * y_rng, stats_txt,
    adj=c(0,0), cex=0.50, font=2, col="#333333", xpd=NA,
    bg="#f7f7f7"
  )
}

cat("Data loaded. Building plots...\n")

png(out("waterfall_5panel.png"), width=26, height=13, units="in", res=DPI, bg="white")
layout(
  matrix(c(1,2,3,4,5,6,
           7,8,9,10,11,12), nrow=2, byrow=TRUE),
  widths=c(1,1,1,1,1,0.72),
  heights=c(8.5,0.65)
)
par(oma=c(0.5,4.8,5.8,0.5), mar=c(0.4,0.4,1.5,0.4), xaxs="i", yaxs="i")

for (tumor in TUMOR_ORDER) {
  sub <- df_trt[df_trt$tumor == tumor, ]
  sub <- sub[order(sub$pct, decreasing=TRUE), ]
  rownames(sub) <- NULL
  n <- nrow(sub)
  bw <- max(0.55, min(0.85, 38 / n))
  xp <- seq_len(n) - 0.5

  style_ax(n, ymin=-110, ylim_bottom=-116, ylim_top=130, tick_max=60, ylabel=(tumor == TUMOR_ORDER[1]))
  reflines(n, thr_labels=(tumor == tail(TUMOR_ORDER, 1)), fs=0.62)
  draw_bars(sub, xp, bw, color_by="resp")

  vals <- orr_vals(sub)
  rc <- table(factor(sub$resp, levels=c("CR","PR","SD","PD")))
  stats_txt <- sprintf(
    "N = %d   CR: %d  PR: %d  SD: %d  PD: %d   ORR = %.0f%%",
    vals["n"], rc["CR"], rc["PR"], rc["SD"], rc["PD"], vals["pct"]
  )
  panel_header(tumor, stats_txt)
}

par(mar=c(0.4,0.1,1.5,0.1))
plot.new()
legend(
  "left",
  inset=c(0.02, 0.18),
  legend=c("Complete Response (CR)", "Partial Response (PR)", "Stable Disease (SD)", "Progressive Disease (PD)"),
  fill=RESP_COLORS[c("CR","PR","SD","PD")],
  border="#555555",
  title="Best Overall\nResponse",
  cex=0.82,
  bty="o",
  box.col="#aaaaaa",
  xpd=NA
)

par(mar=c(0.05,0.4,0.05,0.4))
for (tumor in TUMOR_ORDER) {
  sub <- df_trt[df_trt$tumor == tumor, ]
  sub <- sub[order(sub$pct, decreasing=TRUE), ]
  rownames(sub) <- NULL
  n <- nrow(sub)
  bw <- max(0.55, min(0.85, 38 / n))
  xp <- seq_len(n) - 0.5
  tumor_strip(sub, xp, bw, label=(tumor == TUMOR_ORDER[1]))
}
plot.new()

mtext("Waterfall Plot \u2014 Best Overall Response by Histology   \u00b7   ONCVIZ-001  \u00b7  Treatment Arm",
      outer=TRUE, side=3, line=3.2, cex=1.15, font=2)
mtext(paste0("RECIST 1.1  \u00b7  Best % change from baseline in sum of longest diameters  \u00b7  ", CUTOFF),
      outer=TRUE, side=3, line=1.8, cex=0.85, font=3, col="#555555")
dev.off()
cat("\u2713 waterfall_5panel.png\n")

df_s <- df_trt[order(df_trt$pct, decreasing=TRUE), ]
rownames(df_s) <- NULL
N2 <- nrow(df_s)
bw2 <- 0.72
xp2 <- seq_len(N2) - 0.5

png(out("waterfall_treatment_arm.png"), width=30, height=13, units="in", res=DPI, bg="white")
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(1,0.13), heights=c(8,0.65))
par(oma=c(0.5,4.8,5.8,0.5), mar=c(0.4,0.4,1.5,0.4), xaxs="i", yaxs="i")

style_ax(N2, ymin=-110, ylim_bottom=-116, ylim_top=80, tick_max=60, ylabel=TRUE)
reflines(N2, thr_labels=FALSE)
text(N2 + 1.5, PD_TH + 1.8, "+20% PD", cex=0.85, col="#666666", adj=c(0,0), font=3, xpd=NA)
text(N2 + 1.5, PR_TH - 1.8, "\u221230% PR", cex=0.85, col="#666666", adj=c(0,1), font=3, xpd=NA)
draw_bars(df_s, xp2, bw2, color_by="resp")
vals2 <- orr_vals(df_s)
legend(
  "topleft",
  legend=sprintf("N = %d   ORR = %d/%d (%.0f%%)", vals2["n"], vals2["norr"], vals2["n"], vals2["pct"]),
  bty="o",
  box.col="#999999",
  bg="white",
  text.font=2,
  cex=0.95
)

plot.new()
rc2 <- table(factor(df_s$resp, levels=c("CR","PR","SD","PD")))
legend(
  "topright",
  inset=c(0.0, 0.20),
  legend=sprintf("%s   n = %d", names(rc2), as.integer(rc2)),
  fill=RESP_COLORS[names(rc2)],
  border="#555555",
  title="Response",
  cex=0.85,
  bty="o",
  box.col="#aaaaaa"
)
legend(
  "bottomright",
  inset=c(0.0, 0.16),
  legend=TUMOR_ORDER,
  fill=TUMOR_COLORS[TUMOR_ORDER],
  border=NA,
  title="Tumor Type",
  cex=0.85,
  bty="o",
  box.col="#aaaaaa"
)

par(mar=c(0.05,0.4,0.05,0.4))
tumor_strip(df_s, xp2, bw2, label=TRUE)
plot.new()

mtext(sprintf("Waterfall Plot \u2014 Treatment Arm   \u00b7   ONCVIZ-001   (N = %d)", N2),
      outer=TRUE, side=3, line=3.2, cex=1.2, font=2)
mtext(paste0("RECIST 1.1  \u00b7  Best % change from baseline in SLD  \u00b7  ", CUTOFF),
      outer=TRUE, side=3, line=1.7, cex=0.9, font=3, col="#555555")
dev.off()
cat("\u2713 waterfall_treatment_arm.png\n")

df_a <- df_all[order(df_all$arm, -df_all$pct), ]
rownames(df_a) <- NULL
NA_all <- nrow(df_a)
bw3 <- 0.68
xp3 <- seq_len(NA_all) - 0.5
n_ctrl <- sum(df_a$arm == "CONTROL")
sep_x <- n_ctrl + 0.5

png(out("waterfall_all_patients.png"), width=38, height=16, units="in", res=DPI, bg="white")
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(1,0.10), heights=c(8,0.65))
par(oma=c(0.5,4.8,5.8,0.5), mar=c(0.4,0.4,1.5,0.4), xaxs="i", yaxs="i")

style_ax(NA_all, ymin=-110, ylim_bottom=-116, ylim_top=80, tick_max=60, ylabel=TRUE)
reflines(NA_all, thr_labels=FALSE)
text(NA_all + 1.2, PD_TH + 1.8, "+20% PD", cex=0.85, col="#666666", adj=c(0,0), font=3, xpd=NA)
text(NA_all + 1.2, PR_TH - 1.8, "\u221230% PR", cex=0.85, col="#666666", adj=c(0,1), font=3, xpd=NA)
draw_bars(df_a, xp3, bw3, color_by="resp")
abline(v=sep_x, col="#444444", lwd=1.4, lty=3)
text(sep_x / 2, -11, sprintf("Control  (n = %d)", n_ctrl), adj=c(0.5,1), cex=0.9, font=2, col="#444444")
text((sep_x + NA_all) / 2, -11, sprintf("Treatment  (n = %d)", NA_all - n_ctrl), adj=c(0.5,1), cex=0.9, font=2, col="#444444")

vals_ctl <- orr_vals(df_a[df_a$arm == "CONTROL", ])
vals_trt <- orr_vals(df_a[df_a$arm == "TREATMENT", ])
legend(sep_x / 2, 80, legend=sprintf("ORR = %.0f%%", vals_ctl["pct"]), bty="o", box.col="#999999", bg="white", text.font=2, cex=0.85, xjust=0.5, yjust=1)
legend((sep_x + NA_all) / 2, 80, legend=sprintf("ORR = %.0f%%", vals_trt["pct"]), bty="o", box.col="#999999", bg="white", text.font=2, cex=0.85, xjust=0.5, yjust=1)

plot.new()
rc3 <- table(factor(df_a$resp, levels=c("CR","PR","SD","PD")))
legend(
  "topright",
  inset=c(0.0, 0.20),
  legend=sprintf("%s   n = %d", names(rc3), as.integer(rc3)),
  fill=RESP_COLORS[names(rc3)],
  border="#555555",
  title="Response",
  cex=0.85,
  bty="o",
  box.col="#aaaaaa"
)
legend(
  "bottomright",
  inset=c(0.0, 0.16),
  legend=TUMOR_ORDER,
  fill=TUMOR_COLORS[TUMOR_ORDER],
  border=NA,
  title="Tumor Type",
  cex=0.85,
  bty="o",
  box.col="#aaaaaa"
)

par(mar=c(0.05,0.4,0.05,0.4))
tumor_strip(df_a, xp3, bw3, label=TRUE)
plot.new()

mtext(sprintf("Waterfall Plot \u2014 All Patients   \u00b7   ONCVIZ-001   (N = %d)", NA_all),
      outer=TRUE, side=3, line=3.2, cex=1.2, font=2)
mtext(paste0("RECIST 1.1  \u00b7  Best % change from baseline in SLD  \u00b7  ", CUTOFF),
      outer=TRUE, side=3, line=1.7, cex=0.9, font=3, col="#555555")
dev.off()
cat("\n\u2713 All 3 waterfall plots complete.\n")
