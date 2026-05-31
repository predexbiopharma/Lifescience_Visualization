library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)

CUTOFF  <- as.Date("2026-05-31")
DAY2MO  <- 30.4375

TUMOR_COLORS <- c(NSCLC="#2166AC", BRCA="#C0392B", HCC="#27AE60",
                  CRC="#8B4513",   PDAC="#762A83")
RESP_COLORS  <- c(CR="#1A9641", PR="#4BACC6", SD="#F6A623",
                  PD="#D73027", NE="#BDBDBD")
RESP_SHAPES  <- c(CR=21, PR=21, SD=24, PD=23, NE=4)
RESP_SIZES   <- c(CR=3.5, PR=3.0, SD=2.7, PD=2.7, NE=2.0)
MILESTONE_MO <- c(6, 12, 18, 24, 36)

# ── 1. Data ───────────────────────────────────────────────────
ADSL <- read.csv("ADSL.csv") |>
  mutate(across(c(TRTSDT, TRTEDT, CUTDTC, OSDTC, PFSDTC), as.Date))

ADRS <- read.csv("ADRS.csv") |>
  mutate(ADT = as.Date(ADT))

trt <- ADSL |>
  filter(ARM == "TREATMENT") |>
  mutate(TRTDUR_MO = TRTDURD / DAY2MO)

# ── 2. Build rows ─────────────────────────────────────────────
build_patient_rows <- function(adsl_sub, adrs, cutoff = CUTOFF,
                               sort_by_bor = FALSE, tumor_order = NULL) {
  adrs_resp <- adrs |> filter(PARAMCD == "OVRLRESP")

  rows <- lapply(seq_len(nrow(adsl_sub)), function(i) {
    s       <- adsl_sub[i, ]
    uid     <- s$USUBJID
    trtsdt  <- as.Date(s$TRTSDT)
    trtedt  <- as.Date(s$TRTEDT)
    ongoing <- s$EOSSTT == "ONGOING"
    dcsreas <- trimws(as.character(s$DCSREAS))
    bor     <- trimws(as.character(ifelse(is.na(s$BESTRSPC), "NE", s$BESTRSPC)))
    tumor   <- trimws(as.character(s$TUMORTYPE))
    oscr    <- as.integer(s$OSCR)

    trt_end_dt <- if (ongoing) cutoff else trtedt
    trt_mo     <- max(as.numeric(trt_end_dt - trtsdt) / DAY2MO, 0.2)

    fu_end_mo <- NA_real_
    if (!ongoing && !is.na(s$OSDTC)) {
      osdtc <- as.Date(s$OSDTC)
      if (!is.na(osdtc) && osdtc > trtedt)
        fu_end_mo <- as.numeric(osdtc - trtsdt) / DAY2MO
    }

    pt_rs <- adrs_resp |> filter(USUBJID == uid) |> arrange(ADT)
    resps <- pt_rs |>
      mutate(mo = as.numeric(ADT - trtsdt) / DAY2MO) |>
      filter(mo >= 0, mo <= trt_mo + 2) |>
      select(mo, resp = AVALC)

    list(uid = uid, tumor = tumor, bor = bor, ongoing = ongoing,
         dcsreas = dcsreas, oscr = oscr,
         trt_mo = trt_mo, fu_end_mo = fu_end_mo, resps = resps)
  })

  df <- bind_rows(lapply(rows, function(r) {
    data.frame(uid = r$uid, tumor = r$tumor, bor = r$bor,
               ongoing = r$ongoing, dcsreas = r$dcsreas, oscr = r$oscr,
               trt_mo = r$trt_mo, fu_end_mo = r$fu_end_mo,
               stringsAsFactors = FALSE)
  }))

  resp_df <- bind_rows(lapply(rows, function(r) {
    if (nrow(r$resps) == 0) return(NULL)
    r$resps |> mutate(uid = r$uid)
  }))

  df <- df |> mutate(
    bor_order = case_when(bor == "CR" ~ 0, bor == "PR" ~ 1, TRUE ~ 2)
  )

  if (sort_by_bor) {
    df <- df |> arrange(bor_order, desc(trt_mo))
  } else if (!is.null(tumor_order)) {
    df <- df |>
      mutate(tumor_order = match(tumor, tumor_order)) |>
      arrange(tumor_order, bor_order, desc(trt_mo))
  } else {
    df <- df |> arrange(desc(trt_mo))
  }

  df <- df |> mutate(y = rev(seq_len(n())))
  list(df = df, resp_df = if (!is.null(resp_df)) resp_df else data.frame())
}

# ── 3. Core plot ──────────────────────────────────────────────
draw_swimmer <- function(df, resp_df, title_str,
                         show_bor_labels = TRUE,
                         tumor_order = NULL) {
  n      <- nrow(df)
  max_mo <- max(c(df$trt_mo, df$fu_end_mo), na.rm = TRUE)
  max_mo <- max(max_mo + 2, 14)
  bor_x  <- max_mo + 3.5
  top_y  <- n + 1.5

  df$uid_short <- sub(".*-", "", df$uid)
  df$uid_short <- factor(df$uid_short, levels = df$uid_short[order(df$y)])

  p <- ggplot() +
    # Alternating rows
    geom_rect(data = df |> filter((y %% 2) == 0),
              aes(xmin = -1.5, xmax = bor_x + 2,
                  ymin = y - 0.5, ymax = y + 0.5),
              fill = "#F5F7FA", alpha = 0.6) +
    # Milestone lines
    geom_vline(xintercept = MILESTONE_MO[MILESTONE_MO <= max_mo],
               linetype = "dashed", color = "#AAAAAA",
               linewidth = 0.7, alpha = 0.65) +
    annotate("text",
             x = MILESTONE_MO[MILESTONE_MO <= max_mo],
             y = top_y + 0.3,
             label = paste0(MILESTONE_MO[MILESTONE_MO <= max_mo], "m"),
             size = 2.8, color = "#555555",
             fontface = "bold.italic", vjust = 0) +
    # Tumor color strip
    geom_rect(data = df,
              aes(xmin = -0.15, xmax = 0,
                  ymin = y - 0.3, ymax = y + 0.3,
                  fill = tumor),
              alpha = 0.95) +
    scale_fill_manual(values = TUMOR_COLORS, guide = "none") +
    # Treatment bars
    geom_rect(data = df,
              aes(xmin = 0.05, xmax = trt_mo,
                  ymin = y - 0.28, ymax = y + 0.28,
                  color = tumor),
              fill = NA, linewidth = 0) +
    {
      tc_df <- df |> mutate(fill_col = TUMOR_COLORS[tumor])
      geom_rect(data = tc_df,
                aes(xmin = 0.05, xmax = trt_mo,
                    ymin = y - 0.28, ymax = y + 0.28),
                fill = tc_df$fill_col, alpha = 0.88)
    } +
    # Follow-up bars
    {
      fu_df <- df |> filter(!is.na(fu_end_mo), fu_end_mo > trt_mo)
      if (nrow(fu_df) > 0)
        geom_rect(data = fu_df,
                  aes(xmin = trt_mo, xmax = fu_end_mo,
                      ymin = y - 0.14, ymax = y + 0.14),
                  fill = "#BBBBBB", alpha = 0.45)
      else geom_blank()
    }

  # Response markers
  if (!is.null(resp_df) && nrow(resp_df) > 0) {
    rdf <- resp_df |>
      inner_join(df |> select(uid, y), by = "uid") |>
      filter(resp %in% names(RESP_COLORS)) |>
      mutate(fill_col  = RESP_COLORS[resp],
             shape_val = RESP_SHAPES[resp],
             size_val  = RESP_SIZES[resp])
    p <- p +
      geom_point(data = rdf,
                 aes(x = mo, y = y),
                 shape = rdf$shape_val,
                 fill  = rdf$fill_col,
                 color = "black",
                 size  = rdf$size_val,
                 stroke = 0.4)
  }

  # End-of-bar symbols
  sym_df <- df |>
    mutate(sx = ifelse(!is.na(fu_end_mo) & fu_end_mo > trt_mo,
                       fu_end_mo + 0.5, trt_mo + 0.5))

  p <- p +
    # Ongoing arrows
    geom_segment(data = sym_df |> filter(ongoing),
                 aes(x = sx, xend = sx + 1.6, y = y, yend = y,
                     color = tumor),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 linewidth = 0.9) +
    scale_color_manual(values = TUMOR_COLORS, guide = "none") +
    # Death
    geom_point(data = sym_df |> filter(!ongoing, (oscr == 1 | dcsreas == "Death")),
               aes(x = sx, y = y), shape = 4, size = 2.5,
               color = "#111111", stroke = 1.5) +
    # Discontinued
    geom_point(data = sym_df |> filter(!ongoing, oscr != 1, dcsreas != "Death"),
               aes(x = sx, y = y), shape = 124, size = 2.5,
               color = "#555555", stroke = 1.2) +
    # BOR separator
    geom_vline(xintercept = max_mo + 2, color = "#DDDDDD",
               linewidth = 0.6, linetype = "dotted") +
    # BOR header
    annotate("text", x = bor_x, y = top_y + 0.3,
             label = "BOR", fontface = "bold", size = 3.5,
             color = "#333333", vjust = 0) +
    # BOR values
    geom_text(data = df,
              aes(x = bor_x, y = y, label = bor,
                  color = bor),
              fontface = "bold", size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(TUMOR_COLORS, RESP_COLORS), guide = "none")

  # BOR group labels
  if (show_bor_labels) {
    for (bv in c("CR", "PR")) {
      grp <- df |> filter(bor == bv)
      if (nrow(grp) == 0) next
      mid_y <- (max(grp$y) + min(grp$y)) / 2
      p <- p + annotate("label", x = bor_x + 1.5, y = mid_y,
                        label = bv, size = 2.8, fontface = "bold",
                        color = RESP_COLORS[bv],
                        fill = "white", label.size = 0.4,
                        label.padding = unit(0.2, "lines"))
    }
    cr_g <- df |> filter(bor == "CR")
    pr_g <- df |> filter(bor == "PR")
    if (nrow(cr_g) > 0 && nrow(pr_g) > 0) {
      sep_y <- (min(cr_g$y) + max(pr_g$y)) / 2
      p <- p + geom_hline(yintercept = sep_y,
                          color = "#AAAAAA", linewidth = 0.6,
                          linetype = "dashed", alpha = 0.6)
    }
  }

  p <- p +
    scale_x_continuous(
      name   = "Months from Treatment Start",
      breaks = seq(0, max_mo + 4, 6),
      limits = c(-1.5, bor_x + 2.5)
    ) +
    scale_y_continuous(
      breaks = df$y,
      labels = df$uid_short,
      limits = c(-0.5, top_y + 1)
    ) +
    labs(title = title_str) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background  = element_rect(fill = "#FAFAFA", color = NA),
      panel.background = element_rect(fill = "#FAFAFA", color = NA),
      panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y   = element_text(size = 7, family = "mono"),
      axis.text.x   = element_text(size = 8),
      axis.title.x  = element_text(size = 10, face = "bold", margin = margin(t = 8)),
      axis.title.y  = element_blank(),
      plot.title    = element_text(size = 12, face = "bold",
                                   color = "#1A1A2E", margin = margin(b = 14)),
      legend.position = "none"
    )

  p
}

# ── 4. Summary box ────────────────────────────────────────────
summary_text <- function(df, subtitle) {
  n  <- nrow(df)
  bc <- table(df$bor)
  cr <- as.integer(bc["CR"]); if (is.na(cr)) cr <- 0
  pr <- as.integer(bc["PR"]); if (is.na(pr)) pr <- 0
  sd <- as.integer(bc["SD"]); if (is.na(sd)) sd <- 0
  pd <- as.integer(bc["PD"]); if (is.na(pd)) pd <- 0
  on <- sum(df$ongoing)
  orr <- (cr + pr) / n * 100
  paste0(subtitle, "\n",
         strrep("\u2500", 24), "\n",
         "N = ", n, "    Ongoing: ", on,
         " (", round(on / n * 100), "%)\n",
         "CR ", cr, "    PR ", pr,
         "    SD ", sd, "    PD ", pd, "\n",
         "ORR = ", round(orr, 1), "%")
}

# ── 5. Legend ─────────────────────────────────────────────────
build_legend_plot <- function(resp_subset = c("CR","PR","SD","PD")) {
  df_leg <- data.frame(x = 1, y = 1)
  p <- ggplot(df_leg) +
    # Tumor patches
    lapply(seq_along(TUMOR_COLORS), function(i) {
      annotate("rect", xmin = 0.05, xmax = 0.35,
               ymin = 20 - i * 1.2, ymax = 20 - i * 1.2 + 0.8,
               fill = TUMOR_COLORS[i], alpha = 0.87)
    }) |> (\(x) Reduce("+", x, ggplot(df_leg)))() +
    xlim(0, 4) + ylim(0, 24) +
    theme_void() +
    theme(plot.background = element_rect(fill = "#FAFAFA", color = NA))
  # Note: for production, use cowplot or ggpubr for proper legends
  p
}

# ── 6. Variants ───────────────────────────────────────────────
# Variant A — Responders (CR + PR)
adsl_A  <- trt |> filter(BESTRSPC %in% c("CR","PR"))
res_A   <- build_patient_rows(adsl_A, ADRS, CUTOFF, sort_by_bor = TRUE)
df_A    <- res_A$df;  rsp_A <- res_A$resp_df

p_A <- draw_swimmer(df_A, rsp_A,
  "Swimmer Plot \u2014 ONCVIZ-001 \u00b7 Vizatinib 300 mg QD\nResponders (CR + PR)  \u00b7  Treatment Arm  \u00b7  Cutoff: 31 May 2026",
  show_bor_labels = TRUE)

h_A <- max(nrow(df_A) * 0.42 + 3, 10)
ggsave("swimmer_A_responders.png", p_A,
       width = 16, height = min(h_A, 28),
       dpi = 300, bg = "#FAFAFA")
cat("Saved: swimmer_A_responders.png\n")


# Variant B — ≥12 months + Responders
adsl_B <- trt |> filter(TRTDUR_MO >= 12, BESTRSPC %in% c("CR","PR"))
res_B  <- build_patient_rows(adsl_B, ADRS, CUTOFF, sort_by_bor = TRUE)
df_B   <- res_B$df;  rsp_B <- res_B$resp_df

p_B <- draw_swimmer(df_B, rsp_B,
  "Swimmer Plot \u2014 ONCVIZ-001 \u00b7 Vizatinib 300 mg QD\nResponders (CR + PR) with \u226512 Months on Treatment  \u00b7  Cutoff: 31 May 2026",
  show_bor_labels = TRUE)

h_B <- max(nrow(df_B) * 0.42 + 3, 10)
ggsave("swimmer_B_12mo_responders.png", p_B,
       width = 16, height = min(h_B, 28),
       dpi = 300, bg = "#FAFAFA")
cat("Saved: swimmer_B_12mo_responders.png\n")


# Variant C — Responders by Tumor Type
TUMOR_ORDER <- c("NSCLC","BRCA","HCC","CRC","PDAC")
adsl_C <- trt |> filter(BESTRSPC %in% c("CR","PR"))
res_C  <- build_patient_rows(adsl_C, ADRS, CUTOFF, tumor_order = TUMOR_ORDER)
df_C   <- res_C$df;  rsp_C <- res_C$resp_df

p_C <- draw_swimmer(df_C, rsp_C,
  "Swimmer Plot \u2014 ONCVIZ-001 \u00b7 Vizatinib 300 mg QD\nResponders (CR + PR) by Tumor Type  \u00b7  Treatment Arm  \u00b7  Cutoff: 31 May 2026",
  show_bor_labels = FALSE, tumor_order = TUMOR_ORDER)

# Add tumor group panel labels (left strip)
tumor_groups <- df_C |>
  group_by(tumor) |>
  summarise(y_top = max(y), y_bot = min(y), .groups = "drop") |>
  mutate(mid_y = (y_top + y_bot) / 2,
         fill_col = TUMOR_COLORS[tumor])

for (i in seq_len(nrow(tumor_groups))) {
  tg <- tumor_groups[i, ]
  p_C <- p_C +
    annotate("rect",
             xmin = -1.4, xmax = -0.2,
             ymin = tg$y_bot - 0.45, ymax = tg$y_top + 0.45,
             fill = tg$fill_col, alpha = 0.12) +
    annotate("text",
             x = -0.9, y = tg$mid_y, label = tg$tumor,
             angle = 90, fontface = "bold", size = 3,
             color = tg$fill_col) +
    annotate("segment",
             x = -1.5, xend = max(df_C$trt_mo, na.rm=TRUE) + 6,
             y = tg$y_bot - 0.55, yend = tg$y_bot - 0.55,
             color = "#CCCCCC", linewidth = 0.6)
}

h_C <- max(nrow(df_C) * 0.42 + 3, 10)
ggsave("swimmer_C_responders_by_histology.png", p_C,
       width = 20, height = min(h_C, 30),
       dpi = 300, bg = "#FAFAFA")
cat("Saved: swimmer_C_responders_by_histology.png\n")

cat("\nAll 3 variants complete.\n")
