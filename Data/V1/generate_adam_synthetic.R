generate_adam_oncviz <- function(output_dir = ".", seed = 42, n_patients = 400L, verbose = TRUE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
  })

  set.seed(seed)

  N           <- as.integer(n_patients)
  CUTOFF      <- as.Date("2026-03-05")
  STUDY_START <- as.Date("2022-06-01")
  STUDY_ID    <- "ONCVIZ-001"
  DRUG_NAME   <- "Vizatinib"

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  log_msg <- function(msg) if (verbose) message(msg)

  clamp     <- function(v, lo, hi) pmax(lo, pmin(hi, v))
  fmt_date  <- function(d) format(as.Date(d), "%Y-%m-%d")
  wc        <- function(opts, wts) sample(opts, 1L, prob = wts)
  rweibull2 <- function(shape, scale) scale * (-log(max(runif(1), 1e-10)))^(1 / shape)
  mk        <- function(p) ifelse(runif(1) < p, "Y", "N")

  TUMOR_TYPES <- c("NSCLC", "CRC", "HCC", "PDAC", "BRCA")
  SITES       <- sprintf("SITE-%02d", 1:20)
  RESP_NAMES  <- c("CR", "PR", "SD", "PD")

  RESPONSE_RATES <- list(
    NSCLC = list(
      TREATMENT = c(CR = 0.07, PR = 0.38, SD = 0.32, PD = 0.23),
      CONTROL   = c(CR = 0.01, PR = 0.18, SD = 0.38, PD = 0.43)
    ),
    HCC = list(
      TREATMENT = c(CR = 0.04, PR = 0.22, SD = 0.38, PD = 0.36),
      CONTROL   = c(CR = 0.01, PR = 0.08, SD = 0.35, PD = 0.56)
    ),
    CRC = list(
      TREATMENT = c(CR = 0.02, PR = 0.10, SD = 0.33, PD = 0.55),
      CONTROL   = c(CR = 0.00, PR = 0.03, SD = 0.27, PD = 0.70)
    ),
    BRCA = list(
      TREATMENT = c(CR = 0.06, PR = 0.25, SD = 0.37, PD = 0.32),
      CONTROL   = c(CR = 0.01, PR = 0.12, SD = 0.40, PD = 0.47)
    ),
    PDAC = list(
      TREATMENT = c(CR = 0.01, PR = 0.09, SD = 0.32, PD = 0.58),
      CONTROL   = c(CR = 0.00, PR = 0.02, SD = 0.25, PD = 0.73)
    )
  )

  SURVIVAL_PARAMS <- list(
    NSCLC = list(os_sh = 0.92, os_sc = 24.8, pfs_sh = 0.78, pfs_sc = 15.2, os_hr = 0.49, pfs_hr = 0.52),
    HCC   = list(os_sh = 0.90, os_sc = 19.2, pfs_sh = 0.80, pfs_sc = 10.6, os_hr = 0.58, pfs_hr = 0.59),
    CRC   = list(os_sh = 0.88, os_sc = 14.0, pfs_sh = 0.82, pfs_sc =  7.8, os_hr = 0.72, pfs_hr = 0.75),
    BRCA  = list(os_sh = 0.85, os_sc = 40.0, pfs_sh = 0.80, pfs_sc = 18.0, os_hr = 0.61, pfs_hr = 0.63),
    PDAC  = list(os_sh = 0.95, os_sc =  8.5, pfs_sh = 0.90, pfs_sc =  4.2, os_hr = 0.76, pfs_hr = 0.78)
  )

  MUTATION_PREV <- list(
    NSCLC = c(0.46, 0.30, 0.15, 0.17, 0.12, 0.03, 0.08, 0.07, 0.06, 0.12, 0.08, 0.02, 0.07, 0.05, 0.04),
    CRC   = c(0.60, 0.45, 0.02, 0.03, 0.05, 0.10, 0.03, 0.20, 0.08, 0.05, 0.05, 0.32, 0.12, 0.06, 0.02),
    HCC   = c(0.28, 0.03, 0.01, 0.02, 0.04, 0.02, 0.05, 0.05, 0.04, 0.08, 0.04, 0.03, 0.18, 0.05, 0.02),
    PDAC  = c(0.72, 0.92, 0.02, 0.04, 0.06, 0.03, 0.04, 0.03, 0.03, 0.30, 0.06, 0.22, 0.08, 0.04, 0.02),
    BRCA  = c(0.38, 0.03, 0.02, 0.05, 0.03, 0.02, 0.02, 0.35, 0.22, 0.15, 0.15, 0.02, 0.14, 0.04, 0.02)
  )

  GENE_LIST <- c("TP53", "KRAS", "EGFR", "STK11", "KEAP1", "BRAF",
                 "MET", "PIK3CA", "PTEN", "CDKN2A", "RB1", "SMAD4",
                 "ARID1A", "NF1", "RET")

  TOXICITY_MODS <- list(
    NSCLC = c(hepatic = 1.0, gi = 1.0, derm = 1.0, hemato = 1.0),
    HCC   = c(hepatic = 2.2, gi = 1.1, derm = 0.9, hemato = 0.9),
    CRC   = c(hepatic = 0.8, gi = 1.8, derm = 1.2, hemato = 1.1),
    BRCA  = c(hepatic = 0.9, gi = 1.2, derm = 1.3, hemato = 1.4),
    PDAC  = c(hepatic = 1.2, gi = 2.0, derm = 0.7, hemato = 1.3)
  )

  DISC_TRT <- c("Adverse Event", "Progressive Disease", "Withdrawal of Consent", "Death", "Protocol Deviation")
  DISC_CTL <- c("Progressive Disease", "Withdrawal of Consent", "Death", "Crossover to Treatment", "Protocol Deviation")

  TRANS <- list(
    CR = c(CR = 0.82, PR = 0.10, SD = 0.05, PD = 0.03),
    PR = c(CR = 0.10, PR = 0.58, SD = 0.20, PD = 0.12),
    SD = c(CR = 0.02, PR = 0.10, SD = 0.52, PD = 0.36),
    PD = c(CR = 0.00, PR = 0.00, SD = 0.00, PD = 1.00)
  )

  log_msg("Generating ADRAND...")

  n_screened    <- as.integer(N * 1.15)
  n_screen_fail <- n_screened - N
  sf_reasons    <- c(
    "Does not meet inclusion criteria",
    "Laboratory values outside protocol limits",
    "Prior treatment exclusion",
    "Patient declined",
    "ECOG performance status >2",
    "Active CNS metastases"
  )

  adrand_rows <- map_dfr(seq_len(n_screened), function(i) {
    site     <- sample(SITES, 1L)
    scr_id   <- sprintf("%s-SCR-%04d", STUDY_ID, i)
    scr_date <- STUDY_START + sample(0:540, 1L)
    if (i <= n_screen_fail) {
      tibble(
        STUDYID = STUDY_ID, SCRSEQ = i, SCRID = scr_id, USUBJID = "",
        SITEID = site, SCRDTC = fmt_date(scr_date),
        RANDFL = "N", SFFL = "Y", SFREASN = sample(sf_reasons, 1L),
        ARM = "", SAFFL = "N", ITTFL = "N"
      )
    } else {
      rand_idx <- i - n_screen_fail
      uid      <- sprintf("%s-%s-%04d", STUDY_ID, substr(site, 6, 7), rand_idx)
      arm      <- wc(c("TREATMENT", "CONTROL"), c(2, 1))
      tibble(
        STUDYID = STUDY_ID, SCRSEQ = i, SCRID = scr_id, USUBJID = uid,
        SITEID = site, SCRDTC = fmt_date(scr_date),
        RANDFL = "Y", SFFL = "N", SFREASN = "",
        ARM = arm, SAFFL = "Y", ITTFL = "Y"
      )
    }
  })

  log_msg(sprintf("  ADRAND: %d records | Screened=%d Randomized=%d Screen-fail=%d",
                  nrow(adrand_rows), n_screened, N, n_screen_fail))

  log_msg("Generating ADSL...")

  rand_subjects <- adrand_rows %>% filter(RANDFL == "Y") %>% arrange(SCRSEQ)

  adsl <- map_dfr(seq_len(N), function(i) {
    rr   <- rand_subjects[i, ]
    site <- rr$SITEID
    uid  <- rr$USUBJID
    arm  <- rr$ARM
    armcd <- ifelse(arm == "TREATMENT", "TRT", "CTL")
    tt   <- sample(TUMOR_TYPES, 1L)

    age     <- clamp(round(abs(rnorm(1, 63, 11))), 28L, 85L)
    sex     <- wc(c("M", "F"), c(0.509, 0.491))
    race    <- wc(c("WHITE", "BLACK", "ASIAN", "OTHER"), c(0.88, 0.098, 0.021, 0.001))
    smoking <- wc(c("Former heavy", "Never", "Former light", "Current heavy"), c(0.488, 0.322, 0.178, 0.012))
    ecog    <- wc(0:2, c(0.38, 0.50, 0.12))
    bmi     <- round(clamp(rnorm(1, 25.8, 4.5), 16, 42), 1)

    t_stage <- wc(c("T1", "T2", "T3", "T4"), c(0.22, 0.38, 0.28, 0.12))
    n_stage <- wc(c("N0", "N1", "N2", "N3"), c(0.28, 0.30, 0.32, 0.10))
    m_stage <- wc(c("M0", "M1"), c(0.65, 0.35))
    stage   <- if (m_stage == "M1") "STAGE IV" else
               if (n_stage %in% c("N2", "N3")) wc(c("STAGE IIIA", "STAGE IIIB"), c(0.65, 0.35)) else
               if (n_stage == "N1") wc(c("STAGE IIA", "STAGE IIB"), c(0.55, 0.45)) else
               wc(c("STAGE IA", "STAGE IB"), c(0.55, 0.45))

    pdl1_score <- round(rbeta(1, 1.2, 2.5) * 100, 1)
    pdl1_grp   <- ifelse(pdl1_score >= 50, "HIGH", ifelse(pdl1_score >= 1, "MED", "NEG"))
    tmb        <- clamp(round(rlnorm(1, 2.1, 0.85), 1), 0.5, 80)
    msi_score  <- clamp(round(rexp(1, 1 / 1.2), 3), 0, 15)
    msi_sts    <- ifelse(msi_score >= 3.5, "MSI-H", "MSS")

    mp <- MUTATION_PREV[[tt]]

    trtsdt  <- STUDY_START + sample(0:540, 1L)
    if (trtsdt > CUTOFF - 60) trtsdt <- CUTOFF - 60
    max_dur <- as.integer(CUTOFF - trtsdt)

    still_on <- runif(1) < ifelse(arm == "TREATMENT", 0.42, 0.22)
    r        <- RESPONSE_RATES[[tt]][[arm]]
    bor      <- sample(RESP_NAMES, 1L, prob = r)

    sp      <- SURVIVAL_PARAMS[[tt]]
    os_sc   <- ifelse(arm == "TREATMENT", sp$os_sc, sp$os_sc * sp$os_hr)
    pfs_sc  <- ifelse(arm == "TREATMENT", sp$pfs_sc, sp$pfs_sc * sp$pfs_hr)
    os_m    <- clamp(rweibull2(sp$os_sh, os_sc), 1, 120)
    pfs_m   <- clamp(rweibull2(sp$pfs_sh, pfs_sc), 0.5, os_m)
    os_d    <- as.integer(os_m * 30.44)
    pfs_d   <- as.integer(pfs_m * 30.44)
    dur     <- ifelse(still_on, max_dur, clamp(os_d, 42L, max_dur))
    trtedt  <- trtsdt + dur

    eosstt  <- ifelse(still_on, "ONGOING", "DISCONTINUED")
    dcsreas <- ifelse(still_on, "", sample(if (arm == "TREATMENT") DISC_TRT else DISC_CTL, 1L))

    os_event <- ifelse(still_on, 0L, ifelse(os_d <= max_dur, 1L, 0L))

    noncancer_death_p <- ifelse(arm == "TREATMENT", 0.04, 0.06)
    noncancer_death   <- (!still_on) && os_event == 1 && runif(1) < noncancer_death_p
    pfs_event <- ifelse(noncancer_death, 0L, ifelse(pfs_d <= dur, 1L, 0L))

    comp_event <- 0L
    comp_type  <- "CENSORED"
    if (noncancer_death) {
      comp_event <- 1L; comp_type <- "DEATH_WITHOUT_PROGRESSION"
    } else if (pfs_event == 1L) {
      comp_event <- 1L; comp_type <- "PROGRESSION"
    } else if (os_event == 1L) {
      comp_event <- 1L; comp_type <- "DEATH_WITHOUT_PROGRESSION"
    }

    tibble(
      STUDYID = STUDY_ID, USUBJID = uid, SITEID = site,
      ARM = arm, ARMCD = armcd,
      TRT01P = ifelse(arm == "TREATMENT", DRUG_NAME, "Placebo"),
      TUMORTYPE = tt,
      AGE = age, AGEGR1 = ifelse(age < 65, "<65", ">=65"),
      SEX = sex, RACE = race, SMOKING = smoking, ECOG = ecog, BMI = bmi,
      T_STAGE = t_stage, N_STAGE = n_stage, M_STAGE = m_stage, STAGE = stage,
      LIVERMETS = mk(0.20), BRAINMETS = mk(0.10),
      PRIORSURG = mk(0.35), PRIORRAD = mk(0.28),
      PRIORLINES = wc(1:4, c(0.35, 0.35, 0.20, 0.10)),
      BASESZ = round(rlnorm(1, 3.4, 0.55), 1),
      PDL1SCORE = pdl1_score, PDL1GRP = pdl1_grp,
      TMB = tmb, TMBHIGH = ifelse(tmb >= 10, "Y", "N"),
      MSI_SCORE = msi_score, MSISTS = msi_sts,
      TP53MUT  = mk(mp[1]),  KRASMUT  = mk(mp[2]),  EGFRMUT  = mk(mp[3]),
      STK11MUT = mk(mp[4]),  KEAP1MUT = mk(mp[5]),  BRAFMUT  = mk(mp[6]),
      METAMP   = mk(mp[7]),  PIK3CAMT = mk(mp[8]),  PTENMUT  = mk(mp[9]),
      CDKN2AMT = mk(mp[10]), RB1MUT   = mk(mp[11]), SMAD4MUT = mk(mp[12]),
      ARID1AMT = mk(mp[13]), NF1MUT   = mk(mp[14]), RETFUS   = mk(mp[15]),
      BESTRSPC = bor,
      TRTSDT = fmt_date(trtsdt), TRTEDT = fmt_date(trtedt),
      TRTDURD = dur, EOSSTT = eosstt, DCSREAS = dcsreas,
      OSDUR = os_d, OSCR = os_event,
      OSDTC = fmt_date(trtsdt + os_d),
      PFSDUR = pfs_d, PFSCR = pfs_event,
      PFSDTC = fmt_date(trtsdt + pfs_d),
      COMPEVENT = comp_event, COMPTYPE = comp_type,
      CUTDTC = fmt_date(CUTOFF),
      SAFFL = "Y", ITTFL = "Y",
      PPROTFL = ifelse(runif(1) > 0.07, "Y", "N")
    )
  })

  trt_n <- sum(adsl$ARM == "TREATMENT")
  ctl_n <- sum(adsl$ARM == "CONTROL")
  log_msg(sprintf("  ADSL: %d patients | TRT=%d CTL=%d", nrow(adsl), trt_n, ctl_n))
  for (tt in TUMOR_TYPES) {
    sub <- adsl[adsl$TUMORTYPE == tt & adsl$ARM == "TREATMENT", ]
    orr <- mean(sub$BESTRSPC %in% c("CR", "PR")) * 100
    log_msg(sprintf("    %s: n=%d ORR=%.0f%%", tt, nrow(sub), orr))
  }

  log_msg("Generating ADRS...")

  adrs <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    trtsdt <- as.Date(s$TRTSDT)
    trtedt <- as.Date(s$TRTEDT)
    dur    <- as.integer(trtedt - trtsdt)
    n_vis  <- max(1L, dur %/% 42L)
    r      <- RESPONSE_RATES[[s$TUMORTYPE]][[s$ARM]]
    prev   <- NULL
    rows   <- list()
    for (k in seq_len(n_vis)) {
      aday <- 42L * k
      if (aday > dur + 21L) break
      adt <- trtsdt + aday
      if (adt > CUTOFF) break
      resp <- if (is.null(prev)) sample(RESP_NAMES, 1L, prob = r) else
              sample(RESP_NAMES, 1L, prob = TRANS[[prev]])
      rows[[k]] <- tibble(
        STUDYID = STUDY_ID, USUBJID = s$USUBJID, ARM = s$ARM, TUMORTYPE = s$TUMORTYPE,
        PARAMCD = "OVRLRESP", PARAM = "Overall Response per RECIST 1.1",
        AVISIT = sprintf("CYCLE %d DAY 1", k), AVISITN = k,
        ADT = fmt_date(adt), ADTN = aday,
        AVALC = resp, AVAL = match(resp, RESP_NAMES),
        BICR_CONF = ifelse(runif(1) > 0.05, "Y", "N"),
        ANL01FL = "Y"
      )
      prev <- resp
      if (resp == "PD") break
    }
    bind_rows(rows)
  })

  log_msg(sprintf("  ADRS: %d records", nrow(adrs)))

  log_msg("Generating ADTR...")

  adtr <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    trtsdt <- as.Date(s$TRTSDT)
    trtedt <- as.Date(s$TRTEDT)
    dur    <- as.integer(trtedt - trtsdt)
    bl     <- s$BASESZ
    bor    <- s$BESTRSPC
    n_vis  <- max(1L, dur %/% 42L)

    final_val <- switch(bor,
      CR = bl * runif(1, 0.01, 0.03),
      PR = bl * runif(1, 0.25, 0.65),
      SD = bl * runif(1, 0.82, 1.19),
      PD = bl * runif(1, 1.20, 2.20)
    )
    nadir <- if (bor == "PR") bl * runif(1, 0.25, 0.65) else NA_real_

    rows <- list()
    rows[[1]] <- tibble(
      STUDYID = STUDY_ID, USUBJID = s$USUBJID, ARM = s$ARM, TUMORTYPE = s$TUMORTYPE,
      PARAMCD = "SUMDIAM", PARAM = "Sum of Longest Diameters (mm)",
      AVISIT = "BASELINE", AVISITN = 0L,
      ADT = fmt_date(trtsdt), ADTN = 0L,
      BASE = bl, AVAL = bl, CHG = 0.0, PCHG = 0.0, ANL01FL = "Y"
    )
    cur <- bl
    for (k in seq_len(n_vis)) {
      aday <- 42L * k
      if (aday > dur + 21L) break
      adt  <- trtsdt + aday
      if (adt > CUTOFF) break
      ph   <- k / max(n_vis, 1)
      cur  <- switch(bor,
        CR = max(0.01, cur + (bl * 0.02 - cur) * 0.45 + rnorm(1, 0, bl * 0.008)),
        PR = {
          tgt <- if (ph < 0.55) nadir else nadir * runif(1, 1.0, 1.15)
          max(0.1, cur + (tgt - cur) * 0.40 + rnorm(1, 0, bl * 0.015))
        },
        SD = max(0.1, cur + (bl * runif(1, 0.82, 1.18) - cur) * 0.25 + rnorm(1, 0, bl * 0.03)),
        PD = max(bl * 1.01, cur + (final_val - cur) * 0.35)
      )
      chg  <- round(cur - bl, 1)
      pchg <- round((cur - bl) / bl * 100, 1)
      rows[[k + 1L]] <- tibble(
        STUDYID = STUDY_ID, USUBJID = s$USUBJID, ARM = s$ARM, TUMORTYPE = s$TUMORTYPE,
        PARAMCD = "SUMDIAM", PARAM = "Sum of Longest Diameters (mm)",
        AVISIT = sprintf("CYCLE %d DAY 1", k), AVISITN = k,
        ADT = fmt_date(adt), ADTN = aday,
        BASE = bl, AVAL = round(cur, 1), CHG = chg, PCHG = pchg, ANL01FL = "Y"
      )
    }
    bind_rows(rows)
  })

  log_msg(sprintf("  ADTR: %d records", nrow(adtr)))

  log_msg("Generating ADAE...")

  AE_PROFILE <- tribble(
    ~SOC,                                      ~PT,                                   ~INC_TRT, ~INC_CTL, ~GW1,  ~GW2,  ~GW3,  ~GW4,  ~GW5,  ~TOX,
    "Gastrointestinal disorders",              "Nausea",                               0.58, 0.32, 0.38, 0.40, 0.18, 0.03, 0.01, "gi",
    "Gastrointestinal disorders",              "Diarrhea",                             0.52, 0.25, 0.32, 0.38, 0.24, 0.05, 0.01, "gi",
    "Gastrointestinal disorders",              "Vomiting",                             0.34, 0.18, 0.38, 0.40, 0.18, 0.03, 0.01, "gi",
    "Gastrointestinal disorders",              "Constipation",                         0.26, 0.20, 0.52, 0.36, 0.10, 0.02, 0.00, "gi",
    "Gastrointestinal disorders",              "Abdominal pain",                       0.20, 0.16, 0.42, 0.38, 0.16, 0.04, 0.00, "gi",
    "Gastrointestinal disorders",              "Mucositis",                            0.16, 0.04, 0.38, 0.40, 0.18, 0.04, 0.00, "gi",
    "Skin and subcutaneous tissue disorders",  "Rash",                                 0.45, 0.12, 0.42, 0.38, 0.16, 0.04, 0.00, "derm",
    "Skin and subcutaneous tissue disorders",  "Pruritus",                             0.22, 0.08, 0.55, 0.35, 0.08, 0.02, 0.00, "derm",
    "Skin and subcutaneous tissue disorders",  "Alopecia",                             0.20, 0.08, 0.60, 0.35, 0.05, 0.00, 0.00, "derm",
    "Skin and subcutaneous tissue disorders",  "Palmar-plantar erythrodysaesthesia",   0.22, 0.06, 0.40, 0.38, 0.18, 0.04, 0.00, "derm",
    "General disorders",                       "Fatigue",                              0.65, 0.48, 0.30, 0.40, 0.22, 0.07, 0.01, "gi",
    "General disorders",                       "Decreased appetite",                   0.42, 0.30, 0.40, 0.38, 0.18, 0.04, 0.00, "gi",
    "General disorders",                       "Pyrexia",                              0.22, 0.16, 0.48, 0.36, 0.13, 0.03, 0.00, "gi",
    "General disorders",                       "Peripheral edema",                     0.18, 0.10, 0.52, 0.36, 0.10, 0.02, 0.00, "gi",
    "Investigations",                          "ALT increased",                        0.32, 0.10, 0.40, 0.35, 0.16, 0.07, 0.02, "hepatic",
    "Investigations",                          "AST increased",                        0.30, 0.08, 0.42, 0.35, 0.15, 0.06, 0.02, "hepatic",
    "Investigations",                          "Blood bilirubin increased",            0.14, 0.05, 0.48, 0.34, 0.12, 0.05, 0.01, "hepatic",
    "Investigations",                          "Platelet count decreased",             0.26, 0.07, 0.40, 0.34, 0.18, 0.07, 0.01, "hemato",
    "Investigations",                          "Neutrophil count decreased",           0.32, 0.08, 0.32, 0.30, 0.24, 0.10, 0.04, "hemato",
    "Investigations",                          "Hemoglobin decreased",                 0.38, 0.22, 0.38, 0.36, 0.18, 0.07, 0.01, "hemato",
    "Nervous system disorders",                "Peripheral neuropathy",                0.26, 0.08, 0.42, 0.38, 0.16, 0.04, 0.00, "gi",
    "Nervous system disorders",                "Headache",                             0.20, 0.16, 0.52, 0.36, 0.10, 0.02, 0.00, "gi",
    "Nervous system disorders",                "Dizziness",                            0.16, 0.12, 0.58, 0.35, 0.07, 0.00, 0.00, "gi",
    "Musculoskeletal disorders",               "Arthralgia",                           0.24, 0.16, 0.44, 0.38, 0.14, 0.04, 0.00, "gi",
    "Musculoskeletal disorders",               "Myalgia",                              0.18, 0.14, 0.48, 0.38, 0.12, 0.02, 0.00, "gi",
    "Respiratory disorders",                   "Dyspnea",                              0.20, 0.16, 0.42, 0.38, 0.14, 0.05, 0.01, "gi",
    "Respiratory disorders",                   "Cough",                                0.18, 0.14, 0.58, 0.36, 0.06, 0.00, 0.00, "gi",
    "Metabolism disorders",                    "Hypokalemia",                          0.16, 0.07, 0.44, 0.36, 0.14, 0.05, 0.01, "gi",
    "Metabolism disorders",                    "Hyperglycemia",                        0.20, 0.08, 0.42, 0.36, 0.16, 0.05, 0.01, "gi",
    "Cardiac disorders",                       "QT prolongation",                      0.10, 0.04, 0.52, 0.34, 0.10, 0.04, 0.00, "gi",
    "Cardiac disorders",                       "Hypertension",                         0.16, 0.08, 0.44, 0.36, 0.14, 0.05, 0.01, "gi",
    "Immune system disorders",                 "Hypothyroidism",                       0.12, 0.03, 0.48, 0.38, 0.12, 0.02, 0.00, "gi",
    "Immune system disorders",                 "Pneumonitis",                          0.07, 0.02, 0.38, 0.32, 0.18, 0.10, 0.02, "gi",
    "Infections",                              "Urinary tract infection",              0.10, 0.10, 0.48, 0.38, 0.12, 0.02, 0.00, "gi",
    "Infections",                              "Upper respiratory infection",          0.12, 0.12, 0.52, 0.36, 0.10, 0.02, 0.00, "gi"
  )

  adae_rows <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    trtsdt <- as.Date(s$TRTSDT)
    trtedt <- as.Date(s$TRTEDT)
    dur    <- max(1L, as.integer(trtedt - trtsdt))
    arm    <- s$ARM
    uid    <- s$USUBJID
    tt     <- s$TUMORTYPE
    mods   <- TOXICITY_MODS[[tt]]
    seq_n  <- 0L
    rows   <- list()

    for (j in seq_len(nrow(AE_PROFILE))) {
      ae      <- AE_PROFILE[j, ]
      base_inc <- ifelse(arm == "TREATMENT", ae$INC_TRT, ae$INC_CTL)
      mod     <- mods[ae$TOX]
      inc     <- clamp(base_inc * mod, 0, 0.95)
      if (runif(1) > inc) next
      n_occ <- sample(1:3, 1L, prob = c(0.72, 0.22, 0.06))
      for (occ in seq_len(n_occ)) {
        seq_n  <- seq_n + 1L
        onset  <- sample(seq_len(dur), 1L)
        grade  <- sample(1:5, 1L, prob = c(ae$GW1, ae$GW2, ae$GW3, ae$GW4, ae$GW5))
        dur_ae <- max(1L, as.integer(rexp(1, 1 / (12 + grade * 7))))
        dur_ae <- min(dur_ae, dur - onset + 30L)
        endt   <- min(trtsdt + onset + dur_ae, CUTOFF)
        serious <- grade >= 3
        related <- runif(1) < ifelse(arm == "TREATMENT", 0.62, 0.28)
        action  <- if (grade >= 2) sample(c("DOSE REDUCED", "DRUG INTERRUPTED", "DRUG WITHDRAWN", "NONE"),
                                          1L, prob = c(0.18, 0.22, 0.05, 0.55)) else "NONE"
        outcome <- sample(c("RECOVERED", "RECOVERING", "NOT RECOVERED", "FATAL", "UNKNOWN"),
                          1L, prob = c(0.68, 0.14, 0.11, 0.03, 0.04))
        rows[[length(rows) + 1L]] <- tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = tt,
          AESEQ = seq_n, AESOC = ae$SOC, AEPT = ae$PT, AETOXGR = grade,
          AESER = ifelse(serious, "Y", "N"),
          AEREL = ifelse(related, "RELATED", "NOT RELATED"),
          AEACN = action, AEOUT = outcome,
          AESTDTC = fmt_date(trtsdt + onset),
          AEENDTC = fmt_date(endt),
          AESTDY = onset, AEENDY = onset + dur_ae, AEDUR = dur_ae,
          CTCAE_V = "5.0"
        )
      }
    }

    if (length(rows) == 0L) {
      return(tibble(
        STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = tt,
        AESEQ = 1L, AESOC = "General disorders", AEPT = "Fatigue", AETOXGR = 1L,
        AESER = "N", AEREL = "NOT RELATED", AEACN = "NONE", AEOUT = "RECOVERED",
        AESTDTC = fmt_date(trtsdt + min(7L, dur)),
        AEENDTC = fmt_date(trtsdt + min(21L, dur)),
        AESTDY = 7L, AEENDY = 21L, AEDUR = 14L, CTCAE_V = "5.0"
      ))
    }

    bind_rows(rows)
  })

  log_msg(sprintf("  ADAE: %d records | Grade>=3: %.1f%%",
                  nrow(adae_rows), mean(adae_rows$AETOXGR >= 3) * 100))

  log_msg("Generating ADLB...")

  LAB_TESTS <- tribble(
    ~PARAMCD, ~PARAM,                           ~AVALU,    ~ULN,  ~LLN,  ~MU,    ~SD,   ~EFFECT, ~ETYPE,
    "ALT",   "Alanine Aminotransferase",         "U/L",     45,    7,     22,     8,     1.80,    "hepatic",
    "AST",   "Aspartate Aminotransferase",       "U/L",     40,    10,    20,     7,     1.60,    "hepatic",
    "BILI",  "Total Bilirubin",                  "mg/dL",    1.2,   0.2,   0.7,   0.2,   1.35,    "hepatic",
    "ALP",   "Alkaline Phosphatase",             "U/L",    120,    35,    70,    20,     1.20,    "increase",
    "GGT",   "Gamma-Glutamyl Transferase",       "U/L",     60,    8,     30,    12,     1.45,    "hepatic",
    "CREAT", "Creatinine",                       "mg/dL",    1.2,   0.5,   0.85,  0.18,  1.06,    "increase",
    "BUN",   "Blood Urea Nitrogen",              "mg/dL",   20,    6,     13,     4,     1.05,    "stable",
    "K",     "Potassium",                        "mEq/L",    5.1,   3.5,   4.0,   0.4,  -0.12,   "suppress",
    "SOD",   "Sodium",                           "mEq/L",  145,   135,   139,     2.5,   0.00,    "stable",
    "HGB",   "Hemoglobin",                       "g/dL",    17,    11,    13.5,   1.2,  -0.18,   "suppress",
    "PLT",   "Platelet Count",                  "10^9/L",  400,   150,   220,    45,    -0.20,   "suppress",
    "ANC",   "Absolute Neutrophil Count",       "10^9/L",    7.5,   1.8,   4.0,   1.0,  -0.25,   "suppress",
    "WBC",   "White Blood Cell Count",          "10^9/L",   10,    4.0,   6.5,   1.5,  -0.14,   "suppress",
    "LYMPH", "Lymphocyte Count",                "10^9/L",    4.8,   1.0,   2.2,   0.6,  -0.10,   "suppress",
    "QTCF",  "QTc (Fridericia)",                 "msec",   450,   350,   410,    22,     9.0,    "increase",
    "HR",    "Heart Rate",                       "bpm",    100,    50,    72,    12,     2.5,    "stable",
    "GLUC",  "Glucose",                          "mg/dL",  100,    70,    88,    14,     7.0,    "increase",
    "TSH",   "Thyroid Stimulating Hormone",      "mIU/L",    4.0,   0.4,   2.0,   0.8,   0.7,    "increase",
    "CHOL",  "Total Cholesterol",                "mg/dL",  200,   100,   175,    30,     5.0,    "stable",
    "ALB",   "Albumin",                          "g/dL",     5.2,   3.5,   4.2,   0.4,  -0.08,   "suppress",
    "PROT",  "Total Protein",                    "g/dL",     8.5,   6.0,   7.2,   0.5,  -0.05,   "stable"
  )

  VISITS_LB  <- c("BASELINE", "CYCLE 2 DAY 1", "CYCLE 4 DAY 1", "CYCLE 6 DAY 1",
                  "CYCLE 9 DAY 1", "CYCLE 12 DAY 1", "CYCLE 18 DAY 1", "EOT")
  VISIT_DAYS <- c(0L, 42L, 126L, 210L, 336L, 462L, 714L, NA_integer_)

  adlb <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s       <- adsl[i, ]
    trtsdt  <- as.Date(s$TRTSDT)
    trtedt  <- as.Date(s$TRTEDT)
    dur     <- as.integer(trtedt - trtsdt)
    arm     <- s$ARM
    tt      <- s$TUMORTYPE
    uid     <- s$USUBJID
    pat_re  <- rnorm(1, 0, 0.10)
    has_hep <- arm == "TREATMENT" && runif(1) < ifelse(tt == "HCC", 0.18, 0.08)
    hep_gr  <- if (has_hep) sample(2:4, 1L, prob = c(0.50, 0.35, 0.15)) else 0L

    map_dfr(seq_len(nrow(LAB_TESTS)), function(j) {
      lb       <- LAB_TESTS[j, ]
      base_val <- max(lb$LLN * 0.8, abs(rnorm(1, lb$MU, lb$SD)))

      map_dfr(seq_along(VISITS_LB), function(vi) {
        vday <- ifelse(is.na(VISIT_DAYS[vi]), dur, VISIT_DAYS[vi])
        if (vday > dur + 14L) return(NULL)
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)

        aval <- if (vi == 1L) base_val else {
          phase <- min(vday / max(dur, 1), 1.0)
          if (arm == "TREATMENT") {
            switch(lb$ETYPE,
              hepatic = {
                shp <- exp(-4 * (phase - 0.30)^2)
                if (has_hep && lb$PARAMCD %in% c("ALT", "AST", "BILI")) {
                  txuln <- c("2" = 4.0, "3" = 8.0, "4" = 15.0)[as.character(hep_gr)]
                  pv    <- lb$ULN * txuln * runif(1, 0.85, 1.15)
                  base_val + (pv - base_val) * shp + rnorm(1, 0, lb$SD * 0.20)
                } else {
                  base_val + lb$EFFECT * shp * lb$MU * (1 + pat_re) + rnorm(1, 0, lb$SD * 0.20)
                }
              },
              suppress = {
                drop <- abs(lb$EFFECT) * exp(-5 * (phase - 0.40)^2)
                base_val * (1 - drop) + rnorm(1, 0, lb$SD * 0.12)
              },
              increase = base_val + lb$EFFECT * phase + rnorm(1, 0, ifelse(lb$PARAMCD == "QTCF", 5, lb$SD * 0.15)),
              base_val + rnorm(1, 0, lb$SD * 0.10)
            )
          } else {
            base_val + rnorm(1, 0, lb$SD * 0.08)
          }
        }

        aval <- max(0.01, aval)
        chg  <- round(aval - base_val, 3)
        pchg <- round((aval - base_val) / base_val * 100, 1)
        xuln <- round(aval / lb$ULN, 3)

        grade <- if (lb$PARAMCD %in% c("ALT", "AST")) {
          ifelse(aval > 20 * lb$ULN, 5L,
          ifelse(aval > 10 * lb$ULN, 4L,
          ifelse(aval >  5 * lb$ULN, 3L,
          ifelse(aval >  3 * lb$ULN, 2L,
          ifelse(aval >      lb$ULN, 1L, 0L)))))
        } else if (lb$PARAMCD == "BILI") {
          ifelse(aval > 10 * lb$ULN, 4L,
          ifelse(aval >  3 * lb$ULN, 3L,
          ifelse(aval >  1.5 * lb$ULN, 2L,
          ifelse(aval >      lb$ULN, 1L, 0L))))
        } else if (lb$PARAMCD == "ANC") {
          ifelse(aval < 0.5, 4L, ifelse(aval < 1.0, 3L, ifelse(aval < 1.5, 2L, ifelse(aval < lb$LLN, 1L, 0L))))
        } else if (lb$PARAMCD == "PLT") {
          ifelse(aval < 25, 4L, ifelse(aval < 50, 3L, ifelse(aval < 75, 2L, ifelse(aval < lb$LLN, 1L, 0L))))
        } else if (lb$PARAMCD == "HGB") {
          ifelse(aval < 8.0, 4L, ifelse(aval < 10.0, 3L, ifelse(aval < lb$LLN, 2L, 0L)))
        } else if (lb$PARAMCD == "CREAT") {
          ifelse(aval > 6 * lb$ULN, 4L, ifelse(aval > 3 * lb$ULN, 3L, ifelse(aval > 1.5 * lb$ULN, 2L, ifelse(aval > lb$ULN, 1L, 0L))))
        } else if (lb$PARAMCD == "QTCF") {
          ifelse(aval > 500, 4L, ifelse(aval > 480, 2L, ifelse(aval > 450, 1L, 0L)))
        } else {
          0L
        }

        tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = tt,
          PARAMCD = lb$PARAMCD, PARAM = lb$PARAM, AVALU = lb$AVALU,
          AVISIT = VISITS_LB[vi], AVISITN = vi - 1L,
          ADT = fmt_date(adt), ADTN = vday,
          BASE = round(base_val, 3), AVAL = round(aval, 3),
          CHG = chg, PCHG = pchg, ULN = lb$ULN, LLN = lb$LLN,
          XULNFL = xuln, ATOXGR = grade,
          ANRIND = ifelse(aval > lb$ULN, "HIGH", ifelse(aval < lb$LLN, "LOW", "NORMAL")),
          ANL01FL = "Y"
        )
      })
    })
  })

  log_msg(sprintf("  ADLB: %d records | %d parameters", nrow(adlb), n_distinct(adlb$PARAMCD)))

  log_msg("Generating ADTTE...")

  SUBGRP_VARS <- c("TUMORTYPE", "AGEGR1", "SEX", "ECOG", "PDL1GRP", "MSISTS",
                   "TMBHIGH", "LIVERMETS", "PRIORLINES", "SMOKING",
                   "EGFRMUT", "KRASMUT", "TP53MUT", "STK11MUT", "STAGE")

  adtte <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    uid    <- s$USUBJID
    arm    <- s$ARM
    armcd  <- s$ARMCD
    trtsdt <- as.Date(s$TRTSDT)
    sg     <- s[SUBGRP_VARS]

    mrow <- function(paramcd, param, aval, cnsr, evnt, adt) {
      bind_cols(
        tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, ARMCD = armcd,
          PARAMCD = paramcd, PARAM = param,
          AVAL = aval, AVALU = "DAYS",
          AVALM = round(aval / 30.44, 2),
          CNSR = cnsr, EVNTDESC = evnt,
          ADT = fmt_date(adt), STARTDT = s$TRTSDT,
          LM6MFL  = ifelse(aval >= 182, "Y", "N"),
          LM12MFL = ifelse(aval >= 365, "Y", "N"),
          LM24MFL = ifelse(aval >= 730, "Y", "N"),
          COMPEVENT = s$COMPEVENT,
          COMPTYPE  = s$COMPTYPE,
          ANL01FL = "Y"
        ),
        sg
      )
    }

    rows <- list(
      mrow("OS", "Overall Survival",
           s$OSDUR, 1L - s$OSCR,
           ifelse(s$OSCR == 1L, "DEATH", "CENSORED"),
           as.Date(s$OSDTC)),
      mrow("PFS", "Progression-Free Survival",
           s$PFSDUR, 1L - s$PFSCR,
           ifelse(s$PFSCR == 1L, "PROGRESSION", "CENSORED"),
           as.Date(s$PFSDTC))
    )

    if (s$BESTRSPC %in% c("CR", "PR")) {
      ttr <- sample(28:84, 1L)
      db  <- list(TREATMENT = c(CR = 14.0, PR = 9.5),
                  CONTROL   = c(CR =  6.0, PR = 4.0))[[arm]][[s$BESTRSPC]]
      dor <- max(42L, min(as.integer(rexp(1, 1 / (db * 30.44))), s$OSDUR - ttr))
      de  <- rbinom(1L, 1L, 0.65)
      rows[[3]] <- mrow("DOR", "Duration of Response",
                        dor, 1L - de,
                        ifelse(de == 1L, "PROGRESSION", "CENSORED"),
                        trtsdt + ttr + dor)
      rows[[4]] <- mrow("TTR", "Time to Response",
                        ttr, 0L, "RESPONSE",
                        trtsdt + ttr)
    }

    bind_rows(rows)
  })

  log_msg(sprintf("  ADTTE: %d records", nrow(adtte)))

  log_msg("Generating ADPK...")

  PK_DOSE   <- 300
  KA        <- 0.8
  CL_F      <- 18.0
  VD_F      <- 320.0
  PK_TIMES  <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)
  PK_VISITS <- c("CYCLE 1 DAY 1", "CYCLE 1 DAY 15", "CYCLE 3 DAY 1")
  PK_DAYS   <- c(0L, 14L, 56L)

  adpk <- map_dfr(seq_len(nrow(adsl[adsl$ARM == "TREATMENT", ])), function(i) {
    s      <- adsl[adsl$ARM == "TREATMENT", ][i, ]
    uid    <- s$USUBJID
    trtsdt <- as.Date(s$TRTSDT)
    cl_i   <- CL_F * exp(rnorm(1, 0, 0.28))
    vd_i   <- VD_F * exp(rnorm(1, 0, 0.25))
    ka_i   <- KA   * exp(rnorm(1, 0, 0.35))
    ke_i   <- cl_i / vd_i

    conc_rows <- map_dfr(seq_along(PK_VISITS), function(vi) {
      bdt <- trtsdt + PK_DAYS[vi]
      if (bdt > CUTOFF) return(NULL)
      map_dfr(PK_TIMES, function(t) {
        conc <- if (t == 0) 0.0 else {
          cv <- (PK_DOSE * 1000 * ka_i / (vd_i * (ka_i - ke_i))) * (exp(-ke_i * t) - exp(-ka_i * t))
          max(0.0, cv) * exp(rnorm(1, 0, 0.12)) + abs(rnorm(1, 0, 2.0))
        }
        tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = "TREATMENT",
          TUMORTYPE = s$TUMORTYPE,
          PARAMCD = "CONC", PARAM = paste(DRUG_NAME, "Plasma Concentration"),
          AVISIT = PK_VISITS[vi], AVISITN = vi, NOMTPT = t,
          ADT = fmt_date(bdt + t / 24), AVAL = round(conc, 2), AVALU = "ng/mL",
          DOSE = PK_DOSE, ROUTE = "ORAL", ANL01FL = "Y",
          AGE = s$AGE, SEX = s$SEX
        )
      })
    })

    cmax  <- max(0.1, (PK_DOSE * 1000 * ka_i / (vd_i * (ka_i - ke_i))) *
                 (exp(-ke_i * 3) - exp(-ka_i * 3)) * exp(rnorm(1, 0, 0.10)))
    auc   <- max(0.1, PK_DOSE * 1000 / cl_i * exp(rnorm(1, 0, 0.10)))
    thalf <- 0.693 * vd_i / cl_i

    param_rows <- tibble(
      STUDYID = STUDY_ID, USUBJID = uid, ARM = "TREATMENT",
      TUMORTYPE = s$TUMORTYPE,
      PARAMCD = c("CMAX", "AUCINF", "TMAX", "THALF"),
      PARAM   = paste(DRUG_NAME, c("Cmax", "AUCinf", "Tmax", "t1/2")),
      AVISIT = "CYCLE 1 DAY 1", AVISITN = 1L, NOMTPT = NA_real_,
      ADT = fmt_date(trtsdt),
      AVAL = round(c(cmax, auc, 3.0 + rnorm(1, 0, 0.5), thalf), 1),
      AVALU = c("ng/mL", "ng*h/mL", "h", "h"),
      DOSE = PK_DOSE, ROUTE = "ORAL", ANL01FL = "Y",
      AGE = s$AGE, SEX = s$SEX
    )

    bind_rows(conc_rows, param_rows)
  })

  log_msg(sprintf("  ADPK: %d records | Cmax median=%.0f ng/mL",
                  nrow(adpk), median(adpk$AVAL[adpk$PARAMCD == "CMAX"])))

  log_msg("Generating ADEX...")

  adex <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s       <- adsl[i, ]
    trtsdt  <- as.Date(s$TRTSDT)
    trtedt  <- as.Date(s$TRTEDT)
    dur     <- max(1L, as.integer(trtedt - trtsdt))
    arm     <- s$ARM
    uid     <- s$USUBJID
    planned <- ifelse(arm == "TREATMENT", 300L, 0L)
    level   <- 1.0
    pre_int <- 1.0
    interrupted <- FALSE
    cum_dose <- 0.0
    rows <- list()

    for (cyc in seq_len(max(1L, dur %/% 21L))) {
      cs <- trtsdt + (cyc - 1L) * 21L
      ce <- cs + 20L
      if (ce > CUTOFF) break
      mod_reason <- ""
      if (cyc > 2L && arm == "TREATMENT") {
        if (interrupted) {
          level       <- pre_int
          interrupted <- FALSE
          mod_reason  <- "AE - Dose Resumption"
        } else if (runif(1) < 0.10) {
          level      <- max(0.33, level - 0.33)
          mod_reason <- sample(c("AE - Hepatotoxicity", "AE - Diarrhea", "AE - Fatigue", "AE - Neuropathy"), 1L)
        } else if (runif(1) < 0.04) {
          pre_int     <- level
          level       <- 0.0
          interrupted <- TRUE
          mod_reason  <- "AE - Dose Interruption"
        } else if (level < 1.0 && runif(1) < 0.15) {
          level      <- min(1.0, level + 0.33)
          mod_reason <- "Dose Re-escalation"
        }
      }
      actual   <- round(planned * level)
      days_on  <- ifelse(level > 0, 21L, sample(0:10, 1L))
      intensity <- ifelse(planned > 0, round(actual * days_on / (planned * 21) * 100, 1), 0.0)
      cum_dose  <- cum_dose + actual * days_on

      rows[[cyc]] <- tibble(
        STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = s$TUMORTYPE,
        EXSEQ = cyc, EXTRT = ifelse(arm == "TREATMENT", DRUG_NAME, "Placebo"),
        EXDOSE = actual, EXDOSU = "mg", EXDOSFRQ = "QD", EXROUTE = "ORAL",
        EXSTDTC = fmt_date(cs), EXENDTC = fmt_date(ce),
        EXCYCLE = cyc, PLANDOSE = planned,
        DOSELEVEL = level, DAYSONDRUG = days_on,
        DOSEINT = intensity, CUMDOSE = round(cum_dose, 0),
        MODFL = ifelse(nchar(mod_reason) > 0, "Y", "N"),
        MODREASN = mod_reason
      )
    }
    bind_rows(rows)
  })

  log_msg(sprintf("  ADEX: %d records", nrow(adex)))

  log_msg("Generating ADBM...")

  BM_PARAMS <- tribble(
    ~PARAMCD, ~PARAM,                            ~AVALU,       ~PATTERN,       ~TRT_CHG,
    "TMB",   "Tumor Mutational Burden",           "mut/Mb",     "baseline",     NA_real_,
    "PDL1",  "PD-L1 TPS Score",                  "%",          "baseline",     NA_real_,
    "MSI",   "MSI Sensor Score",                  "score",      "baseline",     NA_real_,
    "CTDNA", "ctDNA Variant Allele Frequency",    "%",          "longitudinal", -0.55,
    "CEA",   "Carcinoembryonic Antigen",           "ng/mL",     "longitudinal", -0.30,
    "CA125", "CA-125",                             "U/mL",      "longitudinal", -0.28,
    "CD8",   "CD8+ T Cell Density",               "cells/mm2",  "longitudinal",  0.42,
    "CD4",   "CD4+ T Cell Density",               "cells/mm2",  "longitudinal",  0.28,
    "NK",    "NK Cell Density",                   "cells/mm2",  "longitudinal",  0.22,
    "TREG",  "Regulatory T Cell Density",         "cells/mm2",  "longitudinal",  0.18,
    "PDL1T", "PD-L1 on Tumor Cells (H-score)",    "H-score",    "longitudinal", -0.20,
    "IFNg",  "IFN-gamma (serum)",                  "pg/mL",     "longitudinal",  0.35
  )

  BM_VISITS <- c("BASELINE", "CYCLE 2 DAY 1", "CYCLE 4 DAY 1", "CYCLE 6 DAY 1", "EOT")
  BM_VDAYS  <- c(0L, 42L, 126L, 210L, NA_integer_)

  RESP_MOD_MAP <- c(CR = 1.3, PR = 1.0, SD = 0.5, PD = -0.2)

  adbm <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    trtsdt <- as.Date(s$TRTSDT)
    trtedt <- as.Date(s$TRTEDT)
    dur    <- as.integer(trtedt - trtsdt)
    arm    <- s$ARM
    uid    <- s$USUBJID
    bor    <- s$BESTRSPC

    map_dfr(seq_len(nrow(BM_PARAMS)), function(j) {
      bm <- BM_PARAMS[j, ]
      base_val <- switch(bm$PARAMCD,
        TMB  = s$TMB,
        PDL1 = s$PDL1SCORE,
        MSI  = s$MSI_SCORE,
        round(rlnorm(1, 2.5, 0.8), 2)
      )
      visits <- if (bm$PATTERN == "longitudinal") BM_VISITS else "BASELINE"
      vdays  <- if (bm$PATTERN == "longitudinal") BM_VDAYS  else 0L

      map_dfr(seq_along(visits), function(vi) {
        vday <- ifelse(is.na(vdays[vi]), dur, vdays[vi])
        if (vday > dur + 14L) return(NULL)
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)

        aval <- if (visits[vi] == "BASELINE" || bm$PATTERN == "baseline") base_val else {
          phase    <- min(vday / max(dur, 1), 1.0)
          resp_mod <- RESP_MOD_MAP[bor]
          if (!is.na(bm$TRT_CHG) && arm == "TREATMENT")
            base_val * (1 + bm$TRT_CHG * phase * resp_mod) * exp(rnorm(1, 0, 0.15))
          else
            base_val * exp(rnorm(1, 0, 0.08))
        }
        aval <- max(0.01, aval)
        tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = s$TUMORTYPE,
          PARAMCD = bm$PARAMCD, PARAM = bm$PARAM, AVALU = bm$AVALU,
          AVISIT = visits[vi],
          AVISITN = which(BM_VISITS == visits[vi])[1] - 1L,
          ADT = fmt_date(adt), ADTN = vday,
          BASE = round(base_val, 3), AVAL = round(aval, 3),
          CHG  = round(aval - base_val, 3),
          PCHG = round((aval - base_val) / base_val * 100, 1),
          BESTRSPC = bor, ANL01FL = "Y"
        )
      })
    })
  })

  log_msg(sprintf("  ADBM: %d records | %d parameters", nrow(adbm), n_distinct(adbm$PARAMCD)))

  log_msg("Generating ADPR...")

  PRO_SCALES <- tribble(
    ~PARAMCD,   ~PARAM,                     ~SCALETYP,     ~MU,  ~SD,  ~DIR,  ~EFFECT,
    "GLQOL",   "Global Health Status/QoL",  "functional",   62,   18,    1,    4.0,
    "PHFUNC",  "Physical Functioning",       "functional",   68,   16,   -1,    3.5,
    "EMOFUNC", "Emotional Functioning",      "functional",   72,   17,    1,    2.5,
    "COGFUNC", "Cognitive Functioning",      "functional",   78,   14,   -1,    2.0,
    "SOCFUNC", "Social Functioning",         "functional",   70,   18,   -1,    3.0,
    "FATIGUE",  "Fatigue",                   "symptom",      32,   22,    1,    8.0,
    "NAUSEA",   "Nausea and Vomiting",       "symptom",      14,   18,    1,    6.0,
    "PAIN",     "Pain",                      "symptom",      28,   22,    1,    4.0,
    "DYSPNEA",  "Dyspnoea",                  "symptom",      20,   18,    1,    3.5,
    "APPETIT",  "Appetite Loss",             "symptom",      22,   20,    1,    5.0,
    "DIARRH",   "Diarrhoea",                 "symptom",      10,   14,    1,    7.0
  )

  PR_VISITS <- c("BASELINE", "CYCLE 2", "CYCLE 4", "CYCLE 6", "CYCLE 9", "CYCLE 12", "EOT")
  PR_VDAYS  <- c(0L, 42L, 126L, 210L, 336L, 462L, NA_integer_)
  MID       <- 10.0

  adpr <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s      <- adsl[i, ]
    trtsdt <- as.Date(s$TRTSDT)
    trtedt <- as.Date(s$TRTEDT)
    dur    <- as.integer(trtedt - trtsdt)
    arm    <- s$ARM
    uid    <- s$USUBJID
    comp   <- rbeta(1, 8, 2)

    map_dfr(seq_len(nrow(PRO_SCALES)), function(j) {
      sc       <- PRO_SCALES[j, ]
      base_val <- pmin(100, pmax(0, rnorm(1, sc$MU, sc$SD)))

      map_dfr(seq_along(PR_VISITS), function(vi) {
        if (vi > 1L && runif(1) > comp) return(NULL)
        vday <- ifelse(is.na(PR_VDAYS[vi]), dur, PR_VDAYS[vi])
        if (vday > dur + 21L) return(NULL)
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)

        aval <- if (vi == 1L) base_val else {
          phase <- min(vday / max(dur, 1), 1.0)
          peak  <- 0.25
          tox   <- sc$EFFECT * sc$DIR * exp(-5 * (phase - peak)^2)
          rec   <- if (phase > peak) sc$EFFECT * sc$DIR * 0.35 * (phase - peak) else 0
          if (arm == "TREATMENT")
            pmin(100, pmax(0, base_val + tox + rec + rnorm(1, 0, sc$SD * 0.18)))
          else
            pmin(100, pmax(0, base_val - 3 * phase + rnorm(1, 0, sc$SD * 0.12)))
        }
        chg <- round(aval - base_val, 1)
        tibble(
          STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = s$TUMORTYPE,
          PARAMCD = sc$PARAMCD, PARAM = sc$PARAM, SCALETYP = sc$SCALETYP,
          AVISIT = PR_VISITS[vi], AVISITN = vi - 1L,
          ADT = fmt_date(adt), ADTN = vday,
          BASE = round(base_val, 1), AVAL = round(aval, 1),
          CHG = chg, PCHG = round(chg / base_val * 100, 1),
          MID = MID,
          MIDRESP = ifelse(abs(chg) >= MID, "Y", "N"),
          DETRFL  = ifelse(
            (sc$SCALETYP == "functional" && chg <= -MID) ||
            (sc$SCALETYP == "symptom"    && chg >= MID), "Y", "N"),
          ANL01FL = "Y"
        )
      })
    })
  })

  log_msg(sprintf("  ADPR: %d records | %d scales", nrow(adpr), n_distinct(adpr$PARAMCD)))

  log_msg("Generating ADMUT...")

  GENE_PROFILE <- list(
    TP53   = list(chr = "17p13.1", types = c("Missense","Nonsense","Frameshift","Splice"),       wts = c(0.65,0.18,0.10,0.07), domains = c("DNA-binding domain","Tetramerization domain"), hotspots = c(175,245,248,249,273,282)),
    KRAS   = list(chr = "12p12.1", types = c("Missense","Amplification"),                        wts = c(0.92,0.08),           domains = c("G-domain"),                                    hotspots = c(12,13,61)),
    EGFR   = list(chr = "7p11.2",  types = c("Missense","In-frame deletion","Amplification"),   wts = c(0.45,0.35,0.20),      domains = c("Kinase domain","Extracellular domain"),         hotspots = c(746,747,748,790,858)),
    STK11  = list(chr = "19p13.3", types = c("Missense","Nonsense","Frameshift","Splice"),       wts = c(0.45,0.28,0.18,0.09), domains = c("Kinase domain"),                               hotspots = c(67,176,281)),
    KEAP1  = list(chr = "19p13.2", types = c("Missense","Nonsense","Frameshift"),                wts = c(0.55,0.25,0.20),      domains = c("Kelch domain","BTB domain"),                   hotspots = c(321,334,414)),
    BRAF   = list(chr = "7q34",    types = c("Missense"),                                        wts = c(1.00),                domains = c("Kinase domain"),                               hotspots = c(600)),
    MET    = list(chr = "7q31.2",  types = c("Amplification","Missense","Splice"),               wts = c(0.55,0.30,0.15),      domains = c("Kinase domain","Sema domain"),                 hotspots = c(14)),
    PIK3CA = list(chr = "3q26.32", types = c("Missense"),                                        wts = c(1.00),                domains = c("Kinase domain","Helical domain"),              hotspots = c(542,545,1047)),
    PTEN   = list(chr = "10q23.31",types = c("Missense","Nonsense","Frameshift","Splice"),       wts = c(0.40,0.28,0.22,0.10), domains = c("Phosphatase domain","C2 domain"),              hotspots = c(130,233)),
    CDKN2A = list(chr = "9p21.3",  types = c("Missense","Nonsense","Frameshift","Deletion"),    wts = c(0.40,0.25,0.20,0.15), domains = c("ANK repeat domain"),                          hotspots = c(58,83,122)),
    RB1    = list(chr = "13q14.2", types = c("Nonsense","Frameshift","Splice","Missense"),       wts = c(0.35,0.30,0.20,0.15), domains = c("Pocket domain A","Pocket domain B"),           hotspots = c(579,661)),
    SMAD4  = list(chr = "18q21.2", types = c("Missense","Nonsense","Frameshift"),                wts = c(0.50,0.30,0.20),      domains = c("MH1 domain","MH2 domain"),                    hotspots = c(361,406)),
    ARID1A = list(chr = "1p36.11", types = c("Nonsense","Frameshift","Missense"),                wts = c(0.40,0.35,0.25),      domains = c("ARID domain"),                                hotspots = c(1721,1760)),
    NF1    = list(chr = "17q11.2", types = c("Nonsense","Frameshift","Missense","Splice"),       wts = c(0.35,0.30,0.22,0.13), domains = c("GRD domain","CSRD domain"),                   hotspots = c(843,1423)),
    RET    = list(chr = "10q11.21",types = c("Fusion","Missense"),                               wts = c(0.70,0.30),           domains = c("Kinase domain"),                              hotspots = c(918))
  )

  MUT_CLASS_MAP <- c(
    Missense = "Missense_Mutation", Nonsense = "Nonsense_Mutation",
    Frameshift = "Frame_Shift_Del", Splice = "Splice_Site",
    "In-frame deletion" = "In_Frame_Del", Amplification = "Amplification",
    Deletion = "Deletion", Fusion = "Fusion"
  )
  AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

  admut <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s   <- adsl[i, ]
    uid <- s$USUBJID
    arm <- s$ARM
    tt  <- s$TUMORTYPE
    mp  <- MUTATION_PREV[[tt]]

    map_dfr(seq_along(GENE_LIST), function(gi) {
      if (runif(1) > mp[gi]) return(NULL)
      gene      <- GENE_LIST[gi]
      gp        <- GENE_PROFILE[[gene]]
      mut_type  <- sample(gp$types, 1L, prob = gp$wts)
      mut_class <- MUT_CLASS_MAP[mut_type]
      pos       <- if (runif(1) < 0.70 && length(gp$hotspots) > 0)
                     max(1L, sample(gp$hotspots, 1L) + sample(-3:3, 1L))
                   else
                     sample(1:1200, 1L)
      ref_aa <- sample(AA, 1L)
      alt_aa <- sample(AA[AA != ref_aa], 1L)
      hgvsp  <- if (mut_type == "Missense")       sprintf("p.%s%d%s", ref_aa, pos, alt_aa) else
                if (mut_type == "Nonsense")        sprintf("p.%s%d*", ref_aa, pos) else
                if (grepl("Frameshift", mut_type)) sprintf("p.%s%dfs", ref_aa, pos) else
                sprintf("%s %s", gene, mut_type)

      clonal  <- runif(1) < 0.60
      vaf_lo  <- ifelse(clonal, 0.35, 0.05)
      vaf_hi  <- ifelse(clonal, 0.55, 0.25)
      vaf     <- clamp(runif(1, vaf_lo, vaf_hi) + rnorm(1, 0, 0.03), 0.01, 0.85)
      depth   <- sample(80:400, 1L)
      alt_ct  <- round(vaf * depth)
      impact  <- if (mut_type == "Missense" && pos %in% gp$hotspots) "HIGH" else
                 if (mut_type %in% c("Nonsense", "Frameshift", "Splice"))  "HIGH" else
                 if (mut_type == "Amplification") "MODERATE" else
                 sample(c("HIGH", "MODERATE", "LOW"), 1L, prob = c(0.30, 0.45, 0.25))

      tibble(
        STUDYID = STUDY_ID, USUBJID = uid, ARM = arm, TUMORTYPE = tt,
        HUGO_SYMBOL = gene, CHROMOSOME = gp$chr,
        VARIANT_CLASS = mut_class, VARIANT_TYPE = mut_type,
        HGVSP = hgvsp, PROTEIN_POS = pos,
        PROTEIN_DOMAIN = sample(gp$domains, 1L),
        REF_ALLELE = ref_aa, ALT_ALLELE = alt_aa,
        VAF = round(vaf, 3), T_DEPTH = depth,
        T_ALT_COUNT = alt_ct, T_REF_COUNT = depth - alt_ct,
        CLONAL = ifelse(clonal, "Y", "N"), IMPACT = impact,
        BESTRSPC = s$BESTRSPC, TMB = s$TMB, TMBHIGH = s$TMBHIGH,
        ANL01FL = "Y"
      )
    })
  })

  log_msg(sprintf("  ADMUT: %d records | %d genes | %d patients",
                  nrow(admut), n_distinct(admut$HUGO_SYMBOL), n_distinct(admut$USUBJID)))

  log_msg("Generating ADSIG...")

  SBS_PROFILES <- list(
    NSCLC = list(
      list(sig = "SBS4",  desc = "Tobacco smoking",          w = 0.35),
      list(sig = "SBS2",  desc = "APOBEC activity",          w = 0.20),
      list(sig = "SBS13", desc = "APOBEC activity",          w = 0.15),
      list(sig = "SBS1",  desc = "Spontaneous deamination",  w = 0.18),
      list(sig = "SBS5",  desc = "Unknown clock-like",       w = 0.12)
    ),
    CRC = list(
      list(sig = "SBS1",  desc = "Spontaneous deamination",  w = 0.30),
      list(sig = "SBS5",  desc = "Unknown clock-like",       w = 0.25),
      list(sig = "SBS15", desc = "Defective MMR",            w = 0.20),
      list(sig = "SBS6",  desc = "Defective MMR",            w = 0.15),
      list(sig = "SBS44", desc = "Defective MMR",            w = 0.10)
    ),
    HCC = list(
      list(sig = "SBS4",  desc = "Tobacco smoking",          w = 0.25),
      list(sig = "SBS22", desc = "Aristolochic acid",        w = 0.20),
      list(sig = "SBS24", desc = "Aflatoxin exposure",       w = 0.20),
      list(sig = "SBS1",  desc = "Spontaneous deamination",  w = 0.20),
      list(sig = "SBS5",  desc = "Unknown clock-like",       w = 0.15)
    ),
    PDAC = list(
      list(sig = "SBS1",  desc = "Spontaneous deamination",  w = 0.35),
      list(sig = "SBS5",  desc = "Unknown clock-like",       w = 0.30),
      list(sig = "SBS2",  desc = "APOBEC activity",          w = 0.15),
      list(sig = "SBS3",  desc = "HR deficiency",            w = 0.12),
      list(sig = "SBS18", desc = "ROS damage",               w = 0.08)
    ),
    BRCA = list(
      list(sig = "SBS3",  desc = "HR deficiency / BRCA1/2",  w = 0.30),
      list(sig = "SBS2",  desc = "APOBEC activity",          w = 0.25),
      list(sig = "SBS13", desc = "APOBEC activity",          w = 0.20),
      list(sig = "SBS1",  desc = "Spontaneous deamination",  w = 0.15),
      list(sig = "SBS5",  desc = "Unknown clock-like",       w = 0.10)
    )
  )

  adsig <- map_dfr(seq_len(nrow(adsl)), function(i) {
    s    <- adsl[i, ]
    uid  <- s$USUBJID
    tt   <- s$TUMORTYPE
    sigs <- SBS_PROFILES[[tt]]
    total_muts <- max(10L, as.integer(s$TMB * 50))

    raw_wts  <- sapply(sigs, function(sg) max(0.01, sg$w + rnorm(1, 0, 0.05)))
    norm_wts <- raw_wts / sum(raw_wts)

    map_dfr(seq_along(sigs), function(si) {
      sg <- sigs[[si]]
      tibble(
        STUDYID   = STUDY_ID, USUBJID = uid,
        ARM       = s$ARM, TUMORTYPE = tt,
        SIG_NAME  = sg$sig, SIG_DESC = sg$desc,
        SIG_WEIGHT = round(norm_wts[si], 4),
        N_MUTS    = round(total_muts * norm_wts[si]),
        TMB       = s$TMB, TMBHIGH = s$TMBHIGH,
        MSISTS    = s$MSISTS, BESTRSPC = s$BESTRSPC,
        ANL01FL   = "Y"
      )
    })
  })

  log_msg(sprintf("  ADSIG: %d records | %d signatures", nrow(adsig), n_distinct(adsig$SIG_NAME)))

  log_msg("\nSaving datasets...")

  datasets <- list(
    ADSL  = adsl,    ADRS  = adrs,  ADTR  = adtr,
    ADAE  = adae_rows, ADLB = adlb, ADTTE = adtte,
    ADPK  = adpk,   ADEX  = adex,  ADBM  = adbm,
    ADPR  = adpr,   ADMUT = admut, ADRAND = adrand_rows,
    ADSIG = adsig
  )

  for (nm in names(datasets)) {
    path <- file.path(output_dir, paste0(nm, ".csv"))
    write.csv(datasets[[nm]], path, row.names = FALSE)
    log_msg(sprintf("  %s: %s rows x %s cols", nm, nrow(datasets[[nm]]), ncol(datasets[[nm]])))
  }

  total <- sum(sapply(datasets, nrow))
  log_msg(sprintf("\nTotal records : %s", format(total, big.mark = ",")))
  log_msg(sprintf("Patients      : %d", N))
  log_msg(sprintf("Study         : %s | %s | Phase II/III Basket Trial", STUDY_ID, DRUG_NAME))
  log_msg(sprintf("Seed          : %d", seed))
  log_msg(sprintf("Output dir    : %s", normalizePath(output_dir, mustWork = FALSE)))

  invisible(datasets)
}

datasets <- generate_adam_oncviz(output_dir = "./data")
