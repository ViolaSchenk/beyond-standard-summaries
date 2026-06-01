# ============================================================
# MASTER THESIS — ANALYTICAL PIPELINE
# Adverse Event Reporting and Information Loss in Clinical Trials
# Boehringer Ingelheim / Otto-Friedrich Universität Bamberg
# Author : Viola Schenk
# Date   : XXXXX
# Purpose:
#   This script implements the full quantitative pipeline for the thesis.
#   It processes CDISC ADaM safety data (ADSL/ADAE), derives treatment-emergent
#   adverse events (TEAEs), and applies a series of reporting thresholds to
#   quantify how information is lost as thresholds increase.
# Pipeline stages:
#   0  Utility functions (shared helpers used throughout)
#   1  Data loading      (pharmaverseadam preferred; admiral.test fallback)
#   2  ADSL cleaning     (treatment arm, first-dose date)
#   3  ADAE cleaning     (PT/SOC/severity standardization, date resolution)
#   3a Severity mapping  (AESEV string → CTCAE-like integer grade)
#   4  TEAE derivation   (merge with ADSL; flag treatment-emergent events)
#   5  Final dataset     (study-day calculation, deduplication)
#   6  CTCAE ≥3 baseline (grade flags; conditional + unconditional denominators)
#   7  Threshold rules   (pure any-arm; keep-≥G3 safeguard)
#   8  Information-loss metrics
#        Shannon entropy (H_event, H_subj)
#        Jensen–Shannon divergence (JSD_event, JSD_subj)
#        Coverage (events, subjects, PTs) normalized to 1% baseline
#        Grade ≥3 clinical visibility (unconditional denominator)
#   9  Multivariate methods
#        Correspondence Analysis (CA) — PT × arm contingency table
#        Principal Component Analysis (PCA) — binary subject × PT matrix
#  10  Survey visualizations
#        Butterfly plot, Volcano plot, Temporal cumulative-incidence plot,
#        Jaccard co-occurrence heatmap, CA strip (one-axis biplot)
# Data:
#   pharmaverseadam::adsl / adae  (CDISC synthetic ADaM; 3 arms:
#   Placebo, Xanomeline Low Dose, Xanomeline High Dose)
# Key outputs:
#   summary_table   — information-loss metrics across thresholds × variants
#   ge3_comparison  — conditional vs unconditional ≥G3 visibility table
#   Figures         — saved via ggsave() or rendered in the RStudio viewer
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)       # data manipulation
  library(tidyr)       # pivoting and reshaping
  library(purrr)       # functional iteration (map_dfr etc.)
  library(stringr)     # string cleaning and wrapping
  library(lubridate)   # date parsing
  library(janitor)     # clean_names() for snake_case columns
  library(ggplot2)     # plotting
  library(scales)      # axis formatting (percent, comma)
  library(forcats)     # factor reordering
  library(ggh4x)       # extended ggplot2 faceting (butterfly plot)
  library(ca)          # correspondence analysis
  library(ggrepel)     # label repulsion (temporal, volcano, CA plots)
})

# ------------------------------------------------------------
# 0) Utility functions
#    These are shared helpers used across all pipeline stages.
#    Defining them here avoids repetition and ensures consistent
#    behavior throughout (e.g. same string normalization everywhere).
# ------------------------------------------------------------

# Infix null-coalescing operator: returns 'a' if non-NULL, else 'b'.
# Used throughout to provide sensible defaults without overwriting valid values.
# Note: unlike dplyr::coalesce(), this only triggers on NULL — not on NA.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Print a clearly visible section header to the R console.
# Useful for tracking pipeline progress in long output logs and
# for legibility when sourcing the script in a thesis appendix.
print_header <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n",
      "= ", txt, " =\n",
      paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
}

# Return the first variable name from 'choices' that is present in 'nms'.
# Used to handle variable name differences across ADaM dataset versions
# (e.g. TRT01A vs TRTA for treatment arm). Returns NA_character_ if none found.
first_present <- function(nms, choices) {
  hit <- intersect(choices, nms)
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

# Normalize a character vector: trim whitespace, convert to upper case,
# replace blank strings with NA. Applied to PT, SOC, and severity strings
# to prevent spurious duplicate categories from casing or spacing differences.
norm_str <- function(x) {
  x %>% as.character() %>% str_squish() %>% toupper() %>% na_if("")
}

# Resolve a date from multiple candidate columns into a single Date vector.
# Tries columns in order and takes the first non-missing value per row.
# Handles Date, POSIXt, and character inputs; attempts ymd_hms then ymd parsing.
# This is necessary because ADaM date variables differ between dataset versions
# (e.g. ASTDT vs ASTDTM vs AESTDTC).
coalesce_date <- function(df, candidates) {
  out <- rep(NA_Date_, nrow(df))
  for (nm in candidates) {
    if (!nm %in% names(df)) next
    x <- df[[nm]]
    y <- if (inherits(x, "Date")) x else if (inherits(x, "POSIXt")) as_date(x) else as.character(x)
    y1 <- suppressWarnings(lubridate::ymd_hms(y, quiet = TRUE))
    y2 <- suppressWarnings(lubridate::ymd(y, quiet = TRUE))
    y_fin <- dplyr::coalesce(as_date(y1), as_date(y2))
    idx <- is.na(out) & !is.na(y_fin)
    out[idx] <- y_fin[idx]
  }
  out
}

# Compute subject-level denominators (n per arm) from an AE table.
# Returns one row per treatment arm with column n_trt.
# NOTE: This counts subjects who APPEAR in the AE table, not all randomized
# subjects. For unconditional denominators (all treated), use denom_all_treated.
denom_from_ae <- function(ae_tbl, trt_col = "trt") {
  ae_tbl %>%
    distinct(.data[[trt_col]], usubjid) %>%
    count(.data[[trt_col]], name = "n_trt") %>%
    rename(trt = .data[[trt_col]])
}

# ------------------------------------------------------------
# 1) Load ADSL / ADAE
#    Data source: pharmaverseadam (preferred) — a CDISC-compliant synthetic
#    ADaM dataset distributed via the pharmaverse R ecosystem. It contains
#    three treatment arms: Placebo, Xanomeline Low Dose, Xanomeline High Dose.
#    Fallback: admiral.test — an earlier synthetic ADaM dataset with the same
#    general structure but fewer variables. Used only when pharmaverseadam is
#    unavailable, e.g. in restricted computing environments.
#    Both datasets are in-memory R objects; no external file I/O is required.
# ------------------------------------------------------------
if (requireNamespace("pharmaverseadam", quietly = TRUE)) {
  adsl <- pharmaverseadam::adsl   # subject-level: demographics, treatment assignment
  adae <- pharmaverseadam::adae   # event-level: one row per AE per subject
  message("✅ Using pharmaverseadam data")
} else if (requireNamespace("admiral.test", quietly = TRUE)) {
  adsl <- admiral.test::adsl
  adae <- admiral.test::adae
  message("✅ Using admiral.test data (fallback)")
} else {
  stop("❌ Install 'pharmaverseadam' (preferred) or 'admiral.test' (fallback).")
}

print_header("ROW COUNTS")
tibble(n_adsl = nrow(adsl), n_adae = nrow(adae)) %>% print(n = Inf)

# ------------------------------------------------------------
# 2) Clean ADSL: extract treatment arm and first dose date
#    ADSL provides one row per subject. We extract:
#      trt    — actual treatment arm (from TRT01A / TRTA / TRT01P)
#      trtsdt — first treatment start date, used as TEAE anchor
#    Different ADaM versions use different variable names for treatment arm;
#    first_present() selects whichever is available.
# ------------------------------------------------------------

# Treatment arm variable name differs across ADaM versions; pick the first available.
trt_var <- first_present(names(adsl), c("TRT01A", "TRTA", "TRT01P"))
stopifnot(!is.na(trt_var))

# Resolve first-dose date from multiple candidate columns (handles POSIX, Date, character).
# TRTSDT is the critical anchor: TEAEs are defined as AEs starting on or after this date.
adsl$TRTSDT <- coalesce_date(adsl, c("TRTSDT", "TRTSDTM", "TRSDTM", "TRSDT"))

adsl1 <- adsl %>%
  transmute(
    USUBJID,
    TRT = as.character(.data[[trt_var]]),
    TRTSDT
  ) %>%
  clean_names()   # converts to snake_case: USUBJID → usubjid, etc.

print_header("ADSL SUMMARY")
adsl1 %>%
  summarise(
    n_subj           = n_distinct(usubjid),
    n_missing_trtsdt = sum(is.na(trtsdt))
  ) %>% print(n = Inf)

# ------------------------------------------------------------
# 3) Clean ADAE: standardize PT, SOC, severity, and dates
#    ADAE provides one row per adverse event per subject.
#    Key derivations:
#      AEDECOD  — Preferred Term (MedDRA), normalized to upper case
#      AEBODSYS — System Organ Class (MedDRA), normalized
#      AESEV    — severity string (MILD/MODERATE/SEVERE etc.), normalized
#      ASTDT    — AE start date (resolved from multiple candidate columns)
#      AEENDT   — AE end date (resolved similarly)
#      TEAEFL   — treatment-emergent flag (used if present; derived later otherwise)
#      AESEQ    — sequence number (for deduplication / traceability)
# ------------------------------------------------------------

# Identify which variable names are present in this version of ADAE
aedecod_var  <- first_present(names(adae), c("AEDECOD"))
aebodsys_var <- first_present(names(adae), c("AEBODSYS"))
aesev_var    <- first_present(names(adae), c("AESEV"))
teae_var     <- first_present(names(adae), c("TEAEFL"))

# Resolve AE start and end dates from multiple possible column formats
adae$ASTDT  <- coalesce_date(adae, c("ASTDT", "ASTDTM", "AESTDTC"))
adae$AEENDT <- coalesce_date(adae, c("AEENDT", "AENDTM", "AEENDTC"))

# Use AESEQ where available for deduplication and record traceability
aseq <- if ("AESEQ" %in% names(adae)) adae$AESEQ else if ("ASEQ" %in% names(adae)) adae$ASEQ else NA_real_

adae1 <- adae %>%
  transmute(
    USUBJID,
    AEDECOD   = if (!is.na(aedecod_var))  .data[[aedecod_var]]  else NA_character_,
    AEBODSYS  = if (!is.na(aebodsys_var)) .data[[aebodsys_var]] else NA_character_,
    AESEV_RAW = if (!is.na(aesev_var))    .data[[aesev_var]]    else NA_character_,
    ASTDT, AEENDT,
    TEAEFL    = if (!is.na(teae_var)) .data[[teae_var]] else NA_character_,
    AESEQ     = aseq
  ) %>%
  mutate(
    # Normalize text fields to prevent spurious category splits
    # (e.g. "Mild" and "MILD" would otherwise be treated as separate levels)
    AEDECOD   = norm_str(AEDECOD),
    AEBODSYS  = norm_str(AEBODSYS),
    AESEV_RAW = norm_str(AESEV_RAW)
  ) %>%
  clean_names()

print_header("UNIQUE SEVERITY STRINGS (NORMALIZED)")
adae1 %>% count(aesev_raw, sort = TRUE) %>% print(n = Inf)

# ------------------------------------------------------------
# 3a) Map AESEV string → CTCAE-like integer grade
#    The pharmaverseadam dataset stores severity as a free-text string
#    (MILD, MODERATE, SEVERE) rather than a numeric CTCAE grade.
#    We map these to integers 1–5 following CTCAE v5.0 conventions:
#      1 = Mild       2 = Moderate     3 = Severe
#      4 = Life-threatening             5 = Fatal
#    Unmapped strings are left as NA (conservative: not assumed mild).
#    This grade is used as a CTCAE-like proxy throughout the pipeline,
#    particularly for the clinical visibility metric (Grade ≥3 retention).
# ------------------------------------------------------------

# Lookup table: integer grade → vector of accepted severity strings (normalized)
sev_map <- list(
  `1` = c("MILD", "GRADE 1", "1", "LOW"),
  `2` = c("MODERATE", "MOD", "GRADE 2", "2", "MEDIUM"),
  `3` = c("SEVERE", "SEV", "GRADE 3", "3"),
  `4` = c("LIFE-THREATENING", "LIFE THREATENING", "VERY SEVERE", "GRADE 4", "4"),
  `5` = c("DEATH", "FATAL", "GRADE 5", "5")
)

# Vectorised mapper: returns integer grade for each severity string, NA if unrecognized
sev_to_grade <- function(x) {
  x <- norm_str(x)
  out <- rep(NA_integer_, length(x))
  for (g in names(sev_map)) out[x %in% sev_map[[g]]] <- as.integer(g)
  out
}

adae2 <- adae1 %>% mutate(grade = sev_to_grade(aesev_raw))

print_header("MAPPING CHECK (AESEV_RAW -> grade)")
adae2 %>% count(aesev_raw, grade, sort = TRUE) %>% print(n = Inf)

# ------------------------------------------------------------
# 4) Merge ADAE with dosed ADSL; derive TEAE flag
#    A treatment-emergent adverse event (TEAE) is an AE that started
#    on or after the date of first treatment administration (TRTSDT).
#    Priority:
#      1. Use TEAEFL from ADAE if present (ADaM-derived value is authoritative)
#      2. Derive from dates: TEAE if ASTDT >= TRTSDT
#    Subjects without a recorded first-dose date are excluded to avoid
#    ambiguous TEAE classification (n_missing_trtsdt from step 2 applies here).
# ------------------------------------------------------------

# Restrict to subjects with a valid first-dose date
adsl_dosed <- adsl1 %>% filter(!is.na(trtsdt))

ae0 <- adae2 %>%
  inner_join(adsl_dosed, by = "usubjid") %>%
  mutate(
    teaefl = case_when(
      !is.na(teaefl)                              ~ teaefl,          # prefer ADaM flag
      !is.na(astdt) & !is.na(trtsdt) & astdt >= trtsdt ~ "Y",       # date-based derivation
      TRUE                                        ~ "N"
    )
  )

# Keep only treatment-emergent events
ae_teae <- ae0 %>% filter(teaefl == "Y")

# ------------------------------------------------------------
# 5) Final TEAE dataset: study-day calculation and deduplication
#    Study Day convention: Day 1 = TRTSDT (i.e. onset_day = ASTDT - TRTSDT + 1).
#    This is standard in clinical trial reporting and ensures Day 1 is not Day 0.
#    Negative onset_day and duration_days values are set to NA.
#    These arise from partial or imputed dates and would distort temporal plots.
#    Deduplication: distinct() across all key columns removes duplicate AE rows
#    that arise from ADAE records being repeated across analysis windows.
#    This is critical to avoid inflating subject-level incidence counts.
# ------------------------------------------------------------

ae_clean <- ae_teae %>%
  mutate(
    # Study Day convention: Day 1 = TRTSDT (so ASTDT - TRTSDT + 1)
    onset_day = if_else(!is.na(astdt) & !is.na(trtsdt),
                        as.integer(astdt - trtsdt) + 1L, NA_integer_),
    duration_days = if_else(!is.na(astdt) & !is.na(aeendt),
                            as.integer(aeendt - astdt) + 1L, NA_integer_)
  ) %>%
  mutate(
    # Guardrails: negative values often mean partial/imputed dates; set to NA to avoid misleading plots
    duration_days = if_else(!is.na(duration_days) & duration_days < 0L, NA_integer_, duration_days),
    onset_day     = if_else(!is.na(onset_day) & onset_day < 0L, NA_integer_, onset_day)
  ) %>%
  mutate(
    # Convert to factors after cleaning/normalization (stable ordering + smaller memory footprint)
    trt      = factor(trt),
    aedecod  = factor(aedecod),
    aebodsys = factor(aebodsys)
  ) %>%
  # Dedupe: ensures downstream incidence counts aren’t inflated by repeated identical records
  distinct(across(any_of(c(
    "usubjid","trt","aedecod","aebodsys","astdt","aeendt",
    "grade","onset_day","duration_days","aesev_raw","teaefl","aeseq"
  )))) %>%
  arrange(trt, usubjid, aedecod, astdt) %>%
  relocate(usubjid, trt, aedecod, aebodsys, astdt, aeendt,
           grade, onset_day, duration_days)

print_header("AE_CLEAN OVERVIEW")
ae_clean %>%
  summarise(
    n          = n(),
    n_subj     = n_distinct(usubjid),
    n_terms    = n_distinct(aedecod),
    miss_grade = sum(is.na(grade)),
    miss_onset = sum(is.na(onset_day))
  ) %>% print(n = Inf)

# ------------------------------------------------------------
# 6) CTCAE ≥3 baseline and denominator definitions
#    ctcae_ge3 is the primary clinical visibility flag used throughout
#    the pipeline. It marks any TEAE with a derived grade of 3 or higher
#    (Severe, Life-threatening, or Fatal).
#    Two denominator definitions are used, reflecting different
#    analytical questions (see thesis Section 3.5.4):
#      denom_teae        — subjects with at least one TEAE record.
#                          Used for CONDITIONAL ≥G3 rates: "given that a
#                          subject has any TEAE, what fraction have ≥G3?"
#      denom_all_treated — all subjects with a valid first-dose date.
#                          Used for UNCONDITIONAL ≥G3 rates: "among all
#                          treated subjects, what fraction have ≥G3?"
#                          This is the preferred denominator for the
#                          clinical visibility metric (H2) as it avoids
#                          the selection bias of conditioning on TEAE presence.
# ------------------------------------------------------------

ae_ctcae <- ae_clean %>%
  mutate(
    # Use derived grade as CTCAE-like grade; define the key safety visibility flag ≥3
    ctcae_grade = grade,
    ctcae_ge3   = !is.na(ctcae_grade) & ctcae_grade >= 3
  )

print_header("ANY GRADE ≥3?")
ae_ctcae %>%
  distinct(usubjid, ctcae_ge3) %>%
  summarise(any_ge3 = any(ctcae_ge3)) %>%
  print(n = Inf)

# Denominators:
# - denom_teae: among subjects with at least one TEAE record
# - denom_all_treated: among all treated subjects with TRTSDT
denom_teae <- ae_ctcae %>% denom_from_ae()
denom_all_treated <- adsl1 %>% filter(!is.na(trtsdt)) %>% count(trt, name = "n_all_treated")

# ≥G3 among TEAE subjects (conditional on having TEAE data)
ge3_by_arm <- ae_ctcae %>%
  filter(ctcae_ge3) %>%
  distinct(trt, usubjid) %>%
  count(trt, name = "n_subj_ge3") %>%
  right_join(denom_teae %>% rename(n_subj_total = n_trt), by = "trt") %>%
  mutate(
    n_subj_ge3  = replace_na(n_subj_ge3, 0L),
    inc_ge3     = n_subj_ge3 / n_subj_total,
    inc_ge3_pct = round(100 * inc_ge3, 1)
  ) %>%
  arrange(desc(inc_ge3))

print_header("≥G3 BY ARM (TEAE subjects)")
ge3_by_arm %>% print(n = Inf)

# Unconditional ≥G3 among all treated (preferred for “visibility” interpretation)
ge3_uncond <- ae_ctcae %>%
  filter(ctcae_ge3) %>%
  distinct(trt, usubjid) %>%
  count(trt, name = "n_subj_ge3") %>%
  right_join(denom_all_treated, by = "trt") %>%
  mutate(
    n_subj_ge3 = replace_na(n_subj_ge3, 0L),
    rate_uncond = 100 * n_subj_ge3 / n_all_treated
  )

print_header("≥G3 BY ARM (UNCONDITIONAL, ALL TREATED)")
ge3_uncond %>% print(n = Inf)

# ------------------------------------------------------------
# 7) Threshold infrastructure
#    Two reporting rules are implemented and compared throughout:
#    (A) PURE any-arm threshold
#        A PT is retained if its subject-level incidence in ANY treatment
#        arm meets or exceeds the cutoff (1%, 2%, 5%, or 10%).
#        This mirrors the most common industry practice.
#    (B) SAFEGUARD: keep-≥G3 rule
#        All PTs that ever reach Grade ≥3 are retained regardless of
#        incidence, in addition to the frequency-retained PTs from (A).
#        This models a clinically conservative reporting strategy that
#        preserves severe-event visibility.
#    Comparing (A) and (B) directly tests Hypothesis H2.
#    THRESHOLDS: 1%, 2%, 5%, 10% — covering the common practical range.
#    The 1% threshold serves as the JSD and coverage baseline.
# ------------------------------------------------------------

# Compute subject-level incidence per PT or SOC, per arm.
# Returns: trt, term, n_subj, n_trt, incidence (proportion 0–1)
incidence_by_term <- function(ae, key = c("aedecod", "aebodsys"), denom_tbl = NULL) {
  key <- match.arg(key)
  denom_tbl <- denom_tbl %||% denom_from_ae(ae)
  ae %>%
    distinct(trt, usubjid, .data[[key]]) %>%
    count(trt, .data[[key]], name = "n_subj") %>%
    left_join(denom_tbl, by = "trt") %>%
    mutate(incidence = ifelse(n_trt > 0, n_subj / n_trt, 0)) %>%
    rename(term = .data[[key]])
}

# (A) Pure any-arm threshold filter.
# Retains all records for PTs whose maximum incidence across arms >= cutoff.
# ge3_only = TRUE restricts the incidence calculation to Grade ≥3 records only
# (used for constructing severe-only termsets in sensitivity analysis).
apply_threshold_any_arm <- function(ae, cutoff, key = c("aedecod", "aebodsys"), ge3_only = FALSE) {
  key <- match.arg(key)
  denom_tbl <- denom_from_ae(ae)
  x <- if (ge3_only) ae %>% filter(ctcae_ge3) else ae

  inc <- incidence_by_term(x, key = key, denom_tbl = denom_tbl)

  keep_terms <- inc %>%
    group_by(term) %>%
    summarise(max_inc = max(incidence, na.rm = TRUE), .groups = "drop") %>%
    filter(max_inc >= cutoff) %>%
    pull(term) %>%
    as.character()

  ae %>%
    filter(as.character(.data[[key]]) %in% keep_terms)
}

# Safeguard: pure threshold UNION all Grade ≥G3 PTs.
# Any PT that ever appears as Grade ≥3 is retained regardless of frequency.
# union() ensures no ≥3 PT is ever suppressed by the frequency rule.
apply_threshold_with_safeguard <- function(ae, cutoff) {
  pt_ge3 <- ae %>%
    filter(ctcae_ge3) %>%
    distinct(aedecod) %>%
    pull(aedecod) %>%
    as.character()

  ae_pure <- apply_threshold_any_arm(ae, cutoff = cutoff, key = "aedecod", ge3_only = FALSE)

  keep_freq <- ae_pure %>%
    distinct(aedecod) %>%
    pull(aedecod) %>%
    as.character()

  keep <- union(keep_freq, pt_ge3)
  ae %>% filter(as.character(aedecod) %in% keep)
}

# ≥G3 Unconditional ≥3 incidence: denominator = all treated subjects.
# This is the primary clinical visibility metric for H2.
ge3_uncond_from_ae <- function(ae_tbl, denom_all_treated_tbl = denom_all_treated) {
  ae_tbl %>%
    filter(ctcae_ge3) %>%
    distinct(trt, usubjid) %>%
    count(trt, name = "n_subj_ge3") %>%
    right_join(denom_all_treated_tbl, by = "trt") %>%
    mutate(
      n_subj_ge3         = replace_na(n_subj_ge3, 0L),
      inc_ge3_uncond     = n_subj_ge3 / n_all_treated,
      inc_ge3_uncond_pct = round(100 * inc_ge3_uncond, 1)
    ) %>%
    arrange(desc(inc_ge3_uncond))
}

# Conditional ≥3 rate (denominator = subjects in the thresholded AE table).
# Used in sensitivity analysis to show the contrast with unconditional rates.
check_ge3_at_conditional <- function(ae, cutoff) {
  ae_thr <- apply_threshold_any_arm(ae, cutoff = cutoff, key = "aedecod", ge3_only = FALSE)
  denom <- ae_thr %>% denom_from_ae() %>% rename(n_subj_total = n_trt)

  ae_thr %>%
    filter(ctcae_ge3) %>%
    distinct(trt, usubjid) %>%
    count(trt, name = "n_subj_ge3") %>%
    right_join(denom, by = "trt") %>%
    mutate(
      n_subj_ge3  = replace_na(n_subj_ge3, 0L),
      inc_ge3     = ifelse(n_subj_total > 0, n_subj_ge3 / n_subj_total, 0),
      inc_ge3_pct = round(100 * inc_ge3, 1),
      threshold   = paste0(round(cutoff * 100), "%")
    ) %>%
    arrange(desc(inc_ge3))
}

# Unconditional ≥3 rate (denominator = all treated). Primary metric for threshold comparison.
check_ge3_at_alltreated <- function(ae, cutoff) {
  ae_thr <- apply_threshold_any_arm(ae, cutoff = cutoff, key = "aedecod", ge3_only = FALSE)
  ge3_uncond_from_ae(ae_thr) %>%
    mutate(threshold = paste0(round(cutoff * 100), "%")) %>%
    select(trt, threshold, n_all_treated, n_subj_ge3, inc_ge3_uncond, inc_ge3_uncond_pct)
}

# Study thresholds: 1% is the most permissive baseline and JSD reference point.
# Results at 1%, 2%, 5%, 10% characterize the monotonic information-loss pattern (H1).
THRESHOLDS <- c(0.01, 0.02, 0.05, 0.10)
LABELS <- c("1%", "2%", "5%", "10%")
threshold_levels <- LABELS

ge3_sensitivity <- map_dfr(THRESHOLDS, ~check_ge3_at_conditional(ae_ctcae, .x))
print_header("≥G3 BY ARM AFTER THRESHOLDS – CONDITIONAL")
ge3_sensitivity %>% arrange(threshold, desc(inc_ge3)) %>% print(n = Inf)

ge3_sensitivity_alltreated <- map_dfr(THRESHOLDS, ~check_ge3_at_alltreated(ae_ctcae, .x)) %>%
  arrange(trt, match(threshold, threshold_levels))
print_header("≥G3 BY ARM AFTER THRESHOLDS – UNCONDITIONAL (ALL TREATED)")
ge3_sensitivity_alltreated %>% print(n = Inf)

ge3_sensitivity_alltreated <- ge3_sensitivity_alltreated %>%
  mutate(threshold = factor(threshold, levels = threshold_levels))

# Side-by-side comparison: highlights how denominators shift interpretation
ge3_comparison <- ge3_sensitivity %>%
  rename(
    n_subj_ge3_cond   = n_subj_ge3,
    n_subj_total_cond = n_subj_total,
    inc_ge3_cond      = inc_ge3,
    inc_ge3_cond_pct  = inc_ge3_pct
  ) %>%
  select(trt, threshold, n_subj_ge3_cond, n_subj_total_cond, inc_ge3_cond, inc_ge3_cond_pct) %>%
  left_join(
    ge3_sensitivity_alltreated %>%
      rename(
        n_subj_ge3_uncond      = n_subj_ge3,
        inc_ge3_uncond_pct_all = inc_ge3_uncond_pct
      ) %>%
      select(trt, threshold, n_all_treated, n_subj_ge3_uncond, inc_ge3_uncond, inc_ge3_uncond_pct_all),
    by = c("trt", "threshold")
  ) %>%
  arrange(trt, match(threshold, threshold_levels))

print_header("≥G3 CONDITIONAL vs UNCONDITIONAL – SIDE BY SIDE")
ge3_comparison %>% print(n = Inf)

# Drop 1% -> 10% (unconditional, percentage points)
ge3_drop_1_to_10 <- ge3_sensitivity_alltreated %>%
  select(trt, threshold, inc_ge3_uncond) %>%
  pivot_wider(names_from = threshold, values_from = inc_ge3_uncond) %>%
  mutate(drop_1_to_10_pp = round((`1%` - `10%`) * 100, 1))

print_header("DROP IN ≥G3 (1% → 10%) – UNCONDITIONAL")
ge3_drop_1_to_10 %>% print(n = Inf)

# Plot: visibility of ≥G3 subjects as thresholds increase
ggplot(ge3_sensitivity_alltreated,
       aes(x = threshold, y = inc_ge3_uncond, group = trt, color = trt)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Reporting threshold (any arm)",
    y = "Subjects with ≥G3 TEAE (unconditional % of all treated)",
    title = "Visibility of Grade ≥3 TEAEs under increasing thresholds",
    color = "Treatment"
  ) +
  theme_minimal()

# Identify ≥G3 PTs that disappear at higher thresholds (pure rule)
severe_terms_lost_at <- function(ae, cutoff) {
  kept <- apply_threshold_any_arm(ae, cutoff, key = "aedecod") %>%
    filter(ctcae_ge3) %>% distinct(aedecod) %>% pull(aedecod) %>% as.character()

  all <- ae %>%
    filter(ctcae_ge3) %>% distinct(aedecod) %>% pull(aedecod) %>% as.character()

  setdiff(all, kept)
}

dropped_severe_terms_5  <- severe_terms_lost_at(ae_ctcae, 0.05)
dropped_severe_terms_10 <- severe_terms_lost_at(ae_ctcae, 0.10)

print_header("≥G3 PTs LOST AT 5%")
print(dropped_severe_terms_5)

print_header("≥G3 PTs LOST AT 10%")
print(dropped_severe_terms_10)

dropped_soc_5 <- ae_ctcae %>%
  filter(ctcae_ge3, aedecod %in% dropped_severe_terms_5) %>%
  distinct(aebodsys, aedecod) %>%
  arrange(aebodsys, aedecod)

print_header("DROPPED ≥G3 PTs @5% BY SOC")
print(dropped_soc_5, n = Inf)

# How many ≥G3 subjects become “hidden” at 5% threshold (pure)?
ae_ge3_all  <- ae_ctcae %>% filter(ctcae_ge3) %>% distinct(trt, usubjid)
ae_ge3_5pct <- apply_threshold_any_arm(ae_ctcae, 0.05, key = "aedecod") %>%
  filter(ctcae_ge3) %>% distinct(trt, usubjid)

hidden_5pct <- left_join(
  ae_ge3_all  %>% count(trt, name = "n_ge3_all"),
  ae_ge3_5pct %>% count(trt, name = "n_ge3_visible_5pct"),
  by = "trt"
) %>% mutate(
  n_ge3_visible_5pct = replace_na(n_ge3_visible_5pct, 0L),
  n_hidden_5pct = n_ge3_all - n_ge3_visible_5pct
)

print_header("HIDDEN ≥G3 SUBJECTS AT 5% THRESHOLD")
hidden_5pct %>% print(n = Inf)

# Compare pure vs safeguard restoration (unconditional) at a cutoff
compare_threshold_vs_safeguard <- function(ae, cutoff) {
  ae_pure <- apply_threshold_any_arm(ae, cutoff, key = "aedecod")
  pure_tbl <- ge3_uncond_from_ae(ae_pure) %>%
    rename(n_ge3_pure = n_subj_ge3, pct_ge3_pure = inc_ge3_uncond_pct)

  ae_safe <- apply_threshold_with_safeguard(ae, cutoff)
  safe_tbl <- ge3_uncond_from_ae(ae_safe) %>%
    rename(n_ge3_safe = n_subj_ge3, pct_ge3_safe = inc_ge3_uncond_pct)

  out <- left_join(pure_tbl, safe_tbl, by = c("trt", "n_all_treated")) %>%
    mutate(
      n_ge3_restored = n_ge3_safe - n_ge3_pure,
      pp_restored    = round(pct_ge3_safe - pct_ge3_pure, 1)
    ) %>%
    arrange(desc(n_ge3_restored), desc(pp_restored))

  # Which ≥G3 PTs were rescued by safeguard but lost under pure threshold?
  restored_pts <- setdiff(
    ae_safe %>% filter(ctcae_ge3) %>% distinct(aedecod) %>% pull(aedecod) %>% as.character(),
    ae_pure %>% filter(ctcae_ge3) %>% distinct(aedecod) %>% pull(aedecod) %>% as.character()
  ) %>% sort()

  list(summary = out, restored_pts = restored_pts)
}

cmp_5  <- compare_threshold_vs_safeguard(ae_ctcae, 0.05)
cmp_10 <- compare_threshold_vs_safeguard(ae_ctcae, 0.10)

cat("\n=== ≥G3 restoration with safeguard @5% (unconditional, all treated) ===\n")
print(cmp_5$summary, n = Inf)
cat("\nRestored ≥G3 PTs @5%:\n")
print(cmp_5$restored_pts)

cat("\n=== ≥G3 restoration with safeguard @10% (unconditional, all treated) ===\n")
print(cmp_10$summary, n = Inf)
cat("\nRestored ≥G3 PTs @10%:\n")
print(cmp_10$restored_pts)

export_cmp <- bind_rows(
  cmp_5$summary  %>% mutate(threshold = "5%"),
  cmp_10$summary %>% mutate(threshold = "10%")
) %>%
  select(trt, threshold,
         n_all_treated,
         n_ge3_pure, pct_ge3_pure,
         n_ge3_safe, pct_ge3_safe,
         n_ge3_restored, pp_restored)

plot_df <- export_cmp %>%
  pivot_longer(cols = c(pct_ge3_pure, pct_ge3_safe),
               names_to = "rule", values_to = "pct") %>%
  mutate(rule = ifelse(rule == "pct_ge3_pure", "Pure threshold", "Keep-≥G3 safeguard"))

plot_df <- plot_df %>%
  mutate(threshold = factor(threshold, levels = c("5%", "10%")))

ggplot(plot_df, aes(x = trt, y = pct, fill = rule)) +
  geom_col(position = position_dodge()) +
  facet_wrap( ~ threshold) +
  labs(x = NULL, y = "Subjects with ≥G3 TEAE (unconditional %)", title = "Effect of a keep-≥G3 safeguard on thresholded reporting") +
  theme_minimal() + theme(
    panel.border = element_rect(
      color = "grey50",
      fill = NA,
      linetype = "dashed"
    ),
    panel.spacing = unit(1, "lines")
  )

# ============================================================
# 8) Information-loss metrics
#    Six metrics are computed at each threshold × variant combination.
#    All metrics are pooled into 'all_metrics' and normalized into
#    'summary_table' for reporting in Chapter 4.
#    Shannon entropy (H_event, H_subj)
#      Measures distributional diversity of the PT frequency distribution.
#      Higher = more evenly spread across terms; lower = more concentrated.
#      Defined as H = -Σ p_i * log(p_i) (Shannon 1948).
#    Jensen–Shannon divergence (JSD_event, JSD_subj)
#      Symmetric measure of structural deviation from the 1% baseline distribution.
#      JSD = 0 means identical structure; JSD → log(2) means maximal dissimilarity.
#      Computed as 0.5 * KL(p||m) + 0.5 * KL(q||m) where m = (p+q)/2.
#    Coverage (cover_events, cover_subjects, cover_pts)
#      Raw counts of retained events, subjects, and unique PTs.
#      Normalized to the 1% baseline in summary_table for interpretability.
#    Grade ≥3 visibility (ge3_uncond_pct)
#      Unconditional rate of subjects with ≥1 Grade ≥3 TEAE among all treated.
#      Key clinical visibility indicator for H2.
# ============================================================

# Shannon entropy H(p) = -Σ p_i * log(p_i), operating on a normalized probability vector.
# Larger H = more diverse/evenly spread distribution. Zero-probability terms excluded.
# Note: Pielou evenness J = H/log(K) is defined below but excluded from primary reporting.

entropy <- function(p) {
  p <- p[!is.na(p) & p > 0]
  -sum(p * log(p))
}

evenness <- function(p) {
  p <- p[!is.na(p) & p > 0]
  K <- length(p)
  if (K <= 1) return(NA_real_)
  entropy(p) / log(K)
}

kl_div <- function(p, q) {
  idx <- p > 0 & q > 0
  p <- p[idx]; q <- q[idx]
  sum(p * log(p / q))
}

jsd <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
}

# Build normalized probability distributions over PTs, per arm.
#   type = "events"   — p_i = (# AE records for PT i) / (total AE records in arm)
#   type = "subjects" — p_i = (# subjects with PT i) / (total subjects with any AE in arm)
# These are the distributions used as inputs to entropy and JSD.
build_distributions <- function(ae_tbl, type = c("events", "subjects")) {
  type <- match.arg(type)
  base <- ae_tbl %>% mutate(aedecod = as.character(aedecod))
  if (type == "events") {
    base %>%
      count(trt, aedecod, name = "n") %>%
      group_by(trt) %>% mutate(p = n / sum(n)) %>% ungroup()
  } else {
    base %>%
      distinct(trt, usubjid, aedecod) %>%
      count(trt, aedecod, name = "n") %>%
      group_by(trt) %>% mutate(p = n / sum(n)) %>% ungroup()
  }
}

# Compute all information-loss metrics for one threshold x variant combination.
# Returns one row per arm with all metrics. Called via map_dfr() across all combinations.
compute_metrics_for_table <- function(ae_tbl,
                                      threshold_label,
                                      variant_label,
                                      baseline_event_dist,
                                      baseline_subj_dist,
                                      denom_all_treated_tbl = denom_all_treated) {
  # Coverage: how much “mass” remains after thresholding
  cover_events   <- nrow(ae_tbl)
  cover_subjects <- ae_tbl %>% distinct(usubjid) %>% nrow()
  cover_pts      <- ae_tbl %>% distinct(aedecod) %>% nrow()

  dist_ev   <- build_distributions(ae_tbl, "events")
  dist_subj <- build_distributions(ae_tbl, "subjects")

  metrics_ev <- dist_ev %>%
    group_by(trt) %>%
    summarise(H_event = entropy(p), J_event = evenness(p), .groups = "drop")

  metrics_subj <- dist_subj %>%
    group_by(trt) %>%
    summarise(H_subj = entropy(p), J_subj = evenness(p), .groups = "drop")

  # JSD vs 1% baseline: full_join ensures PTs dropped at higher thresholds
  # contribute zero-probability mass (not missing), giving correct divergence.
  jsd_ev <- dist_ev %>%
    rename(p_curr = p) %>% select(trt, aedecod, p_curr) %>%
    full_join(baseline_event_dist %>% rename(p_base = p), by = c("trt", "aedecod")) %>%
    mutate(p_curr = replace_na(p_curr, 0), p_base = replace_na(p_base, 0)) %>%
    group_by(trt) %>% summarise(JSD_event = jsd(p_curr, p_base), .groups = "drop")

  jsd_subj <- dist_subj %>%
    rename(p_curr = p) %>% select(trt, aedecod, p_curr) %>%
    full_join(baseline_subj_dist %>% rename(p_base = p), by = c("trt", "aedecod")) %>%
    mutate(p_curr = replace_na(p_curr, 0), p_base = replace_na(p_base, 0)) %>%
    group_by(trt) %>% summarise(JSD_subj = jsd(p_curr, p_base), .groups = "drop")

  # Clinical visibility: ≥G3 subjects among ALL treated for each arm
  ge3_tbl <- ae_tbl %>%
    filter(ctcae_ge3) %>%
    distinct(trt, usubjid) %>%
    count(trt, name = "n_ge3") %>%
    right_join(denom_all_treated_tbl, by = "trt") %>%
    mutate(
      n_ge3          = replace_na(n_ge3, 0L),
      ge3_uncond_pct = round(100 * n_ge3 / n_all_treated, 1)
    ) %>%
    select(trt, n_all_treated, n_ge3, ge3_uncond_pct)

  metrics_ev %>%
    left_join(metrics_subj, by = "trt") %>%
    left_join(jsd_ev,      by = "trt") %>%
    left_join(jsd_subj,    by = "trt") %>%
    left_join(ge3_tbl,     by = "trt") %>%
    mutate(
      threshold      = threshold_label,
      variant        = variant_label,
      cover_events   = cover_events,
      cover_subjects = cover_subjects,
      cover_pts      = cover_pts
    ) %>%
    relocate(trt, threshold, variant)
}

# JSD baseline: 1% pure threshold is the most permissive rule in this study.
# All higher-threshold distributions are compared against this reference.
# Using 1% as baseline means JSD = 0 at 1% and increases with stricter thresholds.
ae_1pct_pure <- apply_threshold_any_arm(ae_ctcae, 0.01, key = "aedecod")
baseline_event_dist <- build_distributions(ae_1pct_pure, "events")
baseline_subj_dist  <- build_distributions(ae_1pct_pure, "subjects")

all_metrics <- map_dfr(
  THRESHOLDS,
  function(th) {
    lab <- paste0(round(th * 100), "%")

    ae_pure <- apply_threshold_any_arm(ae_ctcae, th, key = "aedecod")
    m_pure <- compute_metrics_for_table(
      ae_tbl = ae_pure, threshold_label = lab, variant_label = "pure",
      baseline_event_dist = baseline_event_dist, baseline_subj_dist = baseline_subj_dist
    )

    ae_safe <- apply_threshold_with_safeguard(ae_ctcae, th)
    m_safe <- compute_metrics_for_table(
      ae_tbl = ae_safe, threshold_label = lab, variant_label = "safeguard",
      baseline_event_dist = baseline_event_dist, baseline_subj_dist = baseline_subj_dist
    )

    bind_rows(m_pure, m_safe)
  }
)

print_header("INFO-LOSS METRICS (head)")
all_metrics %>% head(12) %>% print(n = Inf)

summary_table <- all_metrics %>%
  mutate(threshold = factor(threshold, levels = threshold_levels)) %>%
  group_by(trt, variant) %>%
  mutate(
    # Normalize coverage relative to 1% baseline for interpretability
    ref_events   = cover_events[threshold == "1%"][1],
    ref_subjects = cover_subjects[threshold == "1%"][1],
    ref_pts      = cover_pts[threshold == "1%"][1],
    cover_events   = cover_events   / ref_events,
    cover_subjects = cover_subjects / ref_subjects,
    cover_pts      = cover_pts      / ref_pts
  ) %>%
  ungroup() %>%
  select(
    trt, threshold, variant,
    cover_events, cover_subjects, cover_pts,
    H_event, J_event, H_subj, J_subj,
    JSD_event, JSD_subj,
    n_all_treated, n_ge3, ge3_uncond_pct
  )

print_header("INFO-LOSS SUMMARY")
summary_table %>% arrange(trt, variant, threshold) %>% print(n = Inf)

ggplot(summary_table,
       aes(x = threshold, y = cover_pts, group = variant, color = variant)) +
  geom_line() + geom_point() +
  facet_wrap(~ trt) +
  labs(
    x = "Reporting threshold (any-arm, subject-incidence)",
    y = "Relative PT coverage (vs 1% baseline)",
    color = "Rule",
    title = "Loss of PT coverage under increasing thresholds:\nPure threshold vs keep-≥G3 safeguard"
  ) +
  theme_minimal()

ggplot(summary_table,
       aes(x = threshold, y = H_event, group = variant, color = variant)) +
  geom_line() + geom_point() +
  facet_wrap(~ trt) +
  labs(
    x = "Reporting threshold",
    y = "Entropy of AE event distribution (H)",
    color = "Rule",
    title = "Change in diversity of AE events under thresholding"
  ) +
  theme_minimal()

ggplot(summary_table,
       aes(x = threshold, y = JSD_event, group = variant, color = variant)) +
  geom_line() + geom_point() +
  facet_wrap(~ trt) +
  labs(
    x = "Reporting threshold",
    y = "Jensen–Shannon divergence (events)",
    color = "Rule",
    title = "Structural deviation from 1% baseline under higher thresholds"
  ) +
  theme_minimal()

# ============================================================
# 9) Multivariate methods: Correspondence Analysis and PCA
#    CA (Correspondence Analysis)
#      Applied to a PT × arm contingency table (subject-level counts).
#      CA dimension 1 captures the primary between-arm variation in PT structure.
#      Inertia = proportion of variance explained by each dimension.
#      Implementation: ca::ca() (Nenadic & Greenacre 2007, JSS 20(3)).
#    PCA (Principal Component Analysis)
#      Applied to a binary subject × PT indicator matrix.
#      Each subject is a row; each PT is a column (1 = had PT, 0 = did not).
#      Variance explained by PC1+PC2 captures dimensional richness of the
#      AE structure. Implemented via base R prcomp() with centering and scaling.
#    Both methods are run at 1% and 5% thresholds, pure and safeguard variants,
#    to show how structural information changes with the reporting rule.
# ============================================================


# Build a PT × arm contingency table: cells = number of subjects with that PT in that arm.
# Rows = Preferred Terms; columns = treatment arms.
# This is the input matrix to correspondence analysis.
build_pt_arm_table <- function(ae_tbl) {
  ae_tbl %>%
    distinct(trt, usubjid, aedecod) %>%
    count(aedecod, trt, name = "n_subj") %>%
    pivot_wider(names_from = trt, values_from = n_subj, values_fill = 0) %>%
    tibble::column_to_rownames("aedecod")
}

run_ca_for_threshold <- function(ae, cutoff, variant = c("pure", "safeguard")) {
  variant <- match.arg(variant)
  ae_thr <- if (variant == "pure") apply_threshold_any_arm(ae, cutoff, key = "aedecod") else apply_threshold_with_safeguard(ae, cutoff)
  tab <- build_pt_arm_table(ae_thr)

  # CA requires at least 2 rows and 2 columns; guard to avoid hard failures
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    warning("Too few PTs or arms for CA at cutoff=", cutoff, " variant=", variant)
    return(NULL)
  }
  list(cutoff = cutoff, variant = variant, table = tab, ca_fit = ca::ca(tab))
}

plot_ca_map <- function(ca_obj, title_suffix = "") {
  if (is.null(ca_obj)) return(invisible(NULL))
  ca_fit <- ca_obj$ca_fit
  tab <- ca_obj$table

  row_coords <- as.data.frame(ca_fit$rowcoord) %>% mutate(PT = rownames(tab))
  col_coords <- as.data.frame(ca_fit$colcoord) %>% mutate(TRT = colnames(tab))

  p_rows <- ggplot(row_coords, aes(x = Dim1, y = Dim2)) +
    geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(aes(label = PT), size = 3, max.overlaps = 30) +
    labs(title = paste("CA map of PTs", title_suffix), x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()

  p_cols <- ggplot(col_coords, aes(x = Dim1, y = Dim2, label = TRT)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 4) +
    labs(title = paste("CA map of treatment arms", title_suffix), x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()

  list(rows = p_rows, cols = p_cols)
}

ca_1pct_pure <- run_ca_for_threshold(ae_ctcae, 0.01, "pure")
ca_5pct_pure <- run_ca_for_threshold(ae_ctcae, 0.05, "pure")
ca_5pct_safe <- run_ca_for_threshold(ae_ctcae, 0.05, "safeguard")

print_header("CA – 1% PURE (eigenvalues)")
if (!is.null(ca_1pct_pure)) print(ca_1pct_pure$ca_fit$sv^2)

print_header("CA – 5% PURE (eigenvalues)")
if (!is.null(ca_5pct_pure)) print(ca_5pct_pure$ca_fit$sv^2)

print_header("CA – 5% SAFEGUARD (eigenvalues)")
if (!is.null(ca_5pct_safe)) print(ca_5pct_safe$ca_fit$sv^2)

plot_ca_biplot <- function(ca_obj, title_suffix = "") {
  if (is.null(ca_obj)) return(invisible(NULL))
  
  ca_fit <- ca_obj$ca_fit
  tab <- ca_obj$table
  
  row_coords <- as.data.frame(ca_fit$rowcoord) %>%
    mutate(name = rownames(tab), type = "Preferred Term")
  
  col_coords <- as.data.frame(ca_fit$colcoord) %>%
    mutate(name = colnames(tab), type = "Treatment")
  
  both <- bind_rows(row_coords, col_coords)
  
  ggplot(both, aes(x = Dim1, y = Dim2, color = type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    
    geom_point(alpha = 0.8) +
    
    ggrepel::geom_text_repel(
      aes(label = name),
      size = 3,
      max.overlaps = 30
    ) +
    
    labs(
      title = paste("Correspondence Analysis Biplot", title_suffix),
      x = "Dimension 1",
      y = "Dimension 2",
      color = NULL
    ) +
    theme_minimal()
}

plot_ca_biplot(ca_1pct_pure, "– 1% Pure Threshold")
plot_ca_biplot(ca_5pct_pure, "– 5% Pure Threshold")
plot_ca_biplot(ca_5pct_safe, "– 5% Safeguard Threshold")

# Build a binary subject × PT indicator matrix for PCA.
# Rows = subjects; columns = Preferred Terms; values = 1 (had PT) or 0 (did not).
# pt_filter allows restricting to a specific termset (e.g. the threshold-retained set).
build_subject_pt_matrix <- function(ae_tbl, pt_filter = NULL) {
  base <- ae_tbl %>% distinct(usubjid, trt, aedecod)
  if (!is.null(pt_filter)) base <- base %>% filter(aedecod %in% pt_filter)
  base %>%
    mutate(has_pt = 1L) %>%
    pivot_wider(id_cols = c(usubjid, trt), names_from = aedecod, values_from = has_pt, values_fill = 0L)
}

run_pca_for_threshold <- function(ae, cutoff, variant = c("pure", "safeguard")) {
  variant <- match.arg(variant)
  ae_thr <- if (variant == "pure") apply_threshold_any_arm(ae, cutoff, key = "aedecod") else apply_threshold_with_safeguard(ae, cutoff)

  if (n_distinct(ae_thr$aedecod) < 2) {
    warning("Too few PTs for PCA at cutoff=", cutoff, " variant=", variant)
    return(NULL)
  }

  subj_pt <- build_subject_pt_matrix(ae_thr)
  X <- subj_pt %>% select(-usubjid, -trt) %>% as.matrix()

  if (ncol(X) < 2) {
    warning("Too few PT columns for PCA at cutoff=", cutoff, " variant=", variant)
    return(NULL)
  }

  # Remove PTs with zero variance (all subjects either had or did not have them).
  # Zero-variance columns contribute nothing to PCA and cause numerical warnings.
  keep <- apply(X, 2, sd) > 0
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2) {
    warning("Only 1 non-constant PT after filtering at cutoff=", cutoff, " variant=", variant)
    return(NULL)
  }

  pca_fit <- prcomp(X, center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca_fit$x) %>% mutate(usubjid = subj_pt$usubjid, trt = subj_pt$trt)
  loadings <- as.data.frame(pca_fit$rotation) %>% tibble::rownames_to_column("PT")

  list(cutoff = cutoff, variant = variant, subj_matrix = subj_pt, pca_fit = pca_fit, scores = scores, loadings = loadings)
}

plot_pca_scores <- function(pca_obj, dims = c(1, 2), title_suffix = "") {
  if (is.null(pca_obj)) return(invisible(NULL))
  d1 <- dims[1]; d2 <- dims[2]
  scores <- pca_obj$scores
  pc_names <- colnames(pca_obj$pca_fit$x)

  ggplot(scores, aes(x = .data[[pc_names[d1]]], y = .data[[pc_names[d2]]], color = trt)) +
    geom_point(alpha = 0.7) +
    labs(title = paste0("PCA scores ", title_suffix), x = pc_names[d1], y = pc_names[d2], color = "Treatment arm") +
    theme_minimal()
}

plot_pca_loadings <- function(pca_obj, dims = c(1, 2), title_suffix = "") {
  if (is.null(pca_obj)) return(invisible(NULL))
  d1 <- dims[1]; d2 <- dims[2]
  loadings <- pca_obj$loadings
  pc_names <- colnames(pca_obj$pca_fit$rotation)

  ggplot(loadings, aes(x = .data[[pc_names[d1]]], y = .data[[pc_names[d2]]], label = PT)) +
    geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = paste0("PCA loadings ", title_suffix),
         x = paste0(pc_names[d1], " loading"),
         y = paste0(pc_names[d2], " loading")) +
    theme_minimal()
}

pca_1pct_pure <- run_pca_for_threshold(ae_ctcae, 0.01, "pure")
pca_5pct_pure <- run_pca_for_threshold(ae_ctcae, 0.05, "pure")
pca_5pct_safe <- run_pca_for_threshold(ae_ctcae, 0.05, "safeguard")

print_header("PCA – 1% PURE (variance explained)")
if (!is.null(pca_1pct_pure)) print(summary(pca_1pct_pure$pca_fit))

print_header("PCA – 5% PURE (variance explained)")
if (!is.null(pca_5pct_pure)) print(summary(pca_5pct_pure$pca_fit))

print_header("PCA – 5% SAFEGUARD (variance explained)")
if (!is.null(pca_5pct_safe)) print(summary(pca_5pct_safe$pca_fit))


build_scree_df <- function(pca_obj, label) {
  if (is.null(pca_obj)) return(NULL)
  
  var_explained <- (pca_obj$pca_fit$sdev)^2
  var_explained <- var_explained / sum(var_explained)
  
  data.frame(
    PC = seq_along(var_explained),
    Variance = var_explained,
    Setting = label
  )
}

scree_df <- bind_rows(
  build_scree_df(pca_1pct_pure, "1% Pure"),
  build_scree_df(pca_5pct_pure, "5% Pure"),
  build_scree_df(pca_5pct_safe, "5% Safeguard")
)

scree_df <- scree_df %>% filter(PC <= 20)

ggplot(scree_df, aes(x = PC, y = Variance, color = Setting)) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1, max(scree_df$PC), by = 1)) +
  
  # Match color logic from your example
  scale_color_manual(values = c(
    "1% Pure" = "#F8766D",   # soft red
    "5% Pure" = "#00BA38",   # green
    "5% Safeguard" = "#619CFF"  # blue
  )) +
  
  labs(
    title = "PCA Scree Plot – Variance Explained per Component",
    x = "Principal Component",
    y = "Proportion of Variance Explained",
    color = "Adverse Event Threshold"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================
# 10) Survey visualizations
#     All figures share a single threshold (default 5%) and a common
#     termset helper (termset_anyarm) to ensure consistency across plots.
#     Each visualization targets a different dimension of information loss:
#     Butterfly plot   — incidence comparison, two arms side-by-side
#     Volcano plot     — severity-weighted PT signals with statistical support
#     Temporal plot    — cumulative incidence of first onset over study time
#     Jaccard heatmap  — PT co-occurrence within an arm
#     CA strip plot    — one-axis biplot mapping PTs to placebo vs active
# ============================================================

# Shared color palette used consistently across all figures.
# Consistency is important for the expert survey: participants should
# immediately associate blue with active and grey with placebo across plots.
COL_PLACEBO <- "#7A7A7A"
COL_ACTIVE  <- "#0072B2"
COL_G3      <- "#D55E00"
COL_NONG3   <- "#8C8C8C"

# Returns the set of PT/SOC labels retained under the any-arm threshold.
# All survey figures use this with the same cutoff to ensure a consistent PT universe.
termset_anyarm <- function(ae_tbl, cutoff, key = c("aedecod", "aebodsys"), ge3_only = FALSE) {
  key <- match.arg(key)
  apply_threshold_any_arm(ae_tbl, cutoff = cutoff, key = key, ge3_only = ge3_only) %>%
    distinct(.data[[key]]) %>%
    pull(.data[[key]]) %>%
    as.character()
}

# Compute subject-level incidence per PT or SOC for one arm,
# returning a tidy table used as input to both butterfly and volcano plots.
make_incidence <- function(ae_tbl, key = c("aedecod", "aebodsys"), ge3_only = FALSE, denom_tbl = NULL) {
  key <- match.arg(key)
  denom_tbl <- denom_tbl %||% denom_from_ae(ae_tbl)
  x <- if (ge3_only) ae_tbl %>% filter(ctcae_ge3) else ae_tbl

  x %>%
    distinct(trt, usubjid, .data[[key]]) %>%
    count(trt, .data[[key]], name = "n") %>%
    left_join(denom_tbl, by = "trt") %>%
    mutate(inc = ifelse(n_trt > 0, n / n_trt, 0)) %>%
    transmute(trt, term = as.character(.data[[key]]), n, n_trt, inc)
}

# ------------------------------------------------------------
# Butterfly plot
#
#   Three-panel faceted layout: left arm | PT labels | right arm.
#   PT selection: hybrid — top N by maximum incidence (common events)
#   plus top M by absolute risk difference (differentiating events).
#   Display order: descending active-arm incidence, abs(RD) as tie-breaker.
# ------------------------------------------------------------


prepare_butterfly_faceted <- function(inc_tbl,
                                      arm_left,
                                      arm_right,
                                      top_n = 20,
                                      ordering = c("max_inc", "abs_rd", "hybrid"),
                                      n_common = 10,
                                      n_diff   = 5,
                                      center_label = "Preferred Term") {
  ordering <- match.arg(ordering)

  wide <- inc_tbl %>%
    dplyr::filter(trt %in% c(arm_left, arm_right)) %>%
    dplyr::select(trt, term, inc) %>%
    tidyr::pivot_wider(names_from = trt, values_from = inc, values_fill = 0)

  if (!arm_left %in% names(wide))  wide[[arm_left]]  <- 0
  if (!arm_right %in% names(wide)) wide[[arm_right]] <- 0

  wide <- wide %>%
    dplyr::mutate(
      inc_left  = .data[[arm_left]],
      inc_right = .data[[arm_right]],
      rd        = inc_right - inc_left,
      abs_rd    = abs(rd),
      max_inc   = pmax(inc_left, inc_right)
    )

  # hybrid PT selection: combine top-by-incidence and top-by-risk-difference
  wide_sel <- if (ordering == "hybrid") {

    top_common <- wide %>%
      dplyr::arrange(dplyr::desc(max_inc)) %>%
      dplyr::slice_head(n = n_common)

    top_diff <- wide %>%
      dplyr::arrange(dplyr::desc(abs_rd)) %>%
      dplyr::slice_head(n = n_diff)

    dplyr::bind_rows(top_common, top_diff) %>%
      dplyr::distinct(term, .keep_all = TRUE)

  } else {
    metric <- if (ordering == "abs_rd") "abs_rd" else "max_inc"
    wide %>%
      dplyr::arrange(dplyr::desc(.data[[metric]])) %>%
      dplyr::slice_head(n = top_n)
  }

  # display order: descending active-arm incidence, abs(RD) as tie-breaker
  wide_sel <- wide_sel %>%
    dplyr::arrange(dplyr::desc(inc_right), dplyr::desc(abs_rd)) %>%
    dplyr::mutate(
      term = forcats::fct_reorder(term, inc_right, .desc = TRUE)
    )

  bars <- wide_sel %>%
    dplyr::transmute(
      term,
      facet = arm_left,
      arm = arm_left,
      inc = inc_left,
      inc_plot = -inc_left
    ) %>%
    dplyr::bind_rows(
      wide_sel %>%
        dplyr::transmute(
          term,
          facet = arm_right,
          arm = arm_right,
          inc = inc_right,
          inc_plot = inc_right
        )
    )

  pts <- wide_sel %>%
    dplyr::transmute(
      term,
      facet = center_label,
      x = 0,
      label = as.character(term)
    )

  list(
    bars = bars,
    pts = pts,
    wide = wide_sel,
    center_label = center_label
  )
}


plot_butterfly_faceted <- function(bfly,
                                   arm_left,
                                   arm_right,
                                   title,
                                   subtitle = NULL,
                                   pct_accuracy = 1,
                                   wrap_width = 28,
                                   x_breaks_left  = c(0, -0.05, -0.10, -0.15),
                                   x_breaks_right = c(0,  0.10,  0.20,  0.30),
                                   bottom_vjust = 1.4,
                                   center_halfwidth = 0.08,
                                   panel_spacing_x = 0.06,
                                   header_margin_top = 18,
                                   header_fontsize = 4.2,
                                   bar_width = 0.6,
                                   pt_fontsize = 3.6,
                                   pt_lineheight = 1.05) {

  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Please install ggh4x: install.packages('ggh4x')")
  }

  center_label <- bfly$center_label %||% "Preferred Term"

  bars <- bfly$bars %>%
    dplyr::mutate(term_w = stringr::str_wrap(as.character(term), wrap_width))

  pts <- bfly$pts %>%
    dplyr::mutate(
      term_w = stringr::str_wrap(as.character(term), wrap_width),
      label  = stringr::str_wrap(as.character(label), wrap_width)
    )

  # add one extra dummy factor level used as a header spacer row
  lvl <- unique(bars$term_w)
  header_level <- "\u200B"   # zero-width space as a "blank" label row
  lvl2 <- c(lvl, header_level)

  bars$term_w <- factor(bars$term_w, levels = lvl2)
  pts$term_w  <- factor(pts$term_w,  levels = lvl2)

  facet_levels <- c(arm_left, center_label, arm_right)
  bars$facet <- factor(bars$facet, levels = facet_levels)
  pts$facet  <- factor(pts$facet,  levels = facet_levels)

  is_pbo <- function(x) grepl("PLACEBO|PBO", toupper(x))
  fill_vals <- setNames(
    c(if (is_pbo(arm_left)) COL_PLACEBO else COL_ACTIVE,
      if (is_pbo(arm_right)) COL_PLACEBO else COL_ACTIVE),
    c(arm_left, arm_right)
  )

  # manual percentage labels for bar facets only (center panel has none)
  tick_df <- dplyr::bind_rows(
    tibble::tibble(facet = arm_left,  x = x_breaks_left),
    tibble::tibble(facet = arm_right, x = x_breaks_right)
  ) %>%
    dplyr::mutate(
      facet = factor(facet, levels = facet_levels),
      lab   = scales::percent(abs(x), accuracy = pct_accuracy),
      y     = -Inf
    )

  # ghost points to enforce equal left/right facet widths
  max_lim <- max(bars$inc, na.rm = TRUE)

  lr_blank <- tibble::tibble(
    facet  = factor(c(arm_left, arm_left, arm_right, arm_right), levels = facet_levels),
    term_w = factor(lvl2[1], levels = lvl2),
    x      = c(-max_lim, 0, 0, max_lim)
  )

  # ghost points to constrain center facet width
  center_blank <- tibble::tibble(
    facet  = factor(center_label, levels = facet_levels),
    term_w = factor(lvl2[1], levels = lvl2),
    x      = c(-center_halfwidth, center_halfwidth)
  )

  # header label row: arm names rendered at top of each facet
  header_df <- tibble::tibble(
    facet  = factor(facet_levels, levels = facet_levels),
    x      = c(-0.15, 0, 0.15),
    term_w = factor(header_level, levels = lvl2),
    lab    = facet_levels,
    col    = c(fill_vals[[arm_left]], "black", fill_vals[[arm_right]]),
    hjust  = c(1, 0.5, 0)
  )

  # --- Draw "y-axis" lines as segments that stop BEFORE the header row ---
  y0 <- 0.5
  y1 <- length(lvl) + 0.5

  axis_df <- tibble::tibble(
    facet = factor(c(arm_left, arm_right), levels = facet_levels),
    x    = 0,
    xend = 0,
    y    = y0,
    yend = y1
  )

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = bars,
      ggplot2::aes(x = inc_plot, y = term_w, fill = arm),
      width = bar_width
    ) +
    ggplot2::geom_segment(
      data = axis_df,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 1.1,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = bars,
      ggplot2::aes(
        x = inc_plot,
        y = term_w,
        label = scales::percent(inc, accuracy = pct_accuracy),
        hjust = ifelse(inc_plot < 0, 1.1, -0.1)
      ),
      size = 3
    ) +
    ggplot2::geom_text(
      data = pts,
      ggplot2::aes(x = 0, y = term_w, label = label),
      size = pt_fontsize,
      lineheight = pt_lineheight,
      hjust = 0.5
    ) +
    ggplot2::geom_blank(data = lr_blank,     ggplot2::aes(x = x, y = term_w)) +
    ggplot2::geom_blank(data = center_blank, ggplot2::aes(x = x, y = term_w)) +
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(x = x, y = term_w, label = lab, hjust = hjust),
      inherit.aes = FALSE,
      fontface = "bold",
      size = header_fontsize,
      color = header_df$col,
      vjust = 0.5
    ) +
    ggh4x::facet_grid2(
      . ~ facet,
      scales = "free_x",
      space  = "free_x"
    ) +
    ggplot2::scale_fill_manual(values = fill_vals) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::geom_text(
      data = tick_df,
      ggplot2::aes(x = x, y = y, label = lab),
      inherit.aes = FALSE,
      vjust = bottom_vjust,
      size = 3,
      color = "grey30"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Incidence (% of subjects)",
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),

      legend.position = "none",

      panel.spacing.x = grid::unit(panel_spacing_x, "lines"),
      strip.text       = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(header_margin_top, 15, 28, 15)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}


# --- Example ---
cutoff   <- 0.05
key      <- "aedecod"
ge3_only <- FALSE

ae_use <- apply_threshold_any_arm(
  ae = ae_ctcae,
  cutoff = cutoff,
  key = key,
  ge3_only = ge3_only
)

denom <- denom_from_ae(ae_ctcae)

inc <- make_incidence(
  ae_tbl   = ae_use,
  key      = key,
  ge3_only = ge3_only,
  denom_tbl = denom
)

bf <- prepare_butterfly_faceted(
  inc_tbl   = inc,
  arm_left  = "Placebo",
  arm_right = "Xanomeline High Dose",
  ordering  = "hybrid",
  n_common  = 10,
  n_diff    = 5,
  center_label = "Preferred Term"
)

plot_butterfly_faceted(
  bf,
  arm_left  = "Placebo",
  arm_right = "Xanomeline High Dose",
  title     = "PT incidence (all grades) | 5% any-arm threshold",
  subtitle  = "Placebo vs Xanomeline High Dose",
  x_breaks_left  = c(0, -0.10, -0.20, -0.30),
  x_breaks_right = c(0,  0.10,  0.20,  0.30),
  center_halfwidth = 0.10,
  panel_spacing_x  = 0.04,
  header_margin_top = 18,
  header_fontsize = 4
)


# ------------------------------------------------------------
# Volcano plot
#
#   X axis: difference in severity burden between arms (active − placebo),
#           where burden = sum of max CTCAE grade per subject per PT, normalized by N.
#   Y axis: −log10(Fisher p-value), optionally BH-adjusted.
#   Color: whether any subject had a ≥G3 event for that PT.
# ------------------------------------------------------------

pt_termset <- function(ae_tbl, cutoff = 0.05, ge3_only = FALSE) {
  termset_anyarm(ae_tbl, cutoff = cutoff, key = "aedecod", ge3_only = ge3_only)
}

volcano_severity_df <- function(ae_tbl,
                                arm_a,
                                arm_b,
                                cutoff = 0.05,
                                termset = NULL,
                                p_adjust = c("BH", "none", "holm", "bonferroni", "BY")) {
  p_adjust <- match.arg(p_adjust)

  if (is.null(termset)) termset <- pt_termset(ae_tbl, cutoff, ge3_only = FALSE)

  denom <- denom_from_ae(ae_tbl)

  N_a <- denom %>% dplyr::filter(trt == arm_a) %>% dplyr::pull(n_trt)
  N_b <- denom %>% dplyr::filter(trt == arm_b) %>% dplyr::pull(n_trt)

  if (length(N_a) != 1 || length(N_b) != 1) {
    stop("Could not find unique denominators for both arms in denom_from_ae().")
  }

  ae_use <- ae_tbl %>%
    dplyr::filter(trt %in% c(arm_a, arm_b), aedecod %in% termset)

  # --- Severity burden: max grade per subject per PT; sum within arm; normalize by N ---
  sev <- ae_use %>%
    dplyr::group_by(trt, usubjid, aedecod) %>%
    dplyr::summarise(max_grade = max(grade, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(trt, aedecod) %>%
    dplyr::summarise(sev_sum = sum(max_grade, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = trt, values_from = sev_sum, values_fill = 0)

  if (!arm_a %in% names(sev)) sev[[arm_a]] <- 0
  if (!arm_b %in% names(sev)) sev[[arm_b]] <- 0

  sev <- sev %>%
    dplyr::transmute(
      aedecod,
      delta_burden = (.data[[arm_b]] / N_b) - (.data[[arm_a]] / N_a)
    )

  # --- Incidence for Fisher test (subject-level) ---
  inc <- ae_use %>%
    dplyr::distinct(trt, usubjid, aedecod) %>%
    dplyr::count(trt, aedecod, name = "n") %>%
    tidyr::pivot_wider(names_from = trt, values_from = n, values_fill = 0)

  if (!arm_a %in% names(inc)) inc[[arm_a]] <- 0
  if (!arm_b %in% names(inc)) inc[[arm_b]] <- 0

  inc <- inc %>%
    dplyr::transmute(
      aedecod,
      n_a = .data[[arm_a]],
      n_b = .data[[arm_b]]
    )

  ge3_terms <- ae_tbl %>%
    dplyr::filter(ctcae_ge3) %>%
    dplyr::distinct(aedecod) %>%
    dplyr::pull(aedecod) %>%
    as.character()

  out <- dplyr::left_join(sev, inc, by = "aedecod")

  out$p <- mapply(function(nb, na) {
    fisher.test(matrix(c(nb, N_b - nb, na, N_a - na), nrow = 2, byrow = TRUE))$p.value
  }, out$n_b, out$n_a)

  out %>%
    dplyr::mutate(
      p_adjust_method = p_adjust,
      p_plot = if (p_adjust == "none") p else p.adjust(p, method = p_adjust),
      p_plot = pmax(p_plot, 1e-300),
      neglog10p = -log10(p_plot),
      threshold = paste0(round(cutoff * 100), "% any-arm"),
      # store comparison labels in both directions for flexible subtitle rendering
      comparison_ab = paste(arm_a, "vs", arm_b),
      comparison_ba = paste(arm_b, "vs", arm_a),
      arm_a = arm_a,
      arm_b = arm_b,
      ge3_flag = factor(
        ifelse(as.character(aedecod) %in% ge3_terms, "≥G3 PT", "Grade <3 only"),
        levels = c("Grade <3 only", "≥G3 PT")
      )
    )
}

shared_xlim <- function(df1, df2 = NULL, pad = 0.05, min_width = 1e-4) {
  x <- df1$delta_burden
  if (!is.null(df2)) x <- c(x, df2$delta_burden)

  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-min_width, min_width))

  xmax <- max(abs(x), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- min_width

  xmax <- xmax * (1 + pad)
  c(-xmax, xmax)
}

plot_volcano_severity <- function(
    df,
    x_limits = NULL,
    p_line = 0.05,
    label_n_each_side = 6,
    label_only_sig = TRUE,
    label_all = FALSE,
    title = NULL,
    subtitle_extra = NULL,
    label_size = 3.0,
    wrap_width = 18,
    subtitle_order = c("arm_a_vs_arm_b", "arm_b_vs_arm_a")
) {
  subtitle_order <- match.arg(subtitle_order)

  df <- df %>%
    dplyr::mutate(
      sig = p_plot < p_line,
      label_pretty = stringr::str_to_title(as.character(aedecod)) %>%
        stringr::str_wrap(width = wrap_width)
    )

  # label selection: top N per side by combined score (|delta_burden| * -log10 p)
  if (label_all) {
    lab_df <- df
  } else {
    lab_pool <- if (label_only_sig) df %>% dplyr::filter(sig) else df
    lab_df <- if (nrow(lab_pool) == 0) lab_pool else {
      lab_pool %>%
        dplyr::mutate(
          side = ifelse(delta_burden >= 0, "right", "left"),
          score = abs(delta_burden) * neglog10p
        ) %>%
        dplyr::group_by(side) %>%
        dplyr::slice_max(order_by = score, n = label_n_each_side, with_ties = FALSE) %>%
        dplyr::ungroup()
    }
  }

  if (is.null(x_limits)) {
    xmax <- max(abs(df$delta_burden), na.rm = TRUE)
    if (!is.finite(xmax) || xmax == 0) xmax <- 1e-4
    x_limits <- c(-xmax, xmax)
  }

  y_top <- max(df$neglog10p, na.rm = TRUE)
  if (!is.finite(y_top)) y_top <- -log10(p_line)

  y_annot <- y_top * 1.06
  x_left  <- x_limits[1] + 0.03 * diff(x_limits)
  x_right <- x_limits[2] - 0.03 * diff(x_limits)
  y_line <- -log10(p_line)

  # annotation positions for directional labels (above/below significance line)
  x_y <- x_limits[1] + 0.02 * diff(x_limits)
  y_up   <- y_line + 0.5 * (y_annot - y_line)
  y_down <- 0.5 * y_line


  method <- df$p_adjust_method[1] %||% "none"
  adj_txt <- if (method == "none") "p" else paste0(method, "-adjusted p")

  # --- Subtitle ---
  comp_txt <- if (subtitle_order == "arm_a_vs_arm_b") df$comparison_ab[1] else df$comparison_ba[1]

  ggplot2::ggplot(df, ggplot2::aes(x = delta_burden, y = neglog10p)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(p_line), linetype = 2, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = ge3_flag, alpha = sig), size = 2.1) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.10), guide = "none") +

    ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = label_pretty, color = ge3_flag),
      size = label_size,
      min.segment.length = 0,
      box.padding = if (label_all) 0.15 else 0.35,
      point.padding = if (label_all) 0.10 else 0.25,
      max.overlaps = Inf
    ) +

    ggplot2::scale_color_manual(values = c("Grade <3 only" = COL_NONG3, "≥G3 PT" = COL_G3)) +

    ggplot2::coord_cartesian(xlim = x_limits, ylim = c(0, y_annot)) +

    ggplot2::annotate("text", x = x_left,  y = y_annot,
                      hjust = 0, vjust = 0,
                      label = "More on Placebo \u2190", size = 3.6) +
    ggplot2::annotate("text", x = x_right, y = y_annot,
                      hjust = 1, vjust = 0,
                      label = "\u2192 More on Active", size = 3.6) +
    ggplot2::annotate(
      "text",
      x = x_y, y = y_up,
      label = "More statistical evidence \u2191",
      hjust = 0.5, vjust = 0,
      size = 3.4
    ) +
    ggplot2::annotate(
      "text",
      x = x_y, y = y_down,
      label = "Less statistical evidence \u2193",
      hjust = 0.5, vjust = 1,
      size = 3.4
    ) +


    ggplot2::labs(
      title = title %||% paste0("Severity-weighted PT signals: ", comp_txt),
      subtitle = paste0(
        df$threshold[1], " | ", comp_txt," = ",
        if (!is.null(subtitle_extra)) paste0(" | ", subtitle_extra) else ""
      ),
      x = "Difference in severity burden (Active − Placebo)",
      y = "−log10 p-value",
      color = "PT severity class",
      caption =
        if (label_all) "All PTs labeled."
      else if (label_only_sig) paste0("Labels shown only for PTs with ", adj_txt, " < ", p_line, ".")
      else NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text  = ggplot2::element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )
}

# --- Example ---

termset <- pt_termset(ae_ctcae, 0.05)

vol_hd <- volcano_severity_df(
  ae_ctcae, "Placebo", "Xanomeline High Dose",
  cutoff = 0.05, termset = termset,
  p_adjust = "none"
)

xlims <- shared_xlim(vol_hd)

plot_volcano_severity(
  vol_hd,
  x_limits = xlims,
  p_line = 0.05,
  label_all = FALSE,
  label_only_sig = FALSE,   # Option 2
  label_n_each_side = 6,
  title = "Severity-weighted PT signals: Placebo vs High Dose",
  subtitle_order = "arm_a_vs_arm_b",  # <-- Placebo vs High Dose
  label_size = 3.0,
  wrap_width = 18
)


#indicating broadly similar severity between Xanomeline High Dose and placebo. A small number of itching-related terms, particularly pruritus and application-site pruritus, remain clearly separated, demonstrating a focused and robust increase in severity burden associated with Xanomeline.”


# ------------------------------------------------------------
# Temporal plot (cumulative incidence of first onset)
#
#   Each panel shows a step curve for one PT, comparing first-onset
#   cumulative incidence over study day between two arms.
#   PTs are selected by largest absolute end-of-study incidence difference.
# ------------------------------------------------------------

# --- First onset per subject ---
make_first_onset <- function(ae_tbl, arm_a, arm_b) {
  ae_tbl %>%
    dplyr::filter(trt %in% c(arm_a, arm_b), !is.na(onset_day)) %>%
    dplyr::distinct(trt, usubjid, aedecod, onset_day) %>%
    dplyr::group_by(trt, usubjid, aedecod) %>%
    dplyr::summarise(first_day = min(onset_day, na.rm = TRUE), .groups = "drop")
}

# --- Cumulative incidence ---
make_cuminc <- function(first_onset, ae_tbl, arm_a, arm_b, baseline_day = 0) {

  denom <- denom_from_ae(
    ae_tbl %>% dplyr::filter(trt %in% c(arm_a, arm_b))
  ) %>% dplyr::distinct(trt, .keep_all = TRUE)

  cum <- first_onset %>%
    dplyr::count(trt, aedecod, first_day, name = "n_day") %>%
    dplyr::left_join(denom, by = "trt") %>%
    dplyr::arrange(trt, aedecod, first_day) %>%
    dplyr::group_by(trt, aedecod) %>%
    dplyr::mutate(
      cum_n  = cumsum(n_day),
      cuminc = ifelse(n_trt > 0, cum_n / n_trt, 0)
    ) %>%
    dplyr::ungroup()

  baseline <- cum %>%
    dplyr::distinct(trt, aedecod) %>%
    dplyr::mutate(
      first_day = baseline_day,
      n_day = 0L,
      cum_n = 0L,
      cuminc = 0
    )

  dplyr::bind_rows(baseline, cum) %>%
    dplyr::arrange(trt, aedecod, first_day)
}

# --- End-of-study difference table ---
term_end_diff_table <- function(cum_df, arm_a, arm_b) {
  cum_df %>%
    dplyr::group_by(trt, aedecod) %>%
    dplyr::summarise(end_cuminc = max(cuminc, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = trt, values_from = end_cuminc, values_fill = 0) %>%
    dplyr::mutate(
      end_a = .data[[arm_a]],
      end_b = .data[[arm_b]],
      end_diff = end_b - end_a,
      abs_end_diff = abs(end_diff)
    ) %>%
    dplyr::arrange(dplyr::desc(abs_end_diff))
}

# --- Survey PT selection ---
pick_terms_for_survey <- function(end_tbl, top_n = 9) {
  end_tbl %>% dplyr::slice_head(n = min(top_n, nrow(end_tbl)))
}

# --- Plot ---
plot_temporal <- function(df,
                          title,
                          subtitle,
                          placebo = "Placebo",
                          facet_ncol = 3,
                          wrap_width = 22) {

  if (nrow(df) == 0) stop("No data to plot")

  trt_levels <- unique(df$trt)
  col_vals <- setNames(rep(COL_ACTIVE, length(trt_levels)), trt_levels)
  if (placebo %in% names(col_vals)) col_vals[[placebo]] <- COL_PLACEBO

  df <- df %>%
    dplyr::mutate(
      aedecod_pretty = stringr::str_wrap(
        stringr::str_to_title(as.character(aedecod)),
        wrap_width
      )
    )

  y_top <- max(df$cuminc, na.rm = TRUE) * 1.15

  x_range <- diff(range(df$first_day, na.rm = TRUE))

  end_lab <- df %>%
    dplyr::group_by(trt, aedecod_pretty) %>%
    dplyr::filter(first_day == max(first_day, na.rm = TRUE)) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lab   = scales::percent(cuminc, accuracy = 1),
      x_lab = first_day + x_range * 0.05,   # ← moved further right
      y_lab = pmin(cuminc + y_top * 0.03, y_top * 0.98)
    )


  ggplot2::ggplot(df, ggplot2::aes(first_day, cuminc, color = trt)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::facet_wrap(~ aedecod_pretty, ncol = facet_ncol) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, y_top)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Study day (first onset)",
      y = "Cumulative incidence (% subjects)",
      color = NULL
    ) +
    ggrepel::geom_text_repel(
      data = end_lab,
      ggplot2::aes(x = x_lab, y = y_lab, label = lab, color = trt),
      inherit.aes = FALSE,
      direction = "y",
      hjust = 0,
      size = 3.2,
      segment.color = NA,
      show.legend = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

# --- Build ---
build_temporal_plot <- function(ae_tbl,
                                arm_a = "Placebo",
                                arm_b,
                                cutoff = 0.05,
                                top_n = 9,
                                baseline_day = 0) {

  termset <- pt_termset(ae_tbl, cutoff)
  ae_use  <- ae_tbl %>% dplyr::filter(aedecod %in% termset)

  first <- make_first_onset(ae_use, arm_a, arm_b)
  cum   <- make_cuminc(first, ae_use, arm_a, arm_b, baseline_day)

  end_tbl <- term_end_diff_table(cum, arm_a, arm_b)
  picked  <- pick_terms_for_survey(end_tbl, top_n)

  cum_pick <- cum %>%
    dplyr::filter(aedecod %in% picked$aedecod)

  list(
    picked_table = picked,
    plot = plot_temporal(
      cum_pick,
      title = "Temporal Safety Signals (Survey PTs)",
      subtitle = paste0(
        round(cutoff * 100),
        "% any-arm threshold | ",
        arm_a,
        " vs ",
        arm_b
      ),
      placebo = arm_a
    )
  )
}

# --- Example ---
temp_hd <- build_temporal_plot(
  ae_ctcae,
  arm_a = "Placebo",
  arm_b = "Xanomeline High Dose",
  cutoff = 0.05,
  top_n = 9
)

temp_hd$plot
temp_hd$picked_table


#treatment-related tolerability signal rather than a broad increase in adverse events.”


# ------------------------------------------------------------
# Jaccard co-occurrence heatmap
#
#   Computes pairwise Jaccard similarity for all PT pairs within one arm:
#     Jaccard(A,B) = |A ∩ B| / |A ∪ B|  (subjects with both / subjects with either).
#   Values range from 0 (no shared subjects) to 1 (identical subject sets).
#   The diagonal is excluded; PT labels are abbreviated to single-line form.
# ------------------------------------------------------------

# Build subject × PT matrix for one arm
make_subject_pt_matrix_one_arm <- function(ae_tbl, arm, termset) {
  ae_tbl %>%
    dplyr::filter(trt == arm, aedecod %in% termset) %>%
    dplyr::distinct(usubjid, aedecod) %>%
    dplyr::mutate(value = 1L) %>%
    tidyr::pivot_wider(
      names_from = aedecod,
      values_from = value,
      values_fill = 0L
    )
}

# Convert wide subject × PT table to strict binary matrix
binary_matrix <- function(mat_tbl) {
  stopifnot("usubjid" %in% names(mat_tbl))
  X <- mat_tbl %>% dplyr::select(-usubjid) %>% as.matrix()
  storage.mode(X) <- "integer"
  (X > 0) * 1L
}

# Co-occurrence matrix
cooc_from_binary <- function(X) {
  as.matrix(t(X) %*% X)
}

# Jaccard similarity
jaccard_from_binary <- function(X) {
  C <- t(X) %*% X
  m <- diag(C)
  U <- outer(m, m, "+") - C
  J <- C / ifelse(U == 0, NA, U)
  J[is.na(J)] <- 0
  J
}

# Matrix → long
matrix_to_long <- function(M, value_name) {
  df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(df) <- c("pt1", "pt2", value_name)
  df
}

# Order PTs alphabetically.
order_pts_alphabetical <- function(pt_names) {
  sort(pt_names)
}


apply_order_long <- function(df_long, ord) {
  df_long %>%
    dplyr::mutate(
      pt1 = factor(pt1, levels = ord),
      pt2 = factor(pt2, levels = rev(ord))
    )
}

# --- Single-line PT label helper (abbreviated, no wrapping) ---
pretty_pt_label_one_line <- function(x) {
  x %>%
    stringr::str_to_title() %>%
    stringr::str_replace_all("^Application Site ", "App site ") %>%
    stringr::str_replace_all("^Electrocardiogram ", "ECG ") %>%
    stringr::str_replace_all("^Upper Respiratory Tract Infection$", "URTI") %>%
    stringr::str_replace_all("\\s+", " ")
}

# --- Plot ---
plot_jaccard_heatmap <- function(df_long,
                                 title,
                                 subtitle = NULL,
                                 drop_diag = TRUE,
                                 jaccard_cap = 0.6,
                                 axis_text_size = 9,
                                 tick_length_pt = 3) {

  dfp <- if (drop_diag) df_long %>% dplyr::filter(pt1 != pt2) else df_long
  dfp <- dfp %>% dplyr::mutate(jaccard_plot = pmin(jaccard, jaccard_cap))

  xlabs <- pretty_pt_label_one_line(levels(dfp$pt1))
  names(xlabs) <- levels(dfp$pt1)

  ylabs <- pretty_pt_label_one_line(levels(dfp$pt2))
  names(ylabs) <- levels(dfp$pt2)

  ggplot2::ggplot(dfp, ggplot2::aes(pt1, pt2, fill = jaccard_plot)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(
      name   = "Jaccard\n(similarity)",
      option = "magma",
      limits = c(0, jaccard_cap),
      breaks = seq(0, jaccard_cap, 0.2)
    ) +
    ggplot2::scale_x_discrete(labels = xlabs) +
    ggplot2::scale_y_discrete(labels = ylabs) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        subtitle %||% "",
        ifelse(is.null(subtitle), "", " | "),
        "Color scale capped at ", jaccard_cap
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),

      # X axis
      axis.text.x = ggplot2::element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = axis_text_size
      ),
      axis.ticks.x = ggplot2::element_line(color = "grey50"),
      axis.ticks.length.x = grid::unit(tick_length_pt, "pt"),

      # Y axis
      axis.text.y = ggplot2::element_text(size = axis_text_size),
      axis.ticks.y = ggplot2::element_line(color = "grey50"),
      axis.ticks.length.y = grid::unit(tick_length_pt, "pt"),

      plot.margin = ggplot2::margin(t = 10, r = 10, b = 26, l = 10)
    )
}

# --- Build ---
build_heatmap_data_one_arm <- function(ae_tbl, arm, termset, top_n = NULL) {
  mat_tbl <- make_subject_pt_matrix_one_arm(ae_tbl, arm, termset)
  X <- binary_matrix(mat_tbl)

  C <- cooc_from_binary(X)
  J <- jaccard_from_binary(X)

  colnames(C) <- rownames(C) <- colnames(X)
  colnames(J) <- rownames(J) <- colnames(X)

  all_pts <- colnames(X)

  ord <- order_pts_alphabetical(all_pts)

  if (!is.null(top_n)) {
    ord <- head(ord, top_n)
  }


  list(
    order = ord,
    jac_long = matrix_to_long(J[ord, ord], "jaccard") %>% apply_order_long(ord)
  )
}

# --- Example ---
hm_hd <- build_heatmap_data_one_arm(
  ae_ctcae,
  arm = "Xanomeline High Dose",
  termset = pt_termset(ae_ctcae, cutoff = 0.05),
  top_n = 25
)

p_hd_jac <- plot_jaccard_heatmap(
  hm_hd$jac_long,
  title = "Preferred Term Co-Occurrence Heatmap (High Dose, 5% Pure)",
  subtitle = "Jaccard"
)

p_hd_jac


#suggesting that they are sporadic and not part of a broader syndrome. This pattern supports a focused, mechanistically coherent tolerability signal rather than a diffuse safety concern.”


# ------------------------------------------------------------
# CA strip plot (survey-friendly one-axis biplot)
#
#   With exactly two arms, CA1 captures 100% of between-arm inertia.
#   This collapses the biplot to a single horizontal strip:
#     PTs left of center → relatively more common on placebo.
#     PTs right of center → relatively more common on the active arm.
#   The n_each_side most extreme PTs on each side are labeled.
# ------------------------------------------------------------

pt_pretty_wrap <- function(x, wrap_width = 22) {
  stringr::str_wrap(
    stringr::str_to_title(stringr::str_squish(as.character(x))),
    width = wrap_width
  )
}

make_unique_labels <- function(lbl_vec) {
  make.unique(lbl_vec, sep = " ")
}

build_ca_table_two_arm <- function(ae_tbl, placebo, active, cutoff = 0.05, ge3_only = FALSE) {
  termset <- pt_termset(ae_tbl, cutoff = cutoff, ge3_only = ge3_only)

  tab <- ae_tbl %>%
    dplyr::filter(trt %in% c(placebo, active), aedecod %in% termset) %>%
    dplyr::distinct(trt, usubjid, aedecod) %>%
    dplyr::count(aedecod, trt, name = "n_subj") %>%
    tidyr::pivot_wider(names_from = trt, values_from = n_subj, values_fill = 0)

  M <- tab %>%
    tibble::column_to_rownames("aedecod") %>%
    as.matrix()

  M <- M[rowSums(M) > 0, colSums(M) > 0, drop = FALSE]
  list(termset = termset, M = M)
}

plot_ca_twoarm_survey <- function(ae_tbl,
                                  placebo = "Placebo",
                                  active  = "Xanomeline High Dose",
                                  cutoff = 0.05,
                                  ge3_only = FALSE,
                                  n_each_side = 6,
                                  wrap_width = 22,
                                  title = "Preferred Terms Patterns (CA summary)",
                                  subtitle = NULL,
                                  show_arm_points = TRUE,
                                  show_leaders = TRUE,
                                  pt_jitter = 0.06,
                                  label_x_nudge = 0.10) {

  built <- build_ca_table_two_arm(ae_tbl, placebo, active, cutoff, ge3_only)
  M <- built$M
  if (nrow(M) < 2 || ncol(M) < 2) stop("Too few PTs or arms for CA after filtering.")

  fit <- ca::ca(M)

  # --- PT coordinates (CA1) ---
  pt <- tibble::tibble(
    aedecod = as.character(rownames(fit$rowcoord)),
    ca1     = as.numeric(fit$rowcoord[, 1])
  ) %>%
    dplyr::filter(is.finite(ca1)) %>%
    dplyr::arrange(ca1, aedecod) %>%
    dplyr::mutate(
      # deterministic small y jitter (helps show ties / density)
      y0 = rep(c(-pt_jitter, pt_jitter), length.out = dplyr::n())
    )

  # --- Arm coordinates (optional) ---
  arm <- tibble::tibble(
    trt = as.character(rownames(fit$colcoord)),
    ca1 = as.numeric(fit$colcoord[, 1]),
    y   = 0
  ) %>%
    dplyr::filter(is.finite(ca1))

  inertia <- fit$sv^2
  pct1 <- round(100 * inertia[1] / sum(inertia), 1)

  # --- Label selection: extremes on each side ---
  left  <- pt %>% dplyr::arrange(ca1, aedecod)              %>% dplyr::slice_head(n = n_each_side)
  right <- pt %>% dplyr::arrange(dplyr::desc(ca1), aedecod) %>% dplyr::slice_head(n = n_each_side)

  lab <- dplyr::bind_rows(
    left  %>% dplyr::mutate(side = "left")  %>% dplyr::arrange(ca1, aedecod),
    right %>% dplyr::mutate(side = "right") %>% dplyr::arrange(dplyr::desc(ca1), aedecod)
  ) %>%
    dplyr::mutate(
      aedecod = as.character(aedecod) # coerce to character for safe join
    ) %>%
    dplyr::group_by(side) %>%
    dplyr::mutate(
      y = dplyr::row_number(),
      y = ifelse(side == "left", y, -y),
      pt_lab_raw = pt_pretty_wrap(aedecod, wrap_width),
      pt_lab = make_unique_labels(pt_lab_raw),
      hjust = ifelse(side == "left", 1, 0),
      x_lab = ca1 + ifelse(side == "left", -label_x_nudge, label_x_nudge),
      # default y0 = 0; overwritten by the join below
      y0 = 0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(pt %>% dplyr::select(aedecod, y0_join = y0), by = "aedecod") %>%
    dplyr::mutate(
      y0 = dplyr::coalesce(y0_join, y0)
    ) %>%
    dplyr::select(-y0_join)

  # arm colour mapping
  col_vals <- setNames(c(COL_PLACEBO, COL_ACTIVE), c(placebo, active))

  x_max <- max(abs(pt$ca1), na.rm = TRUE)
  x_pad <- max(0.25, x_max * 0.12)

  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
                      fill = "grey95", color = NA) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
                      fill = "aliceblue", color = NA) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, color = "grey45") +

    # all PT points (jittered y)
    ggplot2::geom_point(
      data = pt,
      ggplot2::aes(x = ca1, y = y0),
      size = 1.6, alpha = 0.22, color = "grey20"
    ) +

    # labeled PT points (use same y0)
    ggplot2::geom_point(
      data = lab,
      ggplot2::aes(x = ca1, y = y0),
      size = 2.6, color = "grey20"
    ) +

    # optional leader lines
    {if (show_leaders)
      ggplot2::geom_segment(
        data = lab,
        ggplot2::aes(x = ca1, xend = x_lab, y = y0, yend = y),
        linewidth = 0.25, color = "grey75"
      )
    } +

    ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = x_lab, y = y, label = pt_lab, hjust = hjust),
      size = 3.4, color = "grey15"
    ) +

    # optional arm points
    {if (show_arm_points)
      ggplot2::geom_point(
        data = arm,
        ggplot2::aes(x = ca1, y = 0, color = trt),
        size = 4.0
      )
    } +
    {if (show_arm_points)
      ggplot2::geom_text(
        data = arm,
        ggplot2::aes(x = ca1, y = 0, label = trt, color = trt),
        fontface = "bold", vjust = -1.2, size = 4.0, show.legend = FALSE
      )
    } +
    {if (show_arm_points)
      ggplot2::scale_color_manual(values = col_vals, breaks = c(placebo, active))
    } +

    ggplot2::annotate(
      "text", x = -Inf, y = max(abs(lab$y)) + 1.2,
      label = paste0("← more on ", placebo),
      hjust = -0.02, size = 3.4, color = "grey35"
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = max(abs(lab$y)) + 1.2,
      label = paste0("more on ", active, " →"),
      hjust = 1.02, size = 3.4, color = "grey35"
    ) +

    ggplot2::coord_cartesian(
      xlim = c(-x_max - x_pad, x_max + x_pad),
      ylim = c(-(n_each_side + 2), (n_each_side + 2))
    ) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle %||% paste0(
        round(cutoff * 100), "% Pure | CA1 explains ", pct1,
        "% of inertia (two-arm CA ⇒ one axis)"
      ),
      x = "CA Dimension 1",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = if (show_arm_points) "bottom" else "none",
      legend.direction = "horizontal"
    )
}

# --- Example ---
p_ca <- plot_ca_twoarm_survey(
  ae_tbl = ae_ctcae,
  placebo = "Placebo",
  active  = "Xanomeline High Dose",
  cutoff = 0.05,
  n_each_side = 6,
  wrap_width = 20,
  show_arm_points = FALSE,
  show_leaders = TRUE
)

p_ca


#the first CA dimension captures the entirety of the between-arm variation, highlighting a focused and internally consistent tolerability pattern.”



