# ============================================================
# SafetyLens — Adverse Event Reporting & Information Loss
# Boehringer Ingelheim / Otto-Friedrich University of Bamberg
# Viola Schenk — Master Thesis 2026
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(janitor)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(ggh4x)
  library(ca)
  library(ggrepel)
  library(DT)
  library(plotly)
})

# ============================================================
# 0) UTILITY FUNCTIONS
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

first_present <- function(nms, choices) {
  hit <- intersect(choices, nms)
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

norm_str <- function(x) {
  x %>% as.character() %>% str_squish() %>% toupper() %>% na_if("")
}

coalesce_date <- function(df, candidates) {
  out <- rep(NA_Date_, nrow(df))
  for (nm in candidates) {
    if (!nm %in% names(df)) next
    x  <- df[[nm]]
    y  <- if (inherits(x, "Date")) x else if (inherits(x, "POSIXt")) as_date(x) else as.character(x)
    y1 <- suppressWarnings(lubridate::ymd_hms(y, quiet = TRUE))
    y2 <- suppressWarnings(lubridate::ymd(y,     quiet = TRUE))
    yf <- dplyr::coalesce(as_date(y1), as_date(y2))
    idx <- is.na(out) & !is.na(yf)
    out[idx] <- yf[idx]
  }
  out
}

denom_from_ae <- function(ae_tbl, trt_col = "trt") {
  ae_tbl %>%
    distinct(.data[[trt_col]], usubjid) %>%
    count(.data[[trt_col]], name = "n_trt") %>%
    rename(trt = .data[[trt_col]])
}

# ============================================================
# 1) DATA LOADING & PREPROCESSING
# ============================================================

load_and_process_data <- function() {
  if (requireNamespace("pharmaverseadam", quietly = TRUE)) {
    adsl <- pharmaverseadam::adsl
    adae <- pharmaverseadam::adae
  } else if (requireNamespace("admiral.test", quietly = TRUE)) {
    adsl <- admiral.test::adsl
    adae <- admiral.test::adae
  } else {
    stop("Install 'pharmaverseadam' (preferred) or 'admiral.test' (fallback).")
  }
  
  trt_var     <- first_present(names(adsl), c("TRT01A", "TRTA", "TRT01P"))
  adsl$TRTSDT <- coalesce_date(adsl, c("TRTSDT", "TRTSDTM", "TRSDTM", "TRSDT"))
  adsl1 <- adsl %>%
    transmute(USUBJID, TRT = as.character(.data[[trt_var]]), TRTSDT) %>%
    clean_names()
  
  aedecod_var  <- first_present(names(adae), c("AEDECOD"))
  aebodsys_var <- first_present(names(adae), c("AEBODSYS"))
  aesev_var    <- first_present(names(adae), c("AESEV"))
  teae_var     <- first_present(names(adae), c("TEAEFL"))
  
  adae$ASTDT  <- coalesce_date(adae, c("ASTDT",  "ASTDTM",  "AESTDTC"))
  adae$AEENDT <- coalesce_date(adae, c("AEENDT", "AENDTM",  "AEENDTC"))
  aseq <- if ("AESEQ" %in% names(adae)) adae$AESEQ else
    if ("ASEQ"  %in% names(adae)) adae$ASEQ  else NA_real_
  
  adae1 <- adae %>%
    transmute(
      USUBJID,
      AEDECOD   = if (!is.na(aedecod_var))  .data[[aedecod_var]]  else NA_character_,
      AEBODSYS  = if (!is.na(aebodsys_var)) .data[[aebodsys_var]] else NA_character_,
      AESEV_RAW = if (!is.na(aesev_var))    .data[[aesev_var]]    else NA_character_,
      ASTDT, AEENDT,
      TEAEFL = if (!is.na(teae_var)) .data[[teae_var]] else NA_character_,
      AESEQ  = aseq
    ) %>%
    mutate(
      AEDECOD   = norm_str(AEDECOD),
      AEBODSYS  = norm_str(AEBODSYS),
      AESEV_RAW = norm_str(AESEV_RAW)
    ) %>%
    clean_names()
  
  sev_map <- list(
    `1` = c("MILD", "GRADE 1", "1", "LOW"),
    `2` = c("MODERATE", "MOD", "GRADE 2", "2", "MEDIUM"),
    `3` = c("SEVERE", "SEV", "GRADE 3", "3"),
    `4` = c("LIFE-THREATENING", "LIFE THREATENING", "VERY SEVERE", "GRADE 4", "4"),
    `5` = c("DEATH", "FATAL", "GRADE 5", "5")
  )
  sev_to_grade <- function(x) {
    x   <- norm_str(x)
    out <- rep(NA_integer_, length(x))
    for (g in names(sev_map)) out[x %in% sev_map[[g]]] <- as.integer(g)
    out
  }
  adae2 <- adae1 %>% mutate(grade = sev_to_grade(aesev_raw))
  
  adsl_dosed <- adsl1 %>% filter(!is.na(trtsdt))
  ae0 <- adae2 %>%
    inner_join(adsl_dosed, by = "usubjid") %>%
    mutate(teaefl = case_when(
      !is.na(teaefl)                                    ~ teaefl,
      !is.na(astdt) & !is.na(trtsdt) & astdt >= trtsdt ~ "Y",
      TRUE                                              ~ "N"
    ))
  ae_teae <- ae0 %>% filter(teaefl == "Y")
  
  ae_clean <- ae_teae %>%
    mutate(
      onset_day     = if_else(!is.na(astdt) & !is.na(trtsdt),
                              as.integer(astdt - trtsdt) + 1L, NA_integer_),
      duration_days = if_else(!is.na(astdt) & !is.na(aeendt),
                              as.integer(aeendt - astdt) + 1L, NA_integer_)
    ) %>%
    mutate(
      duration_days = if_else(!is.na(duration_days) & duration_days < 0L,
                              NA_integer_, duration_days),
      onset_day     = if_else(!is.na(onset_day)     & onset_day     < 0L,
                              NA_integer_, onset_day)
    ) %>%
    mutate(
      trt      = factor(trt),
      aedecod  = factor(aedecod),
      aebodsys = factor(aebodsys)
    ) %>%
    distinct(across(any_of(c(
      "usubjid", "trt", "aedecod", "aebodsys", "astdt", "aeendt",
      "grade", "onset_day", "duration_days", "aesev_raw", "teaefl", "aeseq"
    )))) %>%
    arrange(trt, usubjid, aedecod, astdt)
  
  ae_ctcae <- ae_clean %>%
    mutate(ctcae_grade = grade, ctcae_ge3 = !is.na(grade) & grade >= 3)
  
  denom_all_treated <- adsl1 %>%
    filter(!is.na(trtsdt)) %>%
    count(trt, name = "n_all_treated")
  
  list(ae_ctcae = ae_ctcae, adsl1 = adsl1, denom_all_treated = denom_all_treated)
}

# ============================================================
# 2) THRESHOLD FUNCTIONS
# ============================================================

incidence_by_term <- function(ae, key = "aedecod", denom_tbl = NULL) {
  denom_tbl <- denom_tbl %||% denom_from_ae(ae)
  ae %>%
    distinct(trt, usubjid, .data[[key]]) %>%
    count(trt, .data[[key]], name = "n_subj") %>%
    left_join(denom_tbl, by = "trt") %>%
    mutate(incidence = ifelse(n_trt > 0, n_subj / n_trt, 0)) %>%
    rename(term = .data[[key]])
}

apply_threshold_any_arm <- function(ae, cutoff, key = "aedecod") {
  denom_tbl  <- denom_from_ae(ae)
  inc        <- incidence_by_term(ae, key = key, denom_tbl = denom_tbl)
  keep_terms <- inc %>%
    group_by(term) %>%
    summarise(max_inc = max(incidence, na.rm = TRUE), .groups = "drop") %>%
    filter(max_inc >= cutoff) %>%
    pull(term) %>% as.character()
  ae %>% filter(as.character(.data[[key]]) %in% keep_terms)
}

apply_threshold_with_safeguard <- function(ae, cutoff) {
  pt_ge3    <- ae %>% filter(ctcae_ge3) %>% distinct(aedecod) %>%
    pull(aedecod) %>% as.character()
  keep_freq <- apply_threshold_any_arm(ae, cutoff) %>%
    distinct(aedecod) %>% pull(aedecod) %>% as.character()
  ae %>% filter(as.character(aedecod) %in% union(keep_freq, pt_ge3))
}

apply_threshold <- function(ae, cutoff, variant, key = "aedecod") {
  if (variant == "safeguard") {
    apply_threshold_with_safeguard(ae, cutoff)
  } else {
    apply_threshold_any_arm(ae, cutoff = cutoff, key = key)
  }
}

termset_filtered <- function(ae_tbl, cutoff, variant, key = "aedecod") {
  apply_threshold(ae_tbl, cutoff = cutoff, variant = variant, key = key) %>%
    distinct(.data[[key]]) %>% pull(.data[[key]]) %>% as.character()
}

make_incidence <- function(ae_tbl, key = "aedecod", denom_tbl = NULL) {
  denom_tbl <- denom_tbl %||% denom_from_ae(ae_tbl)
  ae_tbl %>%
    distinct(trt, usubjid, .data[[key]]) %>%
    count(trt, .data[[key]], name = "n") %>%
    left_join(denom_tbl, by = "trt") %>%
    mutate(inc = ifelse(n_trt > 0, n / n_trt, 0)) %>%
    transmute(trt, term = as.character(.data[[key]]), n, n_trt, inc)
}

# ============================================================
# 3) INFORMATION-LOSS METRICS
# ============================================================

entropy  <- function(p) { p <- p[!is.na(p) & p > 0]; -sum(p * log(p)) }
evenness <- function(p) {
  p <- p[!is.na(p) & p > 0]; K <- length(p)
  if (K <= 1) NA_real_ else entropy(p) / log(K)
}
kl_div <- function(p, q) {
  idx <- p > 0 & q > 0
  sum(p[idx] * log(p[idx] / q[idx]))
}
jsd <- function(p, q) {
  p <- p / sum(p); q <- q / sum(q); m <- 0.5 * (p + q)
  0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
}

build_distributions <- function(ae_tbl, type = "events") {
  base <- ae_tbl %>% mutate(aedecod = as.character(aedecod))
  if (type == "events") {
    base %>% count(trt, aedecod, name = "n") %>%
      group_by(trt) %>% mutate(p = n / sum(n)) %>% ungroup()
  } else {
    base %>% distinct(trt, usubjid, aedecod) %>%
      count(trt, aedecod, name = "n") %>%
      group_by(trt) %>% mutate(p = n / sum(n)) %>% ungroup()
  }
}

# ============================================================
# 4) COLOURS & SHARED THEME
# ============================================================

COL_PLACEBO <- "#7A7A7A"
COL_ACTIVE  <- "#0072B2"
COL_G3      <- "#D55E00"
COL_NONG3   <- "#8C8C8C"

theme_thesis <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 1,
                                       colour = "#003865", margin = margin(b = 4)),
      plot.subtitle     = element_text(colour = "grey45", size = base_size - 1,
                                       margin = margin(b = 8)),
      plot.caption      = element_text(colour = "grey55", size = base_size - 2,
                                       hjust = 0),
      plot.background   = element_rect(fill = "white",   colour = NA),
      panel.background  = element_rect(fill = "#fafbfc", colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "#eeeeee", linewidth = 0.4),
      panel.border      = element_rect(colour = "#dddddd", fill = NA, linewidth = 0.4),
      strip.text        = element_text(face = "bold", colour = "#003865"),
      strip.background  = element_rect(fill = "#f0f4f8", colour = NA),
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", colour = NA),
      axis.title        = element_text(size = base_size - 1, colour = "grey30"),
      axis.text         = element_text(size = base_size - 2, colour = "grey40"),
      plot.margin       = margin(14, 14, 14, 14)
    )
}

# ============================================================
# 5) STARTUP DATA LOAD
# ============================================================

message("Loading and processing data...")
data_env <- tryCatch(load_and_process_data(), error = function(e) { message(e); NULL })

if (is.null(data_env)) {
  ae_ctcae          <- data.frame()
  adsl1             <- data.frame()
  denom_all_treated <- data.frame()
  arm_choices       <- character(0)
} else {
  ae_ctcae          <- data_env$ae_ctcae
  adsl1             <- data_env$adsl1
  denom_all_treated <- data_env$denom_all_treated
  arm_choices       <- levels(ae_ctcae$trt)
}
message("Arms: ", paste(arm_choices, collapse = " | "))

# ============================================================
# 6) UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = tags$span(
      "SafetyLens",
      style = "font-size:15px; font-weight:700; letter-spacing:0.03em;"
    ),
    titleWidth = 290
  ),
  
  dashboardSidebar(
    width = 290,
    
    sidebarMenu(
      id = "sidebar_tabs",
      menuItem("Overview",          tabName = "overview",  icon = icon("gauge-high")),
      menuItem("Butterfly Plot",    tabName = "butterfly", icon = icon("chart-bar")),
      menuItem("Volcano Plot",      tabName = "volcano",   icon = icon("mountain")),
      menuItem("Temporal Trends",   tabName = "temporal",  icon = icon("chart-line")),
      menuItem("Co-occurrence",     tabName = "heatmap",   icon = icon("th")),
      menuItem("CA Strip Plot",     tabName = "castrip",   icon = icon("grip-lines")),
      menuItem("Info-Loss Metrics", tabName = "metrics",   icon = icon("table"))
    ),
    
    tags$hr(style = "border-color:#1e3050; margin:6px 0;"),
    
    tags$div(
      style = "padding:4px 15px 2px 15px; color:#7eb8d4; font-size:10px;
               font-weight:700; letter-spacing:0.08em; text-transform:uppercase;",
      "Treatment Arms"
    ),
    div(style = "padding:0 12px;",
        selectInput("arm_placebo", "Reference arm:",
                    choices  = arm_choices,
                    selected = grep("PLACEBO|Placebo", arm_choices,
                                    value = TRUE, ignore.case = TRUE)[1] %||% arm_choices[1]),
        selectInput("arm_active", "Active arm:",
                    choices  = arm_choices,
                    selected = arm_choices[length(arm_choices)])
    ),
    
    tags$hr(style = "border-color:#1e3050; margin:6px 0;"),
    
    tags$div(
      style = "padding:4px 15px 2px 15px; color:#7eb8d4; font-size:10px;
               font-weight:700; letter-spacing:0.08em; text-transform:uppercase;",
      "Reporting Threshold"
    ),
    div(style = "padding:0 12px;",
        sliderInput("threshold", label = NULL,
                    min = 1, max = 15, value = 5, step = 1, post = "%"),
        radioButtons("variant", label = NULL,
                     choices  = c("Pure any-arm"        = "pure",
                                  "Keep \u22653 safeguard" = "safeguard"),
                     selected = "pure",
                     inline   = FALSE)
    ),
    
    conditionalPanel(
      condition = "input.variant == 'safeguard'",
      div(
        style = "margin:6px 12px 10px 12px; padding:8px 12px;
                 background:linear-gradient(135deg,#1a4a2e,#0f3020);
                 border-left:3px solid #2ecc71; border-radius:6px;
                 font-size:11px; color:#aee8c0; line-height:1.5;",
        icon("shield-halved", style = "color:#2ecc71; margin-right:6px;"),
        tags$b("Safeguard active:"), tags$br(),
        "All Grade \u22653 PTs retained regardless of incidence."
      )
    ),
    
    tags$hr(style = "border-color:#1e3050; margin:6px 0;"),
    
    tags$div(
      style = "padding:8px 15px 12px 15px; color:#5a7a9a; font-size:11px; line-height:1.7;",
      tags$div(style = "font-weight:700; color:#7eb8d4; margin-bottom:2px;",
               "Master Thesis 2026"),
      tags$div("Viola Schenk"),
      tags$div("Boehringer Ingelheim"),
      tags$div("Otto-Friedrich Universit\u00e4t Bamberg"),
      tags$div("M.Sc. Survey Statistics and Data Analysis")
    )
  ),
  
  dashboardBody(
    
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo,
      .skin-blue .main-header .navbar        { background-color:#003865 !important; }
      .skin-blue .main-sidebar               { background-color:#0d1929 !important; }
      .skin-blue .sidebar-menu > li.active > a,
      .skin-blue .sidebar-menu > li:hover  > a {
        background-color:#0072B2 !important;
        border-left-color:#41b0e8 !important;
        border-left-width:4px !important;
      }
      .skin-blue .sidebar-menu > li > a     { color:#c8d8e8 !important; font-size:13px; }
      .content-wrapper, .right-side          { background-color:#f4f6f9 !important; }
      .box                                   { border-radius:6px;
                                               box-shadow:0 1px 4px rgba(0,0,0,.08); }
      .box.box-primary                       { border-top-color:#0072B2; }
      .box-header .box-title                 { font-weight:700; font-size:14px; color:#003865; }
      .small-box                             { border-radius:8px;
                                               box-shadow:0 2px 6px rgba(0,0,0,.12); }
      .small-box h3                          { font-size:2rem; font-weight:700; }
      .small-box p                           { font-size:12px; }
      .irs--shiny .irs-bar,
      .irs--shiny .irs-bar-edge              { background:#0072B2 !important;
                                               border-color:#0072B2 !important; }
      .irs--shiny .irs-handle                { border-color:#0072B2 !important; }
      .irs--shiny .irs-single                { background:#0072B2 !important; }
      .btn-dl                                { margin-top:10px; width:100%;
                                               background:#0072B2; color:white;
                                               border:none; border-radius:4px;
                                               padding:6px 12px; font-size:12px; }
      .btn-dl:hover                          { background:#005a8e; color:white; }
      .radio label                           { font-size:12px; color:#c8d8e8; }
      .dataTables_wrapper                    { font-size:12px; }
      table.dataTable thead th               { background:#f0f4f8; color:#003865;
                                               font-weight:700; }
      table.dataTable tbody tr:hover         { background:#edf4fb !important; }
      .tab-content                           { padding-top:6px; }
      .hint-text                             { color:#888; font-size:11px;
                                               line-height:1.5; margin-top:4px; }
    "))),
    
    tabItems(
      
      # ── Overview ──────────────────────────────────────────────────────────
      tabItem(tabName = "overview",
              fluidRow(
                valueBoxOutput("vbox_subjects", width = 3),
                valueBoxOutput("vbox_events",   width = 3),
                valueBoxOutput("vbox_pts",      width = 3),
                valueBoxOutput("vbox_ge3_pct",  width = 3)
              ),
              fluidRow(
                box(title = "Retained AE Data by Arm", width = 12,
                    solidHeader = TRUE, status = "primary",
                    DTOutput("overview_table"))
              ),
              fluidRow(
                box(
                  title = tagList(
                    "Grade \u22653 Clinical Visibility",
                    tags$small(" \u2014 unconditional % of all treated subjects, pure threshold",
                               style = "color:#003865; font-weight:800;")
                  ),
                  width = 12, solidHeader = TRUE, status = "primary",
                  plotlyOutput("overview_ge3_plot", height = "340px")
                )
              )
      ),
      
      # ── Butterfly ─────────────────────────────────────────────────────────
      tabItem(tabName = "butterfly",
              fluidRow(
                box(title = "Options", width = 3, solidHeader = TRUE, status = "primary",
                    numericInput("bf_n_common", "Top N by max incidence:", value = 10, min = 3, max = 30),
                    numericInput("bf_n_diff",   "Top N by |risk diff|:",   value = 5,  min = 0, max = 20),
                    checkboxInput("bf_ge3_only", "Grade \u22653 events only", value = FALSE),
                    tags$hr(),
                    tags$p("Three-panel layout: reference | PT labels | active.", class = "hint-text"),
                    tags$p("PT selection = top-by-incidence \u222a top-by-|RD|.", class = "hint-text"),
                    downloadButton("dl_butterfly", "Download PNG", class = "btn-dl btn-sm")
                ),
                box(title = "Butterfly Plot \u2014 Incidence Comparison",
                    width = 9, solidHeader = TRUE, status = "primary",
                    plotlyOutput("butterfly_plot", height = "660px"))
              )
      ),
      
      # ── Volcano ───────────────────────────────────────────────────────────
      tabItem(tabName = "volcano",
              fluidRow(
                box(title = "Options", width = 3, solidHeader = TRUE, status = "primary",
                    numericInput("vol_label_n", "Max labels per side:", value = 6, min = 2, max = 20),
                    checkboxInput("vol_label_sig_only", "Label significant PTs only", value = FALSE),
                    selectInput("vol_p_adjust", "P-value adjustment:",
                                choices  = c("None" = "none", "BH (FDR)" = "BH",
                                             "Holm" = "holm", "Bonferroni" = "bonferroni"),
                                selected = "none"),
                    numericInput("vol_p_line", "Significance threshold (\u03b1):",
                                 value = 0.05, min = 0.001, max = 0.2, step = 0.005),
                    tags$hr(),
                    tags$p("x = severity burden diff (Active \u2212 Reference).", class = "hint-text"),
                    tags$p("Orange = PT has \u2265 one Grade 3 record.", class = "hint-text"),
                    downloadButton("dl_volcano", "Download PNG", class = "btn-dl btn-sm")
                ),
                box(title = "Volcano Plot \u2014 Severity-Weighted PT Signals",
                    width = 9, solidHeader = TRUE, status = "primary",
                    plotlyOutput("volcano_plot", height = "600px"))
              )
      ),
      
      # ── Temporal ──────────────────────────────────────────────────────────
      tabItem(tabName = "temporal",
              fluidRow(
                box(title = "Options", width = 3, solidHeader = TRUE, status = "primary",
                    numericInput("temp_top_n", "Top N PTs (by |\u0394 end-study|):",
                                 value = 9, min = 3, max = 16),
                    numericInput("temp_ncol", "Facet columns:", value = 3, min = 1, max = 4),
                    tags$hr(),
                    tags$b("Selected PTs:", style = "font-size:11px; color:#555;"),
                    tags$div(style = "margin-top:6px;", DTOutput("temp_pt_table")),
                    tags$br(),
                    downloadButton("dl_temporal", "Download PNG", class = "btn-dl btn-sm")
                ),
                box(title = "Temporal Cumulative Incidence",
                    width = 9, solidHeader = TRUE, status = "primary",
                    plotlyOutput("temporal_plot", height = "660px"))
              )
      ),
      
      # ── Heatmap ───────────────────────────────────────────────────────────
      tabItem(tabName = "heatmap",
              fluidRow(
                box(title = "Options", width = 3, solidHeader = TRUE, status = "primary",
                    selectInput("hm_arm", "Arm:", choices = arm_choices,
                                selected = arm_choices[length(arm_choices)]),
                    numericInput("hm_top_n", "Max PTs shown:",    value = 25, min = 5,   max = 50),
                    numericInput("hm_cap",   "Colour scale cap:", value = 0.6, min = 0.1, max = 1,
                                 step = 0.05),
                    tags$hr(),
                    tags$p("Jaccard(A,B) = |A \u2229 B| / |A \u222a B|.", class = "hint-text"),
                    tags$p("0 = no shared subjects; 1 = identical subject sets.", class = "hint-text"),
                    downloadButton("dl_heatmap", "Download PNG", class = "btn-dl btn-sm")
                ),
                box(title = "PT Co-occurrence Heatmap (Jaccard Similarity)",
                    width = 9, solidHeader = TRUE, status = "primary",
                    plotlyOutput("heatmap_plot", height = "690px"))
              )
      ),
      
      # ── CA Strip ──────────────────────────────────────────────────────────
      tabItem(tabName = "castrip",
              fluidRow(
                box(title = "Options", width = 3, solidHeader = TRUE, status = "primary",
                    numericInput("ca_n_each", "Labels per side:", value = 6, min = 2, max = 15),
                    checkboxInput("ca_ge3_only",     "Grade \u22653 events only", value = FALSE),
                    checkboxInput("ca_show_arms",    "Show arm centroids",        value = FALSE),
                    checkboxInput("ca_show_leaders", "Show leader lines",         value = TRUE),
                    tags$hr(),
                    tags$p("Two-arm CA: Dim 1 captures 100% of between-arm inertia.", class = "hint-text"),
                    tags$p("Left = more associated with reference; right = active arm.", class = "hint-text"),
                    downloadButton("dl_castrip", "Download PNG", class = "btn-dl btn-sm")
                ),
                box(title = "CA Strip Plot \u2014 PT Separation Along CA Dimension 1",
                    width = 9, solidHeader = TRUE, status = "primary",
                    plotlyOutput("castrip_plot", height = "560px"))
              )
      ),
      
      # ── Info-Loss Metrics ─────────────────────────────────────────────────
      tabItem(tabName = "metrics",
              fluidRow(
                box(
                  title = "Information-Loss Metrics across Thresholds",
                  width = 12, solidHeader = TRUE, status = "primary",
                  tags$p(
                    "Metrics computed at 1%, 2%, 5% and 10% for both pure and safeguard variants.
               Coverage is normalized to the 1% pure baseline.
               JSD = 0 means structurally identical to that baseline.",
                    style = "color:#555; font-size:12px; margin-bottom:12px;"
                  ),
                  fluidRow(
                    column(4,
                           plotlyOutput("metrics_coverage_plot", height = "270px")),
                    column(4,
                           plotlyOutput("metrics_entropy_plot", height = "270px")),
                    column(4,
                           plotlyOutput("metrics_jsd_plot", height = "270px"))
                  )
                )
              ),
              fluidRow(
                box(title = "Full Metrics Table", width = 12,
                    solidHeader = TRUE, status = "primary",
                    DTOutput("metrics_table"))
              )
      )
      
    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

# ============================================================
# 7) SERVER
# ============================================================

server <- function(input, output, session) {
  
  ae_filtered <- reactive({
    req(nrow(ae_ctcae) > 0)
    apply_threshold(ae_ctcae,
                    cutoff  = input$threshold / 100,
                    variant = input$variant)
  })
  
  # ── Value boxes ─────────────────────────────────────────────────────────
  output$vbox_subjects <- renderValueBox({
    n <- ae_filtered() %>% distinct(usubjid) %>% nrow()
    valueBox(n, "Subjects in filtered data", icon = icon("users"), color = "blue")
  })
  output$vbox_events <- renderValueBox({
    valueBox(nrow(ae_filtered()), "AE records retained",
             icon = icon("file-medical"), color = "teal")
  })
  output$vbox_pts <- renderValueBox({
    n <- ae_filtered() %>% distinct(aedecod) %>% nrow()
    valueBox(n, "Unique PTs retained", icon = icon("list"), color = "navy")
  })
  output$vbox_ge3_pct <- renderValueBox({
    ge3_n <- ae_filtered() %>% filter(ctcae_ge3) %>% distinct(usubjid) %>% nrow()
    all_n <- if (nrow(denom_all_treated) > 0)
      denom_all_treated %>% summarise(n = sum(n_all_treated)) %>% pull(n) else 0
    pct   <- if (all_n > 0) round(100 * ge3_n / all_n, 1) else 0
    valueBox(paste0(pct, "%"), "Subjects with \u22653 TEAE (all treated)",
             icon = icon("heart-pulse"), color = "red")
  })
  
  # ── Overview table ───────────────────────────────────────────────────────
  output$overview_table <- renderDT({
    df <- ae_filtered() %>%
      group_by(trt) %>%
      summarise(Subjects   = n_distinct(usubjid),
                AE_Records = n(),
                Unique_PTs = n_distinct(aedecod),
                Pct_Ge3    = round(100 * mean(ctcae_ge3, na.rm = TRUE), 1),
                .groups    = "drop")
    datatable(df,
              options  = list(dom = "t", pageLength = 10),
              rownames = FALSE,
              colnames = c("Treatment Arm", "Subjects", "AE Records",
                           "Unique PTs", "% \u22653"))
  })
  
  # ── Overview >=G3 plot ───────────────────────────────────────────────────
  output$overview_ge3_plot <- renderPlotly({
    THRESHOLDS <- c(0.01, 0.02, 0.05, 0.10)
    LABELS     <- c("1%", "2%", "5%", "10%")
    df <- map_dfr(seq_along(THRESHOLDS), function(i) {
      ae_thr <- apply_threshold_any_arm(ae_ctcae, THRESHOLDS[i])
      ae_thr %>%
        filter(ctcae_ge3) %>%
        distinct(trt, usubjid) %>%
        count(trt, name = "n_ge3") %>%
        right_join(denom_all_treated, by = "trt") %>%
        mutate(n_ge3      = replace_na(n_ge3, 0L),
               pct_ge3    = 100 * n_ge3 / n_all_treated,
               threshold  = LABELS[i])
    }) %>%
      mutate(threshold = factor(threshold, levels = LABELS))
    
    p <- ggplot(df, aes(x = threshold, y = pct_ge3, group = trt, colour = trt,
                        text = paste0(trt, "<br>Threshold: ", threshold,
                                      "<br>", round(pct_ge3, 1), "%"))) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 3) +
      scale_y_continuous(labels = percent_format(scale = 1, accuracy = 0.1),
                         expand = expansion(mult = c(0.05, 0.18))) +
      scale_colour_brewer(palette = "Set1") +
      labs(x = "Reporting threshold (pure any-arm)",
           y = "% of all treated with \u22653 TEAE",
           colour = "Treatment arm") +
      theme_thesis(13) +
      theme(legend.position = "top")
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = 5))
  })
  
  # ── Butterfly ─────────────────────────────────────────────────────────────
  make_butterfly_plot <- function() {
    req(nrow(ae_ctcae) > 0)
    cutoff    <- input$threshold / 100
    variant   <- input$variant
    arm_left  <- input$arm_placebo
    arm_right <- input$arm_active
    
    ae_thr <- apply_threshold(ae_ctcae, cutoff, variant)
    ae_use <- if (input$bf_ge3_only) ae_thr %>% filter(ctcae_ge3) else ae_thr
    denom  <- denom_from_ae(ae_ctcae)
    inc    <- make_incidence(ae_use, denom_tbl = denom)
    
    wide <- inc %>%
      filter(trt %in% c(arm_left, arm_right)) %>%
      select(trt, term, inc) %>%
      pivot_wider(names_from = trt, values_from = inc, values_fill = 0)
    if (!arm_left  %in% names(wide)) wide[[arm_left]]  <- 0
    if (!arm_right %in% names(wide)) wide[[arm_right]] <- 0
    wide <- wide %>% mutate(
      inc_left  = .data[[arm_left]],
      inc_right = .data[[arm_right]],
      abs_rd    = abs(inc_right - inc_left),
      max_inc   = pmax(inc_left, inc_right)
    )
    
    sel <- bind_rows(
      wide %>% arrange(desc(max_inc)) %>% slice_head(n = input$bf_n_common),
      wide %>% arrange(desc(abs_rd))  %>% slice_head(n = input$bf_n_diff)
    ) %>%
      distinct(term, .keep_all = TRUE) %>%
      arrange(desc(inc_right), desc(abs_rd)) %>%
      mutate(term = fct_reorder(term, inc_right, .desc = TRUE))
    
    if (nrow(sel) == 0) return(ggplot() +
                                 annotate("text", x = 0.5, y = 0.5,
                                          label = "No PTs retained at this threshold.",
                                          size = 5, colour = "grey50") + theme_void())
    
    wrap_w       <- 26
    center_label <- "Preferred Term"
    facet_levels <- c(arm_left, center_label, arm_right)
    
    bars <- bind_rows(
      sel %>% transmute(term, facet = arm_left,  arm = arm_left,
                        inc = inc_left,  inc_plot = -inc_left),
      sel %>% transmute(term, facet = arm_right, arm = arm_right,
                        inc = inc_right, inc_plot =  inc_right)
    ) %>%
      mutate(term_w = str_wrap(as.character(term), wrap_w),
             facet  = factor(facet, levels = facet_levels))
    
    pts <- sel %>%
      transmute(term, facet = center_label, label = as.character(term)) %>%
      mutate(term_w = str_wrap(as.character(term), wrap_w),
             label  = str_wrap(label, wrap_w),
             facet  = factor(facet, levels = facet_levels))
    
    lvl         <- unique(bars$term_w)
    bars$term_w <- factor(bars$term_w, levels = lvl)
    pts$term_w  <- factor(pts$term_w,  levels = lvl)
    
    is_pbo    <- function(x) grepl("PLACEBO|PBO", toupper(x))
    fill_vals <- setNames(
      c(if (is_pbo(arm_left))  COL_PLACEBO else COL_ACTIVE,
        if (is_pbo(arm_right)) COL_PLACEBO else COL_ACTIVE),
      c(arm_left, arm_right)
    )
    max_lim <- max(sel$max_inc, na.rm = TRUE)
    
    lr_blank <- tibble(
      facet  = factor(c(arm_left, arm_left, arm_right, arm_right), levels = facet_levels),
      term_w = factor(lvl[1], levels = lvl),
      x      = c(-max_lim, 0, 0, max_lim)
    )
    center_blank <- tibble(
      facet  = factor(center_label, levels = facet_levels),
      term_w = factor(lvl[1], levels = lvl),
      x      = c(-0.12, 0.12)
    )
    
    vlbl <- if (variant == "safeguard") " | \u22653 safeguard" else ""
    
    ggplot() +
      geom_col(data = bars, aes(x = inc_plot, y = term_w, fill = arm), width = 0.65) +
      geom_text(data = bars,
                aes(x = inc_plot, y = term_w,
                    label = percent(inc, accuracy = 1),
                    hjust = ifelse(inc_plot < 0, 1.12, -0.12)),
                size = 2.9, colour = "grey25") +
      geom_text(data = pts, aes(x = 0, y = term_w, label = label),
                size = 3.4, hjust = 0.5, colour = "#1a2332") +
      geom_blank(data = lr_blank,     aes(x = x, y = term_w)) +
      geom_blank(data = center_blank, aes(x = x, y = term_w)) +
      ggh4x::facet_grid2(. ~ facet, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = fill_vals) +
      scale_x_continuous(breaks = NULL, labels = NULL,
                         expand = expansion(mult = 0.02)) +
      labs(
        title    = paste0("PT Incidence | ", input$threshold, "% threshold", vlbl,
                          if (input$bf_ge3_only) " | Grade \u22653 only" else " | All grades"),
        subtitle = paste(arm_left, "vs", arm_right),
        x = NULL, y = NULL, fill = NULL
      ) +
      theme_thesis(12) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        legend.position    = "none",
        strip.text         = element_blank(),
        strip.background   = element_blank(),
        panel.spacing.x    = unit(0.04, "lines"),
        plot.margin        = margin(16, 14, 26, 14)
      ) +
      coord_cartesian(clip = "off")
  }
  
  output$butterfly_plot <- renderPlotly({
    req(nrow(ae_ctcae) > 0)
    cutoff    <- input$threshold / 100
    variant   <- input$variant
    arm_left  <- input$arm_placebo
    arm_right <- input$arm_active

    ae_thr <- apply_threshold(ae_ctcae, cutoff, variant)
    ae_use <- if (input$bf_ge3_only) ae_thr %>% filter(ctcae_ge3) else ae_thr
    denom  <- denom_from_ae(ae_ctcae)
    inc    <- make_incidence(ae_use, denom_tbl = denom)

    wide <- inc %>%
      filter(trt %in% c(arm_left, arm_right)) %>%
      select(trt, term, inc) %>%
      pivot_wider(names_from = trt, values_from = inc, values_fill = 0)
    if (!arm_left  %in% names(wide)) wide[[arm_left]]  <- 0
    if (!arm_right %in% names(wide)) wide[[arm_right]] <- 0
    wide <- wide %>% mutate(
      inc_left  = .data[[arm_left]],
      inc_right = .data[[arm_right]],
      abs_rd    = abs(inc_right - inc_left),
      max_inc   = pmax(inc_left, inc_right)
    )

    sel <- bind_rows(
      wide %>% arrange(desc(max_inc)) %>% slice_head(n = input$bf_n_common),
      wide %>% arrange(desc(abs_rd))  %>% slice_head(n = input$bf_n_diff)
    ) %>%
      distinct(term, .keep_all = TRUE) %>%
      arrange(inc_right, abs_rd) %>%
      mutate(term_label = str_wrap(str_to_title(as.character(term)), 26))

    if (nrow(sel) == 0) return(plotly_empty() %>% layout(title = "No PTs retained at this threshold."))

    is_pbo <- function(x) grepl("PLACEBO|PBO", toupper(x))
    col_left  <- if (is_pbo(arm_left))  COL_PLACEBO else COL_ACTIVE
    col_right <- if (is_pbo(arm_right)) COL_PLACEBO else COL_ACTIVE
    vlbl <- if (variant == "safeguard") " | \u22653 safeguard" else ""

    max_inc_val <- max(c(sel$inc_left, sel$inc_right), na.rm = TRUE)

    # Left panel — bars go LEFT from 0
    p_left <- plot_ly(sel, x = ~-inc_left, y = ~term_label,
                      type = "bar", orientation = "h",
                      marker = list(color = col_left),
                      text  = ~paste0(round(inc_left * 100, 1), "%"),
                      textposition = "outside",
                      hovertemplate = paste0(arm_left, ": %{text}<extra></extra>"),
                      showlegend = FALSE) %>%
      layout(xaxis = list(
               title = arm_left,
               range = c(-max_inc_val * 1.35, 0),
               tickvals = list(0, -max_inc_val * 0.25, -max_inc_val * 0.5,
                               -max_inc_val * 0.75, -max_inc_val),
               ticktext = list("0%",
                               paste0(round(max_inc_val * 25), "%"),
                               paste0(round(max_inc_val * 50), "%"),
                               paste0(round(max_inc_val * 75), "%"),
                               paste0(round(max_inc_val * 100), "%")),
               showgrid = FALSE),
             yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE,
                          categoryorder = "array",
                          categoryarray = rev(sel$term_label)))

    # Center label panel
    p_center <- plot_ly(sel, x = rep(0, nrow(sel)), y = ~term_label,
                        type = "scatter", mode = "text",
                        text = ~term_label,
                        textfont = list(size = 10, color = "#1a2332"),
                        hoverinfo = "none",
                        showlegend = FALSE) %>%
      layout(xaxis = list(visible = FALSE, range = c(-1, 1), zeroline = FALSE),
             yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE,
                          categoryorder = "array",
                          categoryarray = rev(sel$term_label)))

    # Right panel — bars go RIGHT from 0
    p_right <- plot_ly(sel, x = ~inc_right, y = ~term_label,
                       type = "bar", orientation = "h",
                       marker = list(color = col_right),
                       text  = ~paste0(round(inc_right * 100, 1), "%"),
                       textposition = "outside",
                       hovertemplate = paste0(arm_right, ": %{text}<extra></extra>"),
                       showlegend = FALSE) %>%
      layout(xaxis = list(
               title = arm_right,
               range = c(0, max_inc_val * 1.35),
               tickformat = ".0%",
               showgrid = FALSE),
             yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE,
                          categoryorder = "array",
                          categoryarray = rev(sel$term_label)))

    subplot(p_left, p_center, p_right,
            nrows = 1, shareY = TRUE,
            widths = c(0.38, 0.24, 0.38),
            titleX = TRUE) %>%
      layout(
        title = list(
          text = paste0("PT Incidence | ", input$threshold, "% threshold", vlbl,
                        if (input$bf_ge3_only) " | Grade \u22653 only" else " | All grades"),
          font = list(color = "#003865", size = 14)
        ),
        margin = list(l = 10, r = 10, t = 50, b = 40),
        plot_bgcolor  = "#fafbfc",
        paper_bgcolor = "white",
        bargap = 0.3
      )
  })
  output$dl_butterfly <- downloadHandler(
    filename = function() paste0("butterfly_", input$threshold, "pct_", input$variant, ".png"),
    content  = function(f) ggsave(f, make_butterfly_plot(), width = 12, height = 9, dpi = 180)
  )
  
  # ── Volcano ───────────────────────────────────────────────────────────────
  make_volcano_plot <- function() {
    req(nrow(ae_ctcae) > 0)
    cutoff  <- input$threshold / 100
    variant <- input$variant
    arm_a   <- input$arm_placebo
    arm_b   <- input$arm_active
    p_line  <- input$vol_p_line
    
    termset <- termset_filtered(ae_ctcae, cutoff, variant)
    denom   <- denom_from_ae(ae_ctcae)
    N_a     <- denom %>% filter(trt == arm_a) %>% pull(n_trt)
    N_b     <- denom %>% filter(trt == arm_b) %>% pull(n_trt)
    if (length(N_a) != 1 || length(N_b) != 1) return(NULL)
    
    ae_use    <- ae_ctcae %>% filter(trt %in% c(arm_a, arm_b), aedecod %in% termset)
    ge3_terms <- ae_ctcae %>% filter(ctcae_ge3) %>% distinct(aedecod) %>%
      pull(aedecod) %>% as.character()
    
    sev <- ae_use %>%
      group_by(trt, usubjid, aedecod) %>%
      summarise(max_grade = max(grade, na.rm = TRUE), .groups = "drop") %>%
      group_by(trt, aedecod) %>%
      summarise(sev_sum = sum(max_grade, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = trt, values_from = sev_sum, values_fill = 0)
    if (!arm_a %in% names(sev)) sev[[arm_a]] <- 0
    if (!arm_b %in% names(sev)) sev[[arm_b]] <- 0
    sev <- sev %>%
      transmute(aedecod, delta_burden = .data[[arm_b]] / N_b - .data[[arm_a]] / N_a)
    
    inc <- ae_use %>%
      distinct(trt, usubjid, aedecod) %>%
      count(trt, aedecod, name = "n") %>%
      pivot_wider(names_from = trt, values_from = n, values_fill = 0)
    if (!arm_a %in% names(inc)) inc[[arm_a]] <- 0
    if (!arm_b %in% names(inc)) inc[[arm_b]] <- 0
    inc <- inc %>% transmute(aedecod, n_a = .data[[arm_a]], n_b = .data[[arm_b]])
    
    out <- left_join(sev, inc, by = "aedecod")
    out$p_val <- mapply(function(nb, na_) {
      fisher.test(matrix(c(nb, N_b - nb, na_, N_a - na_),
                         nrow = 2, byrow = TRUE))$p.value
    }, out$n_b, out$n_a)
    
    p_adj <- input$vol_p_adjust
    out <- out %>% mutate(
      p_plot    = if (p_adj == "none") p_val else p.adjust(p_val, method = p_adj),
      p_plot    = pmax(p_plot, 1e-300),
      neglog10p = -log10(p_plot),
      sig       = p_plot < p_line,
      ge3_flag  = factor(
        ifelse(as.character(aedecod) %in% ge3_terms, "\u22653 PT", "Grade <3 only"),
        levels = c("Grade <3 only", "Grade \u22653 PT")
      ),
      label_pretty = str_to_title(as.character(aedecod)) %>% str_wrap(18)
    )
    
    lab_pool <- if (input$vol_label_sig_only && any(out$sig)) {
      out %>% filter(sig)
    } else {
      out
    }
    
    lab_df <- lab_pool %>%
      mutate(side  = ifelse(delta_burden >= 0, "right", "left"),
             score = abs(delta_burden) * neglog10p) %>%
      group_by(side) %>%
      slice_max(score, n = input$vol_label_n, with_ties = FALSE) %>%
      ungroup()
    
    xmax    <- max(abs(out$delta_burden), na.rm = TRUE)
    if (!is.finite(xmax) || xmax == 0) xmax <- 1e-4
    x_lim   <- c(-xmax * 1.08, xmax * 1.08)
    y_top   <- max(out$neglog10p, na.rm = TRUE)
    if (!is.finite(y_top)) y_top <- 3
    y_annot <- y_top * 1.10
    vlbl    <- if (variant == "safeguard") " | \u22653 safeguard" else ""
    
    out <- out %>% mutate(pt_alpha = ifelse(sig, 0.95, 0.15))
    
    ggplot(out, aes(x = delta_burden, y = neglog10p)) +
      geom_vline(xintercept = 0,              linetype = 2,
                 linewidth = 0.4, colour = "grey55") +
      geom_hline(yintercept = -log10(p_line), linetype = 2,
                 linewidth = 0.4, colour = "grey55") +
      geom_point(aes(colour = ge3_flag, alpha = pt_alpha), size = 2.3) +
      scale_alpha_identity(guide = "none") +
      geom_text_repel(data = lab_df,
                      aes(label = label_pretty, colour = ge3_flag),
                      size = 3, min.segment.length = 0,
                      box.padding = 0.4, max.overlaps = Inf) +
      scale_colour_manual(
        values = c("Grade <3 only" = COL_NONG3, "\u22653 PT" = COL_G3),
        labels = c("Grade <3 only" = "All grades < 3",
                   "\u22653 PT"    = "Has \u2265 Grade 3 event")
      ) +
      coord_cartesian(xlim = x_lim, ylim = c(0, y_annot)) +
      annotate("text", x = x_lim[1] + 0.02 * diff(x_lim), y = y_annot,
               hjust = 0, vjust = 0,
               label = paste0("\u2190 More on ", arm_a),
               size = 3.5, colour = "grey40") +
      annotate("text", x = x_lim[2] - 0.02 * diff(x_lim), y = y_annot,
               hjust = 1, vjust = 0,
               label = paste0("More on ", arm_b, " \u2192"),
               size = 3.5, colour = "grey40") +
      labs(
        title    = paste0("Severity-Weighted PT Signals | ",
                          input$threshold, "% threshold", vlbl),
        subtitle = paste(arm_a, "vs", arm_b,
                         "| x = severity burden diff (sum max-grade / N)"),
        x      = paste0("Severity burden diff (", arm_b, " \u2212 ", arm_a, ")"),
        y      = "-log10 p-value",
        colour = NULL
      ) +
      theme_thesis(12) +
      theme(legend.position = "top")
  }
  
  output$volcano_plot <- renderPlotly({
    p <- make_volcano_plot()
    # get lab_df from the plot environment — recompute for annotations
    req(nrow(ae_ctcae) > 0)
    cutoff  <- input$threshold / 100
    variant <- input$variant
    arm_a   <- input$arm_placebo
    arm_b   <- input$arm_active
    p_line  <- input$vol_p_line

    termset <- termset_filtered(ae_ctcae, cutoff, variant)
    denom   <- denom_from_ae(ae_ctcae)
    N_a <- denom %>% filter(trt == arm_a) %>% pull(n_trt)
    N_b <- denom %>% filter(trt == arm_b) %>% pull(n_trt)
    if (length(N_a) != 1 || length(N_b) != 1) return(ggplotly(p))

    ae_use    <- ae_ctcae %>% filter(trt %in% c(arm_a, arm_b), aedecod %in% termset)
    ge3_terms <- ae_ctcae %>% filter(ctcae_ge3) %>% distinct(aedecod) %>% pull(aedecod) %>% as.character()

    sev <- ae_use %>%
      group_by(trt, usubjid, aedecod) %>%
      summarise(max_grade = max(grade, na.rm = TRUE), .groups = "drop") %>%
      group_by(trt, aedecod) %>%
      summarise(sev_sum = sum(max_grade, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = trt, values_from = sev_sum, values_fill = 0)
    if (!arm_a %in% names(sev)) sev[[arm_a]] <- 0
    if (!arm_b %in% names(sev)) sev[[arm_b]] <- 0
    sev <- sev %>% transmute(aedecod, delta_burden = .data[[arm_b]] / N_b - .data[[arm_a]] / N_a)

    inc <- ae_use %>%
      distinct(trt, usubjid, aedecod) %>%
      count(trt, aedecod, name = "n") %>%
      pivot_wider(names_from = trt, values_from = n, values_fill = 0)
    if (!arm_a %in% names(inc)) inc[[arm_a]] <- 0
    if (!arm_b %in% names(inc)) inc[[arm_b]] <- 0
    inc <- inc %>% transmute(aedecod, n_a = .data[[arm_a]], n_b = .data[[arm_b]])

    out <- left_join(sev, inc, by = "aedecod")
    out$p_val <- mapply(function(nb, na_) {
      fisher.test(matrix(c(nb, N_b - nb, na_, N_a - na_), nrow = 2, byrow = TRUE))$p.value
    }, out$n_b, out$n_a)

    p_adj <- input$vol_p_adjust
    out <- out %>% mutate(
      p_plot    = if (p_adj == "none") p_val else p.adjust(p_val, method = p_adj),
      p_plot    = pmax(p_plot, 1e-300),
      neglog10p = -log10(p_plot),
      sig       = p_plot < p_line,
      ge3_flag  = as.character(aedecod) %in% ge3_terms,
      label_pretty = str_to_title(as.character(aedecod)),
      score     = abs(delta_burden) * neglog10p
    )

    # Top points to always label: top 3 by score per side among sig or overall top scorers
    lab_always <- out %>%
      mutate(side = ifelse(delta_burden >= 0, "right", "left")) %>%
      group_by(side) %>%
      slice_max(score, n = 3, with_ties = FALSE) %>%
      ungroup()

    plt <- ggplotly(p, tooltip = c("x", "y", "text")) %>%
      layout(legend = list(orientation = "h", y = -0.15))

    # Add permanent annotations for top points
    for (i in seq_len(nrow(lab_always))) {
      row <- lab_always[i, ]
      plt <- plt %>% add_annotations(
        x = row$delta_burden,
        y = row$neglog10p,
        text = row$label_pretty,
        showarrow = TRUE,
        arrowhead = 2,
        arrowsize = 0.5,
        arrowcolor = "grey60",
        ax = ifelse(row$delta_burden >= 0, 35, -35),
        ay = -25,
        font = list(size = 10,
                    color = if (row$ge3_flag) COL_G3 else COL_NONG3)
      )
    }
    plt
  })
  output$dl_volcano <- downloadHandler(
    filename = function() paste0("volcano_", input$threshold, "pct_", input$variant, ".png"),
    content  = function(f) {
      p <- ggplotly(make_volcano_plot())
      plotly::save_image(p, f, width = 1100, height = 800)
    }
  )
  
  # ── Temporal ──────────────────────────────────────────────────────────────
  temporal_data <- reactive({
    req(nrow(ae_ctcae) > 0)
    cutoff  <- input$threshold / 100
    variant <- input$variant
    arm_a   <- input$arm_placebo
    arm_b   <- input$arm_active
    
    termset <- termset_filtered(ae_ctcae, cutoff, variant)
    ae_use  <- ae_ctcae %>% filter(aedecod %in% termset)
    
    first <- ae_use %>%
      filter(trt %in% c(arm_a, arm_b), !is.na(onset_day)) %>%
      distinct(trt, usubjid, aedecod, onset_day) %>%
      group_by(trt, usubjid, aedecod) %>%
      summarise(first_day = min(onset_day, na.rm = TRUE), .groups = "drop")
    
    denom2 <- denom_from_ae(ae_use %>% filter(trt %in% c(arm_a, arm_b))) %>%
      distinct(trt, .keep_all = TRUE)
    
    cum <- first %>%
      count(trt, aedecod, first_day, name = "n_day") %>%
      left_join(denom2, by = "trt") %>%
      arrange(trt, aedecod, first_day) %>%
      group_by(trt, aedecod) %>%
      mutate(cum_n  = cumsum(n_day),
             cuminc = ifelse(n_trt > 0, cum_n / n_trt, 0)) %>%
      ungroup()
    
    baseline <- cum %>%
      distinct(trt, aedecod) %>%
      mutate(first_day = 0, n_day = 0L, cum_n = 0L, cuminc = 0)
    
    cum_full <- bind_rows(baseline, cum) %>%
      arrange(trt, aedecod, first_day)
    
    end_tbl <- cum_full %>%
      group_by(trt, aedecod) %>%
      summarise(end_cuminc = max(cuminc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = trt, values_from = end_cuminc, values_fill = 0)
    if (!arm_a %in% names(end_tbl)) end_tbl[[arm_a]] <- 0
    if (!arm_b %in% names(end_tbl)) end_tbl[[arm_b]] <- 0
    end_tbl <- end_tbl %>%
      mutate(abs_end_diff = abs(.data[[arm_b]] - .data[[arm_a]])) %>%
      arrange(desc(abs_end_diff))
    
    n_pick   <- min(input$temp_top_n, nrow(end_tbl))
    picked   <- end_tbl %>% slice_head(n = n_pick)
    cum_pick <- cum_full %>% filter(aedecod %in% picked$aedecod)
    list(cum_pick = cum_pick, picked = picked)
  })
  
  output$temp_pt_table <- renderDT({
    d <- temporal_data()$picked %>%
      select(aedecod, abs_end_diff) %>%
      mutate(aedecod      = str_to_title(as.character(aedecod)),
             abs_end_diff = round(abs_end_diff * 100, 1)) %>%
      rename("PT" = aedecod, "|\u0394| (%)" = abs_end_diff)
    datatable(d, options = list(dom = "tp", pageLength = 9), rownames = FALSE)
  })
  
  make_temporal_plot <- function() {
    d     <- temporal_data()$cum_pick
    arm_a <- input$arm_placebo
    arm_b <- input$arm_active
    if (nrow(d) == 0) return(ggplot() +
                               annotate("text", x = 0.5, y = 0.5,
                                        label = "No temporal data available at this threshold.",
                                        size = 5, colour = "grey50") + theme_void())
    
    trt_levels <- unique(d$trt)
    col_vals   <- setNames(rep(COL_ACTIVE, length(trt_levels)), trt_levels)
    if (arm_a %in% names(col_vals)) col_vals[[arm_a]] <- COL_PLACEBO
    
    d <- d %>% mutate(
      aedecod_pretty = str_wrap(str_to_title(as.character(aedecod)), 80)
    )
    
    n_pts    <- n_distinct(d$aedecod_pretty)
    ncol_use <- min(input$temp_ncol, n_pts)
    
    y_top   <- max(d$cuminc, na.rm = TRUE) * 1.18
    x_range <- diff(range(d$first_day, na.rm = TRUE))
    if (x_range == 0) x_range <- 1
    
    end_lab <- d %>%
      group_by(trt, aedecod_pretty) %>%
      filter(first_day == max(first_day)) %>%
      slice_tail(n = 1) %>% ungroup() %>%
      mutate(lab   = percent(cuminc, accuracy = 1),
             x_lab = first_day + max(x_range * 0.04, 5),
             y_lab = pmin(cuminc + y_top * 0.0, y_top * 0.0))
    
    vlbl <- if (input$variant == "safeguard") " | \u22653 safeguard" else ""
    
    ggplot(d, aes(first_day, cuminc, colour = trt)) +
      geom_step(linewidth = 0.8) +
      facet_wrap(~ aedecod_pretty, ncol = ncol_use) +
      scale_colour_manual(values = col_vals) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.8))) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         limits = c(0, y_top)) +
      geom_text_repel(data = end_lab,
                      aes(x = x_lab, y = y_lab, label = lab, colour = trt),
                      inherit.aes = FALSE, direction = "y", hjust = 0,
                      size = 3.2, segment.colour = NA, show.legend = FALSE, lineheight = 0.2) +
      labs(
        title    = paste0("Temporal Cumulative Incidence | ",
                          input$threshold, "% threshold", vlbl),
        subtitle = paste(arm_a, "vs", arm_b, "\u2014 top",
                         n_pts, "PTs by end-study |\u0394|"),
        x = "Study day (first onset)",
        y = "Cumulative incidence (% subjects)",
        colour = NULL
      ) +
      theme_thesis(12) +
      theme(strip.text = element_text(face = "bold", size = 9), legend.position = "bottom")
  }
  
  output$temporal_plot <- renderPlotly({
    p <- make_temporal_plot()
    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(legend = list(orientation = "h", y = -0.15))
  })
  output$dl_temporal <- downloadHandler(
    filename = function() paste0("temporal_", input$threshold, "pct_", input$variant, ".png"),
    content  = function(f) {
      p <- ggplotly(make_temporal_plot())
      plotly::save_image(p, f, width = 1200, height = 900)
    }
  )
  
  # ── Heatmap ───────────────────────────────────────────────────────────────
  make_heatmap_plot <- function() {
    req(nrow(ae_ctcae) > 0)
    cutoff  <- input$threshold / 100
    variant <- input$variant
    arm     <- input$hm_arm
    top_n   <- input$hm_top_n
    
    termset <- termset_filtered(ae_ctcae, cutoff, variant)
    
    mat_tbl <- ae_ctcae %>%
      filter(trt == arm, aedecod %in% termset) %>%
      distinct(usubjid, aedecod) %>%
      mutate(value = 1L) %>%
      pivot_wider(names_from  = aedecod,
                  values_from = value,
                  values_fill = 0L)
    
    if (ncol(mat_tbl) < 3) return(ggplot() +
                                    annotate("text", x = 0.5, y = 0.5,
                                             label = "Not enough PTs retained at this threshold.",
                                             size = 5, colour = "grey50") + theme_void())
    
    X <- mat_tbl %>% select(-usubjid) %>% as.matrix()
    storage.mode(X) <- "integer"
    X <- (X > 0) * 1L
    
    C  <- t(X) %*% X
    m  <- diag(C)
    U  <- outer(m, m, "+") - C
    Jm <- C / ifelse(U == 0, NA_real_, U)
    Jm[is.na(Jm)] <- 0
    colnames(Jm) <- rownames(Jm) <- colnames(X)
    
    top_n_use <- min(top_n, ncol(X))
    ord       <- head(sort(colnames(X)), top_n_use)
    J_sub     <- Jm[ord, ord]
    
    jac_long <- as.data.frame(as.table(J_sub), stringsAsFactors = FALSE) %>%
      setNames(c("pt1", "pt2", "jaccard")) %>%
      filter(pt1 != pt2) %>%
      mutate(
        pt1          = factor(pt1, levels = ord),
        pt2          = factor(pt2, levels = rev(ord)),
        jaccard_plot = pmin(jaccard, input$hm_cap)
      )
    
    pretty_lbl <- function(x) x %>%
      str_to_title() %>%
      str_replace_all("^Application Site ", "App site ") %>%
      str_replace_all("^Electrocardiogram ", "ECG ")
    xlabs <- pretty_lbl(levels(jac_long$pt1)); names(xlabs) <- levels(jac_long$pt1)
    ylabs <- pretty_lbl(levels(jac_long$pt2)); names(ylabs) <- levels(jac_long$pt2)
    
    vlbl <- if (variant == "safeguard") " | \u22653 safeguard" else ""
    
    ggplot(jac_long, aes(pt1, pt2, fill = jaccard_plot)) +
      geom_tile(colour = "white", linewidth = 0.15) +
      scale_fill_viridis_c(name   = "Jaccard\nsimilarity",
                           option = "magma",
                           limits = c(0, input$hm_cap),
                           breaks = seq(0, input$hm_cap, 0.2)) +
      scale_x_discrete(labels = xlabs) +
      scale_y_discrete(labels = ylabs) +
      coord_equal() +
      labs(
        title    = paste0("PT Co-occurrence Heatmap \u2014 ", arm),
        subtitle = paste0(input$threshold, "% threshold", vlbl,
                          " | top ", top_n_use, " PTs | colour capped at ", input$hm_cap),
        x = NULL, y = NULL
      ) +
      theme_thesis(11) +
      theme(
        panel.grid      = element_blank(),
        axis.text.x     = element_text(angle = 55, hjust = 1, vjust = 1, size = 8),
        axis.text.y     = element_text(size = 8),
        legend.position = "right",
        plot.margin     = margin(10, 10, 26, 10)
      )
  }
  
  output$heatmap_plot <- renderPlotly({
    p <- make_heatmap_plot()
    ggplotly(p, tooltip = c("x", "y", "fill")) %>%
      layout(yaxis = list(tickfont = list(size = 9)),
             xaxis = list(tickfont = list(size = 9)))
  })
  output$dl_heatmap <- downloadHandler(
    filename = function() paste0("heatmap_", input$hm_arm, "_", input$threshold, "pct.png"),
    content  = function(f) {
      p <- ggplotly(make_heatmap_plot())
      plotly::save_image(p, f, width = 1100, height = 1000)
    }
  )
  
  # ── CA Strip ──────────────────────────────────────────────────────────────
  make_castrip_plot <- function() {
    req(nrow(ae_ctcae) > 0)
    cutoff   <- input$threshold / 100
    variant  <- input$variant
    placebo  <- input$arm_placebo
    active   <- input$arm_active
    ge3_only <- input$ca_ge3_only
    
    ae_base <- if (ge3_only) ae_ctcae %>% filter(ctcae_ge3) else ae_ctcae
    termset <- apply_threshold(ae_base, cutoff, variant) %>%
      distinct(aedecod) %>% pull(aedecod) %>% as.character()
    
    ae_use <- ae_ctcae %>%
      filter(trt %in% c(placebo, active), aedecod %in% termset)
    if (ge3_only) ae_use <- ae_use %>% filter(ctcae_ge3)
    
    tab <- ae_use %>%
      distinct(trt, usubjid, aedecod) %>%
      count(aedecod, trt, name = "n_subj") %>%
      pivot_wider(names_from = trt, values_from = n_subj, values_fill = 0)
    
    if (nrow(tab) < 2) return(ggplot() +
                                annotate("text", x = 0.5, y = 0.5,
                                         label = "Not enough PTs at this threshold.",
                                         size = 5, colour = "grey50") + theme_void())
    
    M <- tab %>% tibble::column_to_rownames("aedecod") %>% as.matrix()
    M <- M[rowSums(M) > 0, colSums(M) > 0, drop = FALSE]
    if (nrow(M) < 2 || ncol(M) < 2) return(ggplot() +
                                             annotate("text", x = 0.5, y = 0.5,
                                                      label = "Not enough PTs or arms for CA.",
                                                      size = 5, colour = "grey50") + theme_void())
    
    fit     <- ca::ca(M)
    inertia <- fit$sv^2
    pct1    <- round(100 * inertia[1] / sum(inertia), 1)
    
    pt <- tibble(
      aedecod = rownames(M),
      ca1     = as.numeric(fit$rowcoord[, 1])
    ) %>% filter(is.finite(ca1)) %>% arrange(ca1) %>%
      mutate(y0 = rep(c(-0.06, 0.06), length.out = n()))
    
    arm_df <- tibble(
      trt = colnames(M),
      ca1 = as.numeric(fit$colcoord[, 1]),
      y   = 0
    ) %>% filter(is.finite(ca1))
    
    n_side <- input$ca_n_each
    nudge  <- max(0.10, max(abs(pt$ca1), na.rm = TRUE) * 0.08)
    
    left  <- pt %>% arrange(ca1)       %>% slice_head(n = n_side)
    right <- pt %>% arrange(desc(ca1)) %>% slice_head(n = n_side)
    
    pt_wrap <- function(x, w = 22) {
      str_wrap(str_to_title(str_squish(as.character(x))), w)
    }
    
    lab <- bind_rows(
      left  %>% mutate(side = "left"),
      right %>% mutate(side = "right")
    ) %>%
      group_by(side) %>%
      mutate(
        y_rank = row_number(),
        y_rank = ifelse(side == "left", y_rank, -y_rank),
        pt_lab = make.unique(pt_wrap(aedecod), sep = " "),
        hjust  = ifelse(side == "left", 1, 0),
        x_lab  = ca1 + ifelse(side == "left", -nudge, nudge)
      ) %>% ungroup()
    
    col_vals <- setNames(c(COL_PLACEBO, COL_ACTIVE), c(placebo, active))
    x_max    <- max(abs(pt$ca1), na.rm = TRUE)
    x_pad    <- max(0.30, x_max * 0.15)
    vlbl     <- if (variant == "safeguard") " | \u22653 safeguard" else ""
    
    p <- ggplot() +
      annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
               fill = "grey96", colour = NA) +
      annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
               fill = "#e8f4fb", colour = NA) +
      geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "grey55") +
      geom_point(data = pt,  aes(x = ca1, y = y0,
                                 text = paste0(str_to_title(str_squish(as.character(aedecod))),
                                               "<br>CA1: ", round(ca1, 3))),
                 size = 1.8, alpha = 0.25, colour = "grey30") +
      geom_point(data = lab, aes(x = ca1, y = y0,
                                 text = paste0(str_to_title(str_squish(as.character(aedecod))),
                                               "<br>CA1: ", round(ca1, 3))),
                 size = 2.8, colour = "grey20")
    
    if (input$ca_show_leaders) {
      p <- p + geom_segment(
        data = lab,
        aes(x = ca1, xend = x_lab, y = y0, yend = y_rank),
        linewidth = 0.25, colour = "grey75"
      )
    }
    
    p <- p + geom_text(
      data = lab,
      aes(x = x_lab, y = y_rank, label = pt_lab, hjust = hjust),
      size = 2.8, colour = "grey15", lineheight = 0.9
    )
    
    if (input$ca_show_arms) {
      arm_vjust <- ifelse(arm_df$ca1 >= 0, -1.4, 2.2)
      p <- p +
        geom_point(data = arm_df,
                   aes(x = ca1, y = 0, colour = trt), size = 5) +
        geom_text(data  = arm_df,
                  aes(x = ca1, y = 0, label = trt, colour = trt),
                  fontface = "bold", vjust = arm_vjust,
                  size = 3.8, show.legend = FALSE) +
        scale_colour_manual(values = col_vals,
                            guide  = guide_legend(title = NULL))
    }
    
    p +
      annotate("text", x = -Inf, y = n_side + 1.8,
               label = paste0("\u2190 more on ", placebo),
               hjust = -0.02, size = 3.5, colour = "grey40") +
      annotate("text", x =  Inf, y = n_side + 1.8,
               label = paste0("more on ", active, " \u2192"),
               hjust = 1.02, size = 3.5, colour = "grey40") +
      coord_cartesian(xlim = c(-x_max - x_pad, x_max + x_pad),
                      ylim = c(-(n_side + 2.5), (n_side + 2.5))) +
      scale_y_continuous(NULL, breaks = NULL) +
      labs(
        title    = paste0("CA Strip Plot | ", input$threshold, "% threshold", vlbl,
                          if (ge3_only) " | Grade \u22653 only" else ""),
        subtitle = paste0(placebo, " vs ", active,
                          " | CA Dim 1 explains ", pct1, "% of inertia"),
        x = "CA Dimension 1", colour = NULL
      ) +
      theme_thesis(12) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        legend.position    = if (input$ca_show_arms) "bottom" else "none"
      )
  }
  
  output$castrip_plot <- renderPlotly({
    p <- make_castrip_plot()
    placebo <- input$arm_placebo
    active  <- input$arm_active
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        showlegend = input$ca_show_arms,
        legend = list(orientation = "h", y = -0.15),
        shapes = list(
          list(type = "rect", layer = "below",
               x0 = -20, x1 = 0, y0 = -100, y1 = 100,
               fillcolor = "rgba(200,210,220,0.30)",
               line = list(width = 0)),
          list(type = "rect", layer = "below",
               x0 = 0, x1 = 20, y0 = -100, y1 = 100,
               fillcolor = "rgba(0,114,178,0.10)",
               line = list(width = 0))
        ),
        annotations = list(
          list(x = 0.01, y = 0.97, xref = "paper", yref = "paper",
               text = paste0("\u2190 more on ", placebo),
               showarrow = FALSE,
               font = list(size = 11, color = "grey50"),
               xanchor = "left"),
          list(x = 0.99, y = 0.97, xref = "paper", yref = "paper",
               text = paste0("more on ", active, " \u2192"),
               showarrow = FALSE,
               font = list(size = 11, color = "#0072B2"),
               xanchor = "right")
        )
      )
  })
  output$dl_castrip <- downloadHandler(
    filename = function() paste0("castrip_", input$threshold, "pct_", input$variant, ".png"),
    content  = function(f) {
      p <- ggplotly(make_castrip_plot())
      plotly::save_image(p, f, width = 1000, height = 700)
    }
  )
  
  # ── Info-Loss Metrics ──────────────────────────────────────────────────────
  metrics_data <- reactive({
    req(nrow(ae_ctcae) > 0)
    THRESHOLDS <- c(0.01, 0.02, 0.05, 0.10)
    LABELS     <- c("1%", "2%", "5%", "10%")
    
    ae_1pct   <- apply_threshold_any_arm(ae_ctcae, 0.01)
    base_ev   <- build_distributions(ae_1pct, "events")
    base_subj <- build_distributions(ae_1pct, "subjects")
    
    compute_one <- function(ae_tbl, label, variant_lbl) {
      dist_ev   <- build_distributions(ae_tbl, "events")
      dist_subj <- build_distributions(ae_tbl, "subjects")
      
      m_ev   <- dist_ev   %>% group_by(trt) %>%
        summarise(H_event = entropy(p), .groups = "drop")
      m_subj <- dist_subj %>% group_by(trt) %>%
        summarise(H_subj  = entropy(p), .groups = "drop")
      
      jsd_ev <- dist_ev %>%
        rename(p_curr = p) %>% select(trt, aedecod, p_curr) %>%
        full_join(base_ev %>% rename(p_base = p), by = c("trt", "aedecod")) %>%
        mutate(p_curr = replace_na(p_curr, 0), p_base = replace_na(p_base, 0)) %>%
        group_by(trt) %>%
        summarise(JSD_event = jsd(p_curr, p_base), .groups = "drop")
      
      jsd_subj <- dist_subj %>%
        rename(p_curr = p) %>% select(trt, aedecod, p_curr) %>%
        full_join(base_subj %>% rename(p_base = p), by = c("trt", "aedecod")) %>%
        mutate(p_curr = replace_na(p_curr, 0), p_base = replace_na(p_base, 0)) %>%
        group_by(trt) %>%
        summarise(JSD_subj = jsd(p_curr, p_base), .groups = "drop")
      
      ge3 <- ae_tbl %>%
        filter(ctcae_ge3) %>% distinct(trt, usubjid) %>%
        count(trt, name = "n_ge3") %>%
        right_join(denom_all_treated, by = "trt") %>%
        mutate(n_ge3          = replace_na(n_ge3, 0L),
               ge3_uncond_pct = round(100 * n_ge3 / n_all_treated, 1))
      
      m_ev %>%
        left_join(m_subj,   "trt") %>%
        left_join(jsd_ev,   "trt") %>%
        left_join(jsd_subj, "trt") %>%
        left_join(ge3,      "trt") %>%
        mutate(threshold      = label,
               variant        = variant_lbl,
               cover_events   = nrow(ae_tbl),
               cover_subjects = n_distinct(ae_tbl$usubjid),
               cover_pts      = n_distinct(ae_tbl$aedecod))
    }
    
    all_m <- map_dfr(seq_along(THRESHOLDS), function(i) {
      lab <- LABELS[i]; th <- THRESHOLDS[i]
      bind_rows(
        compute_one(apply_threshold_any_arm(ae_ctcae, th),        lab, "pure"),
        compute_one(apply_threshold_with_safeguard(ae_ctcae, th), lab, "safeguard")
      )
    })
    
    all_m %>%
      mutate(threshold = factor(threshold, levels = LABELS)) %>%
      group_by(trt, variant) %>%
      mutate(
        ref_ev = cover_events[threshold   == "1%"][1],
        ref_su = cover_subjects[threshold == "1%"][1],
        ref_pt = cover_pts[threshold      == "1%"][1],
        cover_events   = cover_events   / ref_ev,
        cover_subjects = cover_subjects / ref_su,
        cover_pts      = cover_pts      / ref_pt
      ) %>% ungroup()
  })
  
  metric_plot_base <- function(df, y_var, y_label) {
    ggplot(df, aes(x = threshold, y = .data[[y_var]],
                   group    = interaction(trt, variant),
                   colour   = trt,
                   linetype = variant)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 2.5) +
      facet_wrap(~ trt) +
      scale_linetype_manual(values = c("pure" = "solid", "safeguard" = "dashed")) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = "none") +
      labs(x = "Threshold", y = y_label, linetype = "Variant") +
      theme_thesis(11) +
      theme(legend.position = "bottom",
            legend.margin   = margin(0, 0, 0, 0))
  }
  
  output$metrics_coverage_plot <- renderPlotly({
    p <- metric_plot_base(metrics_data(), "cover_pts", "Relative PT coverage") +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      labs(title = "PT Coverage") +
      theme(strip.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5))
    ggplotly(p, tooltip = c("x", "y", "linetype")) %>%
      layout(legend = list(orientation = "h", y = -0.4),
             margin = list(t = 40))
  })
  output$metrics_entropy_plot <- renderPlotly({
    p <- metric_plot_base(metrics_data(), "H_event", "H (nats)") +
      labs(title = "Shannon Entropy H") +
      theme(strip.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5))
    ggplotly(p, tooltip = c("x", "y", "linetype")) %>%
      layout(legend = list(orientation = "h", y = -0.4),
             margin = list(t = 40))
  })
  output$metrics_jsd_plot <- renderPlotly({
    p <- metric_plot_base(metrics_data(), "JSD_event", "JSD") +
      labs(title = "JSD vs 1% baseline") +
      theme(strip.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5))
    ggplotly(p, tooltip = c("x", "y", "linetype")) %>%
      layout(legend = list(orientation = "h", y = -0.4),
             margin = list(t = 40))
  })
  
  output$metrics_table <- renderDT({
    df <- metrics_data() %>%
      select(trt, threshold, variant,
             cover_events, cover_subjects, cover_pts,
             H_event, H_subj, JSD_event, JSD_subj, ge3_uncond_pct) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
    datatable(
      df,
      options  = list(pageLength = 12, scrollX = TRUE,
                      columnDefs = list(
                        list(className = "dt-center", targets = "_all")
                      )),
      rownames = FALSE,
      colnames = c("Arm", "Threshold", "Variant",
                   "Cover Events", "Cover Subj", "Cover PTs",
                   "H (event)", "H (subj)", "JSD (event)", "JSD (subj)",
                   "\u22653 % (uncond)")
    ) %>%
      formatStyle("variant",
                  backgroundColor = styleEqual(
                    c("pure", "safeguard"),
                    c("transparent", "#e8f4e8")
                  ))
  })
  
} # end server

# ============================================================
# RUN
# ============================================================
shinyApp(ui = ui, server = server)