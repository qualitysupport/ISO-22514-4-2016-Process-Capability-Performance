########################################################################
# ISO 22514-4:2016 — Process Capability & Performance Explorer
#
# Implements the procedures described in ISO 22514-4:2016(E):
#   Statistical methods in process management — Capability and
#   performance — Part 4: Process capability estimates and performance
#   measures.
#
# Covers:
#   Clause 3   - Basic concepts: location, dispersion (inherent vs. total),
#                reference limits, and reference interval
#   Clause 4   - Capability: Cp / CPU / CPL / Cpk (normal distribution),
#                Process Capability Fraction (PCF), Cpm, Qk, MSE (4.7),
#                and proportion out-of-specification (4.8)
#   Clause 5   - Performance: Pp / PPU / PPL / Ppk (normal distribution)
#   Clause 6   - Reporting of capability and performance indices
#   Annex A    - Estimating standard deviations (subgroup-range and
#                moving-range methods)
#   Annex C    - Distribution identification for non-normal data
#                (probability-paper method; log-normal and Weibull fits)
#   Annex D    - Confidence intervals (normal-approximation / formula method)
#
# Required packages: shiny, bslib (>= 0.5.0), ggplot2, DT, readxl, MASS
# (ships with base R), gridExtra
#
# Author: Dan Lay Jr. | Calibration Support LLC
# www.calibrationsupport.com | linkedin.com/in/dlayjr
#
# Disclaimer: This tool implements calculation methods described in
# ISO 22514-4:2016 but is not affiliated with, endorsed by, or reviewed
# by ISO. It does not reproduce the standard itself; consult the official
# ISO 22514-4:2016 document (available for purchase from ISO or your
# national standards body) for authoritative guidance. Provided as-is,
# for educational and reference use.
########################################################################

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(readxl)

# ---------------------------------------------------------------------------
# Version checks — this app relies on bslib >= 0.5.0 for the persistent
# sidebar() / page_navbar(sidebar = ...) layout. Older bslib versions don't
# have sidebar() at all, which surfaces as a confusing
# "could not find function 'sidebar'" error instead of a clear one.
# ---------------------------------------------------------------------------
if (utils::packageVersion("bslib") < "0.5.0") {
  stop(
    "This app needs bslib version 0.5.0 or later (you have ", utils::packageVersion("bslib"),
    " installed) for the sidebar()/page_navbar() layout.\n",
    "Please update it and restart R:\n",
    "  install.packages(\"bslib\")\n",
    "Then run the app again.",
    call. = FALSE
  )
}
if (utils::packageVersion("shiny") < "1.7.0") {
  stop(
    "This app needs shiny version 1.7.0 or later (you have ", utils::packageVersion("shiny"),
    " installed).\n",
    "Please update it and restart R:\n",
    "  install.packages(\"shiny\")\n",
    "Then run the app again.",
    call. = FALSE
  )
}


# ============================================================================
# ISO 22514-4:2016 calculation engine
# Implements the formulas from Clauses 3-5 and Annexes A, C, D of the standard
# ============================================================================

## ---- Table A.1 constants (subgroup-size factors) ----
d2_table <- c(`2` = 1.128, `3` = 1.693, `4` = 2.059, `5` = 2.326, `6` = 2.534,
              `7` = 2.704, `8` = 2.847, `9` = 2.970, `10` = 3.078)
c4_table <- c(`2` = 0.7979, `3` = 0.8862, `4` = 0.9213, `5` = 0.9400, `6` = 0.9515,
              `7` = 0.9594, `8` = 0.9650, `9` = 0.9693, `10` = 0.9727)

## ---- 3.2/3.3 location & dispersion helpers ----

# Moment-based skewness (gamma1) and kurtosis (beta2), population-style
# denominator, matching the convention used in Annex B
sample_skewness <- function(x) {
  n <- length(x); m <- mean(x); s <- sqrt(sum((x - m)^2) / n)
  if (s == 0 || n < 3) return(NA_real_)
  (sum((x - m)^3) / n) / s^3
}
sample_kurtosis <- function(x) {
  n <- length(x); m <- mean(x); s <- sqrt(sum((x - m)^2) / n)
  if (s == 0 || n < 4) return(NA_real_)
  (sum((x - m)^4) / n) / s^4
}

# Annex A.2.1 - inherent (short-term) sigma from consecutive subgroup ranges
sigma_from_subgroups <- function(x, n) {
  n <- as.integer(n)
  k <- floor(length(x) / n)
  if (k < 2) {
    return(list(sigma = NA_real_, k = k, Rbar = NA_real_,
                note = "Need at least 2 complete subgroups of this size."))
  }
  x_use <- x[seq_len(k * n)]
  m <- matrix(x_use, ncol = n, byrow = TRUE)
  ranges <- apply(m, 1, function(r) max(r) - min(r))
  Rbar <- mean(ranges)
  key <- as.character(n)
  d2 <- if (key %in% names(d2_table)) unname(d2_table[key]) else NA_real_
  sigma <- if (!is.na(d2)) Rbar / d2 else NA_real_
  list(sigma = sigma, k = k, Rbar = Rbar,
       note = if (is.na(d2)) "No d2 factor available for this subgroup size (supported: 2-10)." else NULL)
}

# Individuals / moving-range estimate of short-term sigma (n = 2, d2 = 1.128)
sigma_from_moving_range <- function(x) {
  mr <- abs(diff(x))
  mrbar <- mean(mr)
  list(sigma = mrbar / 1.128, MRbar = mrbar)
}

## ---- Clause 4.4 / 5.2 index formulas (normal distribution) ----
cp_val  <- function(U, L, sigma) (U - L) / (6 * sigma)
cpu_val <- function(U, mu, sigma) (U - mu) / (3 * sigma)
cpl_val <- function(mu, L, sigma) (mu - L) / (3 * sigma)

# Convenience wrapper returning Cp/CPU/CPL/Cpk (or Pp/PPU/PPL/Ppk) together
capability_block <- function(sigma, mu, L, U) {
  Cp  <- if (!is.na(L) && !is.na(U)) cp_val(U, L, sigma) else NA_real_
  CPU <- if (!is.na(U)) cpu_val(U, mu, sigma) else NA_real_
  CPL <- if (!is.na(L)) cpl_val(mu, L, sigma) else NA_real_
  Cpk <- suppressWarnings(min(c(CPU, CPL), na.rm = TRUE))
  if (!is.finite(Cpk)) Cpk <- NA_real_
  list(Cp = Cp, CPU = CPU, CPL = CPL, Cpk = Cpk)
}

## ---- Clause 4.7 other capability measures ----
pcf_val      <- function(Cp) 100 / Cp
mse_val      <- function(sigma, mu, target) sigma^2 + (mu - target)^2
qk_val       <- function(sigma, mu, target) {
  if (is.na(target) || target == 0) return(NA_real_)
  100 * sqrt(sigma^2 + (mu - target)^2) / target
}
cpm_val      <- function(U, L, mu, sigma, target) (U - L) / (6 * sqrt(sigma^2 + (mu - target)^2))
cpm_star_val <- function(U, L, mu, sigma, target) {
  min(U - target, target - L) / (3 * sqrt(sigma^2 + (mu - target)^2))
}

## ---- Clause 4.8 proportion out-of-specification (normal distribution) ----
p_from_cpk <- function(Cpk) if (is.na(Cpk)) 0 else pnorm(-3 * Cpk)

proportion_block <- function(block) {
  pU <- if (!is.na(block$CPU)) p_from_cpk(block$CPU) else 0
  pL <- if (!is.na(block$CPL)) p_from_cpk(block$CPL) else 0
  list(pU = pU, pL = pL, ptotal = pU + pL, yield_pct = 100 * (1 - (pU + pL)))
}

## ---- Annex D.1.2 confidence intervals (formula / normal-approximation method) ----
ci_cp <- function(Cp, N, conf = 0.95) {
  if (is.na(Cp) || N < 3) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  se <- Cp * sqrt(1 / (2 * (N - 1)))
  c(lower = Cp - z * se, upper = Cp + z * se)
}
ci_cpk_one_sided <- function(Cpk, N, conf = 0.95) {
  if (is.na(Cpk) || N < 3) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  se <- sqrt(1 / (9 * N) + (Cpk^2) / (2 * (N - 1)))
  c(lower = Cpk - z * se, upper = Cpk + z * se)
}

## ---- Clause 4.5.2 / 5.3.2 probability-paper method for non-normal data ----
# Empirical percentile using a plotting-position quantile type as a stand-in
# for reading values off probability paper.
empirical_percentile <- function(x, p, type = 6) {
  as.numeric(quantile(x, probs = p, type = type, names = FALSE))
}

nonnormal_block <- function(x, L, U) {
  Y1 <- empirical_percentile(x, 0.00135)
  Y2 <- empirical_percentile(x, 0.99865)
  Ymed <- median(x)
  Pp  <- if (!is.na(L) && !is.na(U)) (U - L) / (Y2 - Y1) else NA_real_
  PU  <- if (!is.na(U)) (U - Ymed) / (Y2 - Ymed) else NA_real_
  PL  <- if (!is.na(L)) (Ymed - L) / (Ymed - Y1) else NA_real_
  Ppk <- suppressWarnings(min(c(PU, PL), na.rm = TRUE))
  if (!is.finite(Ppk)) Ppk <- NA_real_
  list(Y1 = Y1, Y2 = Y2, Ymed = Ymed, Pp = Pp, PU = PU, PL = PL, Ppk = Ppk)
}

## ---- Annex C.3 / C.5 distribution-identification method (lognormal & Weibull) ----
fit_lognormal <- function(x) {
  if (any(x <= 0)) return(NULL)
  tryCatch(MASS::fitdistr(x, "lognormal"), error = function(e) NULL)
}
fit_weibull <- function(x) {
  if (any(x <= 0)) return(NULL)
  tryCatch(MASS::fitdistr(x, "weibull"), error = function(e) NULL)
}

lognormal_indices <- function(fit, L, U) {
  if (is.null(fit)) return(NULL)
  ml <- unname(fit$estimate["meanlog"]); sl <- unname(fit$estimate["sdlog"])
  Y1 <- qlnorm(0.00135, ml, sl); Y2 <- qlnorm(0.99865, ml, sl); Ymed <- qlnorm(0.5, ml, sl)
  pL <- if (!is.na(L)) plnorm(L, ml, sl) else 0
  pU <- if (!is.na(U)) 1 - plnorm(U, ml, sl) else 0
  Pp <- if (!is.na(L) && !is.na(U)) (U - L) / (Y2 - Y1) else NA_real_
  list(meanlog = ml, sdlog = sl, Y1 = Y1, Y2 = Y2, Ymed = Ymed, Pp = Pp, pL = pL, pU = pU)
}
weibull_indices <- function(fit, L, U) {
  if (is.null(fit)) return(NULL)
  sh <- unname(fit$estimate["shape"]); sc <- unname(fit$estimate["scale"])
  Y1 <- qweibull(0.00135, sh, sc); Y2 <- qweibull(0.99865, sh, sc); Ymed <- qweibull(0.5, sh, sc)
  pL <- if (!is.na(L)) pweibull(L, sh, sc) else 0
  pU <- if (!is.na(U)) 1 - pweibull(U, sh, sc) else 0
  Pp <- if (!is.na(L) && !is.na(U)) (U - L) / (Y2 - Y1) else NA_real_
  list(shape = sh, scale = sc, Y1 = Y1, Y2 = Y2, Ymed = Ymed, Pp = Pp, pL = pL, pU = pU)
}

## ---- interpretation helper (industry rule-of-thumb, NOT part of the ISO text) ----
interpret_index <- function(v) {
  if (is.na(v)) return(list(label = "N/A", class = "secondary"))
  if (v >= 1.33) return(list(label = "Capable", class = "success"))
  if (v >= 1.00) return(list(label = "Marginal", class = "warning"))
  return(list(label = "Not capable", class = "danger"))
}


# ---------------------------------------------------------------------------
# Example dataset (used until the person supplies their own data)
# ---------------------------------------------------------------------------
make_example_data <- function() {
  set.seed(42)
  round(rnorm(200, mean = 10.05, sd = 0.30), 3)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
theme <- bs_theme(version = 5, bootswatch = "flatly",
                   primary = "#2c3e50", success = "#18823c",
                   warning = "#e08e0b", danger = "#c0392b")

data_sidebar <- sidebar(
  title = "Data & specification",
  width = 340,
  open = "open",

  radioButtons("data_mode", "Data source",
               choices = c("Use example data" = "example",
                           "Upload a file" = "upload",
                           "Paste values" = "manual"),
               selected = "example"),

  conditionalPanel(
    "input.data_mode == 'upload'",
    fileInput("file", "CSV or Excel file", accept = c(".csv", ".xlsx", ".xls")),
    uiOutput("col_select_ui")
  ),

  conditionalPanel(
    "input.data_mode == 'manual'",
    textAreaInput("manual_vals", "Paste measurements",
                   rows = 6, placeholder = "One value per line, or separated by commas/spaces")
  ),

  hr(),
  strong("Short-term (inherent) variation"),
  helpText("Used for the capability indices Cp/Cpk. Choose how the data are grouped into subgroups (ISO 22514-4, Annex A)."),
  numericInput("subgroup_n", "Subgroup size (1 = individuals / moving range)",
               value = 1, min = 1, max = 10, step = 1),

  hr(),
  strong("Specification limits"),
  checkboxInput("has_L", "Lower limit (L)", value = TRUE),
  conditionalPanel("input.has_L", numericInput("L_val", NULL, value = 9)),
  checkboxInput("has_U", "Upper limit (U)", value = TRUE),
  conditionalPanel("input.has_U", numericInput("U_val", NULL, value = 11)),
  checkboxInput("has_T", "Target (T) — for Cpm / Qk", value = TRUE),
  conditionalPanel("input.has_T", numericInput("T_val", NULL, value = 10)),

  hr(),
  sliderInput("conf_level", "Confidence level for intervals",
              min = 0.80, max = 0.99, value = 0.95, step = 0.01)
)

ui <- page_navbar(
  title = "ISO 22514-4 Capability & Performance Explorer",
  theme = theme,
  sidebar = data_sidebar,
  fillable = TRUE,

  nav_panel("Overview", icon = icon("book-open"),
    layout_column_wrap(
      width = 1,
      card(
        card_header("What this app does"),
        card_body(
          p("ISO 22514-4:2016 sets out the common measures used to describe how well a process's output stays inside its specification limits.
             This app lets you calculate those measures on your own data, and includes an interactive section that builds
             intuition for the underlying concepts."),
          p(strong("Two families of measures, one key distinction:")),
          layout_column_wrap(
            width = 1/2,
            card(class = "border-primary",
              card_header(strong("Capability (Cp, Cpk)")),
              card_body(
                p("Describes what a process is capable of when it is running in a demonstrated state of statistical control —
                   i.e. only common-cause (random) variation remains, assignable causes have been removed, and a control
                   chart confirms stability."),
                p("Uses the ", strong("short-term / inherent"), " standard deviation, typically estimated from
                   subgroup ranges on a control chart.")
              )
            ),
            card(class = "border-secondary",
              card_header(strong("Performance (Pp, Ppk)")),
              card_body(
                p("Describes what a process actually delivered over an observation period, with no requirement that the
                   process be shown to be in statistical control. It reflects everything that happened — including any
                   drift, shifts, or instability."),
                p("Uses the ", strong("total / long-term"), " standard deviation, calculated directly from all the data.")
              )
            )
          ),
          p("Because performance includes extra sources of variation that capability assumes have been eliminated,
             Pp/Ppk will typically be lower than Cp/Cpk for the same process."),
          tags$hr(),
          p(strong("How to use this app:")),
          tags$ol(
            tags$li("Load your data (or start with the built-in example) and set your specification limits in the sidebar."),
            tags$li("Visit ", strong("Concepts"), " to build intuition for dispersion, reference limits, and the capability/performance distinction."),
            tags$li("Visit ", strong("Capability & Performance"), " for the main Cp/Cpk/Pp/Ppk calculation, plots, and confidence intervals."),
            tags$li("If your data look skewed, visit ", strong("Non-Normal Data"), " for probability-paper and distribution-fitting estimates."),
            tags$li("Visit ", strong("Other Indices"), " for PCF, MSE, Qk and Cpm, and ", strong("Report"), " for a shareable summary.")
          )
        )
      )
    )
  ),

  nav_panel("Concepts", icon = icon("lightbulb"),
    navset_card_tab(
      nav_panel("Dispersion & reference limits",
        layout_sidebar(
          sidebar = sidebar(
            position = "right", width = 280,
            sliderInput("c_mu", "Process mean (μ)", min = -5, max = 5, value = 0, step = 0.1),
            sliderInput("c_sigma", "Process standard deviation (σ)", min = 0.2, max = 3, value = 1, step = 0.1),
            checkboxInput("c_show_limits", "Show reference limits (0.135% / 99.865%)", TRUE),
            checkboxInput("c_show_spec", "Overlay specification limits", FALSE),
            conditionalPanel("input.c_show_spec",
              sliderInput("c_L", "Lower spec limit", min = -8, max = 8, value = -3, step = 0.1),
              sliderInput("c_U", "Upper spec limit", min = -8, max = 8, value = 3, step = 0.1)
            )
          ),
          p("A process that is in statistical control has a predictable spread. The ", strong("reference interval"),
            " is the span that captures 99.73% of individual values, leaving 0.135% in each tail. For a normal
            distribution this interval is exactly six standard deviations wide. Move the sliders to see how the
            reference limits shift with the process mean and spread, and (optionally) how that spread compares
            with a specification window."),
          plotOutput("concept_plot1", height = 380)
        )
      ),
      nav_panel("Short-term vs. total dispersion",
        p("A process can look tighter over a short window than it does over the long run. Each little cluster below
           represents a short-term subgroup (its spread reflects only the ", em("inherent"), " variation). If the
           subgroup means drift over time — a common-cause of longer-term instability — the ", strong("total"),
          " dispersion of all the individual values pooled together is wider than any single subgroup's spread.
           This is exactly why ISO 22514-4 keeps capability (short-term σ) and performance (total σ) as separate,
           clearly labelled measures rather than one number."),
        sliderInput("c_drift", "Amount of mean drift between subgroups", min = 0, max = 3, value = 1.2, step = 0.1, width = "50%"),
        plotOutput("concept_plot2", height = 420)
      ),
      nav_panel("Capability vs. performance conditions",
        layout_column_wrap(
          width = 1/2,
          card(card_header(strong("Conditions expected for a capability study")),
            card_body(tags$ul(
              tags$li("The process is monitored with a control chart and shown to be in statistical control."),
              tags$li("Technical/environmental conditions are documented and held consistent."),
              tags$li("Measurement system uncertainty has been assessed and judged adequate."),
              tags$li("The standard deviation used reflects only short-term, inherent (common-cause) variation.")
            ))
          ),
          card(card_header(strong("Conditions expected for a performance study")),
            card_body(tags$ul(
              tags$li("No requirement for statistical control — historical or as-collected data can be used."),
              tags$li("Technical/environmental conditions are still documented for context."),
              tags$li("Measurement system uncertainty should still be assessed."),
              tags$li("The standard deviation used reflects total, long-term variation, including any instability.")
            ))
          )
        )
      )
    )
  ),

  nav_panel("Capability & Performance", icon = icon("chart-column"),
    layout_column_wrap(
      width = 1,
      uiOutput("main_valueboxes"),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Distribution vs. specification limits"),
             card_body(plotOutput("hist_plot", height = 380))),
        card(card_header("Normal probability (Q-Q) plot"),
             card_body(plotOutput("qq_plot", height = 380),
                       p(class = "text-muted small", "Points following the diagonal line support the normal-distribution
                          assumption behind Cp/Cpk and Pp/Ppk. If they curve away, see the Non-Normal Data tab.")))
      ),
      card(card_header("Detailed results"),
        card_body(
          DTOutput("results_table"),
          tags$hr(),
          uiOutput("sigma_method_note")
        )
      ),
      card(card_header("Confidence intervals (Annex D.1.2, normal-approximation method)"),
        card_body(
          p(class = "text-muted", "Calculated indices are point estimates of a true, unknown value; ISO 22514-4 recommends
             reporting a confidence interval alongside them. This formula method is intended for samples of 50 or more."),
          DTOutput("ci_table")
        )
      )
    )
  ),

  nav_panel("Non-Normal Data", icon = icon("chart-area"),
    layout_column_wrap(
      width = 1,
      card(card_header("Is the normal distribution a reasonable fit?"),
        card_body(
          layout_column_wrap(width = 1/3,
            value_box(title = "Skewness (γ1)", value = textOutput("skew_val"), theme = "secondary"),
            value_box(title = "Kurtosis (β2)", value = textOutput("kurt_val"), theme = "secondary"),
            value_box(title = "Shapiro-Wilk p-value", value = textOutput("shapiro_val"), theme = "secondary")
          ),
          p(class = "text-muted small",
            "A normal distribution has skewness 0 and kurtosis 3. A Shapiro-Wilk p-value below 0.05 is common practice
             evidence against normality (n must be between 3 and 5000 for this test).")
        )
      ),
      card(card_header("Probability-paper method (Clause 4.5.2 / 5.3.2)"),
        card_body(
          p("Instead of assuming normality, this method reads the 0.135 and 99.865 percentiles directly from the
             empirical distribution and substitutes them for the ±3σ points used in the normal formulas."),
          DTOutput("nonnormal_table"),
          p(class = "text-muted small", "Extreme percentiles estimated from small samples can be unstable —
             treat this as indicative rather than precise for n well under 100.")
        )
      ),
      card(card_header("Distribution-identification method (Annex C)"),
        card_body(
          p("As an alternative, a specific distribution family can be fitted to the data and its theoretical
             percentiles used instead. This app fits the two families most common in capability work:
             the log-normal (Annex C.3) and Weibull (Annex C.5) distributions. Both require strictly positive data."),
          DTOutput("fit_table")
        )
      )
    )
  ),

  nav_panel("Other Indices", icon = icon("calculator"),
    layout_column_wrap(
      width = 1,
      card(card_header("Target-based measures (Clause 4.7)"),
        card_body(
          p("These measures incorporate a target value, T, penalising both off-centre location and excess spread
             in a single number — useful when the goal is to sit close to a nominal value, not merely inside the limits."),
          DTOutput("other_indices_table")
        )
      ),
      card(card_header("Process Capability Fraction (PCF)"),
        card_body(
          p("The PCF is simply the inverse of Cp expressed as a percentage: the proportion of the specification
             width that the process's natural spread consumes."),
          textOutput("pcf_text")
        )
      )
    )
  ),

  nav_panel("Report", icon = icon("file-lines"),
    layout_column_wrap(
      width = 1,
      card(card_header("Summary report"),
        card_body(
          DTOutput("report_table"),
          tags$br(),
          downloadButton("download_report", "Download report (CSV)")
        )
      ),
      card(card_header("Raw data used"),
        card_body(DTOutput("raw_data_table")))
    )
  ),

  nav_spacer(),
  nav_item(tags$a(href = "https://www.iso.org/standard/72463.html", target = "_blank",
                   "ISO 22514-4:2016 on iso.org", class = "small text-muted"))
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- Data ingestion ----
  dataset_raw <- reactive({
    if (input$data_mode == "example") {
      return(data.frame(value = make_example_data()))
    } else if (input$data_mode == "upload") {
      req(input$file)
      ext <- tolower(tools::file_ext(input$file$name))
      df <- switch(ext,
        csv = read.csv(input$file$datapath, stringsAsFactors = FALSE),
        xlsx = as.data.frame(readxl::read_excel(input$file$datapath)),
        xls  = as.data.frame(readxl::read_excel(input$file$datapath)),
        validate("Unsupported file type — please upload a .csv, .xlsx or .xls file.")
      )
      df
    } else {
      req(input$manual_vals)
      vals <- suppressWarnings(as.numeric(unlist(strsplit(input$manual_vals, "[,\\s]+"))))
      vals <- vals[!is.na(vals)]
      data.frame(value = vals)
    }
  })

  output$col_select_ui <- renderUI({
    req(input$data_mode == "upload")
    df <- tryCatch(dataset_raw(), error = function(e) NULL)
    req(df)
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    validate(need(length(numeric_cols) > 0, "No numeric columns found in this file."))
    selectInput("meas_col", "Measurement column", choices = numeric_cols)
  })

  x_vals <- reactive({
    df <- dataset_raw()
    if (input$data_mode == "upload") {
      req(input$meas_col)
      v <- suppressWarnings(as.numeric(df[[input$meas_col]]))
    } else {
      v <- df$value
    }
    v <- v[is.finite(v)]
    validate(need(length(v) >= 5, "Please provide at least 5 numeric values."))
    v
  })

  ## ---- Specification limits ----
  spec <- reactive({
    L <- if (isTRUE(input$has_L)) input$L_val else NA_real_
    U <- if (isTRUE(input$has_U)) input$U_val else NA_real_
    Tg <- if (isTRUE(input$has_T)) input$T_val else NA_real_
    validate(need(!is.na(L) || !is.na(U), "Please specify at least one specification limit in the sidebar."))
    if (!is.na(L) && !is.na(U)) validate(need(U > L, "The upper limit must be greater than the lower limit."))
    list(L = L, U = U, T = Tg)
  })

  ## ---- Short-term sigma ----
  sigma_short_info <- reactive({
    x <- x_vals()
    n <- input$subgroup_n
    if (is.null(n) || n <= 1) {
      info <- sigma_from_moving_range(x)
      list(sigma = info$sigma,
           method = sprintf("Moving-range method (individuals): MR-bar = %.4f, d2 = 1.128", info$MRbar))
    } else {
      info <- sigma_from_subgroups(x, n)
      if (is.na(info$sigma)) {
        list(sigma = NA_real_, method = info$note %||% "Unable to estimate short-term sigma with this subgroup size.")
      } else {
        list(sigma = info$sigma,
             method = sprintf("Subgroup range method: %d complete subgroups of size %d, R-bar = %.4f, d2 = %.3f",
                               info$k, n, info$Rbar, unname(d2_table[as.character(n)])))
      }
    }
  })
  `%||%` <- function(a, b) if (is.null(a)) b else a

  ## ---- Core results ----
  results <- reactive({
    x <- x_vals(); sp <- spec()
    mu <- mean(x); med <- median(x)
    sigma_t <- sd(x)
    ss <- sigma_short_info()
    sigma_s <- ss$sigma

    cap  <- if (!is.na(sigma_s)) capability_block(sigma_s, mu, sp$L, sp$U) else list(Cp=NA,CPU=NA,CPL=NA,Cpk=NA)
    perf <- capability_block(sigma_t, mu, sp$L, sp$U)

    list(N = length(x), mu = mu, med = med, sigma_t = sigma_t, sigma_s = sigma_s,
         ss_method = ss$method, cap = cap, perf = perf, spec = sp)
  })

  ## ---- Value boxes ----
  output$main_valueboxes <- renderUI({
    r <- results()
    mk <- function(label, val, sub) {
      interp <- interpret_index(val)
      value_box(title = label,
                 value = if (is.na(val)) "N/A" else sprintf("%.2f", val),
                 showcase = tags$span(class = paste0("badge bg-", interp$class), interp$label),
                 theme = interp$class, p(class = "small", sub))
    }
    layout_column_wrap(
      width = 1/4,
      mk("Cp (capability)", r$cap$Cp, "Spread vs. tolerance, short-term σ"),
      mk("Cpk (capability)", r$cap$Cpk, "Spread + centring, short-term σ"),
      mk("Pp (performance)", r$perf$Cp, "Spread vs. tolerance, total σ"),
      mk("Ppk (performance)", r$perf$Cpk, "Spread + centring, total σ")
    )
  })

  ## ---- Plots ----
  output$hist_plot <- renderPlot({
    x <- x_vals(); sp <- spec(); r <- results()
    df <- data.frame(x = x)
    g <- ggplot(df, aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = min(30, max(10, round(length(x)/5))),
                      fill = "#2c3e50", alpha = 0.75, color = "white") +
      stat_function(fun = dnorm, args = list(mean = r$mu, sd = r$sigma_t),
                    color = "#c0392b", linewidth = 1) +
      labs(x = "Measured value", y = "Density") +
      theme_minimal(base_size = 13)
    if (!is.na(sp$L)) g <- g + geom_vline(xintercept = sp$L, color = "#e08e0b", linewidth = 1, linetype = "dashed")
    if (!is.na(sp$U)) g <- g + geom_vline(xintercept = sp$U, color = "#e08e0b", linewidth = 1, linetype = "dashed")
    if (!is.na(sp$T)) g <- g + geom_vline(xintercept = sp$T, color = "#18823c", linewidth = 1, linetype = "dotted")
    g
  })

  output$qq_plot <- renderPlot({
    x <- x_vals()
    df <- data.frame(x = x)
    ggplot(df, aes(sample = x)) +
      stat_qq(color = "#2c3e50", alpha = 0.7) +
      stat_qq_line(color = "#c0392b", linewidth = 1) +
      labs(x = "Theoretical normal quantiles", y = "Sample quantiles") +
      theme_minimal(base_size = 13)
  })

  ## ---- Results table ----
  output$results_table <- renderDT({
    r <- results()
    df <- data.frame(
      Measure = c("N", "Mean", "Median", "Short-term σ (capability)", "Total σ (performance)",
                  "Cp", "CPU", "CPL", "Cpk",
                  "Pp", "PPU", "PPL", "Ppk"),
      Value = c(r$N, r$mu, r$med, r$sigma_s, r$sigma_t,
                r$cap$Cp, r$cap$CPU, r$cap$CPL, r$cap$Cpk,
                r$perf$Cp, r$perf$CPU, r$perf$CPL, r$perf$Cpk)
    )
    df$Value <- round(df$Value, 4)
    datatable(df, options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  output$sigma_method_note <- renderUI({
    r <- results()
    p(class = "text-muted small", strong("Short-term σ estimation: "), r$ss_method)
  })

  ## ---- Confidence intervals ----
  output$ci_table <- renderDT({
    r <- results(); conf <- input$conf_level
    rows <- list(
      c("Cp",  r$cap$Cp,  ci_cp(r$cap$Cp, r$N, conf)),
      c("Cpk", r$cap$Cpk, ci_cpk_one_sided(r$cap$Cpk, r$N, conf)),
      c("Pp",  r$perf$Cp, ci_cp(r$perf$Cp, r$N, conf)),
      c("Ppk", r$perf$Cpk, ci_cpk_one_sided(r$perf$Cpk, r$N, conf))
    )
    df <- do.call(rbind.data.frame, rows)
    names(df) <- c("Index", "Estimate", "Lower", "Upper")
    df[ , 2:4] <- lapply(df[ , 2:4], function(col) round(as.numeric(col), 3))
    datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ---- Non-normal analysis ----
  output$skew_val <- renderText(sprintf("%.3f", sample_skewness(x_vals())))
  output$kurt_val <- renderText(sprintf("%.3f", sample_kurtosis(x_vals())))
  output$shapiro_val <- renderText({
    x <- x_vals()
    if (length(x) >= 3 && length(x) <= 5000) {
      sprintf("%.4f", shapiro.test(x)$p.value)
    } else "n out of range"
  })

  output$nonnormal_table <- renderDT({
    x <- x_vals(); sp <- spec()
    nb <- nonnormal_block(x, sp$L, sp$U)
    df <- data.frame(
      Quantity = c("Y1 (0.135th percentile)", "Y2 (99.865th percentile)", "Median",
                   "Pp (probability-paper)", "Upper-side index", "Lower-side index", "Ppk (probability-paper)"),
      Value = round(c(nb$Y1, nb$Y2, nb$Ymed, nb$Pp, nb$PU, nb$PL, nb$Ppk), 4)
    )
    datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$fit_table <- renderDT({
    x <- x_vals(); sp <- spec()
    ln <- suppressWarnings(lognormal_indices(fit_lognormal(x), sp$L, sp$U))
    wb <- suppressWarnings(weibull_indices(fit_weibull(x), sp$L, sp$U))
    if (is.null(ln) && is.null(wb)) {
      return(datatable(data.frame(Note = "Data must be strictly positive to fit log-normal or Weibull distributions."),
                        options = list(dom = "t"), rownames = FALSE))
    }
    rows <- list()
    if (!is.null(ln)) rows$Lognormal <- c("Log-normal", round(ln$Pp, 4), round(ln$pL, 5), round(ln$pU, 5))
    if (!is.null(wb)) rows$Weibull <- c("Weibull", round(wb$Pp, 4), round(wb$pL, 5), round(wb$pU, 5))
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(df) <- c("Distribution", "Pp (fitted)", "Prop. below L", "Prop. above U")
    df[ , 2:4] <- lapply(df[ , 2:4], as.numeric)
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- Other indices ----
  output$other_indices_table <- renderDT({
    r <- results(); sp <- spec()
    if (is.na(sp$T)) {
      return(datatable(data.frame(Note = "Set a target value (T) in the sidebar to see these measures."),
                        options = list(dom = "t"), rownames = FALSE))
    }
    sigma <- r$sigma_t; mu <- r$mu; Tg <- sp$T
    mse <- mse_val(sigma, mu, Tg)
    qk  <- qk_val(sigma, mu, Tg)
    cpm <- if (!is.na(sp$L) && !is.na(sp$U)) cpm_val(sp$U, sp$L, mu, sigma, Tg) else NA
    cpm_star <- if (!is.na(sp$L) && !is.na(sp$U)) cpm_star_val(sp$U, sp$L, mu, sigma, Tg) else NA
    df <- data.frame(
      Measure = c("Mean square error (MSE)", "Qk index (%)", "Cpm", "Cpm* (non-central target)"),
      Value = round(c(mse, qk, cpm, cpm_star), 4)
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  output$pcf_text <- renderText({
    r <- results()
    if (is.na(r$cap$Cp)) return("PCF requires both spec limits and a valid short-term sigma.")
    sprintf("PCF = %.1f%% of the specification width is used by the process's natural (short-term) spread.",
            pcf_val(r$cap$Cp))
  })

  ## ---- Concepts tab plots ----
  output$concept_plot1 <- renderPlot({
    mu <- input$c_mu; sigma <- input$c_sigma
    xr <- seq(mu - 5 * sigma - 2, mu + 5 * sigma + 2, length.out = 600)
    df <- data.frame(x = xr, y = dnorm(xr, mu, sigma))
    g <- ggplot(df, aes(x, y)) +
      geom_line(color = "#2c3e50", linewidth = 1) +
      labs(x = "Value", y = "Density") +
      theme_minimal(base_size = 13)
    if (isTRUE(input$c_show_limits)) {
      lo <- qnorm(0.00135, mu, sigma); hi <- qnorm(0.99865, mu, sigma)
      shade <- df[df$x >= lo & df$x <= hi, ]
      g <- g + geom_area(data = shade, aes(x, y), fill = "#2c3e50", alpha = 0.15) +
        geom_vline(xintercept = c(lo, hi), color = "#2c3e50", linetype = "dashed") +
        annotate("text", x = lo, y = max(df$y) * 1.05, label = "0.135%", size = 3.5, hjust = 1) +
        annotate("text", x = hi, y = max(df$y) * 1.05, label = "99.865%", size = 3.5, hjust = 0)
    }
    if (isTRUE(input$c_show_spec)) {
      g <- g + geom_vline(xintercept = input$c_L, color = "#e08e0b", linewidth = 1) +
        geom_vline(xintercept = input$c_U, color = "#e08e0b", linewidth = 1)
    }
    g
  })

  output$concept_plot2 <- renderPlot({
    set.seed(11)
    n_sub <- 5; n_per <- 25; sigma_short <- 0.4
    drift <- input$c_drift
    means <- seq(-drift * 2, drift * 2, length.out = n_sub)
    dat <- do.call(rbind, lapply(seq_along(means), function(i) {
      data.frame(subgroup = paste("Subgroup", i),
                 x = rnorm(n_per, means[i], sigma_short))
    }))
    dat$subgroup <- factor(dat$subgroup, levels = unique(dat$subgroup))
    g1 <- ggplot(dat, aes(x, fill = subgroup)) +
      geom_histogram(bins = 12, alpha = 0.85, position = "identity", color = "white") +
      facet_wrap(~subgroup, nrow = 1) +
      labs(title = "Short-term dispersion (within each subgroup)", x = NULL, y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    g2 <- ggplot(dat, aes(x)) +
      geom_histogram(bins = 25, fill = "#c0392b", alpha = 0.8, color = "white") +
      labs(title = "Total dispersion (all subgroups pooled)", x = "Value", y = NULL) +
      theme_minimal(base_size = 12)
    cowplot_stack(g1, g2)
  })

  # simple vertical stacking helper without extra dependency on patchwork/cowplot
  cowplot_stack <- function(p1, p2) {
    gridExtra_available <- requireNamespace("gridExtra", quietly = TRUE)
    if (gridExtra_available) {
      gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1, 1.2))
    } else {
      p2
    }
  }

  ## ---- Report ----
  output$report_table <- renderDT({
    r <- results(); sp <- spec()
    df <- data.frame(
      Item = c("Number of values (N)", "Mean", "Median",
               "Lower specification limit", "Upper specification limit", "Target",
               "Short-term sigma method",
               "Cp", "Cpk", "Pp", "Ppk",
               "Confidence level"),
      Value = c(r$N, round(r$mu, 4), round(r$med, 4),
                ifelse(is.na(sp$L), "-", sp$L), ifelse(is.na(sp$U), "-", sp$U), ifelse(is.na(sp$T), "-", sp$T),
                r$ss_method,
                round(r$cap$Cp, 3), round(r$cap$Cpk, 3), round(r$perf$Cp, 3), round(r$perf$Cpk, 3),
                paste0(input$conf_level * 100, "%"))
    )
    datatable(df, options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  output$download_report <- downloadHandler(
    filename = function() paste0("iso22514-4_report_", Sys.Date(), ".csv"),
    content = function(file) {
      r <- results(); sp <- spec()
      df <- data.frame(
        Item = c("N", "Mean", "Median", "L", "U", "Target", "Short_term_sigma_method",
                 "Cp", "Cpk", "Pp", "Ppk", "Confidence_level"),
        Value = c(r$N, r$mu, r$med, sp$L, sp$U, sp$T, r$ss_method,
                  r$cap$Cp, r$cap$Cpk, r$perf$Cp, r$perf$Cpk, input$conf_level)
      )
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$raw_data_table <- renderDT({
    datatable(data.frame(Value = x_vals()), options = list(pageLength = 10))
  })
}

shinyApp(ui, server)
