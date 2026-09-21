################################################################################
## 6b_regional_state_model.R
##
## Gulf Sturgeon CJFAS revision
## Fit the separate four-region multistate model used for Table A5.
## Production workflow: preserves source rows/frequency coding without aggregating histories.
##
## Purpose:
##   1. Start from the SAME adult annual encounter histories used by Script 2.
##   2. Remove fish first entering the adult analysis in 2022 (no subsequent
##      annual survival interval).
##   3. Collapse the seven river states to four regional states:
##
##        W = West             = Pearl + Pascagoula
##        P = Pensacola Bay    = Escambia + Yellow
##        C = Choctawhatchee
##        E = East             = Apalachicola + Suwannee
##
##   4. Fit a separate four-state multistate model:
##
##        S(region:time) p(region:time) Psi(region-to-region)
##
##      Apparent survival and detection vary by region and year.
##      Transition probabilities are constant through time and directional
##      among regions. Direct West <-> East transitions are fixed to zero,
##      consistent with the river-scale transition constraints.
##
##   5. Save the regional model object used by Script 7 and print
##      regional movement/fidelity QA summaries.
##
## IMPORTANT:
##   This regional-state model has a DIFFERENT state space and likelihood from
##   the seven-river candidate models. Its AICc must NOT be compared with
##   Table 3.
##
################################################################################

set.seed(7)

# ==============================================================================
# 1. PROJECT / MARK SETUP
# ==============================================================================

source(file.path("code", "_project_paths.R"), local = TRUE)
results_dir <- regional_results_dir
configure_mark(require_mark = TRUE)

suppressPackageStartupMessages({
  library(RMark)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# Run MARK in its work directory, but use absolute paths for all project files.
setwd(mark_work)

# ==============================================================================
# 2. LOAD THE ADULT ENCOUNTER HISTORIES USED BY THE SEVEN-RIVER MODELS
# ==============================================================================

adult_file <- file.path(
  data_dir,
  "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
)

if (!file.exists(adult_file)) {
  stop("Missing adult encounter-history file:\n", adult_file)
}

dat <- readRDS(adult_file)

required_cols <- c("ch", "freq")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Adult RDS is missing required column(s): ",
       paste(missing_cols, collapse = ", "))
}

if (length(unique(nchar(dat$ch))) != 1L) {
  stop("Encounter histories do not all have the same number of annual occasions.")
}

n_occasions <- unique(nchar(dat$ch))
message("Rows loaded: ", nrow(dat))
message("Annual occasions per encounter history: ", n_occasions)

if (n_occasions != 13L) {
  warning("Expected 13 annual occasions (2010-2022), but found ", n_occasions, ".")
}

# The adult-history RDS should contain only adult river codes plus 0 and .
allowed_chars <- c("A", "B", "C", "E", "F", "G", "H", "0", ".")
all_chars <- unique(unlist(strsplit(dat$ch, "", fixed = TRUE)))
unexpected_chars <- setdiff(all_chars, allowed_chars)

if (length(unexpected_chars) > 0) {
  stop("Unexpected character(s) in adult encounter histories: ",
       paste(unexpected_chars, collapse = ", "),
       "\nThis script expects the adult-only RDS produced by Script 1.")
}

# ==============================================================================
# 3. REMOVE 2022-ONLY ADULT HISTORIES
# ==============================================================================
#
# Match Script 2 exactly: fish first entering as adults in 2022 cannot inform
# a subsequent annual survival interval.

ch2022 <- paste0("000000000000", c("A", "B", "C", "E", "F", "G", "H"))
keep <- !dat$ch %in% ch2022

message("2022-only adult histories removed: ", sum(!keep))

dat <- dat[keep, , drop = FALSE]

message("Rows retained for regional-state model: ", nrow(dat))
if (nrow(dat) != 1017L) {
  warning(
    "Expected 1,017 adult-analysis fish after removing 2022-only histories, ",
    "but retained ", nrow(dat), ". Verify against Script 2."
  )
}

# ==============================================================================
# 4. COLLAPSE SEVEN RIVER STATES TO FOUR REGIONAL STATES
# ==============================================================================
#
# River-state key from the seven-river analysis:
#   A = Apalachicola
#   B = Choctawhatchee
#   C = Escambia
#   E = Pascagoula
#   F = Pearl
#   G = Suwannee
#   H = Yellow
#
# Regional-state key:
#   W = West            (F, E)
#   P = Pensacola Bay   (C, H)
#   C = Choctawhatchee  (B)
#   E = East            (A, G)
#
# IMPORTANT: perform the mapping character-by-character. Sequential gsub()
# replacements would be unsafe because C and E are both old river codes and
# new regional codes.

river_to_region <- c(
  A = "E",  # Apalachicola -> East
  B = "C",  # Choctawhatchee -> Choctawhatchee
  C = "P",  # Escambia -> Pensacola Bay
  E = "W",  # Pascagoula -> West
  F = "W",  # Pearl -> West
  G = "E",  # Suwannee -> East
  H = "P"   # Yellow -> Pensacola Bay
)

collapse_history <- function(ch) {
  x <- strsplit(ch, "", fixed = TRUE)[[1]]

  is_state <- x %in% names(river_to_region)
  x[is_state] <- unname(river_to_region[x[is_state]])

  paste0(x, collapse = "")
}

regional_dat <- dat %>%
  transmute(
    ch = vapply(ch, collapse_history, character(1)),
    freq = as.numeric(freq)
  )

# IMPORTANT: preserve the original row-level frequency coding exactly.
# Do NOT aggregate duplicate encounter histories here. The source MARK input
# can contain frequency values other than +1, and preserving the original
# rows/frequencies is the safest way to reproduce the historical model object
# and its reported likelihood/AICc.

message("Regional encounter-history rows retained: ", nrow(regional_dat))
message("Unique regional encounter histories: ", n_distinct(regional_dat$ch))
message("Sum of source frequencies: ", sum(regional_dat$freq, na.rm = TRUE))
message("Frequency distribution:")
print(table(regional_dat$freq, useNA = "ifany"))

if (nrow(regional_dat) != nrow(dat)) {
  stop("Regional state collapse unexpectedly changed the number of input rows.")
}

regional_chars <- unique(unlist(strsplit(regional_dat$ch, "", fixed = TRUE)))
unexpected_regional <- setdiff(regional_chars, c("W", "P", "C", "E", "0", "."))

if (length(unexpected_regional) > 0) {
  stop("Unexpected regional state code(s): ",
       paste(unexpected_regional, collapse = ", "))
}

message("Regional states present: ",
        paste(sort(intersect(c("W", "P", "C", "E"), regional_chars)),
              collapse = ", "))

# Save the regional encounter histories as a QA artifact.
regional_ch_file <- file.path(results_dir, "regional_state_encounter_histories.csv")
write.csv(regional_dat, regional_ch_file, row.names = FALSE)
message("Saved regional encounter-history audit file: ", regional_ch_file)

# ==============================================================================
# 5. PROCESS THE FOUR-STATE MULTISTRATA DATA
# ==============================================================================

pd_reg <- process.data(regional_dat, model = "Multistrata")

dd_reg <- make.design.data(
  pd_reg,
  parameters = list(
    S   = list(pim.type = "time"),
    p   = list(pim.type = "time"),
    Psi = list(pim.type = "constant")
  )
)

message("Regional S strata: ",
        paste(sort(unique(as.character(dd_reg$S$stratum))), collapse = ", "))
message("Regional p strata: ",
        paste(sort(unique(as.character(dd_reg$p$stratum))), collapse = ", "))
message("Regional Psi origin strata: ",
        paste(sort(unique(as.character(dd_reg$Psi$stratum))), collapse = ", "))

# ==============================================================================
# 6. FIX BIOLOGICALLY IMPLAUSIBLE DIRECT REGIONAL TRANSITIONS
# ==============================================================================
#
# Derived from the seven-river movement constraints:
#
#   West can transition to Pensacola Bay or Choctawhatchee, but not directly East.
#   East can transition to Pensacola Bay or Choctawhatchee, but not directly West.
#
# All other between-region transitions are allowed.
#
# As in Script 2, MARK also carries same-state Psi rows internally; these are not
# treated as impossible transitions here.

psi_we <- as.numeric(
  row.names(dd_reg$Psi[
    dd_reg$Psi$stratum == "W" & dd_reg$Psi$tostratum == "E",
    , drop = FALSE
  ])
)

psi_ew <- as.numeric(
  row.names(dd_reg$Psi[
    dd_reg$Psi$stratum == "E" & dd_reg$Psi$tostratum == "W",
    , drop = FALSE
  ])
)

psi_fixed_idx <- unique(c(psi_we, psi_ew))
psi_fixed_val <- rep(0, length(psi_fixed_idx))

if (length(psi_fixed_idx) == 0) {
  stop("No W<->E Psi rows were found to fix. Inspect dd_reg$Psi before fitting.")
}

Psi.region_to_region <- list(
  formula = ~ -1 + stratum:tostratum,
  fixed = list(
    index = psi_fixed_idx,
    value = psi_fixed_val
  )
)

message("Fixed direct regional transitions to zero: W->E and E->W")
message("Number of fixed Psi design rows: ", length(psi_fixed_idx))

# ==============================================================================
# 7. FIT THE REGIONAL-STATE MODEL
# ==============================================================================

message("\nFitting regional-state model...")

regional_mod <- mark(
  pd_reg,
  dd_reg,
  model = "Multistrata",
  silent = TRUE,
  output = FALSE,
  model.parameters = list(
    S = list(
      formula = ~ -1 + stratum:time,
      link = "sin"
    ),
    p = list(
      formula = ~ -1 + stratum:time,
      link = "sin"
    ),
    Psi = Psi.region_to_region
  ),
  model.name = "regional_S(region:time)p(region:time)Psi(region_to_region)"
)

authoritative_file <- file.path(results_dir, "regional_state_model.RDS")

saveRDS(regional_mod, authoritative_file)

message("Saved regional model: ", authoritative_file)
message("Model K = ", regional_mod$results$npar)
message("Model AICc = ", round(regional_mod$results$AICc, 6))
message("Model -logL = ", round(regional_mod$results$lnl, 6))

if (!is.null(regional_mod$results$singular) &&
    any(regional_mod$results$singular)) {
  warning("Regional model has singular parameter(s). Inspect before use.")
} else {
  message("Singular-parameter check: none flagged.")
}

# ==============================================================================
# 8. EXTRACT REGIONAL MOVEMENT AND FIDELITY
# ==============================================================================

region_state_map <- c(
  W = "West",
  P = "Pensacola Bay",
  C = "Choctawhatchee",
  E = "East"
)

region_order <- c("West", "Pensacola Bay", "Choctawhatchee", "East")

get_real <- function(mod, parameter) {
  rr <- as.data.frame(mod$results$real)
  rr$MARK_row <- rownames(mod$results$real)

  rr %>%
    filter(str_starts(MARK_row, paste0(parameter, " ")))
}

parse_regional_psi <- function(mod) {
  psi <- get_real(mod, "Psi")

  # Remove representative rows for structurally fixed transitions when MARK
  # identifies them as fixed.
  if ("fixed" %in% names(psi)) {
    psi <- psi[
      is.na(psi$fixed) | as.character(psi$fixed) != "Fixed",
      , drop = FALSE
    ]
  }

  psi <- psi %>%
    mutate(
      from_code = str_match(MARK_row, "^Psi s([A-Z])")[, 2],
      to_code   = str_match(MARK_row, " to([A-Z])")[, 2]
    )

  if (any(is.na(psi$from_code)) || any(is.na(psi$to_code))) {
    stop("Could not parse regional Psi row names. Inspect mod$results$real.")
  }

  psi %>%
    mutate(
      From_Region = unname(region_state_map[from_code]),
      To_Region   = unname(region_state_map[to_code])
    ) %>%
    transmute(
      From_Region,
      To_Region,
      Estimate = estimate,
      SE = se,
      LCL_95 = lcl,
      UCL_95 = ucl
    )
}

psi_regional <- parse_regional_psi(regional_mod)

fidelity_regional <- psi_regional %>%
  group_by(From_Region) %>%
  summarise(
    Probability_Leave = sum(Estimate, na.rm = TRUE),
    Fidelity = 1 - Probability_Leave,
    .groups = "drop"
  ) %>%
  mutate(From_Region = factor(From_Region, levels = region_order)) %>%
  arrange(From_Region)

message("\nRegional movement estimates:")
print(
  psi_regional %>%
    mutate(
      From_Region = factor(From_Region, levels = region_order),
      To_Region = factor(To_Region, levels = region_order)
    ) %>%
    arrange(From_Region, To_Region) %>%
    mutate(across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 4)))
)

message("\nRegional fidelity:")
print(
  fidelity_regional %>%
    mutate(
      Probability_Leave = round(Probability_Leave, 4),
      Fidelity = round(Fidelity, 4)
    )
)

# Save a compact QA summary.
qa_summary_file <- file.path(
  results_dir,
  "regional_state_model_movement_QA.csv"
)

write.csv(
  psi_regional %>%
    left_join(
      fidelity_regional %>%
        mutate(From_Region = as.character(From_Region)),
      by = "From_Region"
    ),
  qa_summary_file,
  row.names = FALSE
)

message("Saved regional movement QA file: ", qa_summary_file)

# ==============================================================================
# 9. BENCHMARK AGAINST THE TABLE A5 VALUES WE EXPECT TO RECOVER
# ==============================================================================
#
# These are approximate benchmarks from the existing manuscript workflow.
# They are printed as a diagnostic only; they are NOT used in model fitting.

expected_fidelity <- c(
  "West" = 0.9815,
  "Pensacola Bay" = 0.8989,
  "Choctawhatchee" = 0.9162,
  "East" = 0.9792
)

fidelity_check <- fidelity_regional %>%
  mutate(
    From_Region = as.character(From_Region),
    Expected = unname(expected_fidelity[From_Region]),
    Difference = Fidelity - Expected
  )

message("\nApproximate fidelity benchmark check:")
print(
  fidelity_check %>%
    mutate(across(c(Fidelity, Expected, Difference), ~round(.x, 4)))
)

key_transition_check <- psi_regional %>%
  filter(
    (From_Region == "Pensacola Bay" & To_Region == "Choctawhatchee") |
    (From_Region == "Choctawhatchee" & To_Region == "Pensacola Bay")
  ) %>%
  mutate(
    Expected = case_when(
      From_Region == "Pensacola Bay" &
        To_Region == "Choctawhatchee" ~ 0.0873,
      From_Region == "Choctawhatchee" &
        To_Region == "Pensacola Bay" ~ 0.0709,
      TRUE ~ NA_real_
    ),
    Difference = Estimate - Expected
  )

message("\nKey transition benchmark check:")
print(
  key_transition_check %>%
    mutate(across(c(Estimate, Expected, Difference), ~round(.x, 4)))
)

# ==============================================================================
# 10. SESSION / REPRODUCIBILITY INFORMATION
# ==============================================================================

session_file <- file.path(
  results_dir,
  "regional_state_model_sessionInfo.txt"
)

writeLines(capture.output(sessionInfo()), session_file)
message("Saved sessionInfo(): ", session_file)

message("\nScript 6 complete.")
