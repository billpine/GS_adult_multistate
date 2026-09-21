################################################################################
# Gulf Sturgeon CJFAS Revision 1
# Script 2: Fit seven-river multistate candidate models
# FINAL production version, 2026-09-18
#
# Purpose
#   1. Load the validated adult encounter histories from Script 1.
#   2. Remove fish first entering the adult analysis in 2022 because no
#      subsequent annual survival interval is available.
#   3. Build the Multistrata design data.
#   4. Define the four regional survival groups.
#   5. Fix geographically impossible river-to-river transitions to zero.
#   6. Constrain final-occasion p to equal the preceding occasion within river
#      for all p(river:time) models, resolving terminal S-p confounding.
#   7. Fit/load the 25-model candidate set in a NEW results directory.
#   8. Calculate AICc model-selection quantities from UNROUNDED AICc values.
#
#   The final detection occasion is constrained to equal the preceding
#   occasion within river for p(river:time) models. This resolves the known
#   terminal survival-detection confounding while leaving nonterminal
#   survival estimates essentially unchanged in the validation comparison.
#
#   Profile-CI versions of the top model (mod1) and the event-diagnostic
#   river-by-year model (mod7) are fit/loaded at the end.
#
# Important accounting
#   - Script 1 contains 1,070 adult encounter-history rows.
#   - 53 2022-only histories are removed here.
#   - 1,017 fish therefore inform the adult multistate analysis.
#   - freq = -1 denotes loss on capture/censoring in MARK. sum(freq) is NOT N.
#
# Candidate set
#   Eight survival structures:
#     S(.), S(region), S(river), S(time),
#     S(region+time), S(river+time),
#     S(region:time), S(river:time)
#
#   Each is crossed with:
#     p(river:time), p(river), p(.)
#
#   For p(river:time), the final detection occasion is constrained to equal
#   the preceding occasion within each river. This removes the terminal
#   survival-detection confounding identified in the original fit.
#
#   This gives 24 seven-river models with Psi(stratum:tostratum), plus the
#   fully constant null S(.)p(.)Psi(.), for 25 candidate models total.
#
# Validated final ranking under terminal-p constraint
#   mod1: S(region:time)p(river:time)Psi(strattostrat)
#         K = 148, AICc = 6133.295, weight = 0.980688
#   mod2: Delta AICc = 8.715
#   mod3: Delta AICc = 10.250
#
# Only mod17 and mod19 retained singular S:time12 parameters; both are
# noncompetitive models with p(.).
#
# Downstream products
#   Script 7 is the authoritative manuscript-table script.
#   Script 8 is the authoritative manuscript-figure script.
#   Obsolete appendix-table printing and legacy pairwise contrasts have been
#   removed from this script.
################################################################################

rm(list = ls())
gc()
set.seed(7)

# ---- project / MARK setup ----------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)
results_dir <- seven_results_dir
configure_mark(require_mark = TRUE)

suppressPackageStartupMessages({
  library(RMark)
  library(dplyr)
})

message("Project: ", proj)
message("Model results directory: ", results_dir)
message("MARK working directory: ", mark_work)

# RMark/MARK temporary files are written in the configured MARK work directory.
setwd(mark_work)

# ---- 1. Load validated adult encounter histories -----------------------------

message("\n1. Loading adult encounter histories...")

adult_path <- file.path(
  data_dir,
  "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
)

pre2010_path <- file.path(
  data_dir,
  "pre2010.RDS"
)

at_inf_adults <- readRDS(adult_path)
pre2010       <- readRDS(pre2010_path)

at_inf_adults$FishID <- as.numeric(at_inf_adults$FishID)
pre2010$FishID       <- as.numeric(pre2010$FishID)

message("  Rows loaded: ", nrow(at_inf_adults),
        " [validated expectation = 1070]")
message("  Unique FishID loaded: ", dplyr::n_distinct(at_inf_adults$FishID))
message("  Analysis fish with pre-2010 tagging history: ",
        sum(at_inf_adults$FishID %in% pre2010$FishID))

if (nrow(at_inf_adults) != 1070 ||
    n_distinct(at_inf_adults$FishID) != 1070) {
  stop(
    "Adult encounter-history accounting changed. Expected 1,070 rows / ",
    "1,070 unique fish before removing 2022-only histories."
  )
}


# The prior expression at_inf_adults$first$first.river was invalid.
first_river <- as.character(at_inf_adults$first)

message("  First adult river distribution before 2022-only removal:")
print(table(first_river, useNA = "ifany"))

# Remove fish whose only adult encounter is in the final 2022 occasion.
adult_states <- c("A", "B", "C", "E", "F", "G", "H")
ch2022 <- paste0("000000000000", adult_states)

keep <- !at_inf_adults$ch %in% ch2022

n_removed_2022 <- sum(!keep)

at_inf_adults <- at_inf_adults[keep, , drop = FALSE]
first_river   <- first_river[keep]

message("  2022-only histories removed: ", n_removed_2022,
        " [expected = 53]")
message("  Adult analysis rows retained: ", nrow(at_inf_adults),
        " [expected = 1017]")
message("  Adult analysis unique FishID: ",
        n_distinct(at_inf_adults$FishID))

if (n_removed_2022 != 53 ||
    nrow(at_inf_adults) != 1017 ||
    n_distinct(at_inf_adults$FishID) != 1017) {
  stop(
    "Adult-analysis sample accounting changed. Expected 53 2022-only histories ",
    "removed and 1,017 unique fish retained."
  )
}

message("  First adult river distribution in analysis population:")
print(table(first_river, useNA = "ifany"))

# Frequency audit. Negative frequencies are intentional MARK loss-on-capture
# records and must not be interpreted as negative fish.
message("  MARK frequency distribution:")
print(table(at_inf_adults$freq, useNA = "ifany"))
message("  NOTE: sum(freq) is not the analysis sample size.")

# ---- 2. Process seven-river Multistrata data ---------------------------------

message("\n2. Processing Multistrata data...")

at_inf_adults2 <- at_inf_adults %>%
  select(ch, freq)

pd <- process.data(
  at_inf_adults2,
  model = "Multistrata"
)

dd <- make.design.data(
  pd,
  parameters = list(
    S   = list(pim.type = "time"),
    p   = list(pim.type = "time"),
    Psi = list(pim.type = "constant")
  )
)

observed_states <- sort(unique(as.character(dd$S$stratum)))
expected_states <- sort(adult_states)

message("  MARK states: ", paste(observed_states, collapse = ", "))

if (!identical(observed_states, expected_states)) {
  stop(
    "Unexpected Multistrata state set. Expected ",
    paste(expected_states, collapse = ", "),
    "; observed ",
    paste(observed_states, collapse = ", ")
  )
}

# ---- 3. Resolve terminal survival-detection confounding ----------------------

message("\n3. Pooling final detection occasion with preceding occasion...")

# In time-varying live-recapture models, the final survival interval and final
# detection probability can be confounded. The original candidate set showed
# this directly: S:time12 was singular in the additive-time survival models,
# while final-occasion survival estimates were at/near 1.0.
#
# To make the terminal survival interval estimable, constrain p in the final
# occasion to equal p in the preceding occasion, separately for each river.
#
# We preserve the original dd$p$time column for diagnostics and create p_time
# only for the p(river:time) model formulas.

p_time_num <- suppressWarnings(as.numeric(as.character(dd$p$time)))

if (anyNA(p_time_num)) {
  stop("dd$p$time could not be converted cleanly to numeric time indices.")
}

p_times <- sort(unique(p_time_num))

if (length(p_times) < 2) {
  stop("Need at least two p time indices to pool the terminal occasion.")
}

terminal_p_time <- max(p_times)
previous_p_time <- p_times[length(p_times) - 1]

# S has one fewer time interval than p has encounter occasions.
s_time_num <- suppressWarnings(as.numeric(as.character(dd$S$time)))
if (anyNA(s_time_num)) {
  stop("dd$S$time could not be converted cleanly to numeric time indices.")
}
terminal_s_time <- max(s_time_num)

dd$p$p_time_num <- p_time_num
dd$p$p_time_num[
  dd$p$p_time_num == terminal_p_time
] <- previous_p_time

dd$p$p_time <- factor(
  dd$p$p_time_num,
  levels = sort(unique(dd$p$p_time_num))
)

n_p_rows_pooled <- sum(p_time_num == terminal_p_time)

message("  Original p time indices: ",
        paste(p_times, collapse = ", "))
message("  Final p time index pooled: ", terminal_p_time,
        " -> ", previous_p_time)
message("  p design rows affected: ", n_p_rows_pooled,
        " [one terminal row per river/state design combination]")
message("  p_time levels after pooling: ",
        paste(levels(dd$p$p_time), collapse = ", "))

message("  Original p time -> pooled p_time cross-tab:")
print(table(original_time = p_time_num, pooled_time = dd$p$p_time_num))

# ---- 4. Define regional survival groups --------------------------------------

message("\n3. Defining regional survival groups...")

# Regional grouping used only for S structures in the seven-river model:
#   West          = Pearl (F) + Pascagoula (E)
#   Pensacola Bay = Escambia (C) + Yellow (H)
#   Choctawhatchee= Choctawhatchee (B)
#   East          = Apalachicola (A) + Suwannee (G)
#
# Internal factor labels are retained for exact reproducibility with the
# previously fitted MARK objects:
#   west, ebay, choc, east

dd$S$regs <- factor(
  case_when(
    dd$S$stratum %in% c("F", "E") ~ "west",
    dd$S$stratum %in% c("C", "H") ~ "ebay",
    dd$S$stratum == "B"           ~ "choc",
    dd$S$stratum %in% c("A", "G") ~ "east",
    TRUE                           ~ NA_character_
  ),
  levels = c("choc", "east", "ebay", "west")
)

if (any(is.na(dd$S$regs))) {
  stop("At least one S design-data row could not be assigned to a region.")
}

message("  Regional S design rows:")
print(table(dd$S$stratum, dd$S$regs))

# ---- 5. Fix geographically impossible river transitions ---------------------

message("\n5. Defining allowed river-to-river transitions...")

# Each source river is paired with the destination states allowed by the
# original analysis. Transitions not listed are structurally fixed to zero.
allowed_to <- list(
  F = c("F", "E"),
  E = c("F", "E", "C", "H", "B"),
  C = c("E", "C", "H", "B", "A"),
  H = c("E", "C", "H", "B", "A"),
  B = c("E", "C", "H", "B", "A"),
  A = c("C", "H", "B", "A", "G"),
  G = c("H", "A", "G")
)

get_fixed_psi_indices <- function(source_state, allowed_states) {

  rows <- dd$Psi[
    dd$Psi$stratum == source_state &
      !dd$Psi$tostratum %in% allowed_states,
    ,
    drop = FALSE
  ]

  as.numeric(row.names(rows))
}

Psi.indicesC <- unique(
  unlist(
    Map(
      get_fixed_psi_indices,
      names(allowed_to),
      allowed_to
    )
  )
)

Psi.values <- rep(0, length(Psi.indicesC))

Psi.strattostrat <- list(
  formula = ~ -1 + stratum:tostratum,
  fixed = list(
    index = Psi.indicesC,
    value = Psi.values
  )
)

message("  Psi design rows fixed to zero: ", length(Psi.indicesC))

# Save a human-readable audit of fixed source -> destination routes.
psi_fixed_audit <- dd$Psi[
  as.character(row.names(dd$Psi)) %in% as.character(Psi.indicesC),
  c("stratum", "tostratum"),
  drop = FALSE
]

psi_fixed_audit <- unique(psi_fixed_audit)
psi_fixed_audit <- psi_fixed_audit[
  order(psi_fixed_audit$stratum, psi_fixed_audit$tostratum),
  ,
  drop = FALSE
]

message("  Unique structural-zero routes:")
print(psi_fixed_audit, row.names = FALSE)

# ---- 6. Cache helpers --------------------------------------------------------

rds_path <- function(filename) {
  file.path(results_dir, filename)
}

cached <- function(filename) {

  p <- rds_path(filename)

  if (file.exists(p)) {
    message("  cached: ", filename)
    return(readRDS(p))
  }

  NULL
}

saved <- function(model_object, filename) {

  saveRDS(
    model_object,
    rds_path(filename)
  )

  if (!is.null(model_object$results$singular) &&
      length(model_object$results$singular) > 0) {
    warning(
      "Singular parameter index/indices in ",
      model_object$model.name,
      ": ",
      paste(model_object$results$singular, collapse = ", ")
    )
  }

  model_object
}

fit_or_load <- function(filename, model_call) {

  existing <- cached(filename)

  if (!is.null(existing)) {
    return(existing)
  }

  message("  fitting: ", filename)

  model_object <- eval.parent(substitute(model_call))

  saved(model_object, filename)
}

# ---- 7. Fit/load the 25 candidate models ------------------------------------

message("\n7. Loading/fitting the 25-model candidate set...")

# p(river:time) ---------------------------------------------------------------
# NOTE: model.name retains the manuscript shorthand p(river:time), but the
# actual formula uses stratum:p_time, where the final detection occasion is
# pooled with the preceding occasion within each river.

mod1 <- fit_or_load(
  "mod1.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + regs:time, link = "sin"),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region:time)p(river:time)Psi(strattostrat)"
  )
)

mod2 <- fit_or_load(
  "mod2.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs + time),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region+time)p(river:time)Psi(strattostrat)"
  )
)

mod3 <- fit_or_load(
  "mod3.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum + time),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river+time)p(river:time)Psi(strattostrat)"
  )
)

mod5 <- fit_or_load(
  "mod5.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + stratum:time, link = "sin"),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river:time)p(river:time)Psi(strattostrat)"
  )
)

mod6 <- fit_or_load(
  "mod6.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + time, link = "sin"),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(time)p(river:time)Psi(strattostrat)"
  )
)

mod11 <- fit_or_load(
  "mod11.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river)p(river:time)Psi(strattostrat)"
  )
)

mod12 <- fit_or_load(
  "mod12.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region)p(river:time)Psi(strattostrat)"
  )
)

mod13 <- fit_or_load(
  "mod13.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ 1),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "S(.)p(river:time)Psi(strattostrat)"
  )
)

# p(river) --------------------------------------------------------------------

mod4 <- fit_or_load(
  "mod4.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + regs:time, link = "sin"),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region:time)p(river)Psi(strattostrat)"
  )
)

mod7 <- fit_or_load(
  "mod7.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + stratum:time, link = "sin"),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river:time)p(river)Psi(strattostrat)"
  )
)

mod8 <- fit_or_load(
  "mod8.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum + time),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river+time)p(river)Psi(strattostrat)"
  )
)

mod9 <- fit_or_load(
  "mod9.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs + time),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region+time)p(river)Psi(strattostrat)"
  )
)

mod10 <- fit_or_load(
  "mod10.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + time, link = "sin"),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(time)p(river)Psi(strattostrat)"
  )
)

mod16 <- fit_or_load(
  "mod16.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river)p(river)Psi(strattostrat)"
  )
)

mod18 <- fit_or_load(
  "mod18.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region)p(river)Psi(strattostrat)"
  )
)

mod20 <- fit_or_load(
  "mod20.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ 1),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "S(.)p(river)Psi(strattostrat)"
  )
)

# p(.) ------------------------------------------------------------------------

mod14 <- fit_or_load(
  "mod14.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + regs:time, link = "sin"),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region:time)p(.)Psi(strattostrat)"
  )
)

mod15 <- fit_or_load(
  "mod15.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + stratum:time, link = "sin"),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river:time)p(.)Psi(strattostrat)"
  )
)

mod17 <- fit_or_load(
  "mod17.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum + time),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river+time)p(.)Psi(strattostrat)"
  )
)

mod19 <- fit_or_load(
  "mod19.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs + time),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region+time)p(.)Psi(strattostrat)"
  )
)

mod21 <- fit_or_load(
  "mod21.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ -1 + time, link = "sin"),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(time)p(.)Psi(strattostrat)"
  )
)

mod22 <- fit_or_load(
  "mod22.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ stratum),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(river)p(.)Psi(strattostrat)"
  )
)

mod23 <- fit_or_load(
  "mod23.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ regs),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(region)p(.)Psi(strattostrat)"
  )
)

mod24 <- fit_or_load(
  "mod24.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ 1),
      p = list(formula = ~ 1),
      Psi = Psi.strattostrat
    ),
    model.name = "S(.)p(.)Psi(strattostrat)"
  )
)

# Fully constant null ----------------------------------------------------------

mod25 <- fit_or_load(
  "mod25.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    model.parameters = list(
      S = list(formula = ~ 1),
      p = list(formula = ~ 1),
      Psi = list(formula = ~ 1)
    ),
    model.name = "S(.)p(.)Psi(.)"
  )
)

message("  All 25 candidate model objects available.")

# ---- 8. Candidate-set and AICc audit -----------------------------------------

message("\n8. Auditing candidate model set and AICc ranking...")

all_mods <- list(
  mod1 = mod1, mod2 = mod2, mod3 = mod3, mod4 = mod4, mod5 = mod5,
  mod6 = mod6, mod7 = mod7, mod8 = mod8, mod9 = mod9, mod10 = mod10,
  mod11 = mod11, mod12 = mod12, mod13 = mod13, mod14 = mod14,
  mod15 = mod15, mod16 = mod16, mod17 = mod17, mod18 = mod18,
  mod19 = mod19, mod20 = mod20, mod21 = mod21, mod22 = mod22,
  mod23 = mod23, mod24 = mod24, mod25 = mod25
)

model_names <- vapply(
  all_mods,
  function(m) m$model.name,
  character(1)
)

if (length(all_mods) != 25 ||
    length(unique(model_names)) != 25) {
  stop("Candidate-model audit failed: expected 25 unique model names.")
}

expected_model_names <- c(
  "S(region:time)p(river:time)Psi(strattostrat)",
  "S(region+time)p(river:time)Psi(strattostrat)",
  "S(river+time)p(river:time)Psi(strattostrat)",
  "S(region:time)p(river)Psi(strattostrat)",
  "S(river:time)p(river:time)Psi(strattostrat)",
  "S(time)p(river:time)Psi(strattostrat)",
  "S(river:time)p(river)Psi(strattostrat)",
  "S(river+time)p(river)Psi(strattostrat)",
  "S(region+time)p(river)Psi(strattostrat)",
  "S(time)p(river)Psi(strattostrat)",
  "S(river)p(river:time)Psi(strattostrat)",
  "S(region)p(river:time)Psi(strattostrat)",
  "S(.)p(river:time)Psi(strattostrat)",
  "S(region:time)p(.)Psi(strattostrat)",
  "S(river:time)p(.)Psi(strattostrat)",
  "S(river)p(river)Psi(strattostrat)",
  "S(river+time)p(.)Psi(strattostrat)",
  "S(region)p(river)Psi(strattostrat)",
  "S(region+time)p(.)Psi(strattostrat)",
  "S(.)p(river)Psi(strattostrat)",
  "S(time)p(.)Psi(strattostrat)",
  "S(river)p(.)Psi(strattostrat)",
  "S(region)p(.)Psi(strattostrat)",
  "S(.)p(.)Psi(strattostrat)",
  "S(.)p(.)Psi(.)"
)

if (!setequal(model_names, expected_model_names)) {
  stop("Candidate-model names do not match the intended 24 + null model set.")
}

aicc_values <- vapply(
  all_mods,
  function(m) m$results$AICc,
  numeric(1)
)

delta_aicc <- aicc_values - min(aicc_values)

weights <- exp(-0.5 * delta_aicc) /
  sum(exp(-0.5 * delta_aicc))

table3_audit <- data.frame(
  Model_ID = names(all_mods),
  Model = model_names,
  K = vapply(all_mods, function(m) m$results$npar, numeric(1)),
  AICc_unrounded = aicc_values,
  DeltaAICc_unrounded = delta_aicc,
  Weight_unrounded = weights,
  NegLogLik = vapply(all_mods, function(m) m$results$lnl, numeric(1)),
  stringsAsFactors = FALSE
)

table3_audit <- table3_audit[
  order(table3_audit$AICc_unrounded),
  ,
  drop = FALSE
]

print(
  transform(
    table3_audit,
    AICc = round(AICc_unrounded, 2),
    DeltaAICc = round(DeltaAICc_unrounded, 2),
    Weight = round(Weight_unrounded, 3)
  )[
    ,
    c("Model_ID", "Model", "K", "AICc", "DeltaAICc", "Weight", "NegLogLik")
  ],
  row.names = FALSE
)

write.csv(
  table3_audit,
  file.path(results_dir, "Script2_candidate_model_audit.csv"),
  row.names = FALSE
)

top3 <- table3_audit[1:3, ]

message(
  "  Top 3: ",
  paste0(
    top3$Model_ID,
    " [AICc=",
    sprintf("%.3f", top3$AICc_unrounded),
    ", w=",
    sprintf("%.3f", top3$Weight_unrounded),
    "]",
    collapse = "; "
  )
)

message("  Final top 5 models:")
top5_print <- head(
  transform(
    table3_audit,
    AICc = round(AICc_unrounded, 3),
    DeltaAICc = round(DeltaAICc_unrounded, 3),
    Weight = round(Weight_unrounded, 6)
  )[
    ,
    c("Model_ID", "Model", "K", "AICc", "DeltaAICc", "Weight")
  ],
  5
)
print(top5_print, row.names = FALSE)

# 
expected_top_ids <- c("mod1", "mod2", "mod3", "mod5", "mod6")

if (!identical(head(table3_audit$Model_ID, 5), expected_top_ids)) {
  stop(
    "Final top-five model ordering changed. Expected: ",
    paste(expected_top_ids, collapse = ", "),
    "; observed: ",
    paste(head(table3_audit$Model_ID, 5), collapse = ", "),
    ". Reconcile before running downstream manuscript scripts."
  )
}

if (abs(table3_audit$AICc_unrounded[1] - 6133.295) > 0.05) {
  warning(
    "Top-model AICc differs from validated final value. Observed = ",
    table3_audit$AICc_unrounded[1],
    "; expected approximately 6133.295."
  )
}

if (abs(table3_audit$Weight_unrounded[1] - 0.980688) > 0.005) {
  warning(
    "Top-model weight differs from validated final value. Observed = ",
    table3_audit$Weight_unrounded[1],
    "; expected approximately 0.980688."
  )
}

# Report singular status across all candidate models.
singular_indices <- lapply(
  all_mods,
  function(m) {
    s <- m$results$singular
    if (is.null(s) || length(s) == 0) {
      integer(0)
    } else {
      as.integer(s)
    }
  }
)

singular_flag <- vapply(
  singular_indices,
  function(x) length(x) > 0,
  logical(1)
)

message("  Candidate models with singular parameter flags: ",
        sum(singular_flag))

singular_details <- data.frame(
  Model_ID = character(0),
  Model = character(0),
  Beta_Index = integer(0),
  Beta_Name = character(0),
  Estimate = numeric(0),
  SE = numeric(0),
  stringsAsFactors = FALSE
)

if (any(singular_flag)) {

  for (id in names(all_mods)[singular_flag]) {

    m <- all_mods[[id]]
    idx <- singular_indices[[id]]

    beta <- m$results$beta

    for (j in idx) {

      singular_details <- rbind(
        singular_details,
        data.frame(
          Model_ID = id,
          Model = m$model.name,
          Beta_Index = j,
          Beta_Name = if (j <= nrow(beta)) rownames(beta)[j] else NA_character_,
          Estimate = if (j <= nrow(beta)) beta$estimate[j] else NA_real_,
          SE = if (j <= nrow(beta)) beta$se[j] else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  print(singular_details, row.names = FALSE)

  write.csv(
    singular_details,
    file.path(results_dir, "Script2_singular_parameter_audit.csv"),
    row.names = FALSE
  )
}

expected_singular_models <- c("mod17", "mod19")

if (!identical(sort(names(all_mods)[singular_flag]),
               sort(expected_singular_models))) {
  warning(
    "Singular-model set differs from validated final run. Observed: ",
    paste(names(all_mods)[singular_flag], collapse = ", "),
    "; expected: mod17, mod19."
  )
}

# ---- 9. Terminal-interval real-parameter audit -------------------------------

message("\n9. Auditing final-interval survival and detection...")

extract_terminal_real <- function(m, model_id) {

  real <- m$results$real
  rn <- rownames(real)

  s_rows <- real[
    grepl("^S ", rn) &
      grepl(paste0(" t", terminal_s_time, "$"), rn),
    1:4,
    drop = FALSE
  ]

  p_rows <- real[
    grepl("^p ", rn) &
      grepl(paste0(" t", terminal_p_time, "$"), rn),
    1:4,
    drop = FALSE
  ]

  list(
    model_id = model_id,
    S = s_rows,
    p = p_rows
  )
}

terminal_audit_ids <- head(table3_audit$Model_ID, 5)

terminal_audit <- lapply(
  terminal_audit_ids,
  function(id) extract_terminal_real(all_mods[[id]], id)
)

for (x in terminal_audit) {

  cat("\n", x$model_id, " terminal S rows:\n", sep = "")
  print(x$S)

  cat(x$model_id, " terminal p rows:\n", sep = "")
  print(x$p)
}

# ---- 10. Profile CIs for the final top model --------------------------------

message("\n10. Loading/fitting profile-CI version of final top model...")

mod1_prof <- fit_or_load(
  "mod1_prof.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    profile.int = TRUE,
    model.parameters = list(
      S = list(formula = ~ -1 + regs:time, link = "sin"),
      p = list(formula = ~ -1 + stratum:p_time, link = "sin"),
      Psi = Psi.strattostrat
    ),
    model.name = "prof_S(region:time)p(river:time)"
  )
)

mod1_prof_aicc_diff <- abs(mod1$results$AICc - mod1_prof$results$AICc)
mod1_prof_lnl_diff  <- abs(mod1$results$lnl  - mod1_prof$results$lnl)

message(
  "  mod1 profile audit: |AICc diff| = ",
  format(mod1_prof_aicc_diff, scientific = TRUE),
  "; |NLL diff| = ",
  format(mod1_prof_lnl_diff, scientific = TRUE)
)

if (mod1_prof_aicc_diff > 0.01 || mod1_prof_lnl_diff > 0.01) {
  warning("mod1_prof does not closely match the final mod1 candidate model.")
}

if (!is.null(mod1_prof$results$singular) &&
    length(mod1_prof$results$singular) > 0) {
  warning(
    "mod1_prof has singular parameter index/indices: ",
    paste(mod1_prof$results$singular, collapse = ", ")
  )
}

# ---- 11. Profile CIs for the event-specific mod7 diagnostic ------------------

message("\n11. Loading/fitting profile-CI version of mod7 event diagnostic...")

mod7_prof <- fit_or_load(
  "mod7_prof.RDS",
  mark(
    pd, dd,
    model = "Multistrata",
    silent = TRUE,
    output = FALSE,
    profile.int = TRUE,
    model.parameters = list(
      S = list(formula = ~ -1 + stratum:time, link = "sin"),
      p = list(formula = ~ stratum),
      Psi = Psi.strattostrat
    ),
    model.name = "prof_S(river:time)p(river)Psi(strattostrat)"
  )
)

mod7_prof_aicc_diff <- abs(mod7$results$AICc - mod7_prof$results$AICc)
mod7_prof_lnl_diff  <- abs(mod7$results$lnl  - mod7_prof$results$lnl)

message(
  "  mod7 profile audit: |AICc diff| = ",
  format(mod7_prof_aicc_diff, scientific = TRUE),
  "; |NLL diff| = ",
  format(mod7_prof_lnl_diff, scientific = TRUE)
)

if (mod7_prof_aicc_diff > 0.01 || mod7_prof_lnl_diff > 0.01) {
  warning("mod7_prof does not closely match candidate mod7.")
}

if (!is.null(mod7_prof$results$singular) &&
    length(mod7_prof$results$singular) > 0) {
  warning(
    "mod7_prof has singular parameter index/indices: ",
    paste(mod7_prof$results$singular, collapse = ", ")
  )
}

# ---- final audit -------------------------------------------------------------

message("\n================ SCRIPT 2 FINAL AUDIT ================")
message("Adult histories loaded: 1070")
message("2022-only histories removed: ", n_removed_2022, " [expected 53]")
message("Adult analysis fish: ", nrow(at_inf_adults), " [expected 1017]")
message("Candidate models: ", length(all_mods), " [expected 25]")
message("Unique candidate model names: ", length(unique(model_names)),
        " [expected 25]")
message("Fixed Psi design rows: ", length(Psi.indicesC))
message("Final p time pooled: ", terminal_p_time,
        " -> ", previous_p_time)
message("Terminal S interval index: ", terminal_s_time)
message("Top model: ", top3$Model_ID[1], " | ", top3$Model[1])
message("Top AICc: ", sprintf("%.3f", top3$AICc_unrounded[1]))
message("Top 3 weights: ",
        paste(sprintf("%.6f", top3$Weight_unrounded), collapse = ", "))
message("Candidate singular flags: ", sum(singular_flag))
message("Profile CI objects available: mod1_prof, mod7_prof")
message("Model audit CSV: ",
        file.path(results_dir, "Script2_candidate_model_audit.csv"))
message("======================================================")

message("\nScript 2 FINAL production model workflow complete.")
