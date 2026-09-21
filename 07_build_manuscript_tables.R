################################################################################
## 7_manuscript_tables_revised.R
## Gulf Sturgeon CJFAS manuscript: final/reproducible table generation
## production version 
##
## Key revisions relative to the prior Script 7:
##   1. Table A1 restored as a tag-longevity summary (not constant survival).
##   2. Table A2 summarizes annual river-specific detection from the FINAL top
##      model. In p(river:time) models, final-occasion p is constrained equal
##      to the preceding occasion within river.
##   3. Table A3 retains river x year survival from the alternative
##      S(river+time) model for descriptive appendix use; it is no longer
##      described as a competitive/supported model.
##   4. Table A4 extracts river-to-river Psi from the TOP model and calculates
##      river fidelity = 1 - sum(Psi leaving the river).
##   5. Table A5 reports movement and fidelity from the separate four-region
##      multistate model. Its AICc is NOT compared with the seven-river models.
##   6. Supporting files for top-model regional survival and detection are
##      generated for manuscript Results text.
##   7. AICc weights are calculated from unrounded AICc values.
##   8. Table A1 now summarizes all transmitters that provided active coverage
##      during the adult-analysis period (1,048 tags across 1,017 fish), rather
##      than only pre-2022 deployments on the 985 fish initially tagged as adults.
##   9. Uses the final terminal-p-constrained candidate set in
##      results/models/seven_river/.
##  10. Integrates Table A6 from the river x year / river-detection diagnostic
##      model used for the Hurricane Michael event-specific diagnostic, using
##      likelihood-profile confidence intervals from mod7_prof.RDS.
################################################################################

# ---- project setup -----------------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)

library(tidyverse)
library(RMark)

out_dir     <- table_dir
results_dir <- seven_results_dir
message("Saving manuscript tables to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

# Remove obsolete outputs from earlier versions of this script so the folder
# contains only the current, authoritative table products.
stale_outputs <- c(
  "regional_fidelity.csv",
  "regional_transition_probabilities.csv",
  "TableA1_constant_survival_estimates.csv",
  "TableA2_detection_probability.csv",
  "TableA3_river_year_survival.csv",
  "TableA4_movement_estimates_with_CI.csv",
  "TableA4_movement_matrix_wide.csv",
  "TableA6_river_year_apparent_survival_interactive.csv",
  "Supporting_river_year_apparent_survival_mod3.csv"
)

stale_paths <- file.path(out_dir, stale_outputs)
if (any(file.exists(stale_paths))) {
  message("Removing obsolete outputs from earlier table scripts...")
  unlink(stale_paths[file.exists(stale_paths)])
}

# ---- helper functions --------------------------------------------------------

get_real <- function(mod, parameter) {
  rr <- as.data.frame(mod$results$real)
  rr$MARK_row <- rownames(mod$results$real)

  keep <- intersect(
    c("MARK_row", "estimate", "se", "lcl", "ucl", "fixed", "note"),
    names(rr)
  )

  rr %>%
    filter(str_starts(MARK_row, paste0(parameter, " "))) %>%
    select(all_of(keep))
}

# Parse rows such as "S sA g1 t1" or "p sA g1 t1".
# Time values are converted to sequential study years (2010-2021) by rank,
# rather than assuming that the literal MARK time label always equals 1:12.
parse_stratum_time <- function(tbl, parameter, state_map,
                               first_year = 2010,
                               state_col = "State") {

  parsed <- tbl %>%
    mutate(
      stratum_code = str_match(MARK_row,
                               paste0("^", parameter, " s([A-Z])"))[, 2],
      mark_time = as.numeric(str_match(MARK_row, " t([0-9.]+)")[, 2])
    )

  if (any(is.na(parsed$stratum_code))) {
    stop("Could not parse one or more stratum codes for parameter ", parameter,
         ". Inspect MARK_row values.")
  }
  if (any(is.na(parsed$mark_time))) {
    stop("Could not parse one or more MARK time values for parameter ", parameter,
         ". Inspect MARK_row values.")
  }

  unexpected <- setdiff(unique(parsed$stratum_code), names(state_map))
  if (length(unexpected) > 0) {
    stop("Unexpected stratum code(s): ", paste(unexpected, collapse = ", "))
  }

  time_levels <- sort(unique(parsed$mark_time))
  year_lookup <- setNames(first_year + seq_along(time_levels) - 1,
                          as.character(time_levels))

  parsed$Year <- unname(year_lookup[as.character(parsed$mark_time)])
  parsed[[state_col]] <- unname(state_map[parsed$stratum_code])

  parsed
}

parse_psi <- function(tbl, state_map,
                      from_col = "From", to_col = "To") {

  parsed <- tbl %>%
    mutate(
      from_code = str_match(MARK_row, "^Psi s([A-Z])")[, 2],
      to_code   = str_match(MARK_row, " to([A-Z])")[, 2]
    )

  if (any(is.na(parsed$from_code)) || any(is.na(parsed$to_code))) {
    stop("Could not parse one or more Psi MARK row names. Inspect MARK_row values.")
  }

  unexpected <- setdiff(unique(c(parsed$from_code, parsed$to_code)), names(state_map))
  if (length(unexpected) > 0) {
    stop("Unexpected Psi state code(s): ", paste(unexpected, collapse = ", "))
  }

  parsed[[from_col]] <- unname(state_map[parsed$from_code])
  parsed[[to_col]]   <- unname(state_map[parsed$to_code])

  parsed
}

# ---- labels / mappings -------------------------------------------------------

river_map <- c(
  A = "Apalachicola",
  B = "Choctawhatchee",
  C = "Escambia",
  E = "Pascagoula",
  F = "Pearl",
  G = "Suwannee",
  H = "Yellow"
)

river_alpha <- unname(river_map)

river_order <- c(
  "Pearl", "Pascagoula", "Escambia", "Yellow",
  "Choctawhatchee", "Apalachicola", "Suwannee"
)

# IMPORTANT: these are representative MARK strata for the four survival regions.
# This mapping was verified from the model design/output and fixes the earlier
# East <-> Choctawhatchee label swap in downstream code.
regional_survival_map <- c(
  A = "East",
  B = "Choctawhatchee",
  C = "Pensacola Bay",
  E = "West"
)

regional_state_map <- c(
  W = "West",
  P = "Pensacola Bay",
  C = "Choctawhatchee",
  E = "East"
)

region_order <- c("West", "Pensacola Bay", "Choctawhatchee", "East")

# ---- load models -------------------------------------------------------------

message("Loading river-state candidate models...")

all_mods <- lapply(seq_len(25), function(i) {
  readRDS(file.path(results_dir, paste0("mod", i, ".RDS")))
})

names(all_mods) <- paste0("mod", seq_len(25))

# Profile-CI versions used for manuscript uncertainty intervals.
mod1_prof <- readRDS(file.path(results_dir, "mod1_prof.RDS"))

mod7_prof_path <- file.path(results_dir, "mod7_prof.RDS")
if (!file.exists(mod7_prof_path)) {
  stop(
    "Missing profile-CI diagnostic model: ", mod7_prof_path,
    "\nRun the mod7 profiling step in Script 2 before Script 7."
  )
}
mod7_prof <- readRDS(mod7_prof_path)

# Separate four-region state-space movement model. This model has a different
# state space / likelihood and was independently reconstructed and validated
# during workflow auditing.
# It is stored separately from the seven-river candidate models and its AICc
# must not be compared with the seven-river candidate set.
regional_mod_path <- file.path(regional_results_dir, "regional_state_model.RDS")
if (!file.exists(regional_mod_path)) {
  stop("Missing separate regional model: ", regional_mod_path,
       "\nRun the regional-state model script before Script 7.")
}
regional_mod <- readRDS(regional_mod_path)

# Determine the top model directly from unrounded AICc values.
aicc_values <- vapply(all_mods, function(m) m$results$AICc, numeric(1))
top_idx      <- which.min(aicc_values)
top_mod      <- all_mods[[top_idx]]

message("Top seven-river candidate model: ", names(all_mods)[top_idx],
        " | ", top_mod$model.name,
        " | AICc = ", round(top_mod$results$AICc, 3))

# The current profile object must correspond to the top model. This catches a
# future renumbering/model-set change before tables are silently mislabeled.
if (top_idx != 1) {
  stop("The top model is no longer mod1, but this script expects mod1_prof.RDS ",
       "to be the profile-CI version of the top model. Reconcile model numbering first.")
}


if (abs(mod1_prof$results$AICc - top_mod$results$AICc) > 0.01 ||
    abs(mod1_prof$results$lnl - top_mod$results$lnl) > 0.01) {
  stop(
    "mod1_prof does not match the final mod1 candidate likelihood/AICc. ",
    "Re-run final Script 2 before Script 7."
  )
}

if ((!is.null(top_mod$results$singular) &&
     length(top_mod$results$singular) > 0) ||
    (!is.null(mod1_prof$results$singular) &&
     length(mod1_prof$results$singular) > 0)) {
  stop(
    "Final top model or mod1_prof reports singular parameter indices. ",
    "Inspect Script 2 before generating manuscript tables."
  )
}

# The profile-CI version of mod7 must reproduce the candidate-model likelihood
# and AICc. Profiling changes interval construction, not model selection.
mod7_candidate <- all_mods[["mod7"]]

if (abs(mod7_prof$results$AICc - mod7_candidate$results$AICc) > 0.01 ||
    abs(mod7_prof$results$lnl  - mod7_candidate$results$lnl)  > 0.01) {
  stop(
    "mod7_prof does not match candidate mod7 likelihood/AICc. ",
    "Re-run or reconcile the mod7 profiling step before generating Table A6."
  )
}

if (!is.null(mod7_prof$results$singular) &&
    length(mod7_prof$results$singular) > 0) {
  warning(
    "mod7_prof reports singular parameter index/indices: ",
    paste(mod7_prof$results$singular, collapse = ", "),
    ". Inspect before publication."
  )
}


# Production guardrails for the final terminal-p-constrained candidate set.
delta_full_check  <- aicc_values - min(aicc_values)
weight_full_check <- exp(-0.5 * delta_full_check) /
  sum(exp(-0.5 * delta_full_check))

expected_top5 <- c("mod1", "mod2", "mod3", "mod5", "mod6")
observed_top5 <- names(sort(aicc_values))[1:5]

if (!identical(observed_top5, expected_top5)) {
  stop(
    "Final candidate-model ordering changed. Expected top five: ",
    paste(expected_top5, collapse = ", "),
    "; observed: ",
    paste(observed_top5, collapse = ", "),
    "."
  )
}

if (abs(top_mod$results$AICc - 6133.295) > 0.05 ||
    abs(weight_full_check["mod1"] - 0.980688) > 0.005) {
  warning(
    "Final top-model AICc/weight differs from the validated run. ",
    "Observed AICc = ", round(top_mod$results$AICc, 3),
    "; weight = ", round(weight_full_check["mod1"], 6), "."
  )
}

regional_singular <- !is.null(regional_mod$results$singular) &&
  length(regional_mod$results$singular) > 0

if (regional_singular) {
  warning("regional_state_model.RDS reports singular parameters. Inspect before publication.")
}

# ---- load data ---------------------------------------------------------------

message("Loading adult encounter-history and tag-longevity data...")

at_inf_adults <- readRDS(
  file.path(data_dir, "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS")
)

tags <- readRDS(
  file.path(data_dir, "informed_at_longevity_all_gs_ms_export_20240121.RDS")
)

tagged_as_adults_no_2022 <- readRDS(
  file.path(data_dir, "tagged_as_adults_no_2022_deployments_20240121.RDS")
)

adult_first <- tagged_as_adults_no_2022 %>%
  transmute(FishID = as.numeric(FishID), first_code = as.character(first)) %>%
  distinct(FishID, .keep_all = TRUE)

# Adult-analysis population: remove fish first entering the adult analysis in
# 2022 because there is no subsequent annual survival interval.
adult_states <- names(river_map)
ch2022 <- paste0("000000000000", adult_states)

adult_analysis_fish <- at_inf_adults %>%
  filter(!ch %in% ch2022) %>%
  mutate(FishID = as.numeric(FishID))

if (nrow(adult_analysis_fish) != 1017 ||
    n_distinct(adult_analysis_fish$FishID) != 1017) {
  stop("Adult-analysis accounting changed. Expected 1,017 unique fish after removing 2022-only histories.")
}

# Identify the first adult-analysis state/year for every analysis fish. This is
# needed for Table A1 because 32 fish were initially tagged as subadults but
# later physically recaptured as adults and entered the adult analysis.
get_first_adult_entry <- function(ch) {

  x <- strsplit(ch, "", fixed = TRUE)[[1]]
  idx <- which(x %in% adult_states)

  if (length(idx) == 0) {
    return(c(first_code = NA_character_, adult_entry_year = NA_character_))
  }

  c(
    first_code = x[idx[1]],
    adult_entry_year = as.character(2009 + idx[1])
  )
}

adult_entry_mat <- t(vapply(
  adult_analysis_fish$ch,
  get_first_adult_entry,
  FUN.VALUE = c(first_code = "", adult_entry_year = "")
))

adult_entry_all <- adult_analysis_fish %>%
  transmute(FishID = as.numeric(FishID)) %>%
  bind_cols(as.data.frame(adult_entry_mat, stringsAsFactors = FALSE)) %>%
  mutate(
    adult_entry_year = as.integer(adult_entry_year),
    River = unname(river_map[first_code])
  )

if (any(is.na(adult_entry_all$adult_entry_year)) ||
    any(is.na(adult_entry_all$River))) {
  stop("At least one adult-analysis fish could not be assigned a first adult year/river.")
}

message("Adult-analysis audit: 1,017 unique fish retained; ",
        nrow(at_inf_adults) - nrow(adult_analysis_fish),
        " 2022-only histories excluded.")

# ==============================================================================
# TABLE 1 -- adult acoustic-tag deployments by river and deployment year
# ==============================================================================

message("Building Table 1: adult tag deployments by river and year...")

tags_for_analysis <- tags %>%
  mutate(FishID = as.numeric(FishID)) %>%
  filter(FishID %in% adult_first$FishID,
         start_year <= 2021) %>%
  left_join(adult_first, by = "FishID") %>%
  mutate(River = unname(river_map[first_code]))

if (any(is.na(tags_for_analysis$River))) {
  stop("At least one analysis tag could not be assigned to a river in Table 1.")
}

table1 <- tags_for_analysis %>%
  count(River, start_year, name = "N") %>%
  complete(River = river_order,
           start_year = seq(min(start_year, na.rm = TRUE), 2021),
           fill = list(N = 0)) %>%
  mutate(River = factor(River, levels = river_order)) %>%
  arrange(River, start_year) %>%
  pivot_wider(names_from = start_year, values_from = N) %>%
  mutate(Total = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  arrange(River)

write.csv(
  table1,
  file.path(out_dir, "Table1_adult_deployments_by_river_year.csv"),
  row.names = FALSE,
  na = ""
)

message("  Table 1 saved. Deployment years represented: ",
        min(tags_for_analysis$start_year, na.rm = TRUE), "-2021")
message("  Table 1 audit: ", nrow(tags_for_analysis),
        " pre-2022 deployments on ", n_distinct(tags_for_analysis$FishID),
        " fish tagged as adults.")

if (nrow(tags_for_analysis) != 1015 ||
    n_distinct(tags_for_analysis$FishID) != 985) {
  warning("Table 1 accounting differs from the audited values of 1,015 deployments on 985 adult-tagged fish.")
}

# ==============================================================================
# TABLE 2 -- a priori survival structures used in the candidate model set
# ==============================================================================
#
# These are the eight SURVIVAL structures. Each structure was crossed with
# three detection structures:
#   p(river:time), p(river), and p(.)
# while retaining river-to-river Psi. For p(river:time), detection in the
# final occasion is constrained equal to detection in the preceding occasion
# within each river to resolve terminal survival-detection confounding.
#
# That produces 8 x 3 = 24 river-state models. A separate fully constant null
# model, S(.)p(.)Psi(.), was also fit, for 25 candidate models total.
#
# Table 2 is therefore a hypothesis/notation table, not a table extracted from
# MARK output.

message("Building Table 2: candidate survival structures...")

table2 <- tibble(
  Structure_ID = 1:8,
  Structure_Name = c(
    "Range-wide constant",
    "Region",
    "River",
    "Time",
    "Region + time",
    "River + time",
    "Region x time",
    "River x time"
  ),
  Model_Hypothesis = c(
    "Adult survival is constant across space and time and is represented by a single range-wide rate.",
    "Adult survival differs among regional groups, but survival is constant through time within each region.",
    "Adult survival differs among rivers, but survival is constant through time within each river.",
    "Adult survival is the same among rivers but varies annually across the study period.",
    "Adult survival differs among regions, with a common annual temporal pattern shared among regions.",
    "Adult survival differs among rivers, with a common annual temporal pattern shared among rivers.",
    "Adult survival differs among regions and the annual temporal pattern also differs among regions, yielding a separate estimate for each region-year combination.",
    "Adult survival differs among rivers and the annual temporal pattern also differs among rivers, yielding a separate estimate for each river-year combination."
  ),
  S_Notation = c(
    "S(.)",
    "S(region)",
    "S(river)",
    "S(time)",
    "S(region+time)",
    "S(river+time)",
    "S(region:time)",
    "S(river:time)"
  ),
  Detection_Structures = rep(
    "p(river:time), p(river), p(.)",
    8
  ),
  Detection_Note = rep(
    "For p(river:time), final-occasion detection was constrained equal to the preceding occasion within river.",
    8
  ),
  Psi_Structure = rep(
    "Psi(stratum:tostratum) with geographically impossible transitions fixed to 0",
    8
  )
)

write.csv(
  table2,
  file.path(out_dir, "Table2_survival_model_structures.csv"),
  row.names = FALSE
)

message("  Table 2 saved.")
message("  Candidate-set audit: 8 survival structures x 3 p structures = 24,")
message("  plus S(.)p(.)Psi(.) null = 25 total models.")

# ==============================================================================
# TABLE 3 -- AICc model selection for 25 seven-river candidate models
# ==============================================================================

message("Building Table 3: AICc model selection...")

delta_full  <- aicc_values - min(aicc_values)
weight_full <- exp(-0.5 * delta_full) / sum(exp(-0.5 * delta_full))

table3 <- tibble(
  Model_ID  = names(all_mods),
  Model     = vapply(all_mods, function(m) m$model.name, character(1)),
  K         = vapply(all_mods, function(m) m$results$npar, numeric(1)),
  AICc      = round(aicc_values, 2),
  DeltaAICc = round(delta_full, 2),
  Weight    = round(weight_full, 3),
  NegLogLik = round(vapply(all_mods, function(m) m$results$lnl, numeric(1)), 2)
) %>%
  arrange(AICc)

write.csv(
  table3,
  file.path(out_dir, "Table3_AICc_model_selection.csv"),
  row.names = FALSE
)

message("  Table 3 saved.")
message(
  "  Final model support: mod1 weight = ",
  sprintf("%.6f", weight_full["mod1"]),
  "; mod2 Delta AICc = ",
  sprintf("%.3f", delta_full["mod2"]),
  "; mod3 Delta AICc = ",
  sprintf("%.3f", delta_full["mod3"])
)

# ==============================================================================
# TABLE A1 -- tag longevity for transmitters relevant to the adult analysis
# ==============================================================================

message("Building Table A1: tag longevity for transmitters relevant to the adult analysis...")

# Table 1 describes deployments on fish tagged as adults (985 fish, 1,015
# pre-2022 deployments). Table A1 has a different purpose: document the
# transmitter-longevity information that actually contributed active coverage
# while fish were in the adult multistate analysis.
#
# A transmitter is relevant when:
#   1. it belongs to one of the 1,017 adult-analysis fish;
#   2. it was deployed no later than 2021; and
#   3. its expected active life extends through or beyond that fish's first
#      adult-analysis year.
#
# This intentionally includes tags implanted while a fish was still a subadult
# if the transmitter remained active when that fish later entered the adult
# analysis. It excludes older tags that expired before adult entry and new 2022
# deployments that cannot inform a subsequent annual survival interval.

tags_A1 <- tags %>%
  mutate(FishID = as.numeric(FishID)) %>%
  filter(FishID %in% adult_entry_all$FishID) %>%
  left_join(
    adult_entry_all %>%
      select(FishID, adult_entry_year, first_code, River),
    by = "FishID"
  ) %>%
  filter(
    start_year <= 2021,
    end_year >= adult_entry_year
  ) %>%
  mutate(
    tag_life_years = if ("tag_life_years" %in% names(.)) {
      as.numeric(tag_life_years)
    } else {
      ifelse(
        is.na(Mark_Longevity_Days),
        5,
        floor(as.numeric(Mark_Longevity_Days) / 365)
      )
    },

    # Supporting audit quantity: years of transmitter coverage overlapping the
    # period from adult entry through the end of the 2010-2022 study window.
    adult_coverage_start_year = pmax(start_year, adult_entry_year, 2010),
    adult_coverage_end_year   = pmin(end_year, 2022),
    adult_coverage_years = pmax(
      1,
      adult_coverage_end_year - adult_coverage_start_year
    )
  )

if (any(is.na(tags_A1$River))) {
  stop("At least one Table A1 transmitter could not be assigned to the fish's first adult-analysis river.")
}

# Audited values from the full workflow:
#   1,048 relevant transmitter deployments
#   1,017 unique adult-analysis fish
if (nrow(tags_A1) != 1048 ||
    n_distinct(tags_A1$FishID) != 1017) {
  warning(
    "Table A1 accounting differs from the audited values of ",
    "1,048 relevant tags across 1,017 adult-analysis fish. ",
    "Observed: ", nrow(tags_A1), " tags across ",
    n_distinct(tags_A1$FishID), " fish."
  )
}

# Main A1: counts of relevant tags by assigned manufacturer/assumed longevity,
# with simple summary columns.
A1_counts <- tags_A1 %>%
  count(River, tag_life_years, name = "N") %>%
  mutate(tag_life_years = paste0("Life_", tag_life_years, "yr")) %>%
  pivot_wider(names_from = tag_life_years, values_from = N, values_fill = 0)

A1_summary <- tags_A1 %>%
  group_by(River) %>%
  summarise(
    N_Tags = n(),
    N_Fish = n_distinct(FishID),
    Mean_Tag_Life_Years = mean(tag_life_years, na.rm = TRUE),
    Median_Tag_Life_Years = median(tag_life_years, na.rm = TRUE),
    Min_Tag_Life_Years = min(tag_life_years, na.rm = TRUE),
    Max_Tag_Life_Years = max(tag_life_years, na.rm = TRUE),
    .groups = "drop"
  )

tableA1 <- A1_counts %>%
  full_join(A1_summary, by = "River") %>%
  mutate(River = factor(River, levels = river_order)) %>%
  arrange(River) %>%
  mutate(
    Mean_Tag_Life_Years = round(Mean_Tag_Life_Years, 2),
    Median_Tag_Life_Years = round(Median_Tag_Life_Years, 2)
  )

# Put longevity classes in chronological order regardless of the order in
# which they happen to occur in the data.
life_cols <- grep("^Life_[0-9]+yr$", names(tableA1), value = TRUE)
life_cols <- life_cols[
  order(as.numeric(str_extract(life_cols, "[0-9]+")))
]

tableA1 <- tableA1 %>%
  select(
    River,
    all_of(life_cols),
    N_Tags,
    N_Fish,
    Mean_Tag_Life_Years,
    Median_Tag_Life_Years,
    Min_Tag_Life_Years,
    Max_Tag_Life_Years
  )

message("  Table A1 audit: ", sum(tableA1$N_Tags),
        " relevant transmitter deployments across ",
        n_distinct(tags_A1$FishID),
        " adult-analysis fish.")

# Additional audit: distinguish the 985 fish initially tagged as adults from
# the 32 fish initially tagged as subadults but later entering the adult analysis.
adult_tagged_ids <- unique(as.numeric(tagged_as_adults_no_2022$FishID))

A1_group_audit <- tags_A1 %>%
  mutate(
    Fish_Group = if_else(
      FishID %in% adult_tagged_ids,
      "Initially tagged as adult",
      "Initially tagged as subadult; later entered adult analysis"
    )
  ) %>%
  group_by(Fish_Group) %>%
  summarise(
    N_Tags = n(),
    N_Fish = n_distinct(FishID),
    .groups = "drop"
  )

message("  Table A1 population breakdown:")
print(A1_group_audit)

write.csv(
  tableA1,
  file.path(out_dir, "TableA1_tag_longevity_by_river.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  A1_group_audit,
  file.path(out_dir, "Supporting_A1_population_breakdown.csv"),
  row.names = FALSE
)

# Supporting audit: effective transmitter coverage during the adult-analysis
# portion of each fish's history.
A1_adult_coverage <- tags_A1 %>%
  count(River, adult_coverage_years, name = "N") %>%
  arrange(factor(River, levels = river_order), adult_coverage_years)

write.csv(
  A1_adult_coverage,
  file.path(out_dir, "Supporting_A1_effective_adult_coverage_by_river.csv"),
  row.names = FALSE
)

message("  Table A1 saved.")

# ==============================================================================
# TABLE A2 -- detection probability from TOP MODEL p(river:time), 2011-2022
# ==============================================================================

message("Building Table A2: detection from top model...")

# With the terminal-p identifying constraint, the final two detection occasions
# share the same design-matrix parameter within each river. RMark therefore may
# return only the 77 unique real p rows (7 rivers x 11 estimable p-time levels)
# rather than 84 occasion-specific rows. For manuscript reporting, expand the
# pooled terminal estimate back to both 2021 and 2022. The duplicated 2022 row
# is not an additional estimate; it is the same constrained estimate.

p_top_raw <- get_real(mod1_prof, "p") %>%
  parse_stratum_time(
    parameter = "p",
    state_map = river_map,
    first_year = 2011,
    state_col = "River"
  ) %>%
  transmute(
    River,
    Year,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

message(
  "  Raw top-model p output: ",
  nrow(p_top_raw), " rows; ",
  n_distinct(p_top_raw$River), " rivers; ",
  n_distinct(p_top_raw$Year), " unique p-time levels."
)

if (n_distinct(p_top_raw$River) != 7) {
  stop(
    "Unexpected number of rivers in top-model p output. Expected 7; observed ",
    n_distinct(p_top_raw$River), "."
  )
}

if (nrow(p_top_raw) == 77 &&
    n_distinct(p_top_raw$Year) == 11) {

  # Under the final constraint, the last unique p-time level represents both
  # the penultimate and final detection occasions. With first_year = 2011,
  # the 11 unique levels map to 2011-2021; duplicate 2021 as 2022.
  terminal_p_report <- p_top_raw %>%
    filter(Year == max(Year)) %>%
    mutate(Year = 2022)

  if (nrow(terminal_p_report) != 7) {
    stop(
      "Could not reconstruct the pooled terminal detection occasion. ",
      "Expected 7 terminal river rows; observed ", nrow(terminal_p_report), "."
    )
  }

  p_top <- bind_rows(
    p_top_raw,
    terminal_p_report
  )

  message(
    "  Expanded pooled terminal p estimate to 2022: ",
    "77 unique real rows -> 84 river-year reporting rows."
  )

} else if (nrow(p_top_raw) == 84 &&
           n_distinct(p_top_raw$Year) == 12) {

  # Keep this branch for robustness in case a future RMark version returns one
  # real-output row for each occasion even when the final two share a parameter.
  p_top <- p_top_raw

  message(
    "  RMark returned 84 occasion-specific p rows directly; no expansion needed."
  )

} else {

  stop(
    "Unexpected p(river:time) output after terminal pooling. Observed ",
    nrow(p_top_raw), " rows, ",
    n_distinct(p_top_raw$River), " rivers, and ",
    n_distinct(p_top_raw$Year), " unique p-time levels. ",
    "Expected either 77 rows / 11 levels or 84 rows / 12 levels."
  )
}

p_top <- p_top %>%
  mutate(River = factor(River, levels = river_order)) %>%
  arrange(River, Year)

p_years <- sort(unique(as.integer(p_top$Year)))

if (nrow(p_top) != 84 ||
    n_distinct(p_top$River) != 7 ||
    n_distinct(p_top$Year) != 12 ||
    !identical(p_years, 2011:2022)) {
  message("  p reporting rows: ", nrow(p_top))
  message("  p reporting rivers: ", n_distinct(p_top$River))
  message("  p reporting years: ", paste(p_years, collapse = ", "))
  stop(
    "Expanded top-model p table failed final dimensions/year audit. ",
    "Expected 84 rows = 7 rivers x 12 years, 2011-2022."
  )
}

# Verify the identifying constraint on the reporting scale: p_2022 = p_2021
# within each river, including SE and interval limits.
p_terminal_check <- p_top %>%
  filter(Year %in% c(2021, 2022)) %>%
  select(River, Year, Estimate, SE, LCL_95, UCL_95) %>%
  pivot_wider(
    names_from = Year,
    values_from = c(Estimate, SE, LCL_95, UCL_95)
  )

terminal_diffs <- c(
  abs(p_terminal_check$Estimate_2022 - p_terminal_check$Estimate_2021),
  abs(p_terminal_check$SE_2022       - p_terminal_check$SE_2021),
  abs(p_terminal_check$LCL_95_2022   - p_terminal_check$LCL_95_2021),
  abs(p_terminal_check$UCL_95_2022   - p_terminal_check$UCL_95_2021)
)

if (nrow(p_terminal_check) != 7 ||
    any(terminal_diffs > 1e-8, na.rm = TRUE)) {
  stop(
    "Terminal p reporting rows do not reproduce the intended constraint ",
    "p_2022 = p_2021 within river."
  )
}

# Full 84-row output for reproducibility / supporting appendix detail.
write.csv(
  p_top %>% mutate(across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 4))),
  file.path(out_dir, "Supporting_detection_by_river_year_top_model.csv"),
  row.names = FALSE
)

# Compact manuscript Table A2: summaries of the annual estimates from the top
# model, not estimates from an unsupported p(river) model.
tableA2 <- p_top %>%
  group_by(River) %>%
  summarise(
    Mean_p = mean(Estimate, na.rm = TRUE),
    Median_p = median(Estimate, na.rm = TRUE),
    Min_p = min(Estimate, na.rm = TRUE),
    Min_Year = Year[which.min(Estimate)],
    Max_p = max(Estimate, na.rm = TRUE),
    Max_Year = Year[which.max(Estimate)],
    .groups = "drop"
  ) %>%
  mutate(
    across(c(Mean_p, Median_p, Min_p, Max_p), ~round(.x, 3)),
    River = factor(River, levels = river_order)
  ) %>%
  arrange(River)

write.csv(
  tableA2,
  file.path(out_dir, "TableA2_detection_summary_top_model.csv"),
  row.names = FALSE
)

message("  Table A2 saved (+ full river x year supporting file).")

# ==============================================================================
# TABLE A3 -- river x year apparent survival from alternative S(river+time) model
# ==============================================================================
#
# The final candidate set strongly supports mod1. mod3 is retained here only to
# provide transparent river-specific descriptive estimates in the appendix.
# It should not be described as a competitive or equally supported model.

message("Building Table A3: descriptive river x year apparent survival from mod3...")

mod3 <- all_mods[["mod3"]]

s_river <- get_real(mod3, "S") %>%
  parse_stratum_time(
    parameter = "S",
    state_map = river_map,
    first_year = 2010,
    state_col = "River"
  ) %>%
  transmute(
    River,
    Year,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

if (nrow(s_river) != 84 ||
    n_distinct(s_river$River) != 7 ||
    n_distinct(s_river$Year) != 12) {
  stop(
    "Unexpected S(river+time) dimensions. Expected 84 rows = ",
    "7 rivers x 12 years."
  )
}

tableA3 <- s_river %>%
  mutate(
    River = factor(River, levels = river_order),
    across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 3))
  ) %>%
  arrange(River, Year)

write.csv(
  tableA3,
  file.path(out_dir, "TableA3_river_year_apparent_survival.csv"),
  row.names = FALSE
)

write.csv(
  tableA3,
  file.path(out_dir, "Supporting_river_year_apparent_survival_mod3.csv"),
  row.names = FALSE
)

message(
  "  Table A3 saved. mod3 Delta AICc = ",
  sprintf("%.2f", delta_full["mod3"]),
  "; weight = ",
  sprintf("%.4f", weight_full["mod3"]),
  ". Interpret descriptively, not as a competitive model."
)

# ==============================================================================
# TABLE A4 -- river-to-river movement + river fidelity from TOP MODEL
# ==============================================================================

message("Building Table A4: river movement and fidelity from top model...")

psi_river <- get_real(top_mod, "Psi") %>%
  # MARK returns one representative real-parameter row for transitions that
  # were structurally fixed to zero. Do not present that row as if it were an
  # estimated movement route.
  filter(!("fixed" %in% names(.)) | is.na(fixed) | fixed != "Fixed") %>%
  parse_psi(state_map = river_map, from_col = "From_River", to_col = "To_River") %>%
  transmute(
    From_River,
    To_River,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

# Long output retains uncertainty for every estimated transition.
tableA4_long <- psi_river %>%
  mutate(across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 4))) %>%
  arrange(factor(From_River, levels = river_order),
          factor(To_River, levels = river_order))

write.csv(
  tableA4_long,
  file.path(out_dir, "TableA4_river_movement_estimates_with_CI.csv"),
  row.names = FALSE
)

# Fidelity is the probability of remaining in the source river conditional on
# survival: 1 - sum of estimated transitions to other rivers.
river_fidelity <- psi_river %>%
  group_by(From_River) %>%
  summarise(
    Probability_Leave = sum(Estimate, na.rm = TRUE),
    Fidelity = 1 - Probability_Leave,
    .groups = "drop"
  )

# Wide matrix for manuscript presentation. Impossible/fixed transitions appear
# blank; Fidelity is appended as the final column.
movement_matrix <- matrix(
  NA_real_,
  nrow = length(river_order),
  ncol = length(river_order),
  dimnames = list(river_order, river_order)
)

for (i in seq_len(nrow(psi_river))) {
  movement_matrix[psi_river$From_River[i], psi_river$To_River[i]] <-
    psi_river$Estimate[i]
}

diag(movement_matrix) <- NA_real_

tableA4_wide <- data.frame(
  From_River = rownames(movement_matrix),
  round(movement_matrix, 3),
  check.names = FALSE
) %>%
  left_join(
    river_fidelity %>%
      transmute(From_River,
                Probability_Leave = round(Probability_Leave, 3),
                Fidelity = round(Fidelity, 3)),
    by = "From_River"
  )

write.csv(
  tableA4_wide,
  file.path(out_dir, "TableA4_river_movement_matrix_and_fidelity.csv"),
  row.names = FALSE,
  na = ""
)

message("  Table A4 saved (wide matrix + long CI file).")

# ==============================================================================
# TABLE A5 -- separate four-region multistate movement model
# ==============================================================================

message("Building Table A5: regional movement and fidelity...")

psi_region <- get_real(regional_mod, "Psi") %>%
  # As above, omit representative rows for transitions fixed to zero.
  filter(!("fixed" %in% names(.)) | is.na(fixed) | fixed != "Fixed") %>%
  parse_psi(state_map = regional_state_map,
            from_col = "From_Region", to_col = "To_Region") %>%
  transmute(
    From_Region,
    To_Region,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

regional_fidelity <- psi_region %>%
  group_by(From_Region) %>%
  summarise(
    Probability_Leave = sum(Estimate, na.rm = TRUE),
    Fidelity = 1 - Probability_Leave,
    .groups = "drop"
  )

region_matrix <- matrix(
  NA_real_,
  nrow = length(region_order),
  ncol = length(region_order),
  dimnames = list(region_order, region_order)
)

for (i in seq_len(nrow(psi_region))) {
  region_matrix[psi_region$From_Region[i], psi_region$To_Region[i]] <-
    psi_region$Estimate[i]
}

diag(region_matrix) <- NA_real_

tableA5 <- data.frame(
  From_Region = rownames(region_matrix),
  round(region_matrix, 3),
  check.names = FALSE
) %>%
  left_join(
    regional_fidelity %>%
      transmute(From_Region,
                Probability_Leave = round(Probability_Leave, 3),
                Fidelity = round(Fidelity, 3)),
    by = "From_Region"
  )

write.csv(
  tableA5,
  file.path(out_dir, "TableA5_regional_movement_matrix_and_fidelity.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  psi_region %>%
    mutate(across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 4))) %>%
    arrange(factor(From_Region, levels = region_order),
            factor(To_Region, levels = region_order)),
  file.path(out_dir, "Supporting_regional_movement_estimates_with_CI.csv"),
  row.names = FALSE
)

message("  Table A5 saved.")
message("  Regional model diagnostics: K = ", regional_mod$results$npar,
        "; AICc = ", round(regional_mod$results$AICc, 3),
        "; singular = ", regional_singular)
message("  NOTE: regional-model AICc is recorded only; do NOT compare it with Table 3.")

# ==============================================================================
# TABLE A6 -- river x year survival from interactive diagnostic model
# ==============================================================================
#
# mod7 = S(river:time)p(river)Psi(strattostrat)
#
# This model is retained for the specific event-year diagnostic used to inspect
# the 2018 Apalachicola estimate associated with Hurricane Michael. It is not a
# competitive model in the final candidate set and should not be used as the
# primary basis for annual river-specific inference.

message("Building Table A6: river x year survival from mod7 diagnostic model with profile CIs...")

# Keep ordinary candidate mod7 for model-selection quantities, but use the
# profile-CI refit for real-parameter uncertainty intervals in Table A6.
mod7 <- mod7_candidate

s_mod7 <- get_real(mod7_prof, "S") %>%
  parse_stratum_time(
    parameter = "S",
    state_map = river_map,
    first_year = 2010,
    state_col = "River"
  ) %>%
  transmute(
    River,
    Year,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

if (nrow(s_mod7) != 84 ||
    n_distinct(s_mod7$River) != 7 ||
    n_distinct(s_mod7$Year) != 12) {
  stop(
    "Unexpected mod7 S dimensions. Expected 84 rows = 7 rivers x 12 years."
  )
}

# Audit upper-boundary estimates. Profile intervals may be informative even
# when the maximum-likelihood estimate is exactly 1.0.
n_boundary_mod7 <- sum(s_mod7$Estimate >= 0.999, na.rm = TRUE)
n_boundary_profile_informative <- sum(
  s_mod7$Estimate >= 0.999 & s_mod7$LCL_95 < 0.999,
  na.rm = TRUE
)

message(
  "  mod7 profile-CI boundary audit: ",
  n_boundary_mod7, " estimates >= 0.999; ",
  n_boundary_profile_informative,
  " have profile LCL < 0.999."
)

tableA6 <- s_mod7 %>%
  mutate(
    River = factor(River, levels = river_order),
    across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 3))
  ) %>%
  arrange(River, Year)

write.csv(
  tableA6,
  file.path(out_dir, "TableA6_river_year_apparent_survival_diagnostic.csv"),
  row.names = FALSE
)

a6_apal_2018 <- s_mod7 %>%
  filter(River == "Apalachicola", Year == 2018)

message(
  "  Table A6 saved. mod7 Delta AICc = ",
  sprintf("%.2f", delta_full["mod7"]),
  "; weight = ",
  format(weight_full["mod7"], scientific = TRUE, digits = 3)
)

if (nrow(a6_apal_2018) == 1) {
  message(
    "  Apalachicola 2018 diagnostic S = ",
    sprintf("%.3f", a6_apal_2018$Estimate),
    " (95% profile CI ",
    sprintf("%.3f", a6_apal_2018$LCL_95),
    "-",
    sprintf("%.3f", a6_apal_2018$UCL_95),
    ")"
  )
}

# ==============================================================================
# SUPPORTING RESULTS -- top-model regional survival for Results text / Figure 2
# ==============================================================================

message("Building supporting top-model regional survival summaries...")

s_region_top <- get_real(mod1_prof, "S") %>%
  parse_stratum_time(
    parameter = "S",
    state_map = regional_survival_map,
    first_year = 2010,
    state_col = "Region"
  ) %>%
  transmute(
    Region,
    Year,
    Estimate = estimate,
    SE = se,
    LCL_95 = lcl,
    UCL_95 = ucl
  )

if (nrow(s_region_top) != 48 || n_distinct(s_region_top$Region) != 4 ||
    n_distinct(s_region_top$Year) != 12) {
  stop("Unexpected top-model regional survival dimensions. Expected 48 rows = 4 regions x 12 years.")
}

s_region_top <- s_region_top %>%
  mutate(Region = factor(Region, levels = region_order)) %>%
  arrange(Region, Year)

write.csv(
  s_region_top %>% mutate(across(c(Estimate, SE, LCL_95, UCL_95), ~round(.x, 4))),
  file.path(out_dir, "Supporting_regional_survival_by_year_top_model.csv"),
  row.names = FALSE
)

regional_survival_summary <- s_region_top %>%
  group_by(Region) %>%
  summarise(
    Mean = mean(Estimate, na.rm = TRUE),
    Median = median(Estimate, na.rm = TRUE),
    Min = min(Estimate, na.rm = TRUE),
    Min_Year = Year[which.min(Estimate)],
    Max = max(Estimate, na.rm = TRUE),
    Max_Year = Year[which.max(Estimate)],
    .groups = "drop"
  ) %>%
  mutate(across(c(Mean, Median, Min, Max), ~round(.x, 4)))

write.csv(
  regional_survival_summary,
  file.path(out_dir, "Supporting_regional_survival_summary_top_model.csv"),
  row.names = FALSE
)

# Additional concise fidelity summaries for Results drafting.
write.csv(
  river_fidelity %>%
    mutate(across(c(Probability_Leave, Fidelity), ~round(.x, 4))) %>%
    arrange(factor(From_River, levels = river_order)),
  file.path(out_dir, "Supporting_river_fidelity_top_model.csv"),
  row.names = FALSE
)

write.csv(
  regional_fidelity %>%
    mutate(across(c(Probability_Leave, Fidelity), ~round(.x, 4))) %>%
    arrange(factor(From_Region, levels = region_order)),
  file.path(out_dir, "Supporting_regional_fidelity.csv"),
  row.names = FALSE
)

message("  Supporting Results files saved.")

# ==============================================================================
# FINAL AUDIT SUMMARY
# ==============================================================================

message("\n================ SCRIPT 7 FINAL AUDIT ================")
message("Top river-state model: ", top_mod$model.name)
message("Top-model AICc: ", round(top_mod$results$AICc, 3))
message("Top-model weight: ", sprintf("%.6f", weight_full["mod1"]))
message("Candidate-model directory: ", results_dir)
message("Regional-state model directory: ", regional_results_dir)
message("Top-model Psi rows: ", nrow(psi_river))
message("Top-model unique p rows from RMark: ", nrow(p_top_raw),
        " (expected 77 under terminal pooling)")
message("Top-model p reporting rows after terminal expansion: ",
        nrow(p_top), " (expected 84)")
message("Alternative mod3 river S rows: ", nrow(s_river), " (expected 84)")
message("Top regional S rows: ", nrow(s_region_top), " (expected 48)")
message("Diagnostic mod7 profile S rows: ", nrow(s_mod7), " (expected 84)")
message("Diagnostic mod7 boundary estimates >= 0.999: ", n_boundary_mod7)
message("Diagnostic mod7 boundary rows with profile LCL < 0.999: ",
        n_boundary_profile_informative)
message("Table 1 deployments/fish: ", nrow(tags_for_analysis), " / ",
        n_distinct(tags_for_analysis$FishID), " (expected 1015 / 985)")
message("Table A1 relevant tags/fish: ", nrow(tags_A1), " / ",
        n_distinct(tags_A1$FishID), " (expected 1048 / 1017)")
message("Regional-state model singular: ", regional_singular)
message("Regional-state model K: ", regional_mod$results$npar)
message("Regional-state model AICc: ", round(regional_mod$results$AICc, 3),
        " [record only; different likelihood]")

files_written <- list.files(out_dir, pattern = "\\.csv$", full.names = FALSE)
message("\nCSV files currently in manuscript_tables/:")
for (f in sort(files_written)) message("  ", f)

message("\nScript 7 FINAL complete.")
