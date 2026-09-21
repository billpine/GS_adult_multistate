################################################################################
# Gulf Sturgeon CJFAS Revision 1
# Script 4: Annual adult tag-pool counts by river
# production version
#
# Purpose
#   Build annual counts of adult-analysis fish with active transmitter coverage
#   under two river-assignment rules:
#
#     1. First-observed-river rule
#        A fish remains assigned to its first observed adult river for all
#        subsequent active annual occasions.
#
#     2. Most-recent-river rule
#        A fish's assignment updates whenever it is subsequently observed in a
#        different river and remains there until another observed river state.
#
#   These are assignment rules for active tag pools. "First observed river" is
#   NOT assumed to be the fish's natal river.
#
#   The most-recent-river pool is used by Script 8 for Figure A2.
#
# Outputs
#   data/derived/tableA_pool_firstrule.RDS
#   data/derived/tableA_pool_recentrule.RDS
#
# Important
#   - "." denotes an unavailable occasion because the transmitter is not active.
#   - "0" denotes an active occasion with no detection.
#   - A fish is not assigned to an adult river until its first observed adult
#     river state. This preserves the original analytical rule and prevents
#     pre-adult/pre-state occasions from being assigned retrospectively.
#   - 2022-only adult histories are removed, matching Script 2.
#

################################################################################

# ---- project setup -----------------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
})

pool_first_path <- file.path(
  data_dir,
  "tableA_pool_firstrule.RDS"
)

pool_recent_path <- file.path(
  data_dir,
  "tableA_pool_recentrule.RDS"
)

# ---- 1. Load adult encounter histories ---------------------------------------

message("\n1. Loading adult encounter histories...")

adult_path <- file.path(
  data_dir,
  "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
)

informed_at <- readRDS(adult_path)
informed_at$FishID <- as.numeric(informed_at$FishID)

message("  Adult histories loaded: ", nrow(informed_at),
        " [expected 1070]")

if (nrow(informed_at) != 1070 ||
    dplyr::n_distinct(informed_at$FishID) != 1070) {
  stop(
    "Expected 1,070 adult encounter-history rows / unique FishID before ",
    "removing 2022-only histories."
  )
}

# Remove fish whose only adult encounter is in the final 2022 occasion.
possible_states <- c("A", "B", "C", "E", "F", "G", "H")
ch2022 <- paste0("000000000000", possible_states)

keep <- !informed_at$ch %in% ch2022
n_removed_2022 <- sum(!keep)

informed_at <- informed_at[
  keep,
  ,
  drop = FALSE
]

message("  2022-only histories removed: ", n_removed_2022,
        " [expected 53]")
message("  Fish retained for tag-pool accounting: ", nrow(informed_at),
        " [expected 1017]")

if (n_removed_2022 != 53 ||
    nrow(informed_at) != 1017 ||
    dplyr::n_distinct(informed_at$FishID) != 1017) {
  stop(
    "Tag-pool population changed. Expected 53 histories removed and ",
    "1,017 unique fish retained."
  )
}

years <- 2010:2022
n_occasions <- length(years)

# Each adult encounter history should have 13 annual occasions.
ch_lengths <- nchar(as.character(informed_at$ch))

if (any(ch_lengths != n_occasions)) {
  stop(
    "Encounter-history length mismatch. Expected ",
    n_occasions,
    " annual occasions (2010-2022)."
  )
}

# ---- 2. Helper to parse one encounter history --------------------------------

parse_history <- function(ch) {

  states <- strsplit(
    as.character(ch),
    split = "",
    fixed = TRUE
  )[[1]]

  first_positions <- which(states %in% possible_states)

  if (length(first_positions) == 0) {
    stop("Encounter history contains no observed adult river state: ", ch)
  }

  list(
    states = states,
    first_pos = min(first_positions)
  )
}

# ---- 3. First-observed-river pool --------------------------------------------

message("\n2. Building first-observed-river tag pool...")

pool_first <- matrix(
  0L,
  nrow = length(possible_states),
  ncol = n_occasions,
  dimnames = list(
    possible_states,
    as.character(years)
  )
)

for (i in seq_len(nrow(informed_at))) {

  parsed <- parse_history(informed_at$ch[i])
  b <- parsed$states
  first_pos <- parsed$first_pos
  first_river <- b[first_pos]

  # Count only active occasions at or after the first observed adult river.
  # River letters and "0" indicate active transmitter coverage; "." is
  # unavailable because the tag is not active.
  active_cols <- which(
    b %in% c(possible_states, "0")
  )

  active_cols <- active_cols[
    active_cols >= first_pos
  ]

  pool_first[
    first_river,
    active_cols
  ] <- pool_first[
    first_river,
    active_cols
  ] + 1L
}

# ---- 4. Most-recent-river pool -----------------------------------------------

message("\n3. Building most-recent-river tag pool...")

pool_recent <- matrix(
  0L,
  nrow = length(possible_states),
  ncol = n_occasions,
  dimnames = list(
    possible_states,
    as.character(years)
  )
)

for (i in seq_len(nrow(informed_at))) {

  parsed <- parse_history(informed_at$ch[i])
  b <- parsed$states
  first_pos <- parsed$first_pos

  current_state <- NA_character_

  for (j in seq.int(first_pos, length(b))) {

    # "." = transmitter unavailable, so this fish is not in the active pool.
    if (b[j] == ".") {
      next
    }

    # Update state when a river is observed.
    if (b[j] %in% possible_states) {
      current_state <- b[j]
    }

    # At and after first_pos, current_state should always be known for every
    # active occasion.
    if (is.na(current_state)) {
      stop(
        "Most-recent-river assignment failed for FishID ",
        informed_at$FishID[i],
        " at occasion ", j, "."
      )
    }

    pool_recent[
      current_state,
      j
    ] <- pool_recent[
      current_state,
      j
    ] + 1L
  }
}

# ---- 5. Order and label river rows -------------------------------------------

river_map <- c(
  F = "Pearl",
  E = "Pascagoula",
  C = "Escambia",
  H = "Yellow",
  B = "Choctawhatchee",
  A = "Apalachicola",
  G = "Suwannee"
)

river_codes_ordered <- names(river_map)

pool_first_ordered <- pool_first[
  river_codes_ordered,
  ,
  drop = FALSE
]

pool_recent_ordered <- pool_recent[
  river_codes_ordered,
  ,
  drop = FALSE
]

rownames(pool_first_ordered)  <- unname(river_map)
rownames(pool_recent_ordered) <- unname(river_map)

# ---- 6. Internal consistency checks ------------------------------------------

message("\n4. Auditing tag-pool matrices...")

if (!identical(dim(pool_first_ordered), c(7L, 13L)) ||
    !identical(dim(pool_recent_ordered), c(7L, 13L))) {
  stop("Expected both tag-pool matrices to be 7 rivers x 13 years.")
}

if (!identical(
  colnames(pool_first_ordered),
  as.character(2010:2022)
) ||
    !identical(
      colnames(pool_recent_ordered),
      as.character(2010:2022)
    )) {
  stop("Unexpected year columns in tag-pool matrices.")
}

# The assignment rule can change the river-specific counts but cannot change
# how many fish are active in a year. Annual totals must therefore match.
first_totals  <- colSums(pool_first_ordered)
recent_totals <- colSums(pool_recent_ordered)

if (!identical(
  as.integer(first_totals),
  as.integer(recent_totals)
)) {
  audit <- data.frame(
    year = years,
    first_rule = as.integer(first_totals),
    recent_rule = as.integer(recent_totals),
    difference = as.integer(recent_totals - first_totals)
  )

  print(audit, row.names = FALSE)

  stop(
    "First-river and most-recent-river pools have different annual totals. ",
    "The two rules should alter assignment, not overall active-pool size."
  )
}

message("  Annual total active fish is identical under both assignment rules.")

# Count the river-year cells whose assignment differs.
n_cells_different <- sum(
  pool_first_ordered != pool_recent_ordered
)

message("  River-year cells differing between assignment rules: ",
        n_cells_different)

# ---- 7. Save validated production objects ------------------------------------

saveRDS(
  pool_first_ordered,
  pool_first_path
)

saveRDS(
  pool_recent_ordered,
  pool_recent_path
)

message("\n5. Saved validated tag-pool objects.")

# ---- 8. Save compact audit tables --------------------------------------------

pool_first_long <- as.data.frame(
  as.table(pool_first_ordered),
  stringsAsFactors = FALSE
)

names(pool_first_long) <- c(
  "river",
  "year",
  "N_first_observed_rule"
)

pool_recent_long <- as.data.frame(
  as.table(pool_recent_ordered),
  stringsAsFactors = FALSE
)

names(pool_recent_long) <- c(
  "river",
  "year",
  "N_most_recent_rule"
)

pool_audit <- full_join(
  pool_first_long,
  pool_recent_long,
  by = c("river", "year")
) %>%
  mutate(
    difference = N_most_recent_rule - N_first_observed_rule
  )

write.csv(
  pool_audit,
  file.path(
    audit_dir,
    "Script4_tag_pool_audit.csv"
  ),
  row.names = FALSE
)

annual_totals <- data.frame(
  year = years,
  first_observed_rule = as.integer(first_totals),
  most_recent_rule = as.integer(recent_totals),
  stringsAsFactors = FALSE
)

write.csv(
  annual_totals,
  file.path(
    audit_dir,
    "Script4_tag_pool_annual_totals.csv"
  ),
  row.names = FALSE
)

# ---- final audit -------------------------------------------------------------

message("\n================ SCRIPT 4 FINAL AUDIT ================")
message("Adult analysis fish: ",
        nrow(informed_at),
        " [expected 1017]")
message("Tag-pool matrix dimensions: 7 x 13")
message("Annual totals identical between assignment rules: TRUE")
message("River-year cells differing between rules: ",
        n_cells_different)
message("Most-recent-river RDS used by Figure A2: ",
        pool_recent_path)
message("Audit CSV: ",
        file.path(audit_dir, "Script4_tag_pool_audit.csv"))
message("======================================================")

message("\nScript 4 FINAL complete.")
