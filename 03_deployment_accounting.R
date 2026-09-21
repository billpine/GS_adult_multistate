################################################################################
# Gulf Sturgeon CJFAS Revision 1
# Script 3: Adult-tag deployment accounting and intermediate RDS files
#
# Purpose
#   1. Identify fish that were tagged as adults, rather than first tagged as
#      subadults and later entering the adult analysis.
#   2. Save the two adult-tagged intermediate RDS files used downstream.
#   3. Audit the adult-tag deployment accounting used for manuscript Table 1.
#   4. Confirm the relationship between the 985 fish tagged as adults and the
#      full 1,017-fish adult analysis population.
#
# IMPORTANT
#   The old juvenile Appendix A2-A4 code has
#   been removed because it is no longer part of the final manuscript workflow
#   and contained legacy/fragile string-manipulation code.
#
# Outputs
#   data/derived/tagged_as_adults_20240121.RDS
#   data/derived/tagged_as_adults_no_2022_deployments_20240121.RDS
#
# Validated accounting
#   - 1,017 fish inform the adult multistate analysis after removing 2022-only
#     histories.
#   - 985 of those fish were tagged as adults.
#   - 32 were first tagged as subadults and later physically recaptured and
#     measured as adults before contributing to the adult analysis.
#   - The 985 adult-tagged fish account for 1,015 pre-2022 transmitter
#     deployments used in manuscript Table 1.
################################################################################

# ---- project setup -----------------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})


# ---- 1. Load validated Script 1 products -------------------------------------

message("\n1. Loading Script 1 products...")

all_stages <- readRDS(
  file.path(data_dir, "informed.all.stages.ch.RDS")
)

adult_histories <- readRDS(
  file.path(
    data_dir,
    "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
  )
)

tags <- readRDS(
  file.path(
    data_dir,
    "informed_at_longevity_all_gs_ms_export_20240121.RDS"
  )
)

all_stages$FishID      <- as.numeric(all_stages$FishID)
adult_histories$FishID <- as.numeric(adult_histories$FishID)
tags$FishID            <- as.numeric(tags$FishID)

message("  All-stage histories: ", nrow(all_stages))
message("  Adult histories from Script 1: ", nrow(adult_histories),
        " [expected 1070]")
message("  Unique tag records: ", nrow(tags),
        " [expected 2610]")

# ---- 2. Identify fish tagged as adults ---------------------------------------

message("\n2. Identifying fish tagged as adults...")

# Lowercase state codes denote subadult observations in the all-stage histories.
juvenile_regex <- "[abcefgh]"

ever_observed_subadult <- grepl(
  juvenile_regex,
  all_stages$ch
)

tagged_as_adults <- all_stages[
  !ever_observed_subadult,
  ,
  drop = FALSE
]

# Retain only fish with at least one adult river-state observation.
adult_states <- c("A", "B", "C", "E", "F", "G", "H")

has_adult_detection <- apply(
  sapply(
    adult_states,
    function(s) grepl(s, tagged_as_adults$ch, fixed = TRUE)
  ),
  1,
  any
)

tagged_as_adults <- tagged_as_adults[
  has_adult_detection,
  ,
  drop = FALSE
]

message("  Fish tagged as adults before final-occasion exclusion: ",
        nrow(tagged_as_adults))

saveRDS(
  tagged_as_adults,
  file.path(data_dir, "tagged_as_adults_20240121.RDS")
)

# ---- 3. Remove 2022-only adult histories -------------------------------------

message("\n3. Removing fish whose only adult encounter occurs in 2022...")

ch2022 <- paste0(
  "000000000000",
  adult_states
)

keep_adult_tagged <- !tagged_as_adults$ch %in% ch2022

tagged_as_adults_no_2022 <- tagged_as_adults[
  keep_adult_tagged,
  ,
  drop = FALSE
]

message("  2022-only adult-tagged histories removed: ",
        sum(!keep_adult_tagged))
message("  Adult-tagged fish retained: ",
        nrow(tagged_as_adults_no_2022),
        " [expected 985]")
message("  Unique FishID retained: ",
        n_distinct(tagged_as_adults_no_2022$FishID),
        " [expected 985]")

if (nrow(tagged_as_adults_no_2022) != 985 ||
    n_distinct(tagged_as_adults_no_2022$FishID) != 985) {
  stop(
    "Adult-tagged fish accounting changed. Expected 985 rows / 985 unique FishID."
  )
}

saveRDS(
  tagged_as_adults_no_2022,
  file.path(
    data_dir,
    "tagged_as_adults_no_2022_deployments_20240121.RDS"
  )
)

# ---- 4. Confirm 985 + 32 = 1,017 adult-analysis fish -------------------------

message("\n4. Auditing adult-analysis population...")

keep_analysis <- !adult_histories$ch %in% ch2022

analysis_adults <- adult_histories[
  keep_analysis,
  ,
  drop = FALSE
]

analysis_ids <- sort(unique(analysis_adults$FishID))
adult_tagged_ids <- sort(unique(tagged_as_adults_no_2022$FishID))

transition_ids <- setdiff(
  analysis_ids,
  adult_tagged_ids
)

unexpected_adult_tagged <- setdiff(
  adult_tagged_ids,
  analysis_ids
)

message("  Adult analysis fish: ", length(analysis_ids),
        " [expected 1017]")
message("  Tagged as adults: ", length(adult_tagged_ids),
        " [expected 985]")
message("  Initially subadult, later entered adult analysis: ",
        length(transition_ids),
        " [expected 32]")
message("  Adult-tagged IDs absent from analysis population: ",
        length(unexpected_adult_tagged),
        " [expected 0]")

if (length(analysis_ids) != 1017 ||
    length(adult_tagged_ids) != 985 ||
    length(transition_ids) != 32 ||
    length(unexpected_adult_tagged) != 0) {
  stop(
    "Adult-analysis composition changed. Expected 1017 total = ",
    "985 tagged-as-adults + 32 subadult-to-adult transitions."
  )
}

# ---- 5. Audit manuscript Table 1 deployment accounting -----------------------

message("\n5. Auditing adult-tag deployments used for manuscript Table 1...")

# Table 1 describes deployments on the 985 fish tagged as adults.
# Exclude 2022 deployments because fish first entering in 2022 cannot inform an
# annual survival interval in the 2010-2022 analysis.

tags_table1 <- tags %>%
  filter(
    FishID %in% adult_tagged_ids,
    start_year <= 2021
  ) %>%
  left_join(
    tagged_as_adults_no_2022 %>%
      transmute(
        FishID = as.numeric(FishID),
        first_code = as.character(first)
      ) %>%
      distinct(FishID, .keep_all = TRUE),
    by = "FishID"
  )

if (any(is.na(tags_table1$first_code))) {
  stop(
    "At least one Table 1 tag record could not be assigned to a first adult river."
  )
}

river_map <- c(
  F = "Pearl",
  E = "Pascagoula",
  C = "Escambia",
  H = "Yellow",
  B = "Choctawhatchee",
  A = "Apalachicola",
  G = "Suwannee"
)

river_order <- unname(river_map)

tags_table1 <- tags_table1 %>%
  mutate(
    river = unname(river_map[first_code]),
    river = factor(river, levels = river_order),
    year_group = ifelse(
      start_year < 2010,
      "Pre-2010",
      as.character(start_year)
    )
  )

if (any(is.na(tags_table1$river))) {
  stop("Unexpected first-river code in Table 1 tag records.")
}

table1_long <- tags_table1 %>%
  count(
    river,
    year_group,
    name = "N_tags"
  )

table1 <- table1_long %>%
  mutate(
    year_group = factor(
      year_group,
      levels = c("Pre-2010", as.character(2010:2021))
    )
  ) %>%
  arrange(river, year_group) %>%
  pivot_wider(
    names_from = year_group,
    values_from = N_tags,
    values_fill = 0
  )

message("\n  TABLE 1 AUDIT: adult-tag deployments by first adult river and year")
print(table1_river_totals)

table1_river_totals <- tags_table1 %>%
  count(
    river,
    name = "N_tags"
  ) %>%
  arrange(river)

message("\n  Table 1 tag totals by river:")
print(table1_river_totals, n = Inf)

n_table1_tags <- nrow(tags_table1)
n_table1_fish <- n_distinct(tags_table1$FishID)

message("  Table 1 transmitter deployments: ",
        n_table1_tags,
        " [expected 1015]")
message("  Table 1 unique fish: ",
        n_table1_fish,
        " [expected 985]")

if (n_table1_tags != 1015 ||
    n_table1_fish != 985) {
  stop(
    "Table 1 accounting changed. Expected 1,015 deployments among 985 fish."
  )
}

# ---- 6. Save a compact audit CSV ---------------------------------------------


write.csv(
  table1_long,
  file.path(audit_dir, "Script3_Table1_deployment_audit.csv"),
  row.names = FALSE
)

# ---- final audit -------------------------------------------------------------

message("\n================ SCRIPT 3 FINAL AUDIT ================")
message("Adult histories before 2022-only exclusion: ",
        nrow(adult_histories),
        " [expected 1070]")
message("Adult analysis fish: ",
        length(analysis_ids),
        " [expected 1017]")
message("Fish tagged as adults: ",
        length(adult_tagged_ids),
        " [expected 985]")
message("Subadult-to-adult transition fish: ",
        length(transition_ids),
        " [expected 32]")
message("Table 1 tag deployments: ",
        n_table1_tags,
        " [expected 1015]")
message("Table 1 unique fish: ",
        n_table1_fish,
        " [expected 985]")
message(
  "Saved: ",
  file.path(data_dir, "tagged_as_adults_20240121.RDS")
)
message(
  "Saved: ",
  file.path(
    data_dir,
    "tagged_as_adults_no_2022_deployments_20240121.RDS"
  )
)
message(
  "Audit CSV: ",
  file.path(audit_dir, "Script3_Table1_deployment_audit.csv")
)
message("======================================================")

message("\nScript 3 FINAL complete.")
