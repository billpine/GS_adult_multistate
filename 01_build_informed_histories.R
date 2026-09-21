################################################################################
# Gulf Sturgeon CJFAS Revision 1
# Script 1: Build informed encounter histories
# Revised 2026-09-18
#
# Purpose
#   Read the raw annual encounter-history input and acoustic-tag deployment
#   records, apply transmitter-longevity censoring, remove subadult-only
#   encounter histories, and save the four RDS objects used downstream.
#
# IMPORTANT
#   This script reads the two archived raw inputs in data/raw/ and writes
#   regenerated intermediate objects to data/derived/.
#
# Outputs
#   data/derived/
#     pre2010.RDS
#     informed.all.stages.ch.RDS
#     informed_at_longevity_all_gs_ms_export_20240121.RDS
#     informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS
#
# Notes
#   - Missing manufacturer tag longevity is assigned 5 years, preserving the
#     original analysis rule.
#   - Tag-life indexing is preserved exactly from the validated original
#     workflow. This script does not alter the annual endpoint convention.
#   - The raw .inp already reflects the upstream "no singles" filtering and
#     annual encounter-state construction.
################################################################################

# ---- project setup -----------------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(RMark)
})


message("Project: ", normalizePath(proj, winslash = "/", mustWork = FALSE))
message("Raw inputs: ", normalizePath(raw_dir, winslash = "/", mustWork = FALSE))
message("Derived-data outputs: ",
        normalizePath(data_dir, winslash = "/", mustWork = FALSE))

# ---- 1. Load raw annual encounter histories ---------------------------------

message("\n1. Loading raw encounter histories...")

raw_inp <- file.path(
  raw_dir,
  "MS_TR_1350mmTL_allATc1_MAX_ANNUAL_2010-2022_7basin_20240121.inp"
)

dat <- convert.inp(
  raw_inp,
  covariates = c("length_cutoff_mm", "TL", "FL"),
  use.comments = TRUE
)

dat$FishID <- as.numeric(rownames(dat))

message("  Raw encounter-history rows: ", nrow(dat))

# BASIN KEY
# A = Apalachicola
# B = Choctawhatchee
# C = Escambia
# D = Ochlockonee (not used in seven-river adult model)
# E = Pascagoula
# F = Pearl
# G = Suwannee
# H = Yellow

# ---- 2. Load and clean acoustic-tag deployment information ------------------

message("\n2. Loading acoustic-tag deployment information...")

tag_csv <- file.path(
  raw_dir,
  "admin_query_export_All_Fish_All_ATags_By_Mark_Status_Min_Date_20240121.csv"
)

all_tag_info <- read.csv(
  tag_csv,
  na.strings = c("", "NA")
)

all_tag_info$Deployment_Date <- as.Date(
  all_tag_info$Deployed,
  "%d-%b-%y"
)

all_tag_info$Read_Date <- as.Date(
  all_tag_info$Read,
  "%d-%b-%y"
)

# For each acoustic tag, define the earliest known deployment/read date.
all_tag_info2 <- all_tag_info %>%
  group_by(Mark_Code_Space) %>%
  mutate(
    TrueMinContactDate = min(
      Deployment_Date,
      Read_Date,
      na.rm = TRUE
    )
  ) %>%
  arrange(FishID) %>%
  ungroup()

# Keep records with a nonmissing acoustic-tag code.
#
all_tag_info3 <- all_tag_info2 %>%
  filter(!is.na(Mark_Code_Space))

message("  Raw tag-export rows: ", nrow(all_tag_info))
message("  Rows with missing Mark_Code_Space: ",
        sum(is.na(all_tag_info$Mark_Code_Space)))
message("  Rows retained with nonmissing Mark_Code_Space: ",
        nrow(all_tag_info3))

# Collapse directly to unique fish-tag-longevity-date combinations.
#

tags <- all_tag_info3 %>%
  transmute(
    FishID = as.numeric(FishID),
    Mark_Code_Space = as.character(Mark_Code_Space),
    Mark_Longevity_Days = as.numeric(Mark_Longevity_Days),
    TrueMinContactDate = as.Date(TrueMinContactDate)
  ) %>%
  distinct()

message("  Unique fish-tag-longevity-date records: ", nrow(tags))

if (nrow(tags) != 2610) {
  warning(
    "Expected 2,610 unique tag records from the validated 2024-01-21 export; ",
    "found ", nrow(tags), ". Confirm that the raw export has not changed."
  )
}

# ---- 3. Compute tag-life indices ---------------------------------------------

message("\n3. Computing transmitter-life indices...")

# Preserve original analysis rule:
#   missing manufacturer longevity -> 5 years
#   otherwise floor manufacturer longevity in days to whole years.
tags$tag_life_years <- ifelse(
  is.na(tags$Mark_Longevity_Days),
  5,
  floor(tags$Mark_Longevity_Days / 365)
)

tags$start_year <- as.numeric(
  format(tags$TrueMinContactDate, "%Y")
)

# Preserve the original annual endpoint convention exactly.
tags$end_year <- tags$start_year + tags$tag_life_years

tags$start_index <- ifelse(
  tags$start_year >= 2010,
  tags$start_year - 2009,
  1
)

tags$end_index <- ifelse(
  tags$end_year >= 2010,
  tags$end_year - 2009,
  0
)

message("  Tags with missing manufacturer longevity assigned 5 years: ",
        sum(is.na(tags$Mark_Longevity_Days)))
message("  Earliest tag start year: ",
        min(tags$start_year, na.rm = TRUE))
message("  Latest tag start year: ",
        max(tags$start_year, na.rm = TRUE))

# Fish tagged before 2010 whose transmitters could contribute during the
# 2010-2022 study window.
pre2010 <- dat[
  dat$FishID %in% tags$FishID[tags$start_year < 2010],
  ,
  drop = FALSE
]

saveRDS(
  pre2010,
  file.path(data_dir, "pre2010.RDS")
)

saveRDS(
  tags,
  file.path(
    data_dir,
    "informed_at_longevity_all_gs_ms_export_20240121.RDS"
  )
)

message("  Saved pre-2010 fish: ", nrow(pre2010))
message("  Saved tag table: ", nrow(tags), " records")

# ---- 4. Apply tag-life censoring to annual encounter histories ---------------

message("\n4. Applying transmitter-life censoring...")

freq <- as.matrix(dat$freq)
eh   <- dat$ch

eh.m <- str_split_fixed(
  eh,
  "",
  nchar(eh[1])
)

# add 12 dummy columns so tags whose expected life extends beyond 2022
# can be handled without truncating the tag-life sequence during censoring.
dummycolumns <- matrix(
  "X",
  nrow = nrow(eh.m),
  ncol = 12
)

eh.m.dummy <- cbind(
  eh.m,
  dummycolumns
)

rownames(eh.m.dummy) <- dat$FishID

tags$FishID <- as.numeric(tags$FishID)

for (i in seq_len(nrow(eh.m.dummy))) {

  fish_id <- as.numeric(rownames(eh.m.dummy)[i])

  temp <- tags %>%
    filter(FishID == fish_id)

  # One transmitter: censor all annual occasions after expected tag life.
  if (nrow(temp) == 1) {

    if (!is.na(temp$end_index) &&
        temp$end_index < ncol(eh.m.dummy)) {

      eh.m.dummy[
        i,
        seq(temp$end_index + 1, ncol(eh.m.dummy))
      ] <- "."
    }
  }

  # Multiple transmitters: construct the union of all active tag-life periods
  # and censor gaps between/after deployments.
  if (nrow(temp) > 1) {

    active_indices <- unlist(
      lapply(seq_len(nrow(temp)), function(j) {
        seq(temp$start_index[j], temp$end_index[j])
      })
    )

    all_indices <- seq(
      min(temp$start_index),
      ncol(eh.m.dummy)
    )

    gaps <- setdiff(
      all_indices,
      active_indices
    )

    eh.m.dummy[i, gaps] <- "."
  }
}

# Remove the 12 dummy columns, returning to the 13 annual occasions, 2010-2022.
eh.final <- eh.m.dummy[
  ,
  -seq(
    ncol(eh.m.dummy) - 11,
    ncol(eh.m.dummy)
  ),
  drop = FALSE
]

# ---- 5. Write informed .inp and reload ---------------------------------------

message("\n5. Writing informed encounter-history input...")

informed_inp <- file.path(
  data_dir,
  "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.inp"
)

sink(informed_inp)

for (i in seq_len(nrow(eh.final))) {
  cat(
    paste(eh.final[i, ], collapse = ""),
    freq[i],
    "\n"
  )
}

sink()

temp <- convert.inp(informed_inp)
dat$ch <- temp$ch

# Remove histories containing no observed river state at any life stage.
possible_states <- c(
  "A", "B", "C", "E", "F", "G", "H",
  "a", "b", "c", "e", "f", "g", "h"
)

has_detection <- apply(
  sapply(
    possible_states,
    function(s) grepl(s, dat$ch, fixed = TRUE)
  ),
  1,
  any
)

dat <- dat[
  has_detection,
  ,
  drop = FALSE
]

message("  Histories with at least one observed river state: ", nrow(dat))

# ---- 6. Identify first observed river and life stage -------------------------

message("\n6. Identifying first observed river and life stage...")

eh.m2 <- str_split_fixed(
  dat$ch,
  "",
  nchar(dat$ch[1])
)

get.first <- function(x) {
  min(which(x != "0" & x != "."))
}

first_index <- apply(
  eh.m2,
  1,
  get.first
)

dat$first.stage <- mapply(
  function(row, col) eh.m2[row, col],
  seq_len(nrow(eh.m2)),
  first_index
)

# Convert lowercase subadult codes to their corresponding adult river code only
# for first-river tracking. The encounter history itself is unchanged here.
juvenile_to_adult <- c(
  a = "A",
  b = "B",
  c = "C",
  e = "E",
  f = "F",
  g = "G",
  h = "H"
)

dat$first <- ifelse(
  dat$first.stage %in% names(juvenile_to_adult),
  unname(juvenile_to_adult[dat$first.stage]),
  dat$first.stage
)

saveRDS(
  dat,
  file.path(data_dir, "informed.all.stages.ch.RDS")
)

message("  Saved all-stage encounter histories: ", nrow(dat))

# ---- 7. Build adult-only encounter histories ---------------------------------

message("\n7. Building adult-only encounter histories...")

adult_states <- c(
  "A", "B", "C", "E", "F", "G", "H"
)

# Replace lowercase subadult detections with nondetections while preserving
# transmitter-censoring "." occasions.
dat$ch <- str_replace_all(
  dat$ch,
  c(
    a = "0",
    b = "0",
    c = "0",
    e = "0",
    f = "0",
    g = "0",
    h = "0"
  )
)

has_adult <- apply(
  sapply(
    adult_states,
    function(s) grepl(s, dat$ch, fixed = TRUE)
  ),
  1,
  any
)

dat <- dat[
  has_adult,
  ,
  drop = FALSE
]

# Re-extract the first observed adult river from the adult-only history.
eh.m3 <- str_split_fixed(
  dat$ch,
  "",
  nchar(dat$ch[1])
)

first_adult_index <- apply(
  eh.m3,
  1,
  get.first
)

dat$first <- mapply(
  function(row, col) eh.m3[row, col],
  seq_len(nrow(eh.m3)),
  first_adult_index
)

message("  Adult encounter-history rows: ", nrow(dat),
        " [validated expectation = 1070]")

if (nrow(dat) != 1070) {
  warning(
    "Expected 1,070 adult histories from the validated 2024-01-21 inputs; ",
    "found ", nrow(dat), ". Confirm that upstream files have not changed."
  )
}

saveRDS(
  dat,
  file.path(
    data_dir,
    "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
  )
)

# ---- final audit -------------------------------------------------------------

message("\n================ SCRIPT 1 AUDIT ================")
message("Unique tag records: ", nrow(tags), " [expected 2610]")
message("Pre-2010 fish: ", nrow(pre2010))
message("All-stage encounter histories: ",
        nrow(readRDS(file.path(data_dir, "informed.all.stages.ch.RDS"))))
message("Adult encounter histories: ", nrow(dat), " [expected 1070]")
message("Adult-history output: ",
        file.path(
          data_dir,
          "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS"
        ))
message("================================================")

message("\nScript 1 complete.")
