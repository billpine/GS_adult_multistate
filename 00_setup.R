################################################################################
# Gulf Sturgeon CJFAS reproducibility archive
# Script 00: Check software, inputs, and directory structure
################################################################################

source(file.path("code", "_project_paths.R"), local = TRUE)

required_packages <- c(
  "RMark", "tidyverse", "dplyr", "stringr", "tidyr",
  "ggplot2", "gridExtra", "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the workflow, for example:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

required_raw_files <- c(
  file.path(
    raw_dir,
    "MS_TR_1350mmTL_allATc1_MAX_ANNUAL_2010-2022_7basin_20240121.inp"
  ),
  file.path(
    raw_dir,
    "admin_query_export_All_Fish_All_ATags_By_Mark_Status_Min_Date_20240121.csv"
  )
)

missing_raw <- required_raw_files[!file.exists(required_raw_files)]
if (length(missing_raw) > 0) {
  stop(
    "Missing required raw data file(s):\n  ",
    paste(missing_raw, collapse = "\n  "),
    "\n\nPlace the archived raw data files in data/raw/ and rerun."
  )
}

mark_path <- configure_mark(require_mark = FALSE)

message("Project root: ", proj)
message("Raw data directory: ", raw_dir)
message("Derived data directory: ", data_dir)
message("Seven-river model directory: ", seven_results_dir)
message("Regional model directory: ", regional_results_dir)
message("Table output directory: ", table_dir)
message("Figure output directory: ", fig_dir)
if (is.null(mark_path)) {
  message("Program MARK: not yet configured. It is required by Scripts 02 and 06.")
} else {
  message("Program MARK: ", mark_path)
}
message("\nSetup checks complete.")
