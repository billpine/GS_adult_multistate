################################################################################
## 8_manuscript_figures_FINAL.R
##
## Gulf Sturgeon CJFAS manuscript: final/reproducible figure generation
## FINAL production version, 2026-09-18
##
## Final figure workflow
##   Figure 1  Study-area map produced separately.
##   Figure 2  Annual regional apparent survival from the FINAL top model:
##             S(region:time)p(river:time)Psi(strattostrat).
##   Figure 3  Annual river-specific apparent survival from the alternative
##             S(river+time)p(river:time)Psi(strattostrat) model. This model
##             has substantially less support than the top regional model and
##             Figure 3 is therefore interpreted descriptively, not as evidence
##             for an equally supported river-scale survival structure.
##   Figure A1 Frequencies of successive OBSERVED river-state pairs.
##   Figure A2 Precision vs annual tag-pool size from the fully interactive
##             river x year diagnostic model, using the FINAL terminal-p-
##             constrained candidate set.
##
## Important model revision
##   For p(river:time) models, final-occasion detection is constrained equal to
##   detection in the preceding occasion within river. Final candidate models
##   are read from results/models/seven_river/.
##
## Figure A2
##   Boundary survival estimates (S >= 0.999) are excluded because their
##   back-transformed real-scale SEs approach zero and are not useful measures
##   of estimation precision.
################################################################################

# ---- project setup -----------------------------------------------------------

source(file.path("code", "_project_paths.R"), local = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

results_dir <- seven_results_dir

message("Saving final manuscript figures to: ", fig_dir)
message("Using candidate models from: ", results_dir)

# ---- shared labels / mappings ------------------------------------------------

river_map <- c(
  A = "Apalachicola",
  B = "Choctawhatchee",
  C = "Escambia",
  E = "Pascagoula",
  F = "Pearl",
  G = "Suwannee",
  H = "Yellow"
)

river_order <- c(
  "Pearl", "Pascagoula", "Escambia", "Yellow",
  "Choctawhatchee", "Apalachicola", "Suwannee"
)

rivers_alpha <- unname(river_map)

# Representative MARK strata for the four survival regions.
regional_survival_map <- c(
  A = "East",
  B = "Choctawhatchee",
  C = "Pensacola Bay",
  E = "West"
)

region_order <- c(
  "East",
  "Choctawhatchee",
  "Pensacola Bay",
  "West"
)

# Okabe-Ito-inspired colorblind-safe palettes.
region_colours <- c(
  "Choctawhatchee" = "#0072B2",
  "East"            = "#009E73",
  "Pensacola Bay"   = "#E69F00",
  "West"            = "#CC79A7"
)

river_colours <- c(
  Pearl          = "#CC79A7",
  Pascagoula     = "#56B4E9",
  Escambia       = "#009E73",
  Yellow         = "#999999",
  Choctawhatchee = "#0072B2",
  Apalachicola   = "#D55E00",
  Suwannee       = "#E69F00"
)

gs_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92"),
    axis.text         = element_text(colour = "black"),
    axis.title        = element_text(colour = "black"),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      colour = "grey70",
      linewidth = 0.3
    ),
    legend.key.size   = unit(0.9, "lines"),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9)
  )

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

# Parse rows such as "S sA g1 t1".
parse_stratum_time <- function(
    tbl,
    parameter,
    state_map,
    first_year = 2010,
    state_col = "group") {

  parsed <- tbl %>%
    mutate(
      stratum_code = str_match(
        MARK_row,
        paste0("^", parameter, " s([A-Z])")
      )[, 2],
      mark_time = as.numeric(
        str_match(MARK_row, " t([0-9.]+)")[, 2]
      )
    )

  if (any(is.na(parsed$stratum_code))) {
    stop(
      "Could not parse one or more ",
      parameter,
      " stratum codes from MARK row names."
    )
  }

  if (any(is.na(parsed$mark_time))) {
    stop(
      "Could not parse one or more ",
      parameter,
      " time values from MARK row names."
    )
  }

  unexpected <- setdiff(
    unique(parsed$stratum_code),
    names(state_map)
  )

  if (length(unexpected) > 0) {
    stop(
      "Unexpected ",
      parameter,
      " stratum code(s): ",
      paste(unexpected, collapse = ", ")
    )
  }

  time_levels <- sort(unique(parsed$mark_time))

  year_lookup <- setNames(
    first_year + seq_along(time_levels) - 1,
    as.character(time_levels)
  )

  parsed$year <- unname(
    year_lookup[as.character(parsed$mark_time)]
  )

  parsed[[state_col]] <- unname(
    state_map[parsed$stratum_code]
  )

  parsed
}

extract_S <- function(
    mod,
    state_map,
    state_col = "group") {

  get_real(mod, "S") %>%
    parse_stratum_time(
      parameter = "S",
      state_map = state_map,
      first_year = 2010,
      state_col = state_col
    ) %>%
    transmute(
      year,
      !!state_col := .data[[state_col]],
      estimate,
      se,
      lcl,
      ucl
    )
}

# ---- 1. Load and audit final candidate set -----------------------------------

message("\n1. Loading final candidate models...")

all_mods <- lapply(
  seq_len(25),
  function(i) {
    readRDS(
      file.path(
        results_dir,
        paste0("mod", i, ".RDS")
      )
    )
  }
)

names(all_mods) <- paste0(
  "mod",
  seq_len(25)
)

aicc_values <- vapply(
  all_mods,
  function(m) m$results$AICc,
  numeric(1)
)

delta_aicc <- aicc_values -
  min(aicc_values)

model_weights <- exp(
  -0.5 * delta_aicc
) /
  sum(
    exp(-0.5 * delta_aicc)
  )

ranking <- names(
  sort(aicc_values)
)

expected_top5 <- c(
  "mod1",
  "mod2",
  "mod3",
  "mod5",
  "mod6"
)

if (!identical(
  ranking[1:5],
  expected_top5
)) {
  stop(
    "Final candidate-model ordering changed. Expected top five: ",
    paste(expected_top5, collapse = ", "),
    "; observed: ",
    paste(ranking[1:5], collapse = ", "),
    "."
  )
}

if (abs(
  aicc_values["mod1"] -
    6133.295
) > 0.05) {
  warning(
    "mod1 AICc differs from validated final value. Observed = ",
    aicc_values["mod1"]
  )
}

if (abs(
  model_weights["mod1"] -
    0.980688
) > 0.005) {
  warning(
    "mod1 weight differs from validated final value. Observed = ",
    model_weights["mod1"]
  )
}

message(
  "  Final weights: mod1 = ",
  sprintf("%.6f", model_weights["mod1"]),
  "; mod2 = ",
  sprintf("%.6f", model_weights["mod2"]),
  "; mod3 = ",
  sprintf("%.6f", model_weights["mod3"])
)

mod1 <- all_mods[["mod1"]]

mod1_prof <- readRDS(
  file.path(
    results_dir,
    "mod1_prof.RDS"
  )
)

if (abs(
  mod1_prof$results$AICc -
    mod1$results$AICc
) > 0.01 ||
    abs(
      mod1_prof$results$lnl -
        mod1$results$lnl
    ) > 0.01) {
  stop(
    "mod1_prof does not match final mod1 likelihood/AICc. ",
    "Run final Script 2 before Script 8."
  )
}

if ((!is.null(mod1$results$singular) &&
     length(mod1$results$singular) > 0) ||
    (!is.null(mod1_prof$results$singular) &&
     length(mod1_prof$results$singular) > 0)) {
  stop(
    "Final mod1 or mod1_prof reports singular parameter indices. ",
    "Inspect Script 2 before generating figures."
  )
}

# ==============================================================================
# FIGURE 1 -- study area map
# ==============================================================================

message("\nFigure 1 (study area map) -- skipped; produced separately.")

# ==============================================================================
# FIGURE 2 -- regional apparent survival from FINAL TOP MODEL
# ==============================================================================

message("\n2. Building Figure 2 from final top model...")

s_mod1 <- extract_S(
  mod1_prof,
  state_map = regional_survival_map,
  state_col = "Region"
)

if (nrow(s_mod1) != 48L ||
    n_distinct(s_mod1$Region) != 4L ||
    n_distinct(s_mod1$year) != 12L) {
  stop(
    "Unexpected final top-model S dimensions. Expected 48 rows = ",
    "4 regions x 12 annual intervals."
  )
}

s_mod1 <- s_mod1 %>%
  mutate(
    Region = factor(
      Region,
      levels = region_order
    )
  ) %>%
  arrange(
    Region,
    year
  )

fig2 <- ggplot(
  s_mod1,
  aes(
    x = year,
    y = estimate,
    colour = Region,
    fill = Region,
    group = Region
  )
) +
  geom_ribbon(
    aes(
      ymin = lcl,
      ymax = ucl
    ),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.2
  ) +
  scale_x_continuous(
    breaks = c(
      2010,
      2013,
      2016,
      2019,
      2021
    ),
    expand = c(
      0.05,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0.5,
      1.02
    ),
    breaks = seq(
      0.5,
      1.0,
      0.1
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_colour_manual(
    values = region_colours
  ) +
  scale_fill_manual(
    values = region_colours
  ) +
  facet_wrap(
    ~ Region,
    ncol = 4
  ) +
  labs(
    x = "Year",
    y = expression(
      paste(
        "Annual apparent survival probability (",
        italic(Phi),
        ")"
      )
    )
  ) +
  gs_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    legend.position = "none"
  )

tiff(
  file.path(
    fig_dir,
    "Figure2_regional_survival.tiff"
  ),
  width = 12,
  height = 5,
  units = "in",
  res = 300
)

print(fig2)
dev.off()

write.csv(
  s_mod1,
  file.path(
    fig_dir,
    "Figure2_regional_survival_data.csv"
  ),
  row.names = FALSE
)

message(
  "  Figure 2 saved from ",
  mod1$model.name,
  " | weight = ",
  sprintf("%.6f", model_weights["mod1"])
)

# ==============================================================================
# FIGURE 3 -- river-specific apparent survival from alternative mod3
# ==============================================================================
#
# mod3 = S(river+time)p(river:time)Psi(strattostrat)
#
# This model has substantially less support than the final top model
# (Delta AICc = 10.25; weight ~0.006). Figure 3 is retained because the
# river-specific patterns are biologically useful for interpretation, but it
# should be described as a descriptive alternative-model result rather than as
# evidence for an equally supported river-scale survival structure.

message("\n3. Building Figure 3 from alternative mod3...")

mod3 <- all_mods[["mod3"]]

s_mod3 <- extract_S(
  mod3,
  state_map = river_map,
  state_col = "River"
)

if (nrow(s_mod3) != 84L ||
    n_distinct(s_mod3$River) != 7L ||
    n_distinct(s_mod3$year) != 12L) {
  stop(
    "Unexpected mod3 S dimensions. Expected 84 rows = 7 rivers x 12 years."
  )
}

s_mod3 <- s_mod3 %>%
  mutate(
    River = factor(
      River,
      levels = river_order
    )
  ) %>%
  arrange(
    River,
    year
  )

fig3 <- ggplot(
  s_mod3,
  aes(
    x = year,
    y = estimate,
    colour = River,
    fill = River,
    group = River
  )
) +
  geom_ribbon(
    aes(
      ymin = lcl,
      ymax = ucl
    ),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.2
  ) +
  scale_x_continuous(
    breaks = c(
      2010,
      2013,
      2016,
      2019,
      2021
    ),
    expand = c(
      0.05,
      0
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      1.0,
      0.2
    ),
    expand = c(
      0,
      0
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      1.02
    )
  ) +
  scale_colour_manual(
    values = river_colours
  ) +
  scale_fill_manual(
    values = river_colours
  ) +
  facet_wrap(
    ~ River,
    ncol = 4
  ) +
  labs(
    x = "Year",
    y = expression(
      paste(
        "Annual apparent survival probability (",
        italic(Phi),
        ")"
      )
    )
  ) +
  gs_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

tiff(
  file.path(
    fig_dir,
    "Figure3_river_survival.tiff"
  ),
  width = 12,
  height = 7,
  units = "in",
  res = 300
)

print(fig3)
dev.off()

write.csv(
  s_mod3,
  file.path(
    fig_dir,
    "Figure3_river_survival_data.csv"
  ),
  row.names = FALSE
)

message(
  "  Figure 3 saved from ",
  mod3$model.name,
  " | Delta AICc = ",
  sprintf("%.3f", delta_aicc["mod3"]),
  " | weight = ",
  sprintf("%.6f", model_weights["mod3"])
)

# ==============================================================================
# FIGURE A1 -- successive OBSERVED river-state pairs
# ==============================================================================

message("\n4. Building Figure A1...")

at_inf_adults <- readRDS(
  file.path(data_dir, "informed_TR_MS_AT_ANNUAL_no_singles_2010-2022_7basin_20240121.RDS")
)

at_inf_adults$FishID <- as.numeric(
  at_inf_adults$FishID
)

ch2022 <- paste0(
  "000000000000",
  names(river_map)
)

at_inf_adults <- at_inf_adults[
  !at_inf_adults$ch %in% ch2022,
  ,
  drop = FALSE
]

if (nrow(at_inf_adults) != 1017L ||
    n_distinct(at_inf_adults$FishID) != 1017L) {
  stop(
    "Figure A1 adult-history accounting changed. Expected 1,017 rows / fish."
  )
}

possible_states <- names(
  river_map
)

trans_list <- vector(
  "list",
  nrow(at_inf_adults)
)

for (i in seq_len(
  nrow(at_inf_adults)
)) {

  b <- strsplit(
    at_inf_adults$ch[i],
    "",
    fixed = TRUE
  )[[1]]

  # Keep observed river states only. This skips both nondetection ("0") and
  # unavailable/censored (".") occasions. Pairs therefore represent successive
  # OBSERVED river states and may span one or more intervening years.
  b_obs <- b[
    b %in% possible_states
  ]

  if (length(b_obs) < 2L) {
    next
  }

  pairs <- paste0(
    b_obs[-length(b_obs)],
    b_obs[-1]
  )

  trans_list[[i]] <- data.frame(
    transition = pairs,
    FishID = at_inf_adults$FishID[i],
    stringsAsFactors = FALSE
  )
}

trans_df <- bind_rows(
  trans_list
)

trans_summary <- trans_df %>%
  mutate(
    from = substr(
      transition,
      1,
      1
    ),
    to = substr(
      transition,
      2,
      2
    ),
    source_river = unname(
      river_map[from]
    ),
    dest_river = unname(
      river_map[to]
    )
  ) %>%
  filter(
    !is.na(source_river),
    !is.na(dest_river)
  ) %>%
  count(
    source_river,
    dest_river,
    name = "n"
  ) %>%
  mutate(
    source_river = factor(
      source_river,
      levels = river_order
    ),
    dest_river = factor(
      dest_river,
      levels = river_order
    )
  ) %>%
  arrange(
    source_river,
    dest_river
  )

if (nrow(trans_summary) == 0L) {
  stop(
    "No successive observed river-state pairs were generated for Figure A1."
  )
}

y_max <- max(
  trans_summary$n
) * 1.14

plot_one_river <- function(
    source,
    show_y = FALSE) {

  d <- filter(
    trans_summary,
    source_river == source
  )

  ggplot(
    d,
    aes(
      x = dest_river,
      y = n,
      fill = dest_river,
      colour = dest_river
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    geom_text(
      aes(
        label = n,
        y = n
      ),
      vjust = -0.35,
      size = 3,
      fontface = "bold",
      colour = "black"
    ) +
    scale_fill_manual(
      values = river_colours,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = river_colours,
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(
        0,
        y_max
      ),
      expand = c(
        0,
        0
      )
    ) +
    labs(
      x = source,
      y = if (show_y) {
        "Successive observed-state frequency"
      } else {
        NULL
      }
    ) +
    gs_theme +
    theme(
      axis.text.x = element_text(
        angle = 40,
        hjust = 1,
        size = 8
      ),
      axis.text.y = if (show_y) {
        element_text()
      } else {
        element_blank()
      },
      axis.ticks.y = if (show_y) {
        element_line()
      } else {
        element_blank()
      },
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
}

plots <- mapply(
  plot_one_river,
  source = river_order,
  show_y = c(
    TRUE,
    rep(
      FALSE,
      6
    )
  ),
  SIMPLIFY = FALSE
)

figA1 <- arrangeGrob(
  grobs = plots,
  nrow = 1,
  widths = c(
    1.4,
    rep(
      1,
      6
    )
  )
)

tiff(
  file.path(
    fig_dir,
    "FigureA1_transition_frequencies.tiff"
  ),
  width = 16,
  height = 5.5,
  units = "in",
  res = 300
)

grid.draw(
  figA1
)

dev.off()

write.csv(
  trans_summary,
  file.path(
    fig_dir,
    "FigureA1_successive_observed_state_frequencies.csv"
  ),
  row.names = FALSE
)

message(
  "  Figure A1 saved."
)

# ==============================================================================
# FIGURE A2 -- precision vs tag pool size from fully interactive diagnostic model
# ==============================================================================

message("\n5. Building Figure A2...")

# This model is used only as a precision/sample-size diagnostic because it
# provides separate river-year survival estimates. It is NOT used as the primary
# biological inference model.
mod5 <- all_mods[["mod5"]]

if (!grepl(
  "S\\(river:time\\).*p\\(river:time\\)",
  mod5$model.name
)) {
  stop(
    "mod5 is not the expected S(river:time)p(river:time) model. ",
    "Observed model name: ",
    mod5$model.name
  )
}

message(
  "  Figure A2 model: ",
  mod5$model.name,
  " | AICc = ",
  round(
    mod5$results$AICc,
    3
  ),
  " | Delta AICc = ",
  round(
    delta_aicc["mod5"],
    3
  )
)

pool_recent <- readRDS(
  file.path(data_dir, "tableA_pool_recentrule.RDS")
)

missing_pool_rivers <- setdiff(
  rivers_alpha,
  rownames(pool_recent)
)

if (length(
  missing_pool_rivers
) > 0) {
  stop(
    "Tag-pool matrix is missing river row(s): ",
    paste(
      missing_pool_rivers,
      collapse = ", "
    )
  )
}

required_pool_years <- as.character(
  2010:2021
)

missing_pool_years <- setdiff(
  required_pool_years,
  colnames(pool_recent)
)

if (length(
  missing_pool_years
) > 0) {
  stop(
    "Tag-pool matrix is missing year column(s): ",
    paste(
      missing_pool_years,
      collapse = ", "
    )
  )
}

pool_df <- as.data.frame(
  pool_recent[
    rivers_alpha,
    required_pool_years,
    drop = FALSE
  ]
) %>%
  rownames_to_column(
    "river"
  ) %>%
  pivot_longer(
    cols = all_of(
      required_pool_years
    ),
    names_to = "year",
    values_to = "pool"
  ) %>%
  mutate(
    year = as.integer(
      year
    ),
    pool = as.numeric(
      pool
    )
  )

s_mod5 <- extract_S(
  mod5,
  state_map = river_map,
  state_col = "river"
)

if (nrow(s_mod5) != 84L ||
    n_distinct(s_mod5$river) != 7L ||
    n_distinct(s_mod5$year) != 12L) {
  stop(
    "Unexpected mod5 S dimensions. Expected 84 rows = 7 rivers x 12 years."
  )
}

figA2_df_all <- s_mod5 %>%
  left_join(
    pool_df,
    by = c(
      "river",
      "year"
    )
  )

if (any(
  is.na(
    figA2_df_all$pool
  )
)) {
  stop(
    "At least one Figure A2 survival estimate could not be matched to a tag-pool size."
  )
}

boundary_n <- sum(
  figA2_df_all$estimate >= 0.999,
  na.rm = TRUE
)

message(
  "  Figure A2 boundary estimates excluded (S >= 0.999): ",
  boundary_n
)

figA2_df <- figA2_df_all %>%
  filter(
    estimate < 0.999
  ) %>%
  mutate(
    river = factor(
      river,
      levels = river_order
    )
  ) %>%
  arrange(
    river,
    year
  )

if (any(
  !is.finite(
    figA2_df$se
  )
)) {
  stop(
    "Nonfinite SE remains among Figure A2 plotted observations."
  )
}

figA2 <- ggplot(
  figA2_df,
  aes(
    x = pool,
    y = se,
    fill = river,
    colour = river
  )
) +
  geom_point(
    pch = 21,
    size = 4.5,
    colour = "white",
    stroke = 0.4,
    aes(fill = river)
  ) +
  geom_point(
    pch = 21,
    size = 4.5,
    fill = NA,
    aes(colour = river),
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = river_colours,
    name = "River"
  ) +
  scale_colour_manual(
    values = river_colours,
    name = "River"
  ) +
  scale_x_continuous(
    limits = c(
      0,
      160
    ),
    breaks = seq(
      0,
      150,
      25
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      0.42
    ),
    breaks = seq(
      0,
      0.4,
      0.05
    ),
    expand = c(
      0,
      0
    )
  ) +
  labs(
    x = "Annual tag pool size",
    y = "Standard error of annual apparent survival estimate"
  ) +
  gs_theme +
  theme(
    legend.position = c(
      0.98,
      0.98
    ),
    legend.justification = c(
      1,
      1
    )
  )

tiff(
  file.path(
    fig_dir,
    "FigureA2_precision_vs_pool_size.tiff"
  ),
  width = 9,
  height = 6.5,
  units = "in",
  res = 300
)

print(
  figA2
)

dev.off()

write.csv(
  figA2_df,
  file.path(
    fig_dir,
    "FigureA2_precision_vs_pool_size_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  figA2_df_all,
  file.path(
    fig_dir,
    "Supporting_FigureA2_all_river_years.csv"
  ),
  row.names = FALSE
)

message(
  "  Figure A2 saved."
)

# ==============================================================================
# FINAL FIGURE AUDIT
# ==============================================================================

message(
  "\n================ SCRIPT 8 FINAL AUDIT ================"
)

message(
  "Candidate-model directory: ",
  results_dir
)

message(
  "mod1 weight: ",
  sprintf(
    "%.6f",
    model_weights["mod1"]
  )
)

message(
  "mod2 weight: ",
  sprintf(
    "%.6f",
    model_weights["mod2"]
  )
)

message(
  "mod3 weight: ",
  sprintf(
    "%.6f",
    model_weights["mod3"]
  )
)

message(
  "Figure 2 top-model S rows: ",
  nrow(s_mod1),
  " [expected 48]"
)

message(
  "Figure 3/mod3 S rows: ",
  nrow(s_mod3),
  " [expected 84]"
)

message(
  "Figure 3/mod3 Delta AICc: ",
  sprintf("%.3f", delta_aicc["mod3"]),
  " | weight = ",
  sprintf("%.6f", model_weights["mod3"]),
  " [descriptive alternative model]"
)

message(
  "Figure A1 analysis fish: ",
  nrow(at_inf_adults),
  " [expected 1017]"
)

message(
  "Figure A1 successive observed-state pair records: ",
  nrow(trans_df)
)

message(
  "Figure A2 mod5 S rows: ",
  nrow(s_mod5),
  " [expected 84]"
)

message(
  "Figure A2 boundary estimates excluded: ",
  boundary_n
)

message(
  "Figure A2 plotted river-year rows: ",
  nrow(figA2_df),
  " [84 - boundary_n]"
)

message(
  "Figure A2 model Delta AICc: ",
  sprintf(
    "%.3f",
    delta_aicc["mod5"]
  ),
  " [diagnostic model]"
)

files <- list.files(
  fig_dir,
  pattern = "\\.(tiff|csv)$",
  full.names = FALSE
)

message(
  "\nCurrent final figure products:"
)

for (f in sort(files)) {
  message(
    "  ",
    f
  )
}

message(
  "\nNote: Figure 1 (study area map) is produced separately."
)

message(
  "======================================================"
)

message(
  "\nScript 8 FINAL complete."
)
