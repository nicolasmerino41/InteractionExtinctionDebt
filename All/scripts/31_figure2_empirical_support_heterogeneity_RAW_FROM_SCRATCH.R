## ------------------------------------------------------------
## Script: All/scripts/31_figure2_empirical_support_heterogeneity.R
## Figure 2: Empirical support heterogeneity + raw per-dataset panel
## ------------------------------------------------------------
## Purpose
##   1. Build the main Figure 2 with:
##      A) a dataset-level distribution of local interaction support
##      B) a dataset-balanced heatmap of co-occurrence support vs interaction support
##   2. Also save an additional raw-data 2 x 5 faceted point plot:
##      full_n = sites where pair is recorded together
##      full_K = sites where interaction is observed
##      one point = one co-occurring consumer-resource pair
##
## Notes
##   - Uses all co-occurring pairs for Panel B and the raw 2 x 5 plot.
##   - Keeps pairs with full_K = 0.
##   - Uses only realised regional interaction links, full_K >= 1, for Panel A.
##   - Does not fit models or run site-removal simulations.
## ------------------------------------------------------------

rm(list = ls())

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr", "tidyr", "ggplot2", "tibble", "patchwork",
  "scales", "forcats", "purrr"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

## ---------------------------
## Output folder
## ---------------------------

out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, filename, width, height, dpi = 320){
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

## ---------------------------
## Dataset order
## ---------------------------

preferred_order <- c(
  "Garraf_PP",
  "Garraf_PP2",
  "Montseny",
  "Gottin_PP",
  "Nahuel",
  "Garraf_HP",
  "Quercus",
  "Olot",
  "Gottin_HP",
  "Salix_Galpar",
  "Galpar"
)

available_datasets <- all_dataset_names

dataset_order <- preferred_order[preferred_order %in% available_datasets]
dataset_order <- c(dataset_order, setdiff(available_datasets, dataset_order))

## Keep exactly the 10 original datasets supplied by the helper object.
## If both Salix_Galpar and Galpar naming variants exist, the helper list decides.
dataset_order <- dataset_order[seq_len(min(10, length(dataset_order)))]

## ---------------------------
## Colours and theme
## ---------------------------

support_levels <- c(
  "Observed in 1 site",
  "Observed in 2 sites",
  "Observed in 3–4 sites",
  "Observed in 5–9 sites",
  "Observed in 10 or more sites"
)

support_cols <- c(
  "Observed in 1 site" = "#D55E00",
  "Observed in 2 sites" = "#E69F00",
  "Observed in 3–4 sites" = "#56B4E9",
  "Observed in 5–9 sites" = "#0072B2",
  "Observed in 10 or more sites" = "#004C6D"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "white", colour = "grey35", linewidth = 0.35),
    strip.text = element_text(face = "plain", size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

## ---------------------------
## Data construction
## ---------------------------

make_pair_support <- function(dataset){
  message("Loading pair support for: ", dataset)

  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    transmute(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    ) %>%
    distinct(site, consumer, resource)

  ints <- site_tables$empirical_site_interactions %>%
    transmute(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    ) %>%
    distinct(site, consumer, resource)

  ## Defensive addition: observed interactions imply operational co-occurrence.
  ## This prevents K > n if the helper tables differ in how they encode co-occurrence.
  cooc <- bind_rows(cooc, ints) %>%
    distinct(site, consumer, resource)

  cooc_counts <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  int_counts <- ints %>%
    group_by(consumer, resource) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  out <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      full_B = full_K > 0,
      interaction_occupancy = full_K / full_n
    ) %>%
    select(dataset, consumer, resource, full_n, full_K, full_B, interaction_occupancy)

  if(any(out$full_n < 1)) stop("full_n < 1 in ", dataset)
  if(any(out$full_K < 0)) stop("full_K < 0 in ", dataset)
  if(any(out$full_K > out$full_n)) stop("full_K > full_n in ", dataset)

  out
}

all_pairs <- bind_rows(lapply(dataset_order, make_pair_support)) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

write.csv2(
  all_pairs,
  file.path(out_dir, "31_all_cooccurrence_pair_support.csv"),
  row.names = FALSE
)

## ---------------------------
## Panel A: support distribution for realised links only
## ---------------------------

classify_support <- function(K){
  case_when(
    K == 1 ~ "Observed in 1 site",
    K == 2 ~ "Observed in 2 sites",
    K >= 3 & K <= 4 ~ "Observed in 3–4 sites",
    K >= 5 & K <= 9 ~ "Observed in 5–9 sites",
    K >= 10 ~ "Observed in 10 or more sites",
    TRUE ~ NA_character_
  )
}

realised_links <- all_pairs %>%
  filter(full_K >= 1) %>%
  mutate(support_class = factor(classify_support(full_K), levels = support_levels))

support_by_dataset <- realised_links %>%
  count(dataset, support_class, name = "number_of_links") %>%
  complete(dataset, support_class, fill = list(number_of_links = 0L)) %>%
  group_by(dataset) %>%
  mutate(
    total_links = sum(number_of_links),
    fraction_of_links = if_else(total_links > 0, number_of_links / total_links, NA_real_)
  ) %>%
  ungroup()

across_landscapes <- support_by_dataset %>%
  group_by(support_class) %>%
  summarise(
    fraction_of_links = mean(fraction_of_links, na.rm = TRUE),
    number_of_links = NA_integer_,
    total_links = NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(dataset = "Across landscapes")

sort_by_one_site <- support_by_dataset %>%
  filter(support_class == "Observed in 1 site") %>%
  arrange(desc(fraction_of_links)) %>%
  pull(dataset) %>%
  as.character()

support_distribution <- bind_rows(support_by_dataset, across_landscapes) %>%
  mutate(
    dataset = factor(dataset, levels = rev(c("Across landscapes", sort_by_one_site))),
    support_class = factor(support_class, levels = support_levels)
  )

write.csv2(
  support_distribution,
  file.path(out_dir, "31_support_distribution_by_dataset.csv"),
  row.names = FALSE
)

panel_a <- ggplot(
  support_distribution,
  aes(x = dataset, y = fraction_of_links, fill = support_class)
) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = support_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  base_theme +
  theme(axis.title.y = element_blank()) +
  ylab("Share of regional interaction links") +
  ggtitle("A. Regional interaction links differ greatly in how many sites support them")

## ---------------------------
## Panel B: dataset-balanced heatmap, all co-occurring pairs
## ---------------------------

bin_full_n <- function(x){
  cut(
    x,
    breaks = c(0, 1, 2, 4, 8, 16, 32, Inf),
    labels = c("1", "2", "3–4", "5–8", "9–16", "17–32", "33+"),
    right = TRUE
  )
}

bin_full_K <- function(x){
  cut(
    x,
    breaks = c(-1, 0, 1, 2, 4, 8, 16, 32, Inf),
    labels = c("0", "1", "2", "3–4", "5–8", "9–16", "17–32", "33+"),
    right = TRUE
  )
}

heatmap_data <- all_pairs %>%
  group_by(dataset) %>%
  mutate(dataset_weight = 1 / n()) %>%
  ungroup() %>%
  mutate(
    full_n_bin = bin_full_n(full_n),
    full_K_bin = bin_full_K(full_K)
  ) %>%
  group_by(full_n_bin, full_K_bin) %>%
  summarise(dataset_balanced_count = sum(dataset_weight), .groups = "drop") %>%
  mutate(dataset_balanced_fraction = dataset_balanced_count / sum(dataset_balanced_count))

write.csv2(
  heatmap_data,
  file.path(out_dir, "31_dataset_balanced_pair_support_heatmap.csv"),
  row.names = FALSE
)

panel_b <- ggplot(heatmap_data, aes(x = full_n_bin, y = full_K_bin, fill = dataset_balanced_fraction)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08519C",
    labels = percent_format(accuracy = 0.1),
    name = "Dataset-balanced\nshare"
  ) +
  base_theme +
  xlab("Sites where pair is recorded together") +
  ylab("Sites where interaction is observed") +
  ggtitle(
    "B. Where species meet does not fully determine where interactions are sustained",
    subtitle = "All co-occurring pairs are included, including pairs never observed interacting."
  )

fig2 <- panel_a / panel_b +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(title = "Empirical interaction links differ sharply in local support")

save_both(fig2, "31_Figure2_empirical_support_heterogeneity", width = 12, height = 10)

## ---------------------------
## EXTRA: raw 2 x 5 point panel like Galiana-style dataset panels
## ---------------------------

## This is the plot the user asked for: no dataset weighting, no heatmap.
## One point is one full-network co-occurring consumer-resource pair.
## It includes pairs with full_K = 0.

raw_point_data <- all_pairs %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    ever_interacted = if_else(full_K > 0, "Observed at least once", "Never observed interacting")
  )

write.csv2(
  raw_point_data,
  file.path(out_dir, "31_raw_pair_support_by_dataset_points.csv"),
  row.names = FALSE
)

library(ggh4x)
raw_point_plot <- ggplot(raw_point_data, aes(x = full_n, y = full_K)) +
  geom_point(size = 0.55, alpha = 0.35, colour = "black") +
  facet_wrap2(
    ~ dataset,
    ncol = 5,
    scales = "free_x",
    axes = "all"
  ) +
  base_theme +
  theme(legend.position = "none") +
  xlab("Sites where pair is recorded together") +
  ylab("Sites where interaction is observed") +
  ggtitle(
    "Raw pair-level support across datasets",
    subtitle = "Each point is one co-occurring consumer–resource pair; pairs with no observed interaction are retained at zero."
  )

save_both(raw_point_plot, "31_Figure2B_raw_pair_support_by_dataset", width = 12, height = 6.5)

## A jittered version is saved separately only to make overlapping integer points visible.
## The raw values remain unchanged; jitter is graphical only.
set.seed(123)
raw_point_plot_jittered <- ggplot(raw_point_data, aes(x = full_n, y = full_K)) +
  geom_jitter(width = 0.12, height = 0.12, size = 0.45, alpha = 0.28, colour = "black") +
  facet_wrap(~ dataset, ncol = 5, scales = "free") +
  base_theme +
  theme(legend.position = "none") +
  xlab("Sites where pair is recorded together") +
  ylab("Sites where interaction is observed") +
  ggtitle(
    "Raw pair-level support across datasets, jittered for visibility",
    subtitle = "Jitter is graphical only; each point is one co-occurring consumer–resource pair."
  )

save_both(raw_point_plot_jittered, "31_Figure2B_raw_pair_support_by_dataset_jittered", width = 12, height = 6.5)

## Optional log-scale raw view for large datasets. This keeps K = 0 visible using log1p transformation.
raw_point_plot_log1p <- ggplot(raw_point_data, aes(x = full_n, y = full_K)) +
  geom_point(size = 0.5, alpha = 0.30, colour = "black") +
  facet_wrap(~ dataset, ncol = 5, scales = "free") +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  base_theme +
  theme(legend.position = "none") +
  xlab("Sites where pair is recorded together") +
  ylab("Sites where interaction is observed") +
  ggtitle(
    "Raw pair-level support across datasets, log1p axes",
    subtitle = "The zero row is retained; axes are transformed only for display."
  )

save_both(raw_point_plot_log1p, "31_Figure2B_raw_pair_support_by_dataset_log1p", width = 12, height = 6.5)

## ---------------------------
## Validation checks
## ---------------------------

checks <- tibble(
  check = c(
    "Figure 2B includes all co-occurring pairs with full_K = 0",
    "Raw 2x5 panel includes full_K = 0 pairs",
    "All pairs have full_n >= 1",
    "All pairs have full_K >= 0",
    "All pairs have full_K <= full_n",
    "Figure 2A uses only realised links with full_K >= 1",
    "Dataset-balanced heatmap gives each dataset total weight 1"
  ),
  pass = c(
    any(all_pairs$full_K == 0),
    any(raw_point_data$full_K == 0),
    all(all_pairs$full_n >= 1),
    all(all_pairs$full_K >= 0),
    all(all_pairs$full_K <= all_pairs$full_n),
    all(realised_links$full_K >= 1),
    all(abs((all_pairs %>% group_by(dataset) %>% summarise(w = sum(1 / n()), .groups = "drop"))$w - 1) < 1e-8)
  )
)

write.csv2(
  checks,
  file.path(out_dir, "31_Figure2_checks.csv"),
  row.names = FALSE
)

print(checks)
message("Saved Figure 2 outputs in: ", out_dir)
message("Extra raw point panels saved as:")
message("  - 31_Figure2B_raw_pair_support_by_dataset.png/pdf")
message("  - 31_Figure2B_raw_pair_support_by_dataset_jittered.png/pdf")
message("  - 31_Figure2B_raw_pair_support_by_dataset_log1p.png/pdf")
