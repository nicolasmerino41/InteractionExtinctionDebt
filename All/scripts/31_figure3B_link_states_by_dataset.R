## ------------------------------------------------------------
## Script: All/scripts/31_figure3B_link_states_by_dataset.R
## Extra Figure 3B: original interaction-link states by dataset.
## ------------------------------------------------------------
source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "patchwork", "scales", "future", "future.apply")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Fast mode is useful for checking layout. Set FALSE for final run.
fast_mode <- FALSE
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- if(fast_mode) 50 else 500
set.seed(123)

n_workers <- max(1, min(7, future::availableCores() - 1))
future::plan(future::multisession, workers = n_workers)

save_both <- function(p, name, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

state_cols <- c(
  "Interaction still present" = "#009E73",
  "Species still recorded together, interaction not observed" = "#CC79A7",
  "Species no longer recorded together" = "#999999"
)

state_levels <- names(state_cols)

make_pair_tables <- function(site_tables){
  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource)

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource)

  ## Defensive consistency: an observed interaction is also a recorded-together cell.
  cooc <- bind_rows(cooc, ints) %>%
    distinct(site, consumer, resource)

  full_pairs <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop") %>%
    left_join(
      ints %>%
        group_by(consumer, resource) %>%
        summarise(full_K = n_distinct(site), .groups = "drop"),
      by = c("consumer", "resource")
    ) %>%
    mutate(
      full_K = replace_na(full_K, 0L),
      full_B = full_K > 0
    )

  original_links <- full_pairs %>%
    filter(full_K >= 1) %>%
    select(consumer, resource, full_n, full_K)

  list(cooc = cooc, ints = ints, full_pairs = full_pairs, original_links = original_links)
}

summarise_one_subset <- function(dataset, pair_tables, sites_keep, removal_fraction, replicate_id){
  original_links <- pair_tables$original_links
  n_original_links <- nrow(original_links)

  retained_n <- pair_tables$cooc %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(retained_n = n_distinct(site), .groups = "drop")

  retained_K <- pair_tables$ints %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(retained_K = n_distinct(site), .groups = "drop")

  link_states <- original_links %>%
    left_join(retained_n, by = c("consumer", "resource")) %>%
    left_join(retained_K, by = c("consumer", "resource")) %>%
    mutate(
      retained_n = replace_na(retained_n, 0L),
      retained_K = replace_na(retained_K, 0L),
      link_state = case_when(
        retained_K > 0 ~ "Interaction still present",
        retained_K == 0 & retained_n > 0 ~ "Species still recorded together, interaction not observed",
        retained_K == 0 & retained_n == 0 ~ "Species no longer recorded together"
      )
    )

  counts <- link_states %>%
    count(link_state, name = "n_links") %>%
    right_join(tibble(link_state = state_levels), by = "link_state") %>%
    mutate(
      n_links = replace_na(n_links, 0L),
      dataset = dataset,
      removal_fraction = removal_fraction,
      replicate = replicate_id,
      n_original_links = n_original_links,
      fraction = n_links / n_original_links
    ) %>%
    select(dataset, removal_fraction, replicate, link_state, n_original_links, n_links, fraction)

  counts
}

run_one_dataset <- function(dataset){
  message("Running extra Figure 3B by dataset: ", dataset)

  site_tables <- get_dataset_site_tables(dataset)
  pair_tables <- make_pair_tables(site_tables)

  all_sites <- sort(unique(pair_tables$cooc$site))
  site_subsets <- make_site_subsets(all_sites, removal_levels, n_site_reps)

  out <- future.apply::future_lapply(seq_len(nrow(site_subsets$index)), function(i){
    idx <- site_subsets$index[i, ]
    summarise_one_subset(
      dataset = dataset,
      pair_tables = pair_tables,
      sites_keep = site_subsets$subsets[[i]],
      removal_fraction = idx$removal_fraction,
      replicate_id = idx$site_rep
    )
  }, future.seed = TRUE)

  bind_rows(out)
}

all_states <- bind_rows(future.apply::future_lapply(all_dataset_names, run_one_dataset, future.seed = TRUE)) %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    link_state = factor(link_state, levels = state_levels)
  )

write.csv2(
  all_states,
  file.path(out_dir, "31_Figure3B_original_link_states_by_dataset_replicates.csv"),
  row.names = FALSE
)

plot_states <- all_states %>%
  group_by(dataset, removal_fraction, link_state) %>%
  summarise(
    fraction = mean(fraction, na.rm = TRUE),
    q025 = quantile(fraction, 0.025, na.rm = TRUE),
    q975 = quantile(fraction, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write.csv2(
  plot_states,
  file.path(out_dir, "31_Figure3B_original_link_states_by_dataset_summary.csv"),
  row.names = FALSE
)

## Stacked-area panel for each dataset.
## Uses mean fractions across replicates, so each panel sums to 1 at every removal level.
p_states_by_dataset <- ggplot(plot_states, aes(
  x = removal_fraction,
  y = fraction,
  fill = link_state
)) +
  geom_area(colour = "white", linewidth = 0.15, alpha = 0.95) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_x_continuous(breaks = removal_levels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  xlab("Proportion of sites removed") +
  ylab("Share of original interactions") +
  ggtitle(
    "What happens to original interactions in each dataset?",
    subtitle = "Each panel shows the mean fraction of original regional interaction links in each state after random site removal."
  )

save_both(
  p_states_by_dataset,
  "31_Figure3B_link_states_by_dataset",
  width = 13,
  height = 7.5
)

checks <- all_states %>%
  group_by(dataset, removal_fraction, replicate) %>%
  summarise(
    state_sum = sum(fraction, na.rm = TRUE),
    n_states = n_distinct(link_state),
    .groups = "drop"
  ) %>%
  summarise(
    all_state_sums_equal_1 = all(abs(state_sum - 1) < 1e-10),
    all_three_states_present_in_table = all(n_states == 3),
    n_dataset_removal_replicates = n(),
    .groups = "drop"
  )

write.csv2(
  checks,
  file.path(out_dir, "31_Figure3B_link_states_by_dataset_checks.csv"),
  row.names = FALSE
)

print(checks)
message("Saved extra Figure 3B by-dataset outputs in: ", out_dir)

future::plan(future::sequential)
