## ------------------------------------------------------------
## Script: All/scripts/08_repeatability_heterogeneity_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script8_repeatability_heterogeneity")
sep_out <- dirs$separated
combined_out <- dirs$combined

pair_level_for_subset <- function(dataset, site_tables, sites_keep, subset_row){
  cooc_counts <- site_tables$cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")
  
  int_counts <- site_tables$empirical_site_interactions %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop")
  
  pair_df <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      n_interacting_sites = replace_na(n_interacting_sites, 0L),
      pair_interaction_rate = n_interacting_sites / n_cooccurring_sites,
      realised_link = n_interacting_sites > 0,
      repeatability = ifelse(realised_link, pair_interaction_rate, NA_real_),
      dataset = dataset,
      subset_id = subset_row$subset_id,
      removal_fraction = subset_row$removal_fraction,
      site_rep = subset_row$site_rep,
      n_sites_kept = subset_row$n_sites_kept
    )
  
  realised <- pair_df %>% filter(realised_link)
  
  summary <- data.frame(
    dataset = dataset,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_kept = subset_row$n_sites_kept,
    n_cooc_pairs = nrow(pair_df),
    n_realised_links = nrow(realised),
    mean_repeatability = mean(realised$repeatability, na.rm = TRUE),
    median_repeatability = median(realised$repeatability, na.rm = TRUE),
    q95_repeatability = quantile(realised$repeatability, 0.95, na.rm = TRUE),
    mean_interacting_sites_realised = mean(realised$n_interacting_sites, na.rm = TRUE),
    median_interacting_sites_realised = median(realised$n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites_realised = quantile(realised$n_interacting_sites, 0.95, na.rm = TRUE),
    proportion_one_site_realised_links = mean(realised$n_interacting_sites == 1, na.rm = TRUE),
    gini_interacting_sites_realised = gini_coefficient(realised$n_interacting_sites),
    variance_pair_interaction_rate = var(pair_df$pair_interaction_rate, na.rm = TRUE),
    cv_pair_interaction_rate = sd(pair_df$pair_interaction_rate, na.rm = TRUE) / mean(pair_df$pair_interaction_rate, na.rm = TRUE),
    proportion_zero_interaction_pairs = mean(pair_df$pair_interaction_rate == 0, na.rm = TRUE),
    proportion_high_interaction_pairs = mean(pair_df$pair_interaction_rate > 0.5, na.rm = TRUE)
  )
  
  list(pair_level = pair_df, summary = summary)
}

run_one_dataset <- function(dataset){
  message("Running script 8: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  site_tables <- get_dataset_site_tables(dataset)
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  
  outs <- lapply(seq_len(nrow(subset_object$index)), function(i){
    pair_level_for_subset(dataset, site_tables, subset_object$subsets[[i]], subset_object$index[i,])
  })
  
  pair_level <- bind_rows(lapply(outs, `[[`, "pair_level"))
  summary <- bind_rows(lapply(outs, `[[`, "summary"))
  
  summary_by_removal <- summary %>%
    group_by(dataset, removal_fraction) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  write.csv2(pair_level, file.path(out_dir, paste0(dataset, "_script8_pair_level_repeatability.csv")), row.names = FALSE)
  write.csv2(summary, file.path(out_dir, paste0(dataset, "_script8_repeatability_heterogeneity_summary.csv")), row.names = FALSE)
  write.csv2(summary_by_removal, file.path(out_dir, paste0(dataset, "_script8_repeatability_heterogeneity_summary_by_removal.csv")), row.names = FALSE)
  list(pair_level = pair_level, summary = summary, summary_by_removal = summary_by_removal)
}

outs <- lapply(all_dataset_names, run_one_dataset)
names(outs) <- all_dataset_names
combined_pair <- bind_rows(lapply(outs, `[[`, "pair_level"))
combined <- bind_rows(lapply(outs, `[[`, "summary"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary_by_removal"))

write.csv2(combined_pair, file.path(combined_out, "script8_pair_level_repeatability_combined.csv"), row.names = FALSE)
write.csv2(combined, file.path(combined_out, "script8_repeatability_heterogeneity_combined.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "script8_repeatability_heterogeneity_summary_combined.csv"), row.names = FALSE)

combined_summary$dataset <- factor(combined_summary$dataset, levels = all_dataset_names)

plot_metric <- function(metric, ylab, filename){
  p <- ggplot(combined_summary, aes(removal_fraction, .data[[metric]])) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~dataset, ncol = 5, scales = "free_y") +
    theme_classic(base_size = 11) +
    xlab("Fraction of sites removed") + ylab(ylab) +
    ggtitle(ylab)
  ggsave(file.path(combined_out, filename), p, width = 14, height = 7, dpi = 300)
}

plot_metric("mean_repeatability", "Mean repeatability", "08_mean_repeatability_vs_site_removal.png")
plot_metric("proportion_one_site_realised_links", "Proportion of one-site realised links", "08_one_site_links_vs_site_removal.png")
plot_metric("variance_pair_interaction_rate", "Pair-level interaction-rate variance", "08_interaction_rate_variance_vs_site_removal.png")
plot_metric("proportion_zero_interaction_pairs", "Proportion of zero-interaction co-occurring pairs", "08_zero_interaction_pairs_vs_site_removal.png")

message("Finished script 8.")
