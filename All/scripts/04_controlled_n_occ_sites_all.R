
## ------------------------------------------------------------
## Script: All/scripts/04_controlled_n_occ_sites_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

n_model_reps <- 1000
min_empirical_links <- 5

dirs <- make_output_dirs("controlled_n_occ_sites")
sep_out <- dirs$separated
combined_out <- dirs$combined

summarise_repeatability <- function(df, group_var){
  df %>%
    group_by({{ group_var }}) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      .groups = "drop"
    )
}

simulate_model_pairs <- function(cooc_counts, p_fixed, model_rep){
  simulated <- cooc_counts

  simulated$n_interacting_sites <- rbinom(
    n = nrow(simulated),
    size = simulated$n_cooccurring_sites,
    prob = p_fixed
  )

  simulated %>%
    filter(n_interacting_sites > 0) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      model_rep = model_rep,
      cooc_bin = make_nocc_bin(n_cooccurring_sites)
    )
}

make_envelope <- function(model_summary, group_var){
  model_summary %>%
    group_by({{ group_var }}) %>%
    summarise(
      model_n_links_median = median(n_links, na.rm = TRUE),

      model_mean_interacting_sites_q025 = quantile(mean_interacting_sites, 0.025, na.rm = TRUE),
      model_mean_interacting_sites_q500 = quantile(mean_interacting_sites, 0.500, na.rm = TRUE),
      model_mean_interacting_sites_q975 = quantile(mean_interacting_sites, 0.975, na.rm = TRUE),

      model_mean_repeatability_q025 = quantile(mean_repeatability, 0.025, na.rm = TRUE),
      model_mean_repeatability_q500 = quantile(mean_repeatability, 0.500, na.rm = TRUE),
      model_mean_repeatability_q975 = quantile(mean_repeatability, 0.975, na.rm = TRUE),

      model_prop_single_q025 = quantile(proportion_single_site_links, 0.025, na.rm = TRUE),
      model_prop_single_q500 = quantile(proportion_single_site_links, 0.500, na.rm = TRUE),
      model_prop_single_q975 = quantile(proportion_single_site_links, 0.975, na.rm = TRUE),

      .groups = "drop"
    )
}

run_one_controlled <- function(dataset){
  message("Running controlled n-occ-sites: ", dataset)

  p_fixed <- get_p_fixed(dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions

  cooc_counts <- cooc_triples %>%
    group_by(consumer, resource) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

  empirical_pairs <- empirical_site_interactions %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      cooc_bin = make_nocc_bin(n_cooccurring_sites),
      dataset = dataset
    )

  empirical_exact <- summarise_repeatability(empirical_pairs, n_cooccurring_sites) %>%
    rename(
      empirical_n_links = n_links,
      empirical_mean_interacting_sites = mean_interacting_sites,
      empirical_median_interacting_sites = median_interacting_sites,
      empirical_q95_interacting_sites = q95_interacting_sites,
      empirical_mean_repeatability = mean_repeatability,
      empirical_median_repeatability = median_repeatability,
      empirical_proportion_single_site_links = proportion_single_site_links
    )

  empirical_binned <- summarise_repeatability(empirical_pairs, cooc_bin) %>%
    rename(
      empirical_n_links = n_links,
      empirical_mean_interacting_sites = mean_interacting_sites,
      empirical_median_interacting_sites = median_interacting_sites,
      empirical_q95_interacting_sites = q95_interacting_sites,
      empirical_mean_repeatability = mean_repeatability,
      empirical_median_repeatability = median_repeatability,
      empirical_proportion_single_site_links = proportion_single_site_links
    )

  model_pairs_list <- vector("list", n_model_reps)

  for(r in seq_len(n_model_reps)){
    message(dataset, ": controlled simulation ", r, " / ", n_model_reps)

    model_pairs_list[[r]] <- simulate_model_pairs(
      cooc_counts = cooc_counts,
      p_fixed = p_fixed,
      model_rep = r
    ) %>%
      mutate(dataset = dataset)
  }

  model_pairs <- bind_rows(model_pairs_list)

  model_exact_by_rep <- model_pairs %>%
    group_by(model_rep, n_cooccurring_sites) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      .groups = "drop"
    )

  model_binned_by_rep <- model_pairs %>%
    group_by(model_rep, cooc_bin) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      .groups = "drop"
    )

  model_exact_envelope <- make_envelope(model_exact_by_rep, n_cooccurring_sites)
  model_binned_envelope <- make_envelope(model_binned_by_rep, cooc_bin)

  comparison_exact <- empirical_exact %>%
    left_join(model_exact_envelope, by = "n_cooccurring_sites") %>%
    mutate(
      dataset = dataset,
      enough_empirical_links = empirical_n_links >= min_empirical_links,
      empirical_repeatability_above_model = empirical_mean_repeatability > model_mean_repeatability_q975,
      empirical_repeatability_below_model = empirical_mean_repeatability < model_mean_repeatability_q025,
      empirical_interacting_sites_above_model = empirical_mean_interacting_sites > model_mean_interacting_sites_q975,
      empirical_interacting_sites_below_model = empirical_mean_interacting_sites < model_mean_interacting_sites_q025,
      empirical_single_site_below_model = empirical_proportion_single_site_links < model_prop_single_q025,
      empirical_single_site_above_model = empirical_proportion_single_site_links > model_prop_single_q975
    )

  comparison_binned <- empirical_binned %>%
    left_join(model_binned_envelope, by = "cooc_bin") %>%
    mutate(
      dataset = dataset,
      enough_empirical_links = empirical_n_links >= min_empirical_links,
      empirical_repeatability_above_model = empirical_mean_repeatability > model_mean_repeatability_q975,
      empirical_repeatability_below_model = empirical_mean_repeatability < model_mean_repeatability_q025,
      empirical_interacting_sites_above_model = empirical_mean_interacting_sites > model_mean_interacting_sites_q975,
      empirical_interacting_sites_below_model = empirical_mean_interacting_sites < model_mean_interacting_sites_q025,
      empirical_single_site_below_model = empirical_proportion_single_site_links < model_prop_single_q025,
      empirical_single_site_above_model = empirical_proportion_single_site_links > model_prop_single_q975
    )

  write.csv2(empirical_pairs,
            file.path(out_dir, paste0(dataset, "_empirical_pair_repeatability.csv")),
            row.names = FALSE)

  write.csv2(model_pairs,
            file.path(out_dir, paste0(dataset, "_model_generated_pair_repeatability.csv")),
            row.names = FALSE)

  write.csv2(comparison_exact,
            file.path(out_dir, paste0(dataset, "_controlled_exact_n_cooccurring_sites.csv")),
            row.names = FALSE)

  write.csv2(comparison_binned,
            file.path(out_dir, paste0(dataset, "_controlled_binned_n_cooccurring_sites.csv")),
            row.names = FALSE)

  plot_exact_data <- comparison_exact %>%
    filter(enough_empirical_links)

  p_exact_repeatability <- ggplot(plot_exact_data, aes(x = n_cooccurring_sites)) +
    geom_ribbon(aes(ymin = model_mean_repeatability_q025,
                    ymax = model_mean_repeatability_q975),
                fill = "grey80") +
    geom_line(aes(y = model_mean_repeatability_q500),
              linetype = 2, linewidth = 1) +
    geom_point(aes(y = empirical_mean_repeatability,
                   size = empirical_n_links)) +
    geom_line(aes(y = empirical_mean_repeatability), linewidth = 1.1) +
    theme_classic(base_size = 14) +
    xlab("Number of co-occurring sites") +
    ylab("Mean repeatability") +
    ggtitle(paste0(dataset, ": repeatability controlled for co-occurring sites"),
            subtitle = "Ribbon = 95% model envelope; points/solid line = empirical")

  ggsave(file.path(out_dir, paste0(dataset, "_exact_controlled_mean_repeatability.png")),
         p_exact_repeatability, width = 8, height = 5.5, dpi = 300)

  p_binned_repeatability <- ggplot(comparison_binned, aes(x = cooc_bin)) +
    geom_errorbar(aes(ymin = model_mean_repeatability_q025,
                      ymax = model_mean_repeatability_q975),
                  width = 0.15, linewidth = 1) +
    geom_point(aes(y = model_mean_repeatability_q500), shape = 1, size = 3) +
    geom_point(aes(y = empirical_mean_repeatability,
                   size = empirical_n_links), size = 3) +
    theme_classic(base_size = 14) +
    xlab("Number of co-occurring sites, binned") +
    ylab("Mean repeatability") +
    ggtitle(paste0(dataset, ": binned control for co-occurring sites"),
            subtitle = "Error bars = 95% model envelope; filled points = empirical")

  ggsave(file.path(out_dir, paste0(dataset, "_binned_controlled_mean_repeatability.png")),
         p_binned_repeatability, width = 8, height = 5.5, dpi = 300)

  list(exact = comparison_exact, binned = comparison_binned)
}

all_outputs <- lapply(all_dataset_names, run_one_controlled)
names(all_outputs) <- all_dataset_names

combined_exact <- bind_rows(lapply(all_outputs, `[[`, "exact"))
combined_binned <- bind_rows(lapply(all_outputs, `[[`, "binned"))

write.csv2(combined_exact,
          file.path(combined_out, "controlled_exact_n_cooccurring_sites_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_binned,
          file.path(combined_out, "controlled_binned_n_cooccurring_sites_all_datasets.csv"),
          row.names = FALSE)

combined_binned$dataset <- factor(combined_binned$dataset, levels = all_dataset_names)

p_combined_binned <- ggplot(combined_binned, aes(x = cooc_bin)) +
  geom_errorbar(aes(ymin = model_mean_repeatability_q025,
                    ymax = model_mean_repeatability_q975),
                width = 0.15) +
  geom_point(aes(y = model_mean_repeatability_q500), shape = 1, size = 2) +
  geom_point(aes(y = empirical_mean_repeatability), size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Co-occurring sites, binned") +
  ylab("Mean repeatability") +
  ggtitle("Repeatability controlled for co-occurring sites across all datasets",
          subtitle = "Error bars = 95% model envelope; open points = model median; filled points = empirical")

ggsave(file.path(combined_out, "combined_binned_controlled_repeatability.png"),
       p_combined_binned, width = 14, height = 7, dpi = 300)

message("Finished controlled n-occurring-sites analysis.")
