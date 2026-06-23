
## ------------------------------------------------------------
## Script: All/scripts/02_model_envelope_site_removal_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100
n_model_reps <- 100

dirs <- make_output_dirs("model_envelope_site_removal")
sep_out <- dirs$separated
combined_out <- dirs$combined

evaluate_network_against_subsets <- function(site_interactions, subset_object, subset_expected, source, model_rep = NA_integer_){
  out <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    exp_row <- subset_expected[subset_expected$subset_id == subset_id, ]

    observed_links <- observed_links_from_subset(site_interactions, sites_keep)

    out[[i]] <- data.frame(
      source = source,
      model_rep = model_rep,
      subset_id = subset_id,
      removal_fraction = subset_object$index$removal_fraction[i],
      site_rep = subset_object$index$site_rep[i],
      n_sites_kept = subset_object$index$n_sites_kept[i],
      n_observed_links = observed_links,
      n_cooc_links = exp_row$n_cooc_links,
      expected_links_raw = exp_row$expected_links_raw,
      expected_links_conditioned = exp_row$expected_links_conditioned,
      divergence_raw = observed_links - exp_row$expected_links_raw,
      relative_divergence_raw = safe_relative_difference(observed_links, exp_row$expected_links_raw),
      divergence_conditioned = observed_links - exp_row$expected_links_conditioned,
      relative_divergence_conditioned = safe_relative_difference(observed_links, exp_row$expected_links_conditioned)
    )
  }

  bind_rows(out)
}

run_one_envelope <- function(dataset){
  message("Running envelope: ", dataset)

  p_fixed <- get_p_fixed(dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions
  all_sites <- sort(unique(cooc_triples$site))

  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)

  subset_expected <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]

    cur_exp <- expected_links_from_subset(cooc_triples, sites_keep, p_fixed)
    cur_exp$subset_id <- subset_id
    cur_exp
  }))

  empirical_results <- evaluate_network_against_subsets(
    site_interactions = empirical_site_interactions,
    subset_object = subset_object,
    subset_expected = subset_expected,
    source = "empirical"
  ) %>%
    mutate(dataset = dataset)

  model_results_list <- vector("list", n_model_reps)

  for(m in seq_len(n_model_reps)){
    message(dataset, ": model replicate ", m, " / ", n_model_reps)

    synthetic_site_interactions <- simulate_model_site_interactions(cooc_triples, p_fixed)

    model_results_list[[m]] <- evaluate_network_against_subsets(
      site_interactions = synthetic_site_interactions,
      subset_object = subset_object,
      subset_expected = subset_expected,
      source = "model_generated",
      model_rep = m
    ) %>%
      mutate(dataset = dataset)
  }

  model_results <- bind_rows(model_results_list)

  empirical_summary <- empirical_results %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      empirical_mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
      empirical_sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
      empirical_mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
      empirical_mean_observed_links = mean(n_observed_links, na.rm = TRUE),
      empirical_mean_expected_links_raw = mean(expected_links_raw, na.rm = TRUE),
      .groups = "drop"
    )

  model_mean_curves <- model_results %>%
    group_by(dataset, model_rep, removal_fraction) %>%
    summarise(
      model_mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
      model_mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
      .groups = "drop"
    )

  model_envelope <- model_mean_curves %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      model_q025_relative_divergence_raw = quantile(model_mean_relative_divergence_raw, 0.025, na.rm = TRUE),
      model_q500_relative_divergence_raw = quantile(model_mean_relative_divergence_raw, 0.500, na.rm = TRUE),
      model_q975_relative_divergence_raw = quantile(model_mean_relative_divergence_raw, 0.975, na.rm = TRUE),
      model_mean_relative_divergence_raw = mean(model_mean_relative_divergence_raw, na.rm = TRUE),
      .groups = "drop"
    )

  comparison_summary <- left_join(empirical_summary, model_envelope,
                                  by = c("dataset", "removal_fraction")) %>%
    mutate(
      empirical_above_model_975 = empirical_mean_relative_divergence_raw > model_q975_relative_divergence_raw,
      empirical_below_model_025 = empirical_mean_relative_divergence_raw < model_q025_relative_divergence_raw
    )

  write.csv2(empirical_results,
            file.path(out_dir, paste0(dataset, "_empirical_divergence_by_subset.csv")),
            row.names = FALSE)

  write.csv2(model_results,
            file.path(out_dir, paste0(dataset, "_model_generated_divergence_by_subset.csv")),
            row.names = FALSE)

  write.csv2(comparison_summary,
            file.path(out_dir, paste0(dataset, "_empirical_vs_model_envelope_summary.csv")),
            row.names = FALSE)

  p <- ggplot() +
    geom_ribbon(data = model_envelope,
                aes(x = removal_fraction,
                    ymin = model_q025_relative_divergence_raw,
                    ymax = model_q975_relative_divergence_raw),
                fill = "grey80") +
    geom_line(data = model_envelope,
              aes(x = removal_fraction, y = model_q500_relative_divergence_raw),
              linewidth = 1, linetype = 2) +
    geom_line(data = empirical_summary,
              aes(x = removal_fraction, y = empirical_mean_relative_divergence_raw),
              linewidth = 1.2) +
    geom_point(data = empirical_summary,
               aes(x = removal_fraction, y = empirical_mean_relative_divergence_raw),
               size = 3) +
    geom_hline(yintercept = 0, linetype = 3) +
    theme_classic(base_size = 14) +
    xlab("Fraction of sites removed") +
    ylab("Mean relative divergence") +
    ggtitle(paste0(dataset, ": empirical divergence vs model-generated envelope"),
            subtitle = "Ribbon = 95% envelope of model-generated mean curves; solid = empirical")

  ggsave(file.path(out_dir, paste0(dataset, "_empirical_vs_model_envelope.png")),
         p, width = 8, height = 5.5, dpi = 300)

  list(empirical = empirical_results,
       model = model_results,
       comparison = comparison_summary,
       envelope = model_envelope)
}

all_outputs <- lapply(all_dataset_names, run_one_envelope)
names(all_outputs) <- all_dataset_names

combined_comparison <- bind_rows(lapply(all_outputs, `[[`, "comparison"))
combined_model_envelope <- bind_rows(lapply(all_outputs, `[[`, "envelope"))

write.csv2(combined_comparison,
          file.path(combined_out, "model_envelope_summary_all_datasets.csv"),
          row.names = FALSE)

combined_comparison$dataset <- factor(combined_comparison$dataset, levels = all_dataset_names)

p_combined <- ggplot(combined_comparison, aes(x = removal_fraction)) +
  geom_ribbon(aes(ymin = model_q025_relative_divergence_raw,
                  ymax = model_q975_relative_divergence_raw),
              fill = "grey80") +
  geom_line(aes(y = model_q500_relative_divergence_raw), linetype = 2, linewidth = 0.9) +
  geom_line(aes(y = empirical_mean_relative_divergence_raw), linewidth = 1) +
  geom_point(aes(y = empirical_mean_relative_divergence_raw), size = 1.8) +
  geom_hline(yintercept = 0, linetype = 3) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 12) +
  xlab("Fraction of sites removed") +
  ylab("Mean relative divergence") +
  ggtitle("Empirical divergence vs model-generated envelope across all datasets",
          subtitle = "Raw expected links; ribbon = 95% model envelope; solid = empirical")

ggsave(file.path(combined_out, "combined_model_envelope_relative_divergence.png"),
       p_combined, width = 14, height = 7, dpi = 300)

message("Finished model envelope analysis.")
