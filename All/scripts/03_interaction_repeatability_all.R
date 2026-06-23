
## ------------------------------------------------------------
## Script: All/scripts/03_interaction_repeatability_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

n_model_reps <- 1000

dirs <- make_output_dirs("interaction_repeatability")
sep_out <- dirs$separated
combined_out <- dirs$combined

simulate_model_repeatability <- function(cooc_triples, cooc_counts, p_fixed, model_rep){
  sim <- cooc_triples
  sim$interaction <- rbinom(nrow(sim), size = 1, prob = p_fixed)

  sim %>%
    filter(interaction == 1) %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      source = "model_generated",
      model_rep = model_rep
    )
}

run_one_repeatability <- function(dataset){
  message("Running repeatability: ", dataset)

  p_fixed <- get_p_fixed(dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions

  cooc_counts <- cooc_triples %>%
    group_by(consumer, resource) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

  empirical_repeatability <- empirical_site_interactions %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      source = "empirical",
      model_rep = NA_integer_,
      dataset = dataset
    )

  model_repeatability_list <- vector("list", n_model_reps)

  for(r in seq_len(n_model_reps)){
    message(dataset, ": model replicate ", r, " / ", n_model_reps)

    model_repeatability_list[[r]] <- simulate_model_repeatability(
      cooc_triples = cooc_triples,
      cooc_counts = cooc_counts,
      p_fixed = p_fixed,
      model_rep = r
    ) %>%
      mutate(dataset = dataset)
  }

  model_repeatability <- bind_rows(model_repeatability_list)
  all_repeatability <- bind_rows(empirical_repeatability, model_repeatability)

  summary_by_source <- all_repeatability %>%
    group_by(dataset, source) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      .groups = "drop"
    )

  model_summary_by_rep <- model_repeatability %>%
    group_by(dataset, model_rep) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      .groups = "drop"
    )

  empirical_summary <- empirical_repeatability %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE)
    )

  comparison_to_model <- data.frame(
    dataset = dataset,
    metric = c(
      "mean_interacting_sites", "median_interacting_sites",
      "q95_interacting_sites", "max_interacting_sites",
      "proportion_single_site_links", "mean_repeatability",
      "median_repeatability"
    ),
    empirical_value = c(
      empirical_summary$mean_interacting_sites,
      empirical_summary$median_interacting_sites,
      empirical_summary$q95_interacting_sites,
      empirical_summary$max_interacting_sites,
      empirical_summary$proportion_single_site_links,
      empirical_summary$mean_repeatability,
      empirical_summary$median_repeatability
    ),
    model_q025 = c(
      quantile(model_summary_by_rep$mean_interacting_sites, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$median_interacting_sites, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$q95_interacting_sites, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$max_interacting_sites, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$proportion_single_site_links, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$mean_repeatability, 0.025, na.rm = TRUE),
      quantile(model_summary_by_rep$median_repeatability, 0.025, na.rm = TRUE)
    ),
    model_q500 = c(
      quantile(model_summary_by_rep$mean_interacting_sites, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$median_interacting_sites, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$q95_interacting_sites, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$max_interacting_sites, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$proportion_single_site_links, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$mean_repeatability, 0.500, na.rm = TRUE),
      quantile(model_summary_by_rep$median_repeatability, 0.500, na.rm = TRUE)
    ),
    model_q975 = c(
      quantile(model_summary_by_rep$mean_interacting_sites, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$median_interacting_sites, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$q95_interacting_sites, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$max_interacting_sites, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$proportion_single_site_links, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$mean_repeatability, 0.975, na.rm = TRUE),
      quantile(model_summary_by_rep$median_repeatability, 0.975, na.rm = TRUE)
    )
  ) %>%
    mutate(
      empirical_above_model_975 = empirical_value > model_q975,
      empirical_below_model_025 = empirical_value < model_q025
    )

  write.csv2(empirical_repeatability,
            file.path(out_dir, paste0(dataset, "_empirical_interaction_repeatability.csv")),
            row.names = FALSE)

  write.csv2(model_repeatability,
            file.path(out_dir, paste0(dataset, "_model_generated_interaction_repeatability.csv")),
            row.names = FALSE)

  write.csv2(summary_by_source,
            file.path(out_dir, paste0(dataset, "_repeatability_summary_by_source.csv")),
            row.names = FALSE)

  write.csv2(model_summary_by_rep,
            file.path(out_dir, paste0(dataset, "_model_summary_by_rep.csv")),
            row.names = FALSE)

  write.csv2(comparison_to_model,
            file.path(out_dir, paste0(dataset, "_empirical_vs_model_repeatability_summary.csv")),
            row.names = FALSE)

  plot_data <- all_repeatability %>%
    mutate(
      source = factor(source, levels = c("model_generated", "empirical")),
      n_interacting_sites_capped = pmin(n_interacting_sites, 20)
    )

  p_hist <- ggplot(plot_data, aes(x = n_interacting_sites_capped, fill = source)) +
    geom_histogram(position = "identity", alpha = 0.45, bins = 20) +
    theme_classic(base_size = 14) +
    xlab("Number of sites where interaction is observed, capped at 20") +
    ylab("Number of pairwise interactions") +
    ggtitle(paste0(dataset, ": spatial repetition of realised interactions"),
            subtitle = "Empirical vs model-generated links")

  ggsave(file.path(out_dir, paste0(dataset, "_hist_interacting_sites_empirical_vs_model.png")),
         p_hist, width = 8, height = 5.5, dpi = 300)

  p_density <- ggplot(plot_data, aes(x = repeatability, color = source)) +
    geom_density(linewidth = 1.2, na.rm = TRUE) +
    theme_classic(base_size = 14) +
    xlab("Interaction repeatability: interacting sites / co-occurring sites") +
    ylab("Density") +
    ggtitle(paste0(dataset, ": interaction repeatability"),
            subtitle = "Empirical vs model-generated links")

  ggsave(file.path(out_dir, paste0(dataset, "_density_repeatability_empirical_vs_model.png")),
         p_density, width = 8, height = 5.5, dpi = 300)

  list(empirical = empirical_repeatability,
       model = model_repeatability,
       comparison = comparison_to_model,
       model_summary_by_rep = model_summary_by_rep)
}

all_outputs <- lapply(all_dataset_names, run_one_repeatability)
names(all_outputs) <- all_dataset_names

combined_comparison <- bind_rows(lapply(all_outputs, `[[`, "comparison"))
combined_model_summary <- bind_rows(lapply(all_outputs, `[[`, "model_summary_by_rep"))

write.csv2(combined_comparison,
          file.path(combined_out, "interaction_repeatability_summary_all_datasets.csv"),
          row.names = FALSE)

combined_comparison$dataset <- factor(combined_comparison$dataset, levels = all_dataset_names)

p_combined <- ggplot(combined_comparison, aes(x = metric, y = empirical_value)) +
  geom_errorbar(aes(ymin = model_q025, ymax = model_q975), width = 0.15) +
  geom_point(aes(y = model_q500), shape = 1, size = 2.2) +
  geom_point(size = 2.2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Metric value") +
  ggtitle("Interaction repeatability metrics across all datasets",
          subtitle = "Error bars = 95% model envelope; open points = model median; filled points = empirical")

ggsave(file.path(combined_out, "combined_interaction_repeatability_metrics.png"),
       p_combined, width = 16, height = 8, dpi = 300)

message("Finished interaction repeatability analysis.")
