
## ------------------------------------------------------------
## Script: All/scripts/05_degree_site_removal_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

parallel_packages <- c("future", "future.apply")
for(pkg in parallel_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 80#100
n_model_reps <- 80#100
selected_removal_levels <- c(0, 0.4, 0.8)

## Parallel settings.
## On Windows, multisession is safer than multicore.
## Leave one core free by default.
n_workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers for model replicates.")

dirs <- make_output_dirs("degree_site_removal")
sep_out <- dirs$separated
combined_out <- dirs$combined

evaluate_network_for_subsets <- function(site_interactions, subset_object, source, model_rep = NA_integer_){
  metric_list <- vector("list", nrow(subset_object$index))
  freq_list <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]

    degrees <- degree_table_from_interactions(site_interactions, sites_keep)

    cur_metrics <- summarise_degrees(degrees) %>%
      mutate(
        source = source,
        model_rep = model_rep,
        subset_id = subset_id,
        removal_fraction = subset_object$index$removal_fraction[i],
        site_rep = subset_object$index$site_rep[i],
        n_sites_kept = subset_object$index$n_sites_kept[i]
      ) %>%
      select(source, model_rep, subset_id, removal_fraction, site_rep, n_sites_kept,
             trophic_level, everything())

    cur_freq <- frequency_degrees(degrees) %>%
      mutate(
        source = source,
        model_rep = model_rep,
        subset_id = subset_id,
        removal_fraction = subset_object$index$removal_fraction[i],
        site_rep = subset_object$index$site_rep[i],
        n_sites_kept = subset_object$index$n_sites_kept[i]
      ) %>%
      select(source, model_rep, subset_id, removal_fraction, site_rep, n_sites_kept,
             trophic_level, degree, n_species)

    metric_list[[i]] <- cur_metrics
    freq_list[[i]] <- cur_freq
  }

  list(metrics = bind_rows(metric_list), frequency = bind_rows(freq_list))
}

make_metric_envelope <- function(model_metrics, empirical_metrics){
  metric_cols <- c(
    "mean_degree", "median_degree", "variance_degree", "maximum_degree",
    "proportion_degree_1", "n_species_degree_gt0", "gini_degree"
  )

  model_long <- model_metrics %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value")

  model_mean_curves <- model_long %>%
    group_by(model_rep, trophic_level, removal_fraction, metric) %>%
    summarise(model_mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  model_envelope <- model_mean_curves %>%
    group_by(trophic_level, removal_fraction, metric) %>%
    summarise(
      model_q025 = quantile(model_mean_value, 0.025, na.rm = TRUE),
      model_q500 = quantile(model_mean_value, 0.500, na.rm = TRUE),
      model_q975 = quantile(model_mean_value, 0.975, na.rm = TRUE),
      model_mean = mean(model_mean_value, na.rm = TRUE),
      .groups = "drop"
    )

  empirical_long <- empirical_metrics %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "empirical_value") %>%
    group_by(trophic_level, removal_fraction, metric) %>%
    summarise(
      empirical_mean = mean(empirical_value, na.rm = TRUE),
      empirical_sd = sd(empirical_value, na.rm = TRUE),
      .groups = "drop"
    )

  empirical_long %>%
    left_join(model_envelope, by = c("trophic_level", "removal_fraction", "metric")) %>%
    mutate(
      empirical_above_model_975 = empirical_mean > model_q975,
      empirical_below_model_025 = empirical_mean < model_q025
    )
}

make_frequency_envelope <- function(model_frequency, empirical_frequency){
  model_summary <- model_frequency %>%
    group_by(model_rep, removal_fraction, trophic_level, degree) %>%
    summarise(n_species = mean(n_species, na.rm = TRUE), .groups = "drop")

  model_envelope <- model_summary %>%
    group_by(removal_fraction, trophic_level, degree) %>%
    summarise(
      model_q025 = quantile(n_species, 0.025, na.rm = TRUE),
      model_q500 = quantile(n_species, 0.500, na.rm = TRUE),
      model_q975 = quantile(n_species, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  empirical_summary <- empirical_frequency %>%
    group_by(removal_fraction, trophic_level, degree) %>%
    summarise(empirical_mean = mean(n_species, na.rm = TRUE), .groups = "drop")

  degrees_all <- sort(unique(c(model_envelope$degree, empirical_summary$degree)))
  if(length(degrees_all) == 0) degrees_all <- 0

  full_grid <- expand.grid(
    removal_fraction = selected_removal_levels,
    trophic_level = c("consumer", "resource"),
    degree = degrees_all,
    stringsAsFactors = FALSE
  )

  full_grid %>%
    left_join(model_envelope, by = c("removal_fraction", "trophic_level", "degree")) %>%
    left_join(empirical_summary, by = c("removal_fraction", "trophic_level", "degree")) %>%
    mutate(across(c(model_q025, model_q500, model_q975, empirical_mean), ~replace_na(.x, 0)))
}

plot_metric_envelope <- function(comparison, dataset, metric_name, y_label, output_file){
  dat <- comparison %>% filter(metric == metric_name)

  p <- ggplot(dat, aes(x = removal_fraction)) +
    geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") +
    geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) +
    geom_line(aes(y = empirical_mean), linewidth = 1.2) +
    geom_point(aes(y = empirical_mean), size = 2.5) +
    facet_wrap(~ trophic_level, scales = "free_y") +
    theme_classic(base_size = 14) +
    xlab("Fraction of sites removed") +
    ylab(y_label) +
    ggtitle(paste0(dataset, ": ", y_label, " under random site removal"),
            subtitle = "Ribbon = 95% model-generated envelope; dashed = model median; solid = empirical")

  ggsave(output_file, p, width = 8.5, height = 5.5, dpi = 300)
}

plot_degree_frequency <- function(freq_comparison, dataset, output_file){
  dat <- freq_comparison %>%
    filter(removal_fraction %in% selected_removal_levels)

  p <- ggplot(dat, aes(x = degree)) +
    geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") +
    geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) +
    geom_line(aes(y = empirical_mean), linewidth = 1) +
    geom_point(aes(y = empirical_mean), size = 2) +
    facet_grid(trophic_level ~ removal_fraction, scales = "free_y") +
    theme_classic(base_size = 13) +
    xlab("Degree") +
    ylab("Number of species") +
    ggtitle(paste0(dataset, ": degree-frequency distributions"),
            subtitle = "Columns = removal fraction; ribbon = model envelope; solid = empirical")

  ggsave(output_file, p, width = 10, height = 6.5, dpi = 300)
}

run_one_degree <- function(dataset){
  message("Running degree site removal: ", dataset)

  p_fixed <- get_p_fixed(dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions

  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)

  empirical_out <- evaluate_network_for_subsets(
    site_interactions = empirical_site_interactions,
    subset_object = subset_object,
    source = "empirical",
    model_rep = NA_integer_
  )

  ## Parallel model-generated replicates.
  ## Each worker simulates one full site-level model network and evaluates it
  ## against the same site-removal subsets.
  model_outputs <- future.apply::future_lapply(
    seq_len(n_model_reps),
    function(m){
      suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(tibble)
      })

      message(dataset, ": model replicate ", m, " / ", n_model_reps)

      sim_interactions <- simulate_model_site_interactions(cooc_triples, p_fixed)

      evaluate_network_for_subsets(
        site_interactions = sim_interactions,
        subset_object = subset_object,
        source = "model_generated",
        model_rep = m
      )
    },
    future.seed = TRUE
  )

  model_metrics <- bind_rows(lapply(model_outputs, `[[`, "metrics"))
  model_frequency <- bind_rows(lapply(model_outputs, `[[`, "frequency"))

  comparison <- make_metric_envelope(model_metrics, empirical_out$metrics)

  freq_comparison <- make_frequency_envelope(
    model_frequency %>% filter(removal_fraction %in% selected_removal_levels),
    empirical_out$frequency %>% filter(removal_fraction %in% selected_removal_levels)
  )

  empirical_metrics_out <- empirical_out$metrics %>% mutate(dataset = dataset)
  model_metrics_out <- model_metrics %>% mutate(dataset = dataset)
  comparison_out <- comparison %>% mutate(dataset = dataset)
  empirical_freq_out <- empirical_out$frequency %>% mutate(dataset = dataset)
  model_freq_out <- model_frequency %>% mutate(dataset = dataset)
  freq_comparison_out <- freq_comparison %>% mutate(dataset = dataset)

  write.csv2(empirical_metrics_out,
            file.path(out_dir, paste0(dataset, "_degree_metrics_empirical.csv")),
            row.names = FALSE)

  write.csv2(model_metrics_out,
            file.path(out_dir, paste0(dataset, "_degree_metrics_model_generated.csv")),
            row.names = FALSE)

  write.csv2(comparison_out,
            file.path(out_dir, paste0(dataset, "_degree_metrics_empirical_vs_model_envelope.csv")),
            row.names = FALSE)

  write.csv2(empirical_freq_out,
            file.path(out_dir, paste0(dataset, "_degree_frequency_empirical.csv")),
            row.names = FALSE)

  write.csv2(model_freq_out,
            file.path(out_dir, paste0(dataset, "_degree_frequency_model_generated.csv")),
            row.names = FALSE)

  write.csv2(freq_comparison_out,
            file.path(out_dir, paste0(dataset, "_degree_frequency_empirical_vs_model_envelope.csv")),
            row.names = FALSE)

  plot_metric_envelope(comparison, dataset, "mean_degree", "Mean degree",
                       file.path(out_dir, paste0(dataset, "_mean_degree_empirical_vs_model_envelope.png")))

  plot_metric_envelope(comparison, dataset, "maximum_degree", "Maximum degree",
                       file.path(out_dir, paste0(dataset, "_maximum_degree_empirical_vs_model_envelope.png")))

  plot_metric_envelope(comparison, dataset, "proportion_degree_1", "Proportion of degree-1 species",
                       file.path(out_dir, paste0(dataset, "_proportion_degree1_empirical_vs_model_envelope.png")))

  plot_metric_envelope(comparison, dataset, "gini_degree", "Degree Gini coefficient",
                       file.path(out_dir, paste0(dataset, "_degree_gini_empirical_vs_model_envelope.png")))

  plot_degree_frequency(freq_comparison, dataset,
                        file.path(out_dir, paste0(dataset, "_degree_frequency_selected_removals_empirical_vs_model_envelope.png")))

  list(
    empirical_metrics = empirical_metrics_out,
    model_metrics = model_metrics_out,
    comparison = comparison_out,
    empirical_frequency = empirical_freq_out,
    model_frequency = model_freq_out,
    frequency_comparison = freq_comparison_out
  )
}

all_outputs <- lapply(all_dataset_names, run_one_degree)
names(all_outputs) <- all_dataset_names

combined_empirical_metrics <- bind_rows(lapply(all_outputs, `[[`, "empirical_metrics"))
combined_model_metrics <- bind_rows(lapply(all_outputs, `[[`, "model_metrics"))
combined_comparison <- bind_rows(lapply(all_outputs, `[[`, "comparison"))
combined_empirical_frequency <- bind_rows(lapply(all_outputs, `[[`, "empirical_frequency"))
combined_model_frequency <- bind_rows(lapply(all_outputs, `[[`, "model_frequency"))
combined_frequency_comparison <- bind_rows(lapply(all_outputs, `[[`, "frequency_comparison"))

write.csv2(combined_empirical_metrics,
          file.path(combined_out, "degree_metrics_empirical_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_model_metrics,
          file.path(combined_out, "degree_metrics_model_generated_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_comparison,
          file.path(combined_out, "degree_metrics_empirical_vs_model_envelope_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_empirical_frequency,
          file.path(combined_out, "degree_frequency_empirical_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_model_frequency,
          file.path(combined_out, "degree_frequency_model_generated_all_datasets.csv"),
          row.names = FALSE)

write.csv2(combined_frequency_comparison,
          file.path(combined_out, "degree_frequency_empirical_vs_model_envelope_all_datasets.csv"),
          row.names = FALSE)

combined_comparison$dataset <- factor(combined_comparison$dataset, levels = all_dataset_names)

plot_combined_metric <- function(metric_name, y_label, file_name){
  dat <- combined_comparison %>% filter(metric == metric_name)

  p <- ggplot(dat, aes(x = removal_fraction)) +
    geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") +
    geom_line(aes(y = model_q500), linetype = 2, linewidth = 0.8) +
    geom_line(aes(y = empirical_mean), linewidth = 0.9) +
    geom_point(aes(y = empirical_mean), size = 1.6) +
    facet_grid(trophic_level ~ dataset, scales = "free_y") +
    theme_classic(base_size = 10) +
    xlab("Fraction of sites removed") +
    ylab(y_label) +
    ggtitle(paste0(y_label, " under site removal across all datasets"),
            subtitle = "Ribbon = 95% model envelope; dashed = model median; solid = empirical")

  ggsave(file.path(combined_out, file_name), p, width = 16, height = 7, dpi = 300)
}

plot_combined_metric("mean_degree", "Mean degree", "combined_degree_mean_all_datasets.png")
plot_combined_metric("maximum_degree", "Maximum degree", "combined_degree_maximum_all_datasets.png")
plot_combined_metric("proportion_degree_1", "Proportion of degree-1 species", "combined_degree1_proportion_all_datasets.png")
plot_combined_metric("gini_degree", "Degree Gini coefficient", "combined_degree_gini_all_datasets.png")

future::plan(future::sequential)
message("Finished degree site-removal analysis.")
