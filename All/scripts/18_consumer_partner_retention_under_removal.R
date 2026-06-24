## ------------------------------------------------------------
## Script: All/scripts/18_consumer_partner_retention_under_removal.R
##
## Purpose:
## Consumer-level diagnostic linking full-network partner excess and
## repeatability excess to consumer partner-retention error under random
## site removal, using the original conditioned Galiana IR model only.
##
## Self-contained; does not alter earlier scripts.
## Uses only the original 10 datasets from helper script.
## Consumers only. No beta-binomial models. No resource-side analysis.
##
## Replications deliberately set to 100 model simulations per dataset.
## Site-removal subsets use the same removal levels and 100 random
## replicates per nonzero removal level used in previous scripts.
## Parallelised across datasets.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("dplyr", "tidyr", "tibble", "ggplot2", "future", "future.apply", "parallelly")
for(pkg in packages_extra){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

## ---------------------------
## Settings
## ---------------------------

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100
n_model_sims <- 100
min_sites_retained <- 2

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script assumes the original 10 Galiana datasets only.")
}

if(exists("make_output_dirs", mode = "function")){
  dirs <- make_output_dirs("script18_consumer_partner_retention_under_removal")
  sep_out <- dirs$separated
  combined_out <- dirs$combined
} else {
  out_root <- file.path("All", "outputs", "script18_consumer_partner_retention_under_removal")
  sep_out <- file.path(out_root, "separated")
  combined_out <- file.path(out_root, "combined")
  dir.create(sep_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(combined_out, recursive = TRUE, showWarnings = FALSE)
}

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")
message("Using ", n_model_sims, " homogeneous-model simulations per dataset.")

## ---------------------------
## Helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

consumer_conditioned_link_probability <- function(n_ai, N_alpha, p){
  denom <- 1 - (1 - p)^N_alpha
  ifelse(is.finite(denom) & denom > 0,
         (1 - (1 - p)^n_ai) / denom,
         NA_real_)
}

expected_total_links_conditioned <- function(pair_table, p){
  sum(consumer_conditioned_link_probability(
    n_ai = pair_table$n_cooccurring_sites,
    N_alpha = pair_table$N_alpha,
    p = p
  ), na.rm = TRUE)
}

fit_conditioned_p <- function(pair_table, observed_links){
  if(observed_links <= 0) return(0)
  if(observed_links >= nrow(pair_table)) return(1)

  f <- function(p){
    expected_total_links_conditioned(pair_table, p) - observed_links
  }

  f_low <- f(1e-12)
  f_high <- f(1 - 1e-12)

  if(!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
  if(abs(f_low) < 1e-8) return(1e-12)
  if(abs(f_high) < 1e-8) return(1 - 1e-12)
  if(f_low * f_high > 0){
    warning("Could not bracket conditioned p. Returning NA.")
    return(NA_real_)
  }

  tryCatch(
    uniroot(f, interval = c(1e-12, 1 - 1e-12), tol = 1e-12)$root,
    error = function(e) NA_real_
  )
}

expected_repeatability_realised <- function(n, p){
  denom <- 1 - (1 - p)^n
  ifelse(is.finite(denom) & denom > 0, p / denom, NA_real_)
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

mc_p_two_sided <- function(empirical, simulated){
  simulated <- simulated[is.finite(simulated)]
  if(!is.finite(empirical) || length(simulated) == 0) return(NA_real_)
  centre <- mean(simulated)
  (sum(abs(simulated - centre) >= abs(empirical - centre)) + 1) / (length(simulated) + 1)
}

make_site_subsets_local <- function(sites){
  sites <- sort(unique(as.character(sites)))
  n_sites <- length(sites)
  subset_list <- list()
  subset_index <- list()
  counter <- 1

  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_sites * (1 - removal)))
    n_keep <- min(n_keep, n_sites)
    reps_here <- ifelse(removal == 0, 1, n_site_reps)

    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) sites else sample(sites, n_keep, replace = FALSE)
      subset_list[[counter]] <- sort(keep)
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        actual_removal_fraction = 1 - length(keep) / n_sites,
        site_rep = r,
        n_sites_kept = length(keep)
      )
      counter <- counter + 1
    }
  }

  list(subsets = subset_list, index = dplyr::bind_rows(subset_index))
}

make_tables <- function(dataset){
  site_tables <- get_dataset_site_tables(dataset)

  cooc_cells <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource),
      cell_id = paste(consumer, resource, site, sep = "___")
    )

  empirical_cells <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource),
      cell_id = paste(consumer, resource, site, sep = "___")
    )

  cooc_counts <- cooc_cells %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

  int_counts <- empirical_cells %>%
    group_by(pair_id) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop")

  pair_table <- cooc_counts %>%
    left_join(int_counts, by = "pair_id") %>%
    mutate(
      dataset = dataset,
      n_interacting_sites = tidyr::replace_na(n_interacting_sites, 0L),
      is_realised_link = n_interacting_sites > 0
    ) %>%
    group_by(consumer) %>%
    mutate(N_alpha = sum(n_cooccurring_sites)) %>%
    ungroup() %>%
    select(dataset, consumer, resource, pair_id, n_cooccurring_sites,
           n_interacting_sites, is_realised_link, N_alpha)

  list(pair_table = pair_table, cooc_cells = cooc_cells, empirical_cells = empirical_cells)
}

consumer_full_metrics <- function(pair_table, p){
  pair_aug <- pair_table %>%
    mutate(
      expected_link_probability = consumer_conditioned_link_probability(
        n_ai = n_cooccurring_sites,
        N_alpha = N_alpha,
        p = p
      ),
      observed_repeatability = ifelse(n_interacting_sites > 0,
                                      n_interacting_sites / n_cooccurring_sites,
                                      NA_real_),
      expected_repeatability = expected_repeatability_realised(n_cooccurring_sites, p),
      repeatability_excess = observed_repeatability - expected_repeatability
    )

  pair_aug %>%
    group_by(dataset, consumer) %>%
    summarise(
      n_cooccurring_resources = n(),
      full_empirical_partner_number = sum(n_interacting_sites > 0),
      expected_full_partner_number = sum(expected_link_probability, na.rm = TRUE),
      `Realised partners above model expectation` =
        full_empirical_partner_number - expected_full_partner_number,
      n_realised_links_n_ge_2 = sum(n_interacting_sites > 0 & n_cooccurring_sites >= 2),
      mean_observed_repeatability = mean(observed_repeatability[n_interacting_sites > 0 & n_cooccurring_sites >= 2], na.rm = TRUE),
      mean_expected_repeatability = mean(expected_repeatability[n_interacting_sites > 0 & n_cooccurring_sites >= 2], na.rm = TRUE),
      `Extra repeatability of realised links` = mean(repeatability_excess[n_interacting_sites > 0 & n_cooccurring_sites >= 2], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean_observed_repeatability = ifelse(is.nan(mean_observed_repeatability), NA_real_, mean_observed_repeatability),
      mean_expected_repeatability = ifelse(is.nan(mean_expected_repeatability), NA_real_, mean_expected_repeatability),
      `Extra repeatability of realised links` = ifelse(is.nan(`Extra repeatability of realised links`), NA_real_, `Extra repeatability of realised links`)
    )
}

simulate_full_network <- function(dataset, pair_table, cooc_cells, p){
  sim_cells <- cooc_cells %>%
    mutate(interaction = FALSE)

  consumers <- sort(unique(cooc_cells$consumer))

  for(cons in consumers){
    idx <- which(sim_cells$consumer == cons)
    if(length(idx) == 0) next

    draw <- rbinom(length(idx), size = 1, prob = p) == 1
    attempts <- 1
    while(!any(draw)){
      draw <- rbinom(length(idx), size = 1, prob = p) == 1
      attempts <- attempts + 1
      if(attempts > 100000){
        stop("Conditioned simulation failed for consumer ", cons,
             ". p may be too small or available cells may be invalid.")
      }
    }
    sim_cells$interaction[idx] <- draw
  }

  sim_int_cells <- sim_cells %>%
    filter(interaction) %>%
    select(site, consumer, resource, pair_id, cell_id)

  sim_counts <- sim_int_cells %>%
    group_by(pair_id) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop")

  sim_pair_table <- pair_table %>%
    select(dataset, consumer, resource, pair_id, n_cooccurring_sites, N_alpha) %>%
    left_join(sim_counts, by = "pair_id") %>%
    mutate(
      n_interacting_sites = tidyr::replace_na(n_interacting_sites, 0L),
      is_realised_link = n_interacting_sites > 0
    )

  list(sim_pair_table = sim_pair_table, sim_int_cells = sim_int_cells)
}

retention_mean_by_consumer <- function(full_links, interaction_cells, subset_object){
  full_partner_counts <- full_links %>%
    filter(n_interacting_sites > 0) %>%
    group_by(consumer) %>%
    summarise(full_partners = n_distinct(resource), .groups = "drop")

  consumers <- full_partner_counts$consumer

  rows <- lapply(seq_len(nrow(subset_object$index)), function(i){
    subset_row <- subset_object$index[i, ]
    if(subset_row$removal_fraction == 0) return(NULL)

    sites_keep <- subset_object$subsets[[i]]

    retained <- interaction_cells %>%
      filter(site %in% sites_keep) %>%
      semi_join(full_links %>% filter(n_interacting_sites > 0) %>% select(consumer, resource),
                by = c("consumer", "resource")) %>%
      group_by(consumer) %>%
      summarise(retained_partners = n_distinct(resource), .groups = "drop")

    full_partner_counts %>%
      left_join(retained, by = "consumer") %>%
      mutate(
        retained_partners = tidyr::replace_na(retained_partners, 0L),
        retention = retained_partners / full_partners,
        removal_fraction = subset_row$removal_fraction,
        subset_id = subset_row$subset_id
      ) %>%
      select(consumer, removal_fraction, subset_id, retention)
  })

  dplyr::bind_rows(rows) %>%
    group_by(consumer) %>%
    summarise(mean_retention = mean(retention, na.rm = TRUE), .groups = "drop")
}

simulation_correlations <- function(sim_id, sim_pair_table, sim_int_cells,
                                    sim_retention_mean, loo_prediction,
                                    p){
  metrics <- consumer_full_metrics(sim_pair_table, p) %>%
    left_join(sim_retention_mean, by = "consumer") %>%
    left_join(loo_prediction, by = "consumer") %>%
    mutate(`Mean error in retained partner breadth` = mean_retention - loo_mean_prediction)

  cor_A <- safe_spearman(metrics$`Realised partners above model expectation`,
                         metrics$`Mean error in retained partner breadth`)

  metrics_B <- metrics %>% filter(n_realised_links_n_ge_2 > 0)
  cor_B <- safe_spearman(metrics_B$`Extra repeatability of realised links`,
                         metrics_B$`Mean error in retained partner breadth`)

  tibble::tibble(
    replicate = sim_id,
    `Partner breadth` = cor_A,
    `Link repeatability` = cor_B,
    n_A = nrow(metrics),
    n_B = nrow(metrics_B)
  ) %>%
    tidyr::pivot_longer(
      cols = c(`Partner breadth`, `Link repeatability`),
      names_to = "test_type",
      values_to = "simulated_spearman_correlation"
    ) %>%
    mutate(
      n_included_consumers = ifelse(test_type == "Partner breadth", n_A, n_B)
    ) %>%
    select(replicate, test_type, simulated_spearman_correlation, n_included_consumers)
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
  })

  message("Running script 18 consumer retention diagnostic: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  built <- make_tables(dataset)
  pair_table <- built$pair_table
  cooc_cells <- built$cooc_cells
  empirical_cells <- built$empirical_cells

  observed_full_links <- sum(pair_table$n_interacting_sites > 0)
  p_fit <- fit_conditioned_p(pair_table, observed_full_links)
  expected_full_links <- expected_total_links_conditioned(pair_table, p_fit)

  all_sites <- sort(unique(cooc_cells$site))
  subset_object <- make_site_subsets_local(all_sites)

  empirical_metrics <- consumer_full_metrics(pair_table, p_fit)
  empirical_retention <- retention_mean_by_consumer(pair_table, empirical_cells, subset_object) %>%
    rename(mean_empirical_partner_retention = mean_retention)

  ## Simulate full networks and immediately reduce each to consumer-level retention.
  sim_outputs <- lapply(seq_len(n_model_sims), function(sim_id){
    sim <- simulate_full_network(dataset, pair_table, cooc_cells, p_fit)
    sim_retention <- retention_mean_by_consumer(sim$sim_pair_table, sim$sim_int_cells, subset_object) %>%
      rename(mean_retention = mean_retention)

    list(
      sim_id = sim_id,
      sim_pair_table = sim$sim_pair_table,
      sim_int_cells = sim$sim_int_cells,
      sim_retention = sim_retention
    )
  })

  sim_retention_all <- bind_rows(lapply(sim_outputs, function(x){
    x$sim_retention %>% mutate(replicate = x$sim_id)
  }))

  model_prediction <- sim_retention_all %>%
    group_by(consumer) %>%
    summarise(mean_model_predicted_partner_retention = mean(mean_retention, na.rm = TRUE),
              .groups = "drop")

  species_table <- empirical_metrics %>%
    left_join(empirical_retention, by = "consumer") %>%
    left_join(model_prediction, by = "consumer") %>%
    mutate(
      `Mean error in retained partner breadth` =
        mean_empirical_partner_retention - mean_model_predicted_partner_retention
    ) %>%
    transmute(
      dataset,
      consumer,
      full_empirical_partner_number,
      expected_full_partner_number,
      `Realised partners above model expectation`,
      n_realised_links_n_ge_2,
      `Extra repeatability of realised links`,
      mean_empirical_partner_retention,
      mean_model_predicted_partner_retention,
      `Mean error in retained partner breadth`
    )

  ## Empirical correlations.
  empirical_cor_A <- safe_spearman(species_table$`Realised partners above model expectation`,
                                   species_table$`Mean error in retained partner breadth`)
  empirical_B_table <- species_table %>% filter(n_realised_links_n_ge_2 > 0)
  empirical_cor_B <- safe_spearman(empirical_B_table$`Extra repeatability of realised links`,
                                   empirical_B_table$`Mean error in retained partner breadth`)

  ## Leave-one-out prediction means for each simulation and consumer.
  sim_sums <- sim_retention_all %>%
    group_by(consumer) %>%
    summarise(total_retention_sum = sum(mean_retention, na.rm = TRUE),
              n_sims_available = sum(is.finite(mean_retention)),
              .groups = "drop")

  null_correlations <- bind_rows(lapply(sim_outputs, function(x){
    loo_prediction <- x$sim_retention %>%
      left_join(sim_sums, by = "consumer") %>%
      mutate(
        loo_mean_prediction = ifelse(n_sims_available > 1,
                                     (total_retention_sum - mean_retention) / (n_sims_available - 1),
                                     NA_real_)
      ) %>%
      select(consumer, loo_mean_prediction)

    simulation_correlations(
      sim_id = x$sim_id,
      sim_pair_table = x$sim_pair_table,
      sim_int_cells = x$sim_int_cells,
      sim_retention_mean = x$sim_retention,
      loo_prediction = loo_prediction,
      p = p_fit
    )
  })) %>%
    mutate(dataset = dataset) %>%
    select(dataset, replicate, test_type, simulated_spearman_correlation, n_included_consumers)

  null_summary <- null_correlations %>%
    group_by(dataset, test_type) %>%
    summarise(
      simulation_mean = mean(simulated_spearman_correlation, na.rm = TRUE),
      simulation_q025 = quantile(simulated_spearman_correlation, 0.025, na.rm = TRUE),
      simulation_q975 = quantile(simulated_spearman_correlation, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      empirical_spearman_correlation = ifelse(test_type == "Partner breadth", empirical_cor_A, empirical_cor_B),
      monte_carlo_p_value = mapply(
        function(emp, tt){
          sims <- null_correlations$simulated_spearman_correlation[null_correlations$test_type == tt]
          mc_p_two_sided(emp, sims)
        },
        empirical_spearman_correlation,
        test_type
      ),
      n_included_consumers = ifelse(test_type == "Partner breadth", nrow(species_table), nrow(empirical_B_table))
    )

  wide_null <- null_summary %>%
    mutate(prefix = ifelse(test_type == "Partner breadth", "test_A_partner_breadth", "test_B_link_repeatability")) %>%
    select(prefix, empirical_spearman_correlation, simulation_mean, simulation_q025,
           simulation_q975, monte_carlo_p_value, n_included_consumers) %>%
    tidyr::pivot_wider(
      names_from = prefix,
      values_from = c(empirical_spearman_correlation, simulation_mean, simulation_q025,
                      simulation_q975, monte_carlo_p_value, n_included_consumers),
      names_glue = "{prefix}_{.value}"
    )

  summary_row <- tibble::tibble(
    dataset = dataset,
    fitted_conditioned_p = p_fit,
    observed_full_network_link_number = observed_full_links,
    expected_full_network_link_number = expected_full_links,
    mean_consumer_retention_error = mean(species_table$`Mean error in retained partner breadth`, na.rm = TRUE),
    median_consumer_retention_error = median(species_table$`Mean error in retained partner breadth`, na.rm = TRUE)
  ) %>%
    bind_cols(wide_null)

  ## Validation checks.
  cooc_cell_ids <- cooc_cells$cell_id
  sim_validation <- bind_rows(lapply(sim_outputs, function(x){
    full_partner_counts <- x$sim_pair_table %>%
      group_by(consumer) %>%
      summarise(full_partners = sum(n_interacting_sites > 0),
                total_interactions = sum(n_interacting_sites),
                .groups = "drop")

    retained_zero_possible <- retention_mean_by_consumer(x$sim_pair_table, x$sim_int_cells, subset_object) %>%
      summarise(any_zero_mean = any(mean_retention == 0, na.rm = TRUE)) %>%
      pull(any_zero_mean)

    tibble::tibble(
      replicate = x$sim_id,
      all_consumers_have_interaction_full_network = all(full_partner_counts$total_interactions > 0),
      all_consumers_have_partner_full_network = all(full_partner_counts$full_partners > 0),
      any_consumer_zero_mean_retention_after_removal = retained_zero_possible,
      interaction_outside_cooccurrence_detected = any(!x$sim_int_cells$cell_id %in% cooc_cell_ids),
      duplicate_site_cell_detected = any(duplicated(x$sim_int_cells$cell_id))
    )
  }))

  validation_row <- tibble::tibble(
    dataset = dataset,
    link_number_calibration_ok = isTRUE(abs(observed_full_links - expected_full_links) < 1e-6),
    max_abs_link_calibration_error = abs(observed_full_links - expected_full_links),
    every_simulated_full_network_consumer_has_interaction = all(sim_validation$all_consumers_have_interaction_full_network),
    every_simulated_full_network_consumer_has_partner = all(sim_validation$all_consumers_have_partner_full_network),
    at_least_one_simulated_consumer_can_have_zero_retention_after_removal = any(sim_validation$any_consumer_zero_mean_retention_after_removal),
    interaction_outside_cooccurrence_detected = any(sim_validation$interaction_outside_cooccurrence_detected),
    duplicate_site_cell_detected = any(sim_validation$duplicate_site_cell_detected),
    p_refit_after_removal = FALSE,
    same_removal_subsets_used_for_empirical_and_simulated = TRUE,
    n_model_sims = n_model_sims,
    n_site_removal_subsets = nrow(subset_object$index),
    nonzero_removal_levels = paste(removal_levels[removal_levels > 0], collapse = ",")
  )

  if(!validation_row$link_number_calibration_ok ||
     !validation_row$every_simulated_full_network_consumer_has_interaction ||
     validation_row$interaction_outside_cooccurrence_detected ||
     validation_row$duplicate_site_cell_detected ||
     validation_row$p_refit_after_removal ||
     !validation_row$same_removal_subsets_used_for_empirical_and_simulated){
    print(validation_row)
    stop("Validation failed for dataset ", dataset)
  }

  write.csv2(species_table, file.path(out_dir, paste0(dataset, "_18_consumer_partner_retention_species.csv")), row.names = FALSE)
  write.csv2(null_correlations, file.path(out_dir, paste0(dataset, "_18_consumer_partner_retention_null_correlations.csv")), row.names = FALSE)
  write.csv2(validation_row, file.path(out_dir, paste0(dataset, "_18_consumer_partner_retention_validation.csv")), row.names = FALSE)

  list(
    species_table = species_table,
    summary_row = summary_row,
    null_correlations = null_correlations,
    validation_row = validation_row
  )
}

## ---------------------------
## Run across datasets
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)

names(all_outputs) <- all_dataset_names

species_all <- bind_rows(lapply(all_outputs, `[[`, "species_table"))
summary_all <- bind_rows(lapply(all_outputs, `[[`, "summary_row"))
null_correlations_all <- bind_rows(lapply(all_outputs, `[[`, "null_correlations"))
validation_all <- bind_rows(lapply(all_outputs, `[[`, "validation_row"))

## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(
  species_all,
  file.path(combined_out, "18_consumer_partner_retention_species.csv"),
  row.names = FALSE
)

write.csv2(
  summary_all,
  file.path(combined_out, "18_consumer_partner_retention_summary.csv"),
  row.names = FALSE
)

write.csv2(
  null_correlations_all,
  file.path(combined_out, "18_consumer_partner_retention_null_correlations.csv"),
  row.names = FALSE
)

write.csv2(
  validation_all,
  file.path(combined_out, "18_consumer_partner_retention_validation_checks.csv"),
  row.names = FALSE
)

message("Validation summary:")
print(validation_all)

message("Dataset summary:")
print(summary_all)

## ---------------------------
## Plots
## ---------------------------

species_all <- species_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

point_colour <- "grey20"
trend_colour <- "#0072B2"
emp_colour <- "#D55E00"
interval_colour <- "#009E73"

## Plot 1: partner excess versus retention error
plot_A_data <- species_all %>%
  mutate(
    `Consumer points` = `Mean error in retained partner breadth`,
    `Visual trend` = `Mean error in retained partner breadth`
  )

p1 <- ggplot(plot_A_data,
             aes(x = `Realised partners above model expectation`,
                 y = `Mean error in retained partner breadth`)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_point(aes(colour = "Consumers"), alpha = 0.75, size = 1.8) +
  geom_smooth(aes(colour = "Visual trend"), method = "lm", se = TRUE,
              linewidth = 0.9, alpha = 0.18) +
  scale_colour_manual(values = c("Consumers" = point_colour,
                                 "Visual trend" = trend_colour)) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 10) +
  xlab("Realised partners above model expectation") +
  ylab("Mean error in retained partner breadth") +
  labs(colour = "")

ggsave(
  file.path(combined_out, "18_consumer_partner_excess_retention_error.png"),
  p1,
  width = 14,
  height = 7,
  dpi = 300
)

## Plot 2: repeatability excess versus retention error
plot_B_data <- species_all %>%
  filter(n_realised_links_n_ge_2 > 0)

p2 <- ggplot(plot_B_data,
             aes(x = `Extra repeatability of realised links`,
                 y = `Mean error in retained partner breadth`)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_point(aes(colour = "Consumers"), alpha = 0.75, size = 1.8) +
  geom_smooth(aes(colour = "Visual trend"), method = "lm", se = TRUE,
              linewidth = 0.9, alpha = 0.18) +
  scale_colour_manual(values = c("Consumers" = point_colour,
                                 "Visual trend" = trend_colour)) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 10) +
  xlab("Extra repeatability of realised links") +
  ylab("Mean error in retained partner breadth") +
  labs(colour = "")

ggsave(
  file.path(combined_out, "18_consumer_repeatability_retention_error.png"),
  p2,
  width = 14,
  height = 7,
  dpi = 300
)

## Plot 3: null correlation summary
null_summary_long <- null_correlations_all %>%
  group_by(dataset, test_type) %>%
  summarise(
    simulation_q025 = quantile(simulated_spearman_correlation, 0.025, na.rm = TRUE),
    simulation_q975 = quantile(simulated_spearman_correlation, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

empirical_long <- summary_all %>%
  transmute(
    dataset,
    `Partner breadth` = test_A_partner_breadth_empirical_spearman_correlation,
    `Link repeatability` = test_B_link_repeatability_empirical_spearman_correlation
  ) %>%
  pivot_longer(cols = c(`Partner breadth`, `Link repeatability`),
               names_to = "test_type",
               values_to = "empirical_spearman_correlation")

summary_plot_data <- null_summary_long %>%
  left_join(empirical_long, by = c("dataset", "test_type")) %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p3 <- ggplot(summary_plot_data,
             aes(x = dataset, colour = test_type)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_linerange(aes(ymin = simulation_q025,
                     ymax = simulation_q975),
                 position = position_dodge(width = 0.45), linewidth = 1.2) +
  geom_point(aes(y = empirical_spearman_correlation),
             position = position_dodge(width = 0.45), size = 2.4) +
  scale_colour_manual(values = c("Partner breadth" = interval_colour,
                                 "Link repeatability" = emp_colour)) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Spearman correlation") +
  labs(colour = "")

ggsave(
  file.path(combined_out, "18_consumer_retention_error_null_test.png"),
  p3,
  width = 12,
  height = 6,
  dpi = 300
)

future::plan(future::sequential)

message("Finished script 18 consumer partner retention under removal.")
