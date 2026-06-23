## ------------------------------------------------------------
## Script: All/scripts/12.2_calibrated_heterogeneous_models_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

parallel_packages <- c("future", "future.apply", "parallelly")
for(pkg in parallel_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100
run_simulated_metrics <- TRUE
n_model_reps <- 50
alpha_species <- 5
alpha_pair <- 2
lambda_lower <- 1e-10
lambda_initial_upper <- 1
lambda_max_upper <- 1e8
calibration_tolerance <- 1e-6

dirs <- make_output_dirs("script12.2_calibrated_heterogeneous_models")
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")

model_types <- c(
  "M0_homogeneous",
  "M1_consumer_specific",
  "M2_resource_specific",
  "M3_consumer_resource_combined",
  "M4_pair_specific_oracle"
)

build_compact_tables <- function(site_tables){
  cooc <- site_tables$cooc_triples %>% distinct(site, consumer, resource)
  ints <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource)
  site_levels <- sort(unique(cooc$site))
  consumer_levels <- sort(unique(cooc$consumer))
  resource_levels <- sort(unique(cooc$resource))
  cooc_compact <- cooc %>%
    mutate(site_id = match(site, site_levels),
           consumer_id = match(consumer, consumer_levels),
           resource_id = match(resource, resource_levels),
           pair_key = paste(consumer_id, resource_id, sep = "___"))
  pair_levels <- sort(unique(cooc_compact$pair_key))
  cooc_compact <- cooc_compact %>%
    mutate(pair_id = match(pair_key, pair_levels)) %>%
    select(site_id, consumer_id, resource_id, pair_id, site, consumer, resource)
  ints_compact <- ints %>%
    mutate(site_id = match(site, site_levels),
           consumer_id = match(consumer, consumer_levels),
           resource_id = match(resource, resource_levels),
           pair_key = paste(consumer_id, resource_id, sep = "___"),
           pair_id = match(pair_key, pair_levels)) %>%
    filter(!is.na(site_id), !is.na(consumer_id), !is.na(resource_id), !is.na(pair_id)) %>%
    select(site_id, consumer_id, resource_id, pair_id, site, consumer, resource)
  pair_lookup <- cooc_compact %>% distinct(pair_id, consumer_id, resource_id, consumer, resource)
  list(cooc = cooc_compact, interactions = ints_compact, pair_lookup = pair_lookup,
       site_levels = site_levels, consumer_levels = consumer_levels, resource_levels = resource_levels)
}

normalise_weight <- function(w){
  m <- mean(w, na.rm = TRUE)
  if(is.na(m) || m <= 0) return(rep(1, length(w)))
  w / m
}

add_model_weights <- function(compact, p_global){
  cooc <- compact$cooc
  ints <- compact$interactions
  n_opportunities <- nrow(cooc)
  n_events <- nrow(ints)
  global_site_rate <- ifelse(n_opportunities > 0, n_events / n_opportunities, NA_real_)
  if(is.na(global_site_rate) || global_site_rate <= 0) global_site_rate <- p_global
  consumer_rates <- cooc %>% count(consumer_id, name = "n_consumer") %>%
    left_join(ints %>% count(consumer_id, name = "x_consumer"), by = "consumer_id") %>%
    mutate(x_consumer = replace_na(x_consumer, 0),
           consumer_rate = (x_consumer + alpha_species * global_site_rate) / (n_consumer + alpha_species),
           consumer_weight_raw = consumer_rate / global_site_rate) %>%
    select(consumer_id, consumer_weight_raw)
  resource_rates <- cooc %>% count(resource_id, name = "n_resource") %>%
    left_join(ints %>% count(resource_id, name = "x_resource"), by = "resource_id") %>%
    mutate(x_resource = replace_na(x_resource, 0),
           resource_rate = (x_resource + alpha_species * global_site_rate) / (n_resource + alpha_species),
           resource_weight_raw = resource_rate / global_site_rate) %>%
    select(resource_id, resource_weight_raw)
  pair_rates <- cooc %>% count(pair_id, name = "n_pair") %>%
    left_join(ints %>% count(pair_id, name = "x_pair"), by = "pair_id") %>%
    mutate(x_pair = replace_na(x_pair, 0),
           pair_rate = (x_pair + alpha_pair * global_site_rate) / (n_pair + alpha_pair),
           pair_weight_raw = pair_rate / global_site_rate) %>%
    select(pair_id, pair_weight_raw)
  cooc <- cooc %>%
    left_join(consumer_rates, by = "consumer_id") %>%
    left_join(resource_rates, by = "resource_id") %>%
    left_join(pair_rates, by = "pair_id") %>%
    mutate(weight_M0_homogeneous = 1,
           weight_M1_consumer_specific = consumer_weight_raw,
           weight_M2_resource_specific = resource_weight_raw,
           weight_M3_consumer_resource_combined = consumer_weight_raw * resource_weight_raw,
           weight_M4_pair_specific_oracle = pair_weight_raw)
  for(m in model_types){
    col <- paste0("weight_", m)
    cooc[[col]] <- normalise_weight(cooc[[col]])
    cooc[[col]][!is.finite(cooc[[col]]) | cooc[[col]] < 0] <- 0
  }
  list(cooc = cooc, interactions = ints, pair_lookup = compact$pair_lookup,
       site_levels = compact$site_levels, global_site_rate = global_site_rate)
}

expected_links_from_weights <- function(cooc_subset, weight_col, lambda){
  if(nrow(cooc_subset) == 0 || is.na(lambda)) return(NA_real_)
  pair_sum <- rowsum(cooc_subset[[weight_col]], cooc_subset$pair_id, reorder = FALSE)
  sum(1 - exp(-lambda * as.numeric(pair_sum)), na.rm = TRUE)
}

calibrate_lambda <- function(cooc_full, weight_col, empirical_full_links){
  if(nrow(cooc_full) == 0 || empirical_full_links <= 0){
    return(list(lambda = NA_real_, expected_full_links = NA_real_, calibration_error = NA_real_,
                calibration_failed = TRUE, calibration_warning = "empty cooc or zero empirical links"))
  }
  max_possible_links <- length(unique(cooc_full$pair_id))
  if(empirical_full_links >= max_possible_links){
    exp_hi <- expected_links_from_weights(cooc_full, weight_col, lambda_max_upper)
    return(list(lambda = lambda_max_upper, expected_full_links = exp_hi,
                calibration_error = exp_hi - empirical_full_links, calibration_failed = FALSE,
                calibration_warning = "empirical links at or above maximum possible links; lambda set very high"))
  }
  f <- function(lambda) expected_links_from_weights(cooc_full, weight_col, lambda) - empirical_full_links
  upper <- lambda_initial_upper
  f_upper <- f(upper)
  while(is.finite(f_upper) && f_upper < 0 && upper < lambda_max_upper){
    upper <- upper * 10
    f_upper <- f(upper)
  }
  if(!is.finite(f_upper) || f_upper < 0){
    return(list(lambda = NA_real_, expected_full_links = NA_real_, calibration_error = NA_real_,
                calibration_failed = TRUE, calibration_warning = "failed to bracket calibration root"))
  }
  fit <- tryCatch(uniroot(f, interval = c(lambda_lower, upper), tol = calibration_tolerance), error = function(e) NULL)
  if(is.null(fit)){
    return(list(lambda = NA_real_, expected_full_links = NA_real_, calibration_error = NA_real_,
                calibration_failed = TRUE, calibration_warning = "uniroot failed"))
  }
  lambda <- fit$root
  expected_full <- expected_links_from_weights(cooc_full, weight_col, lambda)
  list(lambda = lambda, expected_full_links = expected_full,
       calibration_error = expected_full - empirical_full_links,
       calibration_failed = abs(expected_full - empirical_full_links) > 1e-4,
       calibration_warning = ifelse(abs(expected_full - empirical_full_links) > 1e-4,
                                    "calibration error above tolerance", ""))
}

weight_diagnostics <- function(cooc, weight_col){
  w <- cooc[[weight_col]]
  data.frame(mean_weight = mean(w, na.rm = TRUE), sd_weight = sd(w, na.rm = TRUE),
             min_weight = min(w, na.rm = TRUE), max_weight = max(w, na.rm = TRUE),
             q95_weight = quantile(w, 0.95, na.rm = TRUE), q99_weight = quantile(w, 0.99, na.rm = TRUE))
}

make_site_subsets_compact <- function(site_ids, removal_levels, n_site_reps){
  site_subsets <- list(); subset_index <- NULL; counter <- 1
  for(removal in removal_levels){
    n_keep <- max(1, round(length(site_ids) * (1 - removal)))
    reps_here <- ifelse(removal == 0, 1, n_site_reps)
    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) site_ids else sample(site_ids, n_keep, replace = FALSE)
      site_subsets[[counter]] <- keep
      subset_index <- rbind(subset_index, data.frame(subset_id = counter, removal_fraction = removal,
                                                     site_rep = r, n_sites_kept = length(keep)))
      counter <- counter + 1
    }
  }
  list(index = subset_index, subsets = site_subsets)
}

empirical_links_for_sites <- function(interactions, sites_keep){
  interactions %>% filter(site_id %in% sites_keep) %>% distinct(pair_id) %>% nrow()
}

expected_subset_row <- function(dataset, model_type, cooc, interactions, sites_keep, subset_row,
                                lambda, empirical_full_links, expected_full_links){
  weight_col <- paste0("weight_", model_type)
  cooc_subset <- cooc[cooc$site_id %in% sites_keep, , drop = FALSE]
  empirical_retained_links <- empirical_links_for_sites(interactions, sites_keep)
  expected_retained_links <- expected_links_from_weights(cooc_subset, weight_col, lambda)
  data.frame(dataset = dataset, model_type = model_type, subset_id = subset_row$subset_id,
             removal_fraction = subset_row$removal_fraction, site_rep = subset_row$site_rep,
             n_sites_kept = subset_row$n_sites_kept,
             empirical_retained_links = empirical_retained_links,
             expected_retained_links = expected_retained_links,
             link_error = empirical_retained_links - expected_retained_links,
             relative_link_error = safe_relative_difference(empirical_retained_links, expected_retained_links),
             empirical_retained_link_proportion = empirical_retained_links / empirical_full_links,
             expected_retained_link_proportion = expected_retained_links / expected_full_links)
}

degree_metrics_from_pairs <- function(pair_ids, pair_lookup){
  if(length(pair_ids) == 0){
    return(data.frame(trophic_side = c("consumer", "resource"), mean_degree = NA_real_,
                      maximum_degree = NA_real_, degree_gini = NA_real_, proportion_degree_1 = NA_real_))
  }
  pairs <- pair_lookup %>% filter(pair_id %in% pair_ids)
  bind_rows(
    pairs %>% count(species_id = consumer_id, name = "degree") %>% mutate(trophic_side = "consumer"),
    pairs %>% count(species_id = resource_id, name = "degree") %>% mutate(trophic_side = "resource")
  ) %>%
    group_by(trophic_side) %>%
    summarise(mean_degree = mean(degree, na.rm = TRUE), maximum_degree = max(degree, na.rm = TRUE),
              degree_gini = gini_coefficient(degree), proportion_degree_1 = mean(degree == 1, na.rm = TRUE),
              .groups = "drop") %>%
    right_join(data.frame(trophic_side = c("consumer", "resource")), by = "trophic_side")
}

repeatability_metrics_from_site_interactions <- function(sim_interactions){
  if(nrow(sim_interactions) == 0) return(data.frame(mean_repeatability = NA_real_, proportion_one_site_links = NA_real_))
  pair_counts <- sim_interactions %>% group_by(pair_id) %>% summarise(n_interacting_sites = n_distinct(site_id), .groups = "drop")
  data.frame(mean_repeatability = mean(pair_counts$n_interacting_sites, na.rm = TRUE),
             proportion_one_site_links = mean(pair_counts$n_interacting_sites == 1, na.rm = TRUE))
}

simulate_subset_metrics <- function(dataset, model_type, cooc, pair_lookup, sites_keep, subset_row, lambda, model_rep){
  weight_col <- paste0("weight_", model_type)
  cur <- cooc[cooc$site_id %in% sites_keep, , drop = FALSE]
  if(nrow(cur) == 0){
    degree_out <- data.frame(trophic_side = c("consumer", "resource"), mean_degree = NA_real_,
                             maximum_degree = NA_real_, degree_gini = NA_real_, proportion_degree_1 = NA_real_)
    repeat_out <- data.frame(mean_repeatability = NA_real_, proportion_one_site_links = NA_real_)
  } else {
    p_event <- 1 - exp(-lambda * cur[[weight_col]])
    event <- rbinom(nrow(cur), size = 1, prob = p_event)
    sim_interactions <- cur[event == 1, c("site_id", "pair_id", "consumer_id", "resource_id"), drop = FALSE]
    degree_out <- degree_metrics_from_pairs(unique(sim_interactions$pair_id), pair_lookup)
    repeat_out <- repeatability_metrics_from_site_interactions(sim_interactions)
  }
  degree_long <- degree_out %>%
    mutate(dataset = dataset, model_type = model_type, subset_id = subset_row$subset_id,
           removal_fraction = subset_row$removal_fraction, site_rep = subset_row$site_rep,
           n_sites_kept = subset_row$n_sites_kept, model_rep = model_rep, metric_set = "degree") %>%
    select(dataset, model_type, subset_id, removal_fraction, site_rep, n_sites_kept,
           model_rep, metric_set, trophic_side, everything())
  repeat_long <- data.frame(dataset = dataset, model_type = model_type, subset_id = subset_row$subset_id,
                            removal_fraction = subset_row$removal_fraction, site_rep = subset_row$site_rep,
                            n_sites_kept = subset_row$n_sites_kept, model_rep = model_rep,
                            metric_set = "repeatability", trophic_side = "metaweb",
                            mean_repeatability = repeat_out$mean_repeatability,
                            proportion_one_site_links = repeat_out$proportion_one_site_links)
  list(degree = degree_long, repeatability = repeat_long)
}

summarise_simulated_metrics <- function(sim_degree, sim_repeat){
  degree_summary <- sim_degree %>%
    group_by(dataset, model_type, removal_fraction, trophic_side) %>%
    summarise(mean_degree_q025 = quantile(mean_degree, 0.025, na.rm = TRUE),
              mean_degree_q500 = quantile(mean_degree, 0.500, na.rm = TRUE),
              mean_degree_q975 = quantile(mean_degree, 0.975, na.rm = TRUE),
              maximum_degree_q025 = quantile(maximum_degree, 0.025, na.rm = TRUE),
              maximum_degree_q500 = quantile(maximum_degree, 0.500, na.rm = TRUE),
              maximum_degree_q975 = quantile(maximum_degree, 0.975, na.rm = TRUE),
              degree_gini_q025 = quantile(degree_gini, 0.025, na.rm = TRUE),
              degree_gini_q500 = quantile(degree_gini, 0.500, na.rm = TRUE),
              degree_gini_q975 = quantile(degree_gini, 0.975, na.rm = TRUE),
              proportion_degree_1_q025 = quantile(proportion_degree_1, 0.025, na.rm = TRUE),
              proportion_degree_1_q500 = quantile(proportion_degree_1, 0.500, na.rm = TRUE),
              proportion_degree_1_q975 = quantile(proportion_degree_1, 0.975, na.rm = TRUE),
              .groups = "drop")
  repeat_summary <- sim_repeat %>%
    group_by(dataset, model_type, removal_fraction) %>%
    summarise(mean_repeatability_q025 = quantile(mean_repeatability, 0.025, na.rm = TRUE),
              mean_repeatability_q500 = quantile(mean_repeatability, 0.500, na.rm = TRUE),
              mean_repeatability_q975 = quantile(mean_repeatability, 0.975, na.rm = TRUE),
              proportion_one_site_links_q025 = quantile(proportion_one_site_links, 0.025, na.rm = TRUE),
              proportion_one_site_links_q500 = quantile(proportion_one_site_links, 0.500, na.rm = TRUE),
              proportion_one_site_links_q975 = quantile(proportion_one_site_links, 0.975, na.rm = TRUE),
              .groups = "drop") %>% mutate(trophic_side = "metaweb")
  list(degree = degree_summary, repeatability = repeat_summary)
}

run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble) })
  message("Running calibrated heterogeneous models: ", dataset)
  p_global <- get_p_fixed(dataset)
  site_tables <- get_dataset_site_tables(dataset)
  compact <- build_compact_tables(site_tables)
  weighted <- add_model_weights(compact, p_global)
  cooc <- weighted$cooc; interactions <- weighted$interactions; pair_lookup <- weighted$pair_lookup
  site_ids <- sort(unique(cooc$site_id))
  empirical_full_links <- length(unique(interactions$pair_id))
  subsets <- make_site_subsets_compact(site_ids, removal_levels, n_site_reps)
  calibration_list <- list(); expected_list <- list(); sim_degree_list <- list(); sim_repeat_list <- list()
  for(model_type in model_types){
    weight_col <- paste0("weight_", model_type)
    calib <- calibrate_lambda(cooc, weight_col, empirical_full_links)
    diag <- weight_diagnostics(cooc, weight_col)
    calibration_list[[model_type]] <- data.frame(dataset = dataset, model_type = model_type,
      lambda = calib$lambda, equivalent_probability_when_w_equals_1 = 1 - exp(-calib$lambda),
      empirical_full_links = empirical_full_links, expected_full_links = calib$expected_full_links,
      calibration_error = calib$calibration_error, calibration_failed = calib$calibration_failed,
      calibration_warning = calib$calibration_warning, global_site_rate = weighted$global_site_rate,
      p_global_input = p_global, alpha_species = alpha_species, alpha_pair = alpha_pair,
      empirical_full_site_events = nrow(interactions), n_full_opportunities = nrow(cooc), diag)
    expected_list[[model_type]] <- bind_rows(lapply(seq_len(nrow(subsets$index)), function(i){
      expected_subset_row(dataset, model_type, cooc, interactions, subsets$subsets[[i]], subsets$index[i, ],
                          calib$lambda, empirical_full_links, calib$expected_full_links)
    }))
    if(run_simulated_metrics && !is.na(calib$lambda)){
      sim_out_model <- list(); counter <- 1
      for(i in seq_len(nrow(subsets$index))){
        for(r in seq_len(n_model_reps)){
          sim_out_model[[counter]] <- simulate_subset_metrics(dataset, model_type, cooc, pair_lookup,
                                                              subsets$subsets[[i]], subsets$index[i, ],
                                                              calib$lambda, r)
          counter <- counter + 1
        }
      }
      sim_degree_list[[model_type]] <- bind_rows(lapply(sim_out_model, `[[`, "degree"))
      sim_repeat_list[[model_type]] <- bind_rows(lapply(sim_out_model, `[[`, "repeatability"))
    }
  }
  calibration <- bind_rows(calibration_list)
  expected <- bind_rows(expected_list)
  expected_summary <- expected %>% group_by(dataset, model_type, removal_fraction) %>%
    summarise(mean_empirical_retained_links = mean(empirical_retained_links, na.rm = TRUE),
              mean_expected_retained_links = mean(expected_retained_links, na.rm = TRUE),
              mean_link_error = mean(link_error, na.rm = TRUE),
              mean_relative_link_error = mean(relative_link_error, na.rm = TRUE),
              RMSE = sqrt(mean(link_error^2, na.rm = TRUE)),
              MAE = mean(abs(link_error), na.rm = TRUE), bias = mean(link_error, na.rm = TRUE),
              mean_empirical_retained_link_proportion = mean(empirical_retained_link_proportion, na.rm = TRUE),
              mean_expected_retained_link_proportion = mean(expected_retained_link_proportion, na.rm = TRUE),
              .groups = "drop")
  performance <- expected %>% group_by(dataset, model_type) %>%
    summarise(RMSE = sqrt(mean(link_error^2, na.rm = TRUE)), MAE = mean(abs(link_error), na.rm = TRUE),
              bias = mean(link_error, na.rm = TRUE),
              correlation_empirical_expected = suppressWarnings(cor(empirical_retained_links, expected_retained_links, use = "complete.obs")),
              R2 = suppressWarnings(summary(lm(empirical_retained_links ~ expected_retained_links))$r.squared),
              .groups = "drop") %>%
    group_by(dataset) %>%
    mutate(RMSE_M0 = RMSE[model_type == "M0_homogeneous"][1], MAE_M0 = MAE[model_type == "M0_homogeneous"][1],
           RMSE_relative_to_M0 = RMSE / RMSE_M0, MAE_relative_to_M0 = MAE / MAE_M0) %>% ungroup()
  sim_degree <- if(length(sim_degree_list) > 0) bind_rows(sim_degree_list) else NULL
  sim_repeat <- if(length(sim_repeat_list) > 0) bind_rows(sim_repeat_list) else NULL
  sim_summary <- NULL
  if(!is.null(sim_degree) && !is.null(sim_repeat)){
    sim_sum <- summarise_simulated_metrics(sim_degree, sim_repeat)
    sim_summary <- bind_rows(sim_sum$degree %>% mutate(metric_set = "degree"),
                             sim_sum$repeatability %>% mutate(metric_set = "repeatability"))
  }
  list(calibration = calibration, expected = expected, expected_summary = expected_summary,
       performance = performance, sim_summary = sim_summary)
}

all_outputs <- future.apply::future_lapply(all_dataset_names, run_one_dataset, future.seed = TRUE)
names(all_outputs) <- all_dataset_names
calibration_all <- bind_rows(lapply(all_outputs, `[[`, "calibration"))
expected_all <- bind_rows(lapply(all_outputs, `[[`, "expected"))
expected_summary_all <- bind_rows(lapply(all_outputs, `[[`, "expected_summary"))
performance_all <- bind_rows(lapply(all_outputs, `[[`, "performance"))
sim_summary_all <- bind_rows(lapply(all_outputs, `[[`, "sim_summary"))

write.csv2(calibration_all, file.path(combined_out, "script12.2_calibration_summary.csv"), row.names = FALSE)
write.csv2(expected_all, file.path(combined_out, "script12.2_expected_links_by_subset.csv"), row.names = FALSE)
write.csv2(expected_summary_all, file.path(combined_out, "script12.2_expected_links_summary.csv"), row.names = FALSE)
write.csv2(performance_all, file.path(combined_out, "script12.2_model_performance_summary.csv"), row.names = FALSE)
if(nrow(sim_summary_all) > 0) write.csv2(sim_summary_all, file.path(combined_out, "script12.2_simulated_network_metrics_summary.csv"), row.names = FALSE)

expected_summary_all$dataset <- factor(expected_summary_all$dataset, levels = all_dataset_names)
performance_all$dataset <- factor(performance_all$dataset, levels = all_dataset_names)
calibration_all$dataset <- factor(calibration_all$dataset, levels = all_dataset_names)

p_perf <- ggplot(performance_all, aes(x = model_type, y = RMSE_relative_to_M0)) +
  geom_hline(yintercept = 1, linetype = 2) + geom_col() +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") + theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("RMSE relative to M0") +
  ggtitle("Calibrated model performance relative to homogeneous M0")
ggsave(file.path(combined_out, "script12.2_model_performance_relative_to_M0.png"), p_perf, width = 14, height = 7, dpi = 300)

p_error <- ggplot(expected_summary_all, aes(x = removal_fraction, y = mean_link_error, linetype = model_type)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") + theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") + ylab("Empirical retained links - expected retained links") +
  ggtitle("Calibrated model link prediction error")
ggsave(file.path(combined_out, "script12.2_link_prediction_error_by_model.png"), p_error, width = 14, height = 7, dpi = 300)

retention_long <- expected_summary_all %>%
  select(dataset, model_type, removal_fraction, mean_empirical_retained_link_proportion, mean_expected_retained_link_proportion) %>%
  pivot_longer(cols = c(mean_empirical_retained_link_proportion, mean_expected_retained_link_proportion), names_to = "curve_type", values_to = "retained_link_proportion")
p_retention <- ggplot(retention_long, aes(x = removal_fraction, y = retained_link_proportion, linetype = model_type, alpha = curve_type)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
  scale_alpha_manual(values = c(mean_empirical_retained_link_proportion = 1, mean_expected_retained_link_proportion = 0.55)) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") + theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") + ylab("Retained-link proportion") +
  ggtitle("Empirical and calibrated expected retained-link proportion", subtitle = "Darker curves = empirical; lighter curves = model expectations")
ggsave(file.path(combined_out, "script12.2_retained_link_proportion_by_model.png"), p_retention, width = 14, height = 7, dpi = 300)

p_obs_exp <- ggplot(expected_summary_all, aes(x = mean_empirical_retained_links, y = mean_expected_retained_links, shape = model_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free") + theme_classic(base_size = 10) +
  xlab("Empirical retained links") + ylab("Expected retained links") + ggtitle("Observed vs expected retained links by calibrated model")
ggsave(file.path(combined_out, "script12.2_observed_vs_expected_links_by_model.png"), p_obs_exp, width = 14, height = 7, dpi = 300)

calib0 <- expected_all %>% filter(removal_fraction == 0) %>% group_by(dataset, model_type) %>%
  summarise(empirical_full_links = mean(empirical_retained_links, na.rm = TRUE),
            expected_full_links = mean(expected_retained_links, na.rm = TRUE),
            calibration_error = mean(link_error, na.rm = TRUE), .groups = "drop")
p_calib <- ggplot(calib0, aes(x = empirical_full_links, y = expected_full_links, shape = model_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(size = 2.5) +
  facet_wrap(~ dataset, ncol = 5, scales = "free") + theme_classic(base_size = 10) +
  xlab("Empirical full-network links") + ylab("Expected full-network links after calibration") + ggtitle("Calibration check at removal 0")
ggsave(file.path(combined_out, "script12.2_calibration_check_removal0.png"), p_calib, width = 14, height = 7, dpi = 300)

if(nrow(sim_summary_all) > 0){
  sim_summary_all$dataset <- factor(sim_summary_all$dataset, levels = all_dataset_names)
  repeat_plot_data <- sim_summary_all %>% filter(metric_set == "repeatability", trophic_side == "metaweb")
  p_rep <- ggplot(repeat_plot_data, aes(x = removal_fraction, y = mean_repeatability_q500, linetype = model_type)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) + facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
    theme_classic(base_size = 10) + xlab("Fraction of sites removed") + ylab("Model-generated mean interacting sites per link") +
    ggtitle("Repeatability reproduced by calibrated models")
  ggsave(file.path(combined_out, "script12.2_repeatability_by_model.png"), p_rep, width = 14, height = 7, dpi = 300)
  degree_plot_data <- sim_summary_all %>% filter(metric_set == "degree", trophic_side %in% c("consumer", "resource"))
  p_gini <- ggplot(degree_plot_data, aes(x = removal_fraction, y = degree_gini_q500, linetype = model_type)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) + facet_grid(trophic_side ~ dataset, scales = "free_y") +
    theme_classic(base_size = 9) + xlab("Fraction of sites removed") + ylab("Median simulated degree Gini") + ggtitle("Degree Gini by calibrated model")
  ggsave(file.path(combined_out, "script12.2_degree_gini_by_model.png"), p_gini, width = 16, height = 7, dpi = 300)
  p_max <- ggplot(degree_plot_data, aes(x = removal_fraction, y = maximum_degree_q500, linetype = model_type)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) + facet_grid(trophic_side ~ dataset, scales = "free_y") +
    theme_classic(base_size = 9) + xlab("Fraction of sites removed") + ylab("Median simulated maximum degree") + ggtitle("Maximum degree by calibrated model")
  ggsave(file.path(combined_out, "script12.2_max_degree_by_model.png"), p_max, width = 16, height = 7, dpi = 300)
  p_deg1 <- ggplot(degree_plot_data, aes(x = removal_fraction, y = proportion_degree_1_q500, linetype = model_type)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.2) + facet_grid(trophic_side ~ dataset, scales = "free_y") +
    theme_classic(base_size = 9) + xlab("Fraction of sites removed") + ylab("Median simulated proportion degree-1 species") + ggtitle("Degree-1 species by calibrated model")
  ggsave(file.path(combined_out, "script12.2_degree1_by_model.png"), p_deg1, width = 16, height = 7, dpi = 300)
}

future::plan(future::sequential)
message("Finished script 12.2 calibrated heterogeneous models.")
