## ------------------------------------------------------------
## Script: All/scripts/12_minimal_heterogeneous_models_all.R
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
n_site_reps <- 100
n_model_reps <- 25
alpha_shrinkage <- 10

n_workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers for heterogeneous model simulations.")

dirs <- make_output_dirs("script12_minimal_heterogeneous_models")
sep_out <- dirs$separated
combined_out <- dirs$combined

model_levels <- c("empirical", "M0_homogeneous", "M1_consumer_specific", "M2_resource_specific",
                  "M3_consumer_resource_combined", "M4_pair_specific_oracle")

clip01 <- function(x) pmin(1, pmax(0, x))

degree_table_from_pairs <- function(pairs){
  if(nrow(pairs) == 0){
    return(data.frame(trophic_level = character(), species = character(), degree = integer()))
  }
  consumer_degrees <- pairs %>%
    group_by(species = consumer) %>%
    summarise(degree = n_distinct(resource), .groups = "drop") %>%
    mutate(trophic_level = "consumer")
  resource_degrees <- pairs %>%
    group_by(species = resource) %>%
    summarise(degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(trophic_level = "resource")
  bind_rows(consumer_degrees, resource_degrees) %>%
    filter(!is.na(species), species != "", degree > 0) %>%
    select(trophic_level, species, degree)
}

summarise_degrees_for_model <- function(pairs){
  deg <- degree_table_from_pairs(pairs)
  out <- summarise_degrees(deg) %>%
    select(trophic_level, mean_degree, maximum_degree, gini_degree, proportion_degree_1)
  out
}

repeatability_metrics_from_site_interactions <- function(site_interactions, cooc_triples_subset){
  if(nrow(site_interactions) == 0){
    return(data.frame(
      mean_repeatability = NA_real_,
      median_repeatability = NA_real_,
      proportion_one_site_realised_links = NA_real_,
      mean_interacting_sites_realised = NA_real_
    ))
  }
  cooc_counts <- cooc_triples_subset %>%
    group_by(consumer, resource) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")
  int_counts <- site_interactions %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(repeatability = n_interacting_sites / n_cooccurring_sites)
  data.frame(
    mean_repeatability = mean(int_counts$repeatability, na.rm = TRUE),
    median_repeatability = median(int_counts$repeatability, na.rm = TRUE),
    proportion_one_site_realised_links = mean(int_counts$n_interacting_sites == 1, na.rm = TRUE),
    mean_interacting_sites_realised = mean(int_counts$n_interacting_sites, na.rm = TRUE)
  )
}

estimate_heterogeneous_probabilities <- function(site_tables, p_global, alpha = 10){
  cooc <- site_tables$cooc_triples
  ints <- site_tables$empirical_site_interactions
  
  cooc_cons <- cooc %>% count(consumer, name = "n_consumer")
  int_cons <- ints %>% count(consumer, name = "x_consumer")
  p_cons <- cooc_cons %>%
    left_join(int_cons, by = "consumer") %>%
    mutate(x_consumer = replace_na(x_consumer, 0),
           p_consumer = (x_consumer + alpha * p_global) / (n_consumer + alpha)) %>%
    select(consumer, p_consumer)
  
  cooc_res <- cooc %>% count(resource, name = "n_resource")
  int_res <- ints %>% count(resource, name = "x_resource")
  p_res <- cooc_res %>%
    left_join(int_res, by = "resource") %>%
    mutate(x_resource = replace_na(x_resource, 0),
           p_resource = (x_resource + alpha * p_global) / (n_resource + alpha)) %>%
    select(resource, p_resource)
  
  cooc_pair <- cooc %>% count(consumer, resource, name = "n_pair")
  int_pair <- ints %>% count(consumer, resource, name = "x_pair")
  p_pair <- cooc_pair %>%
    left_join(int_pair, by = c("consumer", "resource")) %>%
    mutate(x_pair = replace_na(x_pair, 0),
           p_pair = (x_pair + alpha * p_global) / (n_pair + alpha)) %>%
    select(consumer, resource, p_pair)
  
  list(p_cons = p_cons, p_res = p_res, p_pair = p_pair)
}

add_model_probabilities <- function(cooc_subset, probs, p_global){
  cooc_subset %>%
    left_join(probs$p_cons, by = "consumer") %>%
    left_join(probs$p_res, by = "resource") %>%
    left_join(probs$p_pair, by = c("consumer", "resource")) %>%
    mutate(
      p_consumer = replace_na(p_consumer, p_global),
      p_resource = replace_na(p_resource, p_global),
      p_pair = replace_na(p_pair, p_global),
      p_combined = clip01(p_consumer * p_resource / p_global),
      p_global = p_global
    )
}

simulate_model_interactions <- function(cooc_probs, model_type){
  p_vec <- switch(
    model_type,
    M0_homogeneous = cooc_probs$p_global,
    M1_consumer_specific = cooc_probs$p_consumer,
    M2_resource_specific = cooc_probs$p_resource,
    M3_consumer_resource_combined = cooc_probs$p_combined,
    M4_pair_specific_oracle = cooc_probs$p_pair
  )
  cooc_probs %>%
    mutate(interaction = rbinom(n(), size = 1, prob = p_vec)) %>%
    filter(interaction == 1) %>%
    select(site, consumer, resource) %>%
    distinct(site, consumer, resource)
}

summarise_one_network <- function(dataset, model_type, model_rep, subset_row, site_interactions,
                                  cooc_subset, L_emp, L_full_empirical){
  pairs <- site_interactions %>% distinct(consumer, resource)
  L_model <- nrow(pairs)
  rep_met <- repeatability_metrics_from_site_interactions(site_interactions, cooc_subset)
  deg_met <- summarise_degrees_for_model(pairs)
  
  base <- data.frame(
    dataset = dataset,
    model_type = model_type,
    model_rep = model_rep,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_kept = subset_row$n_sites_kept,
    retained_metaweb_links = L_model,
    retained_link_proportion = L_model / L_full_empirical,
    empirical_retained_links = L_emp,
    link_prediction_error = L_emp - L_model,
    absolute_error = abs(L_emp - L_model),
    relative_error = ifelse(L_emp > 0, (L_emp - L_model) / L_emp, NA_real_),
    rep_met
  )
  
  bind_rows(lapply(seq_len(nrow(deg_met)), function(i){
    cbind(base, deg_met[i, , drop = FALSE])
  }))
}

run_one_dataset <- function(dataset){
  message("Running script 12: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_global <- get_p_fixed(dataset)
  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions
  probs <- estimate_heterogeneous_probabilities(site_tables, p_global, alpha_shrinkage)
  
  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  L_full_empirical <- empirical_site_interactions %>% distinct(consumer, resource) %>% nrow()
  
  empirical_rows <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    sites_keep <- subset_object$subsets[[i]]
    subset_row <- subset_object$index[i, ]
    cooc_subset <- cooc_triples %>% filter(site %in% sites_keep)
    emp_int <- empirical_site_interactions %>% filter(site %in% sites_keep)
    L_emp <- emp_int %>% distinct(consumer, resource) %>% nrow()
    summarise_one_network(dataset, "empirical", NA_integer_, subset_row, emp_int, cooc_subset, L_emp, L_full_empirical)
  }))
  
  model_tasks <- expand.grid(
    subset_i = seq_len(nrow(subset_object$index)),
    model_type = c("M0_homogeneous", "M1_consumer_specific", "M2_resource_specific",
                   "M3_consumer_resource_combined", "M4_pair_specific_oracle"),
    model_rep = seq_len(n_model_reps),
    stringsAsFactors = FALSE
  )
  
  model_rows_list <- future.apply::future_lapply(seq_len(nrow(model_tasks)), function(tt){
    suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble)})
    task <- model_tasks[tt, ]
    sites_keep <- subset_object$subsets[[task$subset_i]]
    subset_row <- subset_object$index[task$subset_i, ]
    cooc_subset <- cooc_triples %>% filter(site %in% sites_keep)
    cooc_probs <- add_model_probabilities(cooc_subset, probs, p_global)
    sim_int <- simulate_model_interactions(cooc_probs, task$model_type)
    L_emp <- empirical_site_interactions %>%
      filter(site %in% sites_keep) %>%
      distinct(consumer, resource) %>%
      nrow()
    summarise_one_network(dataset, task$model_type, task$model_rep, subset_row, sim_int, cooc_subset, L_emp, L_full_empirical)
  }, future.seed = TRUE)
  
  model_rows <- bind_rows(model_rows_list)
  out <- bind_rows(empirical_rows, model_rows) %>%
    mutate(
      p_global = p_global,
      alpha_shrinkage = alpha_shrinkage,
      model_label = case_when(
        model_type == "M0_homogeneous" ~ "M0 homogeneous global p",
        model_type == "M1_consumer_specific" ~ "M1 consumer-specific p",
        model_type == "M2_resource_specific" ~ "M2 resource-specific p",
        model_type == "M3_consumer_resource_combined" ~ "M3 consumer-resource combined p",
        model_type == "M4_pair_specific_oracle" ~ "M4 pair-specific/oracle p",
        TRUE ~ "Empirical"
      )
    )
  
  write.csv2(out, file.path(out_dir, paste0(dataset, "_script12_heterogeneous_model_by_subset.csv")), row.names = FALSE)
  out
}

all_by_subset <- bind_rows(lapply(all_dataset_names, run_one_dataset))

write.csv2(all_by_subset,
           file.path(combined_out, "script12_heterogeneous_model_by_subset_combined.csv"),
           row.names = FALSE)

summary_combined <- all_by_subset %>%
  group_by(dataset, model_type, model_label, removal_fraction, trophic_level) %>%
  summarise(
    retained_metaweb_links = mean(retained_metaweb_links, na.rm = TRUE),
    retained_link_proportion = mean(retained_link_proportion, na.rm = TRUE),
    link_prediction_error = mean(link_prediction_error, na.rm = TRUE),
    absolute_error = mean(absolute_error, na.rm = TRUE),
    relative_error = mean(relative_error, na.rm = TRUE),
    mean_repeatability = mean(mean_repeatability, na.rm = TRUE),
    median_repeatability = mean(median_repeatability, na.rm = TRUE),
    proportion_one_site_realised_links = mean(proportion_one_site_realised_links, na.rm = TRUE),
    mean_degree = mean(mean_degree, na.rm = TRUE),
    maximum_degree = mean(maximum_degree, na.rm = TRUE),
    gini_degree = mean(gini_degree, na.rm = TRUE),
    proportion_degree_1 = mean(proportion_degree_1, na.rm = TRUE),
    .groups = "drop"
  )

write.csv2(summary_combined,
           file.path(combined_out, "script12_heterogeneous_model_summary_combined.csv"),
           row.names = FALSE)

errors_combined <- all_by_subset %>%
  filter(model_type != "empirical") %>%
  group_by(dataset, model_type, model_label) %>%
  summarise(
    RMSE = sqrt(mean(link_prediction_error^2, na.rm = TRUE)),
    MAE = mean(absolute_error, na.rm = TRUE),
    mean_relative_error = mean(relative_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(dataset) %>%
  mutate(
    RMSE_M0 = RMSE[model_type == "M0_homogeneous"][1],
    MAE_M0 = MAE[model_type == "M0_homogeneous"][1],
    RMSE_relative_to_M0 = RMSE / RMSE_M0,
    MAE_relative_to_M0 = MAE / MAE_M0
  ) %>%
  ungroup()

write.csv2(errors_combined,
           file.path(combined_out, "script12_heterogeneous_model_errors_combined.csv"),
           row.names = FALSE)

write.csv2(errors_combined,
           file.path(combined_out, "script12_model_improvement_over_M0.csv"),
           row.names = FALSE)

## ------------------------------------------------------------
## Combined plots only
## ------------------------------------------------------------

summary_combined$dataset <- factor(summary_combined$dataset, levels = all_dataset_names)
summary_combined$model_type <- factor(summary_combined$model_type, levels = model_levels)
all_by_subset$dataset <- factor(all_by_subset$dataset, levels = all_dataset_names)
all_by_subset$model_type <- factor(all_by_subset$model_type, levels = model_levels)
errors_combined$dataset <- factor(errors_combined$dataset, levels = all_dataset_names)
errors_combined$model_type <- factor(errors_combined$model_type, levels = model_levels)

retention_plot_data <- summary_combined %>%
  group_by(dataset, model_type, model_label, removal_fraction) %>%
  summarise(retained_link_proportion = mean(retained_link_proportion, na.rm = TRUE), .groups = "drop")

p_ret <- ggplot(retention_plot_data, aes(removal_fraction, retained_link_proportion, linetype = model_type)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") + ylab("Retained-link proportion") +
  ggtitle("Retained-link proportion by model",
          subtitle = "M4 is an upper-bound/oracle model using pair-level empirical information")

ggsave(file.path(combined_out, "script12_retained_link_proportion_by_model.png"), p_ret, width = 16, height = 7, dpi = 300)

error_plot_data <- summary_combined %>%
  filter(model_type != "empirical") %>%
  group_by(dataset, model_type, model_label, removal_fraction) %>%
  summarise(link_prediction_error = mean(link_prediction_error, na.rm = TRUE), .groups = "drop")

p_err <- ggplot(error_plot_data, aes(removal_fraction, link_prediction_error, linetype = model_type)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") + ylab("Observed retained links - predicted retained links") +
  ggtitle("Link prediction error by model")

ggsave(file.path(combined_out, "script12_link_prediction_error_by_model.png"), p_err, width = 16, height = 7, dpi = 300)

p_perf <- ggplot(errors_combined, aes(model_type, RMSE_relative_to_M0)) +
  geom_hline(yintercept = 1, linetype = 2) + geom_point(size = 2) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("RMSE relative to M0") +
  ggtitle("Overall model performance relative to homogeneous M0",
          subtitle = "Values below 1 improve on M0; M4 is oracle/upper-bound")

ggsave(file.path(combined_out, "script12_model_performance_summary.png"), p_perf, width = 16, height = 7, dpi = 300)

repeat_plot_data <- summary_combined %>%
  group_by(dataset, model_type, model_label, removal_fraction) %>%
  summarise(mean_repeatability = mean(mean_repeatability, na.rm = TRUE), .groups = "drop")

p_rep <- ggplot(repeat_plot_data, aes(removal_fraction, mean_repeatability, linetype = model_type)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") + ylab("Mean repeatability") +
  ggtitle("Repeatability reproduced by each model")

ggsave(file.path(combined_out, "script12_repeatability_by_model.png"), p_rep, width = 16, height = 7, dpi = 300)

degree_plot_data <- summary_combined %>%
  filter(model_type %in% model_levels) %>%
  group_by(dataset, model_type, model_label, trophic_level, removal_fraction) %>%
  summarise(gini_degree = mean(gini_degree, na.rm = TRUE),
            maximum_degree = mean(maximum_degree, na.rm = TRUE), .groups = "drop")

p_gini <- ggplot(degree_plot_data, aes(removal_fraction, gini_degree, linetype = model_type)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1) +
  facet_grid(trophic_level ~ dataset, scales = "free_y") +
  theme_classic(base_size = 9) +
  xlab("Fraction of sites removed") + ylab("Degree Gini") +
  ggtitle("Degree Gini by model")

ggsave(file.path(combined_out, "script12_degree_gini_by_model.png"), p_gini, width = 16, height = 7, dpi = 300)

p_max <- ggplot(degree_plot_data, aes(removal_fraction, maximum_degree, linetype = model_type)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1) +
  facet_grid(trophic_level ~ dataset, scales = "free_y") +
  theme_classic(base_size = 9) +
  xlab("Fraction of sites removed") + ylab("Maximum degree") +
  ggtitle("Maximum degree by model")

ggsave(file.path(combined_out, "script12_max_degree_by_model.png"), p_max, width = 16, height = 7, dpi = 300)

future::plan(future::sequential)
message("Finished script 12.")
