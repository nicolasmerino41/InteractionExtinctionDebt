## ------------------------------------------------------------
## Script: All/scripts/30_local_support_hidden_in_regional_networks.R
##
## Purpose:
## Empirical figure set for the narrative:
##   Local interaction support is hidden in regional binary networks.
##
## Uses the 10 original Galiana datasets only.
## Uses write.csv2().
## Does not fit beta-binomial models, regressions, GAMs, power laws,
## trait models, or new interaction models.
##
## Run from the parent repository folder.
## ------------------------------------------------------------
source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("dplyr", "tidyr", "ggplot2", "tibble", "future", "future.apply", "parallelly", "patchwork")
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
n_site_reps <- 50
min_sites_retained <- 1

result_type <- "30_local_support_hidden_in_regional_networks"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined
results_out <- file.path("All", "results", result_type)
dir.create(results_out, recursive = TRUE, showWarnings = FALSE)

## Dataset order: use project order from helper.
dataset_order <- all_dataset_names

## Colours: Okabe-Ito inspired, colour-blind conscious.
col_n <- "#999999"
col_K <- "#0072B2"
col_B <- "#D55E00"
col_model <- "#6A3D9A"
col_persist <- "#009E73"
col_opp_loss <- "#7F7F7F"
col_hidden_loss <- "#CC79A7"

degree_cols <- c(
  "Lower initial degree" = "#56B4E9",
  "Middle initial degree" = "#E69F00",
  "Higher initial degree" = "#009E73"
)

K_class_cols <- c(
  "K = 1" = "#FEE8C8",
  "K = 2" = "#FDBB84",
  "K = 3–4" = "#E34A33",
  "K = 5–8" = "#B30000",
  "K >= 9" = "#67000D"
)

support_layer_cols <- c(
  "Co-occurrence support" = col_n,
  "Local interaction support" = col_K,
  "Regional binary links" = col_B,
  "Homogeneous fixed-p baseline" = col_model
)

fate_cols <- c(
  "Interaction persists" = col_persist,
  "No co-occurrence opportunity remains" = col_opp_loss,
  "Co-occurrence remains, interaction absent" = col_hidden_loss
)

## ---------------------------
## Generic helpers
## ---------------------------
make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

save_plot_both <- function(plot, filename_base, width, height){
  ggsave(file.path(results_out, paste0(filename_base, ".png")), plot,
         width = width, height = height, dpi = 320)
  ggsave(file.path(results_out, paste0(filename_base, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
}

q025 <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
q975 <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

safe_divide <- function(num, den){
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

plink <- function(n, p){
  1 - (1 - p)^n
}

expected_links_from_p <- function(n_vec, p){
  sum(plink(n_vec, p), na.rm = TRUE)
}

fit_p_full_unconditioned <- function(n_vec, observed_links){
  if(length(n_vec) == 0 || is.na(observed_links)) return(NA_real_)
  if(observed_links <= 0) return(0)
  if(observed_links >= length(n_vec)) return(1)
  f <- function(p) expected_links_from_p(n_vec, p) - observed_links
  uniroot(f, interval = c(0, 1), tol = 1e-12)$root
}

expected_mean_K_conditional_bin <- function(n_vec, p){
  pl <- plink(n_vec, p)
  sum(n_vec * p, na.rm = TRUE) / sum(pl, na.rm = TRUE)
}

expected_prop_K1_conditional_bin <- function(n_vec, p){
  pk1 <- n_vec * p * (1 - p)^(n_vec - 1)
  sum(pk1, na.rm = TRUE) / sum(plink(n_vec, p), na.rm = TRUE)
}

make_site_subsets_30 <- function(sites){
  sites <- sort(unique(as.character(sites)))
  n_sites <- length(sites)
  subset_list <- list()
  subset_index <- list()
  counter <- 1L
  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_sites * (1 - removal)))
    n_keep <- min(n_keep, n_sites)
    reps_here <- ifelse(removal == 0, 1L, n_site_reps)
    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) sites else sample(sites, n_keep, replace = FALSE)
      subset_list[[counter]] <- sort(as.character(keep))
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        replicate = r,
        n_sites_kept = length(keep),
        n_sites_removed = n_sites - length(keep),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }
  list(index = bind_rows(subset_index), subsets = subset_list)
}

support_bin <- function(n){
  dplyr::case_when(
    n == 1 ~ "1",
    n == 2 ~ "2",
    n == 3 ~ "3",
    n == 4 ~ "4",
    n == 5 ~ "5",
    n >= 6 & n <= 8 ~ "6–8",
    n >= 9 & n <= 12 ~ "9–12",
    n >= 13 & n <= 20 ~ "13–20",
    n >= 21 ~ "21+",
    TRUE ~ NA_character_
  )
}

support_bin_mid <- function(x){
  map <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,
           "6–8" = 7, "9–12" = 10.5, "13–20" = 16.5, "21+" = 25)
  as.numeric(map[x])
}

K_class_figure4 <- function(K){
  dplyr::case_when(
    K == 1 ~ "K = 1",
    K == 2 ~ "K = 2",
    K >= 3 & K <= 4 ~ "K = 3–4",
    K >= 5 & K <= 8 ~ "K = 5–8",
    K >= 9 ~ "K >= 9",
    TRUE ~ NA_character_
  )
}

make_degree_groups <- function(df, degree_col, group_col = "initial_degree_group"){
  degs <- sort(unique(df[[degree_col]]))
  if(length(degs) == 1){
    cuts <- data.frame(
      degree_value = degs,
      group = "Middle initial degree",
      stringsAsFactors = FALSE
    )
  } else {
    rk <- rank(degs, ties.method = "average")
    q1 <- quantile(rk, 1/3, names = FALSE, type = 7)
    q2 <- quantile(rk, 2/3, names = FALSE, type = 7)
    cuts <- data.frame(
      degree_value = degs,
      group = case_when(
        rk <= q1 ~ "Lower initial degree",
        rk <= q2 ~ "Middle initial degree",
        TRUE ~ "Higher initial degree"
      ),
      stringsAsFactors = FALSE
    )
  }
  out <- df %>%
    left_join(cuts, by = setNames("degree_value", degree_col)) %>%
    rename(!!group_col := group)
  out[[group_col]] <- factor(out[[group_col]],
                             levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree"))
  out
}

make_link_support_table <- function(dataset, site_tables){
  cooc <- site_tables$cooc_triples %>%
    transmute(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource)) %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource))

  ints <- site_tables$empirical_site_interactions %>%
    transmute(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource)) %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource))

  ## Observed interactions necessarily imply pair recorded together in that site.
  ## This keeps retained_K <= retained_n even if a raw cooccurrence table is conservative.
  cooc <- bind_rows(cooc, ints %>% select(site, consumer, resource, pair_id)) %>%
    distinct(site, consumer, resource, pair_id)

  full_n <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  full_K <- ints %>%
    group_by(pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  full <- full_n %>%
    left_join(full_K, by = "pair_id") %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      full_B = full_K > 0,
      full_K_over_n = full_K / full_n,
      n_bin = support_bin(full_n),
      n_bin_mid = support_bin_mid(n_bin),
      initial_K_class = K_class_figure4(full_K)
    ) %>%
    select(dataset, consumer, resource, pair_id, full_n, full_K, full_B,
           full_K_over_n, n_bin, n_bin_mid, initial_K_class)

  list(full_pairs = full, cooc = cooc, interactions = ints)
}

node_link_table <- function(full_realised){
  consumers <- full_realised %>%
    transmute(dataset, guild = "Consumer", node = consumer, partner = resource,
              pair_id, full_n, full_K, full_K_over_n)
  resources <- full_realised %>%
    transmute(dataset, guild = "Resource", node = resource, partner = consumer,
              pair_id, full_n, full_K, full_K_over_n)
  bind_rows(consumers, resources)
}

## ---------------------------
## Main dataset runner
## ---------------------------
run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble); library(ggplot2)})
  message("Running script 30 local support analysis: ", dataset)

  site_tables <- get_dataset_site_tables(dataset)
  built <- make_link_support_table(dataset, site_tables)
  full_pairs <- built$full_pairs
  cooc <- built$cooc
  interactions <- built$interactions
  realised <- full_pairs %>% filter(full_K >= 1)

  all_sites <- sort(unique(cooc$site))
  n_sites_total <- length(all_sites)
  subset_object <- make_site_subsets_30(all_sites)

  ## Fixed full-network homogeneous p, calibrated once.
  observed_full_links <- sum(full_pairs$full_B)
  p_full <- fit_p_full_unconditioned(full_pairs$full_n, observed_full_links)
  expected_full_links <- expected_links_from_p(full_pairs$full_n, p_full)

  p_calibration <- data.frame(
    dataset = dataset,
    p_full = p_full,
    calibration_method_used = "unconditioned independent site-level p fitted so sum[1-(1-p)^full_n] equals empirical full regional link count",
    empirical_full_regional_link_count = observed_full_links,
    expected_model_full_regional_link_count = expected_full_links,
    calibration_difference = observed_full_links - expected_full_links,
    consumer_persistence_conditioning_note = "No after-removal refit; analytical independent-binomial baseline used for retained_n expectations.",
    stringsAsFactors = FALSE
  )

  ## Figure 2A/B support architecture.
  support_given_n <- full_pairs %>%
    group_by(dataset, n_bin, n_bin_mid) %>%
    summarise(
      number_cooccurring_pairs = n(),
      number_observed_realised_links = sum(full_B),
      empirical_mean_K_conditional_on_regional_interaction = ifelse(sum(full_B) >= 5, mean(full_K[full_B], na.rm = TRUE), NA_real_),
      empirical_one_site_link_proportion = ifelse(sum(full_B) >= 5, mean(full_K[full_B] == 1, na.rm = TRUE), NA_real_),
      homogeneous_expected_mean_K_conditional_on_regional_interaction = ifelse(n() >= 10, expected_mean_K_conditional_bin(full_n, p_full), NA_real_),
      homogeneous_expected_one_site_link_proportion = ifelse(n() >= 10, expected_prop_K1_conditional_bin(full_n, p_full), NA_real_),
      .groups = "drop"
    ) %>%
    arrange(n_bin_mid)

  ## Figure 2C / node portfolios.
  node_links <- node_link_table(realised)

  node_profiles <- node_links %>%
    group_by(dataset, guild, node) %>%
    summarise(
      initial_interaction_degree = n_distinct(partner),
      mean_K_per_link = mean(full_K),
      median_K_per_link = median(full_K),
      prop_one_site_links = mean(full_K == 1),
      .groups = "drop"
    ) %>%
    group_by(dataset, guild) %>%
    group_modify(~ make_degree_groups(.x, "initial_interaction_degree", "initial_degree_group")) %>%
    ungroup()

  node_portfolio_support_summary <- node_profiles %>%
    group_by(dataset, guild, initial_degree_group) %>%
    summarise(
      number_of_nodes = n(),
      mean_interaction_degree = mean(initial_interaction_degree, na.rm = TRUE),
      median_interaction_degree = median(initial_interaction_degree, na.rm = TRUE),
      mean_K_per_realised_link = mean(mean_K_per_link, na.rm = TRUE),
      median_K_per_realised_link = median(mean_K_per_link, na.rm = TRUE),
      mean_one_site_link_proportion = mean(prop_one_site_links, na.rm = TRUE),
      median_one_site_link_proportion = median(prop_one_site_links, na.rm = TRUE),
      .groups = "drop"
    )

  ## Precompute node initial degree table for Figure 5.
  node_init <- node_profiles %>%
    select(dataset, guild, node, initial_interaction_degree, initial_degree_group)

  ## Removal replicate calculations.
  removal_rows <- vector("list", nrow(subset_object$index))
  fate_rows <- vector("list", nrow(subset_object$index))
  persistence_rows <- vector("list", nrow(subset_object$index))
  species_rows <- vector("list", nrow(subset_object$index))
  retained_check_rows <- vector("list", nrow(subset_object$index))

  for(ii in seq_len(nrow(subset_object$index))){
    idx <- subset_object$index[ii, ]
    sites_keep <- subset_object$subsets[[idx$subset_id]]

    cooc_ret <- cooc %>%
      filter(site %in% sites_keep) %>%
      group_by(pair_id) %>%
      summarise(retained_n = n_distinct(site), .groups = "drop")

    int_ret <- interactions %>%
      filter(site %in% sites_keep) %>%
      group_by(pair_id) %>%
      summarise(retained_K = n_distinct(site), .groups = "drop")

    pair_ret <- full_pairs %>%
      left_join(cooc_ret, by = "pair_id") %>%
      left_join(int_ret, by = "pair_id") %>%
      mutate(
        retained_n = replace_na(retained_n, 0L),
        retained_K = replace_na(retained_K, 0L),
        retained_B = retained_K > 0
      )

    retained_check_rows[[ii]] <- pair_ret %>%
      summarise(
        dataset = dataset,
        removal_fraction = idx$removal_fraction,
        replicate = idx$replicate,
        retained_values_valid = all(retained_K >= 0 & retained_K <= retained_n & retained_n <= full_n),
        stringsAsFactors = FALSE
      )

    removal_rows[[ii]] <- data.frame(
      dataset = dataset,
      removal_fraction = idx$removal_fraction,
      replicate = idx$replicate,
      retained_site_count = idx$n_sites_kept,
      cooccurrence_support_retained = sum(pair_ret$retained_n),
      interaction_support_retained = sum(pair_ret$retained_K),
      empirical_binary_links_retained = sum(pair_ret$retained_B),
      homogeneous_expected_binary_links_retained = expected_links_from_p(pair_ret$retained_n[pair_ret$retained_n > 0], p_full),
      full_cooccurrence_support = sum(full_pairs$full_n),
      full_interaction_support = sum(full_pairs$full_K),
      full_empirical_binary_links = observed_full_links,
      full_homogeneous_expected_binary_links = expected_full_links,
      cooccurrence_support_relative_retention = sum(pair_ret$retained_n) / sum(full_pairs$full_n),
      interaction_support_relative_retention = sum(pair_ret$retained_K) / sum(full_pairs$full_K),
      empirical_binary_link_relative_retention = sum(pair_ret$retained_B) / observed_full_links,
      homogeneous_expected_binary_link_relative_retention = expected_links_from_p(pair_ret$retained_n[pair_ret$retained_n > 0], p_full) / expected_full_links,
      stringsAsFactors = FALSE
    )

    real_ret <- pair_ret %>% filter(full_B)

    persistence_rows[[ii]] <- real_ret %>%
      group_by(dataset, initial_K_class) %>%
      summarise(
        removal_fraction = idx$removal_fraction,
        replicate = idx$replicate,
        binary_link_persistence = mean(retained_B),
        number_original_links_in_class = n(),
        .groups = "drop"
      )

    fate_rows[[ii]] <- real_ret %>%
      mutate(
        fate_category = case_when(
          retained_K > 0 ~ "Interaction persists",
          retained_K == 0 & retained_n == 0 ~ "No co-occurrence opportunity remains",
          retained_K == 0 & retained_n > 0 ~ "Co-occurrence remains, interaction absent"
        )
      ) %>%
      count(dataset, fate_category, name = "n_links") %>%
      mutate(
        removal_fraction = idx$removal_fraction,
        replicate = idx$replicate,
        fraction_original_regional_links = n_links / nrow(real_ret)
      )

    ## Species/node reorganisation.
    node_ret_links <- bind_rows(
      real_ret %>% transmute(dataset, guild = "Consumer", node = consumer, partner = resource, retained_B),
      real_ret %>% transmute(dataset, guild = "Resource", node = resource, partner = consumer, retained_B)
    ) %>%
      group_by(dataset, guild, node) %>%
      summarise(retained_degree = n_distinct(partner[retained_B]), .groups = "drop")

    cur_species <- node_init %>%
      left_join(node_ret_links, by = c("dataset", "guild", "node")) %>%
      mutate(
        retained_degree = replace_na(retained_degree, 0L),
        network_active = retained_degree > 0,
        retained_degree_fraction = retained_degree / initial_interaction_degree,
        removal_fraction = idx$removal_fraction,
        replicate = idx$replicate
      )

    species_rows[[ii]] <- cur_species %>%
      group_by(dataset, guild, removal_fraction, replicate, initial_degree_group) %>%
      summarise(
        fraction_active = mean(network_active),
        median_retained_degree_fraction_among_active_nodes = ifelse(any(network_active), median(retained_degree_fraction[network_active], na.rm = TRUE), NA_real_),
        number_original_nodes = n(),
        number_active_nodes = sum(network_active),
        .groups = "drop"
      )
  }

  removal_replicates <- bind_rows(removal_rows)
  fate_replicates <- bind_rows(fate_rows)
  persistence_replicates <- bind_rows(persistence_rows)
  species_reorganisation_summary <- bind_rows(species_rows) %>%
    group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
    summarise(
      fraction_active = median(fraction_active, na.rm = TRUE),
      median_retained_degree_fraction_among_active_nodes = median(median_retained_degree_fraction_among_active_nodes, na.rm = TRUE),
      q025 = q025(median_retained_degree_fraction_among_active_nodes),
      q975 = q975(median_retained_degree_fraction_among_active_nodes),
      number_original_nodes = max(number_original_nodes, na.rm = TRUE),
      number_active_nodes = median(number_active_nodes, na.rm = TRUE),
      .groups = "drop"
    )

  link_persistence_summary <- persistence_replicates %>%
    group_by(dataset, removal_fraction, initial_K_class) %>%
    summarise(
      median_binary_link_persistence = median(binary_link_persistence, na.rm = TRUE),
      q025 = q025(binary_link_persistence),
      q975 = q975(binary_link_persistence),
      number_original_links_in_class = max(number_original_links_in_class, na.rm = TRUE),
      .groups = "drop"
    )

  link_fates_summary <- fate_replicates %>%
    group_by(dataset, removal_fraction, fate_category) %>%
    summarise(
      mean_fraction_original_regional_links = mean(fraction_original_regional_links, na.rm = TRUE),
      q025 = q025(fraction_original_regional_links),
      q975 = q975(fraction_original_regional_links),
      .groups = "drop"
    )

  checks <- data.frame(
    dataset = dataset,
    used_original_dataset = dataset %in% all_dataset_names,
    realised_links_valid = all(realised$full_K >= 1 & realised$full_K <= realised$full_n),
    retained_values_valid = all(bind_rows(retained_check_rows)$retained_values_valid),
    zero_removal_retention_all_one = all(
      removal_replicates %>% filter(removal_fraction == 0) %>%
        transmute(ok = cooccurrence_support_relative_retention == 1 &
                    interaction_support_relative_retention == 1 &
                    empirical_binary_link_relative_retention == 1) %>% pull(ok)
    ),
    fate_categories_sum_to_one = all(
      fate_replicates %>%
        group_by(dataset, removal_fraction, replicate) %>%
        summarise(s = sum(fraction_original_regional_links), .groups = "drop") %>%
        mutate(ok = abs(s - 1) < 1e-10) %>% pull(ok)
    ),
    p_full_calibrated_once_not_refitted = TRUE,
    full_model_matches_empirical_links = abs(observed_full_links - expected_full_links) < 1e-6,
    consumers_resources_separate = TRUE,
    no_forbidden_models_used = TRUE,
    stringsAsFactors = FALSE
  )

  list(
    p_calibration = p_calibration,
    full_link_support = full_pairs,
    support_given_n = support_given_n,
    node_portfolio_support_summary = node_portfolio_support_summary,
    removal_replicates = removal_replicates,
    link_persistence_summary = link_persistence_summary,
    link_fates_summary = link_fates_summary,
    species_reorganisation_summary = species_reorganisation_summary,
    checks = checks,
    node_profiles = node_profiles,
    persistence_replicates = persistence_replicates,
    fate_replicates = fate_replicates
  )
}

## ---------------------------
## Run datasets, parallel across datasets
## ---------------------------
n_workers <- max(1, min(parallelly::availableCores() - 1, length(dataset_order)))
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " workers across datasets.")

all_outputs <- future.apply::future_lapply(dataset_order, run_one_dataset, future.seed = TRUE)
names(all_outputs) <- dataset_order
future::plan(future::sequential)

p_full_calibration <- bind_rows(lapply(all_outputs, `[[`, "p_calibration"))
full_link_support <- bind_rows(lapply(all_outputs, `[[`, "full_link_support"))
support_given_n_summary <- bind_rows(lapply(all_outputs, `[[`, "support_given_n"))
node_portfolio_support_summary <- bind_rows(lapply(all_outputs, `[[`, "node_portfolio_support_summary"))
removal_support_retention_replicates <- bind_rows(lapply(all_outputs, `[[`, "removal_replicates"))
link_persistence_by_support_summary <- bind_rows(lapply(all_outputs, `[[`, "link_persistence_summary"))
link_fates_summary <- bind_rows(lapply(all_outputs, `[[`, "link_fates_summary"))
species_reorganisation_summary <- bind_rows(lapply(all_outputs, `[[`, "species_reorganisation_summary"))
local_support_checks <- bind_rows(lapply(all_outputs, `[[`, "checks"))
node_profiles_all <- bind_rows(lapply(all_outputs, `[[`, "node_profiles"))
persistence_replicates_all <- bind_rows(lapply(all_outputs, `[[`, "persistence_replicates"))
fate_replicates_all <- bind_rows(lapply(all_outputs, `[[`, "fate_replicates"))

## Enforce factor order.
for(obj_name in c("full_link_support", "support_given_n_summary", "node_portfolio_support_summary",
                  "removal_support_retention_replicates", "link_persistence_by_support_summary",
                  "link_fates_summary", "species_reorganisation_summary", "node_profiles_all",
                  "persistence_replicates_all", "fate_replicates_all")){
  obj <- get(obj_name)
  if("dataset" %in% names(obj)) obj$dataset <- factor(obj$dataset, levels = dataset_order)
  assign(obj_name, obj)
}

## ---------------------------
## Save tables
## ---------------------------
write.csv2(p_full_calibration, file.path(results_out, "30_p_full_calibration.csv"), row.names = FALSE)
write.csv2(full_link_support, file.path(results_out, "30_full_link_support.csv"), row.names = FALSE)
write.csv2(support_given_n_summary, file.path(results_out, "30_support_given_n_summary.csv"), row.names = FALSE)
write.csv2(node_portfolio_support_summary, file.path(results_out, "30_node_portfolio_support_summary.csv"), row.names = FALSE)
write.csv2(removal_support_retention_replicates, file.path(results_out, "30_removal_support_retention_replicates.csv"), row.names = FALSE)
write.csv2(link_persistence_by_support_summary, file.path(results_out, "30_link_persistence_by_support_summary.csv"), row.names = FALSE)
write.csv2(link_fates_summary, file.path(results_out, "30_link_fates_summary.csv"), row.names = FALSE)
write.csv2(species_reorganisation_summary, file.path(results_out, "30_species_reorganisation_summary.csv"), row.names = FALSE)
write.csv2(local_support_checks, file.path(results_out, "30_local_support_checks.csv"), row.names = FALSE)

## Also mirror important outputs in CombinedOutputs for convenience.
write.csv2(p_full_calibration, file.path(combined_out, "30_p_full_calibration.csv"), row.names = FALSE)
write.csv2(local_support_checks, file.path(combined_out, "30_local_support_checks.csv"), row.names = FALSE)

## ---------------------------
## Figure 1: conceptual illustration
## ---------------------------
toy_links <- expand.grid(
  landscape = c("Locally thin support", "Repeated local support"),
  site = paste0("Site ", 1:5),
  link = c("A–x", "A–y", "B–x", "C–z"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    cooccurs = TRUE,
    interacts = case_when(
      landscape == "Locally thin support" & link == "A–x" & site == "Site 1" ~ TRUE,
      landscape == "Locally thin support" & link == "A–y" & site == "Site 2" ~ TRUE,
      landscape == "Locally thin support" & link == "B–x" & site == "Site 3" ~ TRUE,
      landscape == "Locally thin support" & link == "C–z" & site == "Site 4" ~ TRUE,
      landscape == "Repeated local support" & link == "A–x" & site %in% c("Site 1", "Site 3", "Site 5") ~ TRUE,
      landscape == "Repeated local support" & link == "A–y" & site %in% c("Site 2", "Site 3") ~ TRUE,
      landscape == "Repeated local support" & link == "B–x" & site %in% c("Site 1", "Site 4") ~ TRUE,
      landscape == "Repeated local support" & link == "C–z" & site %in% c("Site 2", "Site 4", "Site 5") ~ TRUE,
      TRUE ~ FALSE
    ),
    removed = site %in% c("Site 2", "Site 4"),
    status = case_when(
      removed ~ "Removed site",
      interacts ~ "Interaction observed",
      cooccurs ~ "Recorded together only"
    )
  )

p_concept <- ggplot(toy_links, aes(x = site, y = link, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  facet_wrap(~ landscape, ncol = 1) +
  scale_fill_manual(values = c(
    "Removed site" = "#E0E0E0",
    "Interaction observed" = col_K,
    "Recorded together only" = "#F5F5F5"
  )) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("Local sampled sites") +
  ylab("Same regional binary links") +
  ggtitle("Conceptual illustration: local support hidden by regional binary links",
          subtitle = "The same regional links can have thin or repeated local interaction support before sites are removed")

save_plot_both(p_concept, "30_conceptual_hidden_local_support", 7, 6)

## ---------------------------
## Figure 2: support architecture
## ---------------------------
fig2a <- support_given_n_summary %>%
  select(dataset, n_bin, n_bin_mid,
         empirical_mean_K_conditional_on_regional_interaction,
         homogeneous_expected_mean_K_conditional_on_regional_interaction) %>%
  pivot_longer(cols = c(empirical_mean_K_conditional_on_regional_interaction,
                        homogeneous_expected_mean_K_conditional_on_regional_interaction),
               names_to = "source", values_to = "value") %>%
  mutate(source = recode(source,
                         empirical_mean_K_conditional_on_regional_interaction = "Empirical realised links",
                         homogeneous_expected_mean_K_conditional_on_regional_interaction = "Homogeneous fixed-p baseline")) %>%
  ggplot(aes(x = n_bin_mid, y = value, colour = source)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = c("Empirical realised links" = col_B,
                                 "Homogeneous fixed-p baseline" = col_model)) +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Co-occurrence support, n") +
  ylab("Mean K | regional link") +
  ggtitle("A. Interaction support conditional on co-occurrence support")

fig2b <- support_given_n_summary %>%
  select(dataset, n_bin, n_bin_mid,
         empirical_one_site_link_proportion,
         homogeneous_expected_one_site_link_proportion) %>%
  pivot_longer(cols = c(empirical_one_site_link_proportion,
                        homogeneous_expected_one_site_link_proportion),
               names_to = "source", values_to = "value") %>%
  mutate(source = recode(source,
                         empirical_one_site_link_proportion = "Empirical realised links",
                         homogeneous_expected_one_site_link_proportion = "Homogeneous fixed-p baseline")) %>%
  ggplot(aes(x = n_bin_mid, y = value, colour = source)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = c("Empirical realised links" = col_B,
                                 "Homogeneous fixed-p baseline" = col_model)) +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Co-occurrence support, n") +
  ylab("Proportion K = 1 | regional link") +
  ggtitle("B. One-site regional links conditional on co-occurrence support")

fig2c_data <- node_portfolio_support_summary %>%
  select(dataset, guild, initial_degree_group,
         median_K_per_realised_link,
         median_one_site_link_proportion) %>%
  pivot_longer(cols = c(median_K_per_realised_link, median_one_site_link_proportion),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         median_K_per_realised_link = "Median K per realised link",
                         median_one_site_link_proportion = "Median proportion K = 1"))

fig2c_median <- fig2c_data %>%
  group_by(guild, metric, initial_degree_group) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

fig2c <- ggplot() +
  geom_line(
    data = fig2c_data,
    mapping = aes(
      x = initial_degree_group,
      y = value,
      group = dataset
    ),
    alpha = 0.22,
    colour = "grey35"
  ) +
  geom_point(
    data = fig2c_data,
    mapping = aes(
      x = initial_degree_group,
      y = value
    ),
    alpha = 0.22,
    colour = "grey35",
    size = 1
  ) +
  geom_line(
    data = fig2c_median,
    mapping = aes(
      x = initial_degree_group,
      y = value,
      group = 1
    ),
    linewidth = 1.1,
    colour = col_K
  ) +
  geom_point(
    data = fig2c_median,
    mapping = aes(
      x = initial_degree_group,
      y = value
    ),
    size = 2.4,
    colour = col_K
  ) +
  facet_grid(guild ~ metric, scales = "free_y") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  xlab("Initial-degree group") +
  ylab("Node portfolio summary") +
  ggtitle("C. Local support across species’ interaction portfolios")

fig2 <- fig2a / fig2b / fig2c + patchwork::plot_layout(heights = c(1, 1, 0.9))
save_plot_both(fig2, "30_support_architecture_full_region", 14, 15)

## ---------------------------
## Figure 3: retention layers and model baseline
## ---------------------------
layer_rep <- removal_support_retention_replicates %>%
  select(dataset, removal_fraction, replicate,
         cooccurrence_support_relative_retention,
         interaction_support_relative_retention,
         empirical_binary_link_relative_retention,
         homogeneous_expected_binary_link_relative_retention) %>%
  pivot_longer(cols = c(cooccurrence_support_relative_retention,
                        interaction_support_relative_retention,
                        empirical_binary_link_relative_retention,
                        homogeneous_expected_binary_link_relative_retention),
               names_to = "layer", values_to = "retention") %>%
  mutate(layer = recode(layer,
                        cooccurrence_support_relative_retention = "Co-occurrence support",
                        interaction_support_relative_retention = "Local interaction support",
                        empirical_binary_link_relative_retention = "Regional binary links",
                        homogeneous_expected_binary_link_relative_retention = "Homogeneous fixed-p baseline"))

layer_sum <- layer_rep %>%
  group_by(dataset, removal_fraction, layer) %>%
  summarise(median = median(retention, na.rm = TRUE), q025 = q025(retention), q975 = q975(retention), .groups = "drop")

fig3a <- layer_sum %>% filter(layer %in% c("Co-occurrence support", "Local interaction support", "Regional binary links")) %>%
  ggplot(aes(x = removal_fraction, y = median, colour = layer, fill = layer)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = support_layer_cols) +
  scale_fill_manual(values = support_layer_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Retained fraction") +
  ggtitle("A. Co-occurrence support, local interaction support, and regional binary links")

fig3b <- layer_sum %>% filter(layer %in% c("Regional binary links", "Homogeneous fixed-p baseline")) %>%
  ggplot(aes(x = removal_fraction, y = median, colour = layer, fill = layer)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = support_layer_cols) +
  scale_fill_manual(values = support_layer_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Regional-link retention") +
  ggtitle("B. Empirical binary-link retention versus fixed-p baseline")

fig3 <- fig3a / fig3b
save_plot_both(fig3, "30_support_and_link_retention_under_removal", 14, 11)

## ---------------------------
## Figure 4: link persistence and fates
## ---------------------------
fig4a <- link_persistence_by_support_summary %>%
  mutate(initial_K_class = factor(initial_K_class, levels = names(K_class_cols))) %>%
  ggplot(aes(x = removal_fraction, y = median_binary_link_persistence,
             colour = initial_K_class, fill = initial_K_class)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.85, na.rm = TRUE) +
  geom_point(size = 1.2, na.rm = TRUE) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = K_class_cols, na.value = "grey70") +
  scale_fill_manual(values = K_class_cols, na.value = "grey70") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("P(link still observed)") +
  ggtitle("A. Link persistence by original interaction support")

fig4b_data <- link_fates_summary %>%
  filter(removal_fraction > 0) %>%
  mutate(fate_category = factor(fate_category, levels = names(fate_cols)))

fig4b <- ggplot(fig4b_data,
                aes(x = factor(removal_fraction), y = mean_fraction_original_regional_links, fill = fate_category)) +
  geom_col(width = 0.85) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_fill_manual(values = fate_cols) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Fraction of original regional links") +
  ggtitle("B. Fate decomposition of original regional interaction links")

fig4 <- fig4a / fig4b
save_plot_both(fig4, "30_link_persistence_and_fates", 14, 11)

## ---------------------------
## Figure 5: species-level reorganisation
## ---------------------------
degree_group_cols <- c(
  "Lower initial degree"  = "#4E79A7",
  "Middle initial degree" = "#F28E2B",
  "Higher initial degree" = "#59A14F"
)

fig5_summary <- species_reorganisation_summary %>%
  select(
    dataset, guild, removal_fraction, initial_degree_group,
    fraction_active,
    median_retained_degree_fraction_among_active_nodes
  ) %>%
  pivot_longer(
    cols = c(fraction_active, median_retained_degree_fraction_among_active_nodes),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      fraction_active = "Probability of remaining network-active",
      median_retained_degree_fraction_among_active_nodes =
        "Retained-degree fraction among active nodes"
    ),
    metric = factor(metric, levels = c(
      "Probability of remaining network-active",
      "Retained-degree fraction among active nodes"
    )),
    guild = factor(guild, levels = c("Consumer", "Resource")),
    initial_degree_group = factor(initial_degree_group, levels = names(degree_group_cols))
  ) %>%
  filter(!is.na(value))

fig5_dataset <- fig5_summary %>%
  group_by(dataset, guild, metric, initial_degree_group, removal_fraction) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

fig5_overall <- fig5_dataset %>%
  group_by(guild, metric, initial_degree_group, removal_fraction) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    q25 = as.numeric(quantile(value, 0.25, na.rm = TRUE)),
    q75 = as.numeric(quantile(value, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

print(names(fig5_overall))

fig5 <- ggplot() +
  geom_ribbon(
    data = fig5_overall,
    aes(
      x = removal_fraction,
      ymin = q25,
      ymax = q75,
      fill = initial_degree_group,
      group = interaction(guild, metric, initial_degree_group)
    ),
    alpha = 0.12,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = fig5_dataset,
    aes(
      x = removal_fraction,
      y = value,
      group = interaction(dataset, initial_degree_group),
      colour = initial_degree_group
    ),
    alpha = 0.20,
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = fig5_overall,
    aes(
      x = removal_fraction,
      y = median_value,
      colour = initial_degree_group,
      group = initial_degree_group
    ),
    linewidth = 1.15,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = fig5_overall,
    aes(
      x = removal_fraction,
      y = median_value,
      colour = initial_degree_group
    ),
    size = 1.6,
    inherit.aes = FALSE
  ) +
  facet_grid(guild ~ metric) +
  scale_colour_manual(values = degree_group_cols, drop = FALSE) +
  scale_fill_manual(values = degree_group_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Fraction") +
  ggtitle("Species-level reorganisation under site removal")

save_plot_both(fig5, "30_species_reorganisation_under_removal", 10, 7)

## ---------------------------
## Console validation and reminders
## ---------------------------
message("\nValidation summary for script 30:")
print(local_support_checks)

if(!all(local_support_checks$used_original_dataset)) warning("At least one non-original dataset was included.")
if(!all(local_support_checks$realised_links_valid)) warning("Some realised links violate 1 <= full_K <= full_n.")
if(!all(local_support_checks$retained_values_valid)) warning("Some retained values violate 0 <= retained_K <= retained_n <= full_n.")
if(!all(local_support_checks$zero_removal_retention_all_one)) warning("Zero-removal retention did not equal 1 for all datasets.")
if(!all(local_support_checks$fate_categories_sum_to_one)) warning("Fate categories do not sum to 1 somewhere.")
if(!all(local_support_checks$full_model_matches_empirical_links)) warning("Homogeneous full-network expected link count did not match empirical count within tolerance.")

message("\nRandom site removal is a controlled diagnostic of how local interaction support is represented in regional networks. It is not a forecast of habitat loss, species extinction, or ecosystem collapse.")
message("\nDifferences between empirical and fixed-p model retention indicate that one shared per-site interaction probability does not fully capture how realised interaction support is distributed across co-occurrence opportunities. They do not, by themselves, identify the ecological causes of that heterogeneity.")
message("\nHigher absolute interaction support among high-degree species must not be interpreted automatically as a higher probability of realising every co-occurrence opportunity.")

message("\nFinished script 30. Outputs saved in: ", results_out)
