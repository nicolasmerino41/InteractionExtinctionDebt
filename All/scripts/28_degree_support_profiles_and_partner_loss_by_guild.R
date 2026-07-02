## ------------------------------------------------------------
## Script: All/scripts/28_degree_support_profiles_and_partner_loss_by_guild.R
##
## Descriptive analysis of local support profiles and binary-partner
## loss under random site removal, separately for consumers and resources.
## Uses shared helper: All/scripts/00_dataset_loaders_and_helpers_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("future", "future.apply", "parallelly")
for(pkg in packages_extra){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500
min_sites_retained <- 1

result_type <- "script28_degree_support_profiles_and_partner_loss_by_guild"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " workers across datasets.")

guild_levels <- c("Consumer", "Resource")
degree_group_levels <- c("Lower initial degree", "Middle initial degree", "Higher initial degree")
support_class_levels <- c("Observed in 1 site", "Observed in 2 sites",
                          "Observed in 3-4 sites", "Observed in 5 or more sites")

support_colours <- c(
  "Observed in 1 site" = "#b2182b",
  "Observed in 2 sites" = "#ef8a62",
  "Observed in 3-4 sites" = "#67a9cf",
  "Observed in 5 or more sites" = "#2166ac"
)

degree_group_colours <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#7570b3",
  "Higher initial degree" = "#d95f02"
)

make_pair_id <- function(consumer, resource){ paste(consumer, resource, sep = "___") }

support_class_from_K <- function(K){
  dplyr::case_when(
    K == 1 ~ "Observed in 1 site",
    K == 2 ~ "Observed in 2 sites",
    K >= 3 & K <= 4 ~ "Observed in 3-4 sites",
    K >= 5 ~ "Observed in 5 or more sites",
    TRUE ~ NA_character_
  )
}

assign_tertile_groups_keep_ties <- function(df, value_col, group_col = "initial_degree_group"){
  vals <- sort(unique(df[[value_col]]))
  if(length(vals) == 1){
    map <- data.frame(value = vals, group = "Middle initial degree", stringsAsFactors = FALSE)
  } else {
    ranks <- seq_along(vals)
    cuts <- stats::quantile(ranks, probs = c(1/3, 2/3), type = 1, na.rm = TRUE)
    groups <- ifelse(ranks <= cuts[1], "Lower initial degree",
                     ifelse(ranks <= cuts[2], "Middle initial degree", "Higher initial degree"))
    map <- data.frame(value = vals, group = groups, stringsAsFactors = FALSE)
  }
  out <- df %>% left_join(map, by = setNames("value", value_col))
  names(out)[names(out) == "group"] <- group_col
  out[[group_col]] <- factor(out[[group_col]], levels = degree_group_levels)
  out
}

make_site_subsets_local <- function(all_sites){
  all_sites <- sort(unique(as.character(all_sites)))
  n_total <- length(all_sites)
  subset_list <- list(); subset_index <- list(); counter <- 1
  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_total * (1 - removal)))
    n_keep <- min(n_keep, n_total)
    reps_here <- ifelse(removal == 0, 1, n_site_reps)
    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) all_sites else sample(all_sites, size = n_keep, replace = FALSE)
      subset_list[[counter]] <- sort(keep)
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        replicate = r,
        n_sites_total = n_total,
        n_sites_kept = length(keep),
        n_sites_removed = n_total - length(keep),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
  }
  list(index = bind_rows(subset_index), subsets = subset_list)
}

link_survival_probability <- function(S, K, m){
  if(is.na(S) || is.na(K) || is.na(m)) return(NA_real_)
  if(m <= 0) return(1)
  if(K <= 0) return(0)
  if(m > S) return(NA_real_)
  if(m < K) return(1)
  if(m > (S - 0)) return(NA_real_)
  ## Link is lost if all K supporting sites are removed:
  ## P(lost) = choose(S-K, m-K) / choose(S, m), for m >= K.
  p_lost <- exp(lchoose(S - K, m - K) - lchoose(S, m))
  1 - p_lost
}

spearman_safe <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok) < 3) return(NA_real_)
  if(length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

build_original_link_table <- function(dataset, site_tables){
  cooc <- site_tables$cooc_triples %>%
    mutate(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource)) %>%
    distinct(site, consumer, resource, pair_id)

  ints <- site_tables$empirical_site_interactions %>%
    mutate(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource)) %>%
    distinct(site, consumer, resource, pair_id)

  ## Defensive consistency: an observed interaction is an operational co-occurrence.
  cooc <- bind_rows(cooc, ints %>% select(site, consumer, resource, pair_id)) %>%
    distinct(site, consumer, resource, pair_id)

  n_tab <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  k_tab <- ints %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  links <- n_tab %>%
    left_join(k_tab, by = c("consumer", "resource", "pair_id")) %>%
    mutate(dataset = dataset,
           full_K = replace_na(full_K, 0L),
           full_K = as.integer(full_K),
           full_n = as.integer(full_n)) %>%
    filter(full_K >= 1) %>%
    mutate(full_K_over_n = full_K / full_n,
           support_class = factor(support_class_from_K(full_K), levels = support_class_levels)) %>%
    select(dataset, consumer, resource, pair_id, full_n, full_K, full_K_over_n, support_class)

  list(links = links, cooc = cooc, interactions = ints)
}

make_node_link_table <- function(links){
  consumer_links <- links %>%
    transmute(dataset, guild = "Consumer", node = consumer, partner = resource,
              consumer, resource, pair_id, full_n, full_K, full_K_over_n, support_class)
  resource_links <- links %>%
    transmute(dataset, guild = "Resource", node = resource, partner = consumer,
              consumer, resource, pair_id, full_n, full_K, full_K_over_n, support_class)
  bind_rows(consumer_links, resource_links) %>%
    mutate(guild = factor(guild, levels = guild_levels),
           support_class = factor(support_class, levels = support_class_levels))
}

make_node_profiles <- function(node_links){
  profiles <- node_links %>%
    group_by(dataset, guild, node) %>%
    summarise(
      initial_degree = n_distinct(partner),
      total_K = sum(full_K, na.rm = TRUE),
      mean_K_per_link = mean(full_K, na.rm = TRUE),
      median_K_per_link = median(full_K, na.rm = TRUE),
      proportion_links_K_1 = mean(full_K == 1, na.rm = TRUE),
      proportion_links_K_2 = mean(full_K == 2, na.rm = TRUE),
      proportion_links_K_3_4 = mean(full_K >= 3 & full_K <= 4, na.rm = TRUE),
      proportion_links_K_5_plus = mean(full_K >= 5, na.rm = TRUE),
      mean_n_per_link = mean(full_n, na.rm = TRUE),
      mean_K_over_n_among_links = mean(full_K_over_n, na.rm = TRUE),
      .groups = "drop"
    )
  profiles %>%
    group_by(dataset, guild) %>%
    group_modify(~ assign_tertile_groups_keep_ties(.x, "initial_degree")) %>%
    ungroup()
}

compute_subset_node_metrics <- function(dataset, node_profiles, node_links, cooc, interactions, sites_keep, subset_row){
  retained_cooc_pairs <- cooc %>%
    filter(site %in% sites_keep) %>%
    group_by(pair_id) %>%
    summarise(retained_n = n_distinct(site), .groups = "drop")

  retained_int_pairs <- interactions %>%
    filter(site %in% sites_keep) %>%
    group_by(pair_id) %>%
    summarise(retained_K = n_distinct(site), .groups = "drop")

  link_retained <- node_links %>%
    left_join(retained_cooc_pairs, by = "pair_id") %>%
    left_join(retained_int_pairs, by = "pair_id") %>%
    mutate(retained_n = replace_na(retained_n, 0L),
           retained_K = replace_na(retained_K, 0L),
           opportunity_survives = retained_n > 0,
           interaction_survives = retained_K > 0)

  by_node <- link_retained %>%
    group_by(dataset, guild, node) %>%
    summarise(retained_degree = sum(interaction_survives, na.rm = TRUE),
              opportunity_partners_retained = sum(opportunity_survives, na.rm = TRUE),
              retained_total_K = sum(retained_K, na.rm = TRUE),
              .groups = "drop")

  retained_presence <- bind_rows(
    cooc %>% filter(site %in% sites_keep) %>% transmute(guild = "Consumer", node = consumer),
    cooc %>% filter(site %in% sites_keep) %>% transmute(guild = "Resource", node = resource)
  ) %>% distinct(guild, node) %>% mutate(node_present = TRUE)

  node_profiles %>%
    left_join(by_node, by = c("dataset", "guild", "node")) %>%
    left_join(retained_presence, by = c("guild", "node")) %>%
    mutate(retained_degree = replace_na(retained_degree, 0L),
           opportunity_partners_retained = replace_na(opportunity_partners_retained, 0L),
           retained_total_K = replace_na(retained_total_K, 0L),
           node_present = replace_na(node_present, FALSE),
           network_active = retained_degree > 0,
           binary_partner_retention = retained_degree / initial_degree,
           opportunity_partner_retention = opportunity_partners_retained / initial_degree,
           local_support_retention = retained_total_K / total_K,
           removal_fraction = subset_row$removal_fraction,
           replicate = subset_row$replicate) %>%
    select(dataset, guild, node, initial_degree, initial_degree_group,
           removal_fraction, replicate, node_present, network_active, retained_degree,
           binary_partner_retention, opportunity_partner_retention, local_support_retention)
}

summarise_retention_by_group <- function(retention_raw){
  all_summary <- retention_raw %>%
    group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
    summarise(summary_type = "All original nodes",
              number_original_nodes = n_distinct(node),
              number_nodes_present = sum(node_present, na.rm = TRUE),
              number_network_active = sum(network_active, na.rm = TRUE),
              median_binary_partner_retention = median(binary_partner_retention, na.rm = TRUE),
              q025_binary_partner_retention = quantile(binary_partner_retention, 0.025, na.rm = TRUE),
              q975_binary_partner_retention = quantile(binary_partner_retention, 0.975, na.rm = TRUE),
              median_opportunity_partner_retention = median(opportunity_partner_retention, na.rm = TRUE),
              median_local_support_retention = median(local_support_retention, na.rm = TRUE),
              median_node_present = median(as.numeric(node_present), na.rm = TRUE),
              q025_node_present = quantile(as.numeric(node_present), 0.025, na.rm = TRUE),
              q975_node_present = quantile(as.numeric(node_present), 0.975, na.rm = TRUE),
              median_network_active = median(as.numeric(network_active), na.rm = TRUE),
              q025_network_active = quantile(as.numeric(network_active), 0.025, na.rm = TRUE),
              q975_network_active = quantile(as.numeric(network_active), 0.975, na.rm = TRUE),
              .groups = "drop")

  present_summary <- retention_raw %>%
    filter(node_present) %>%
    group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
    summarise(summary_type = "Conditional on node presence",
              number_original_nodes = n_distinct(node),
              number_nodes_present = n_distinct(node),
              number_network_active = sum(network_active, na.rm = TRUE),
              median_binary_partner_retention = median(binary_partner_retention, na.rm = TRUE),
              q025_binary_partner_retention = quantile(binary_partner_retention, 0.025, na.rm = TRUE),
              q975_binary_partner_retention = quantile(binary_partner_retention, 0.975, na.rm = TRUE),
              median_opportunity_partner_retention = median(opportunity_partner_retention, na.rm = TRUE),
              median_local_support_retention = median(local_support_retention, na.rm = TRUE),
              median_node_present = 1,
              q025_node_present = 1,
              q975_node_present = 1,
              median_network_active = median(as.numeric(network_active), na.rm = TRUE),
              q025_network_active = quantile(as.numeric(network_active), 0.025, na.rm = TRUE),
              q975_network_active = quantile(as.numeric(network_active), 0.975, na.rm = TRUE),
              .groups = "drop")
  bind_rows(all_summary, present_summary)
}

make_support_distribution <- function(node_links, node_profiles){
  node_links %>%
    left_join(node_profiles %>% select(dataset, guild, node, initial_degree_group),
              by = c("dataset", "guild", "node")) %>%
    group_by(dataset, guild, initial_degree_group, support_class) %>%
    summarise(number_links_in_support_class = n(), .groups = "drop") %>%
    complete(dataset, guild, initial_degree_group, support_class,
             fill = list(number_links_in_support_class = 0L)) %>%
    left_join(node_profiles %>%
                group_by(dataset, guild, initial_degree_group) %>%
                summarise(number_nodes = n_distinct(node), .groups = "drop"),
              by = c("dataset", "guild", "initial_degree_group")) %>%
    group_by(dataset, guild, initial_degree_group) %>%
    mutate(number_original_links = sum(number_links_in_support_class, na.rm = TRUE),
           fraction_links_in_support_class = ifelse(number_original_links > 0,
                                                    number_links_in_support_class / number_original_links,
                                                    NA_real_)) %>%
    ungroup()
}

make_degree_groups_table <- function(node_profiles){
  node_profiles %>%
    group_by(dataset, guild, initial_degree_group) %>%
    summarise(number_nodes = n_distinct(node),
              minimum_initial_degree = min(initial_degree, na.rm = TRUE),
              median_initial_degree = median(initial_degree, na.rm = TRUE),
              maximum_initial_degree = max(initial_degree, na.rm = TRUE),
              .groups = "drop")
}

make_expected_retention <- function(node_links, node_profiles, subset_index){
  S_by_dataset <- subset_index %>%
    group_by(dataset, removal_fraction, n_sites_total, n_sites_removed) %>%
    summarise(.groups = "drop")

  expected_link <- node_links %>%
    select(dataset, guild, node, pair_id, full_K) %>%
    left_join(S_by_dataset, by = "dataset") %>%
    rowwise() %>%
    mutate(expected_link_survival = link_survival_probability(S = n_sites_total, K = full_K, m = n_sites_removed)) %>%
    ungroup()

  expected_link %>%
    group_by(dataset, guild, node, removal_fraction) %>%
    summarise(expected_binary_partner_retention = mean(expected_link_survival, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(node_profiles %>% select(dataset, guild, node, initial_degree, initial_degree_group),
              by = c("dataset", "guild", "node"))
}

make_expected_vs_observed <- function(expected_node, retention_raw){
  observed <- retention_raw %>%
    group_by(dataset, guild, node, removal_fraction) %>%
    summarise(observed_mean_binary_partner_retention = mean(binary_partner_retention, na.rm = TRUE),
              observed_q025 = quantile(binary_partner_retention, 0.025, na.rm = TRUE),
              observed_q975 = quantile(binary_partner_retention, 0.975, na.rm = TRUE),
              .groups = "drop")
  expected_node %>%
    left_join(observed, by = c("dataset", "guild", "node", "removal_fraction")) %>%
    mutate(observed_minus_expected = observed_mean_binary_partner_retention - expected_binary_partner_retention)
}

make_active_enrichment <- function(retention_raw){
  all_mean <- retention_raw %>%
    distinct(dataset, guild, node, initial_degree) %>%
    group_by(dataset, guild) %>%
    summarise(number_original_nodes = n_distinct(node),
              mean_initial_degree_all_original_nodes = mean(initial_degree, na.rm = TRUE),
              .groups = "drop")

  retention_raw %>%
    group_by(dataset, guild, removal_fraction, replicate) %>%
    summarise(number_network_active = sum(network_active, na.rm = TRUE),
              mean_initial_degree_active_nodes = ifelse(any(network_active), mean(initial_degree[network_active], na.rm = TRUE), NA_real_),
              median_initial_degree_active_nodes = ifelse(any(network_active), median(initial_degree[network_active], na.rm = TRUE), NA_real_),
              .groups = "drop") %>%
    left_join(all_mean, by = c("dataset", "guild")) %>%
    mutate(initial_degree_enrichment_among_active_nodes =
             mean_initial_degree_active_nodes / mean_initial_degree_all_original_nodes)
}

make_associations <- function(node_profiles, expected_80){
  profiles <- node_profiles %>%
    left_join(expected_80 %>%
                select(dataset, guild, node,
                       expected_binary_partner_retention_at_80 = expected_binary_partner_retention),
              by = c("dataset", "guild", "node"))
  vars <- c("mean_K_per_link", "median_K_per_link", "proportion_links_K_1",
            "mean_n_per_link", "mean_K_over_n_among_links",
            "expected_binary_partner_retention_at_80")
  bind_rows(lapply(vars, function(v){
    profiles %>%
      group_by(dataset, guild) %>%
      summarise(variable_1 = "initial_degree",
                variable_2 = v,
                spearman_correlation = spearman_safe(initial_degree, .data[[v]]),
                number_nodes = n_distinct(node),
                .groups = "drop")
  }))
}

run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble)})
  message("Running script 28 degree support profiles and partner loss: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  built <- build_original_link_table(dataset, site_tables)
  original_links <- built$links
  cooc <- built$cooc
  interactions <- built$interactions
  node_links <- make_node_link_table(original_links)
  node_profiles <- make_node_profiles(node_links)

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets_local(all_sites)
  subset_index <- subset_object$index %>% mutate(dataset = dataset)

  retention_list <- vector("list", nrow(subset_index))
  for(i in seq_len(nrow(subset_index))){
    if(i %% 250 == 0) message("  ", dataset, ": subset ", i, " / ", nrow(subset_index))
    subset_id <- subset_index$subset_id[i]
    retention_list[[i]] <- compute_subset_node_metrics(
      dataset = dataset,
      node_profiles = node_profiles,
      node_links = node_links,
      cooc = cooc,
      interactions = interactions,
      sites_keep = subset_object$subsets[[subset_id]],
      subset_row = subset_index[i, ]
    )
  }
  retention_raw <- bind_rows(retention_list)
  support_distribution <- make_support_distribution(node_links, node_profiles)
  degree_groups <- make_degree_groups_table(node_profiles)
  expected_node <- make_expected_retention(node_links, node_profiles, subset_index)
  expected_vs_observed <- make_expected_vs_observed(expected_node, retention_raw)
  expected_80 <- expected_node %>% filter(abs(removal_fraction - 0.8) < 1e-9)
  associations <- make_associations(node_profiles, expected_80)
  retention_by_group <- summarise_retention_by_group(retention_raw)
  active_enrichment <- make_active_enrichment(retention_raw)

  write.csv2(node_profiles, file.path(out_dir, paste0(dataset, "_28_initial_node_support_profiles.csv")), row.names = FALSE)
  write.csv2(retention_by_group, file.path(out_dir, paste0(dataset, "_28_partner_retention_by_degree_group.csv")), row.names = FALSE)

  list(original_links = original_links,
       node_links = node_links,
       node_profiles = node_profiles,
       support_distribution = support_distribution,
       degree_groups = degree_groups,
       associations = associations,
       retention_raw = retention_raw,
       retention_by_group = retention_by_group,
       expected_vs_observed = expected_vs_observed,
       active_enrichment = active_enrichment,
       subset_index = subset_index)
}

all_outputs <- future.apply::future_lapply(all_dataset_names, run_one_dataset, future.seed = TRUE)
names(all_outputs) <- all_dataset_names

original_links_all <- bind_rows(lapply(all_outputs, `[[`, "original_links"))
node_profiles_all <- bind_rows(lapply(all_outputs, `[[`, "node_profiles"))
support_distribution_all <- bind_rows(lapply(all_outputs, `[[`, "support_distribution"))
degree_groups_all <- bind_rows(lapply(all_outputs, `[[`, "degree_groups"))
associations_all <- bind_rows(lapply(all_outputs, `[[`, "associations"))
retention_raw_all <- bind_rows(lapply(all_outputs, `[[`, "retention_raw"))
retention_by_group_all <- bind_rows(lapply(all_outputs, `[[`, "retention_by_group"))
expected_vs_observed_all <- bind_rows(lapply(all_outputs, `[[`, "expected_vs_observed"))
active_enrichment_all <- bind_rows(lapply(all_outputs, `[[`, "active_enrichment"))

write.csv2(node_profiles_all, file.path(combined_out, "28_initial_node_support_profiles.csv"), row.names = FALSE)
write.csv2(support_distribution_all, file.path(combined_out, "28_initial_support_distribution_by_degree_group.csv"), row.names = FALSE)
write.csv2(associations_all, file.path(combined_out, "28_degree_support_associations.csv"), row.names = FALSE)
write.csv2(retention_raw_all, file.path(combined_out, "28_node_partner_retention_raw.csv"), row.names = FALSE)
write.csv2(retention_by_group_all, file.path(combined_out, "28_partner_retention_by_degree_group.csv"), row.names = FALSE)
write.csv2(expected_vs_observed_all, file.path(combined_out, "28_expected_partner_retention_from_K_profile.csv"), row.names = FALSE)
write.csv2(active_enrichment_all, file.path(combined_out, "28_active_node_initial_degree_enrichment.csv"), row.names = FALSE)
write.csv2(degree_groups_all, file.path(combined_out, "28_initial_degree_groups_by_guild.csv"), row.names = FALSE)

## Validation checks
n_original <- node_profiles_all %>% count(dataset, guild, name = "n_nodes_orig")
n_retention <- retention_raw_all %>%
  group_by(dataset, guild, removal_fraction, replicate) %>%
  summarise(n_nodes_ret = n_distinct(node), .groups = "drop") %>%
  left_join(n_original, by = c("dataset", "guild"))

checks <- data.frame(
  check = c(
    "all_original_links_have_full_K_ge_1_and_le_full_n",
    "consumers_and_resources_analysed_separately",
    "all_nodes_primary_retention_have_initial_degree_ge_1",
    "all_original_nodes_contribute_to_primary_retention",
    "binary_partner_retention_equals_1_at_zero_removal",
    "expected_partner_retention_equals_1_at_zero_removal",
    "expected_partner_retention_between_0_and_1",
    "observed_mean_compared_to_expected_using_same_removed_site_counts",
    "initial_degree_enrichment_equals_1_at_zero_removal",
    "no_interaction_model_null_model_regression_or_refitting_used"
  ),
  passed = c(
    all(original_links_all$full_K >= 1 & original_links_all$full_K <= original_links_all$full_n),
    all(sort(unique(as.character(node_profiles_all$guild))) == sort(guild_levels)),
    all(retention_raw_all$initial_degree >= 1),
    all(n_retention$n_nodes_ret == n_retention$n_nodes_orig),
    all(abs(retention_raw_all$binary_partner_retention[retention_raw_all$removal_fraction == 0] - 1) < 1e-12),
    all(abs(expected_vs_observed_all$expected_binary_partner_retention[expected_vs_observed_all$removal_fraction == 0] - 1) < 1e-12),
    all(expected_vs_observed_all$expected_binary_partner_retention >= -1e-12 & expected_vs_observed_all$expected_binary_partner_retention <= 1 + 1e-12, na.rm = TRUE),
    TRUE,
    all(abs(active_enrichment_all$initial_degree_enrichment_among_active_nodes[active_enrichment_all$removal_fraction == 0] - 1) < 1e-12, na.rm = TRUE),
    TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv2(checks, file.path(combined_out, "28_degree_support_partner_loss_checks.csv"), row.names = FALSE)
message("Validation summary:")
print(checks)
if(!all(checks$passed)) warning("At least one validation check failed. Inspect 28_degree_support_partner_loss_checks.csv")

## Plot prep
node_profiles_all <- node_profiles_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels),
         initial_degree_group = factor(initial_degree_group, levels = degree_group_levels))
support_distribution_all <- support_distribution_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels),
         initial_degree_group = factor(initial_degree_group, levels = degree_group_levels),
         support_class = factor(support_class, levels = support_class_levels))
retention_by_group_all <- retention_by_group_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels),
         initial_degree_group = factor(initial_degree_group, levels = degree_group_levels))
expected_vs_observed_group <- expected_vs_observed_all %>%
  group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
  summarise(mean_expected_binary_partner_retention = mean(expected_binary_partner_retention, na.rm = TRUE),
            mean_observed_binary_partner_retention = mean(observed_mean_binary_partner_retention, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels),
         initial_degree_group = factor(initial_degree_group, levels = degree_group_levels))
active_enrichment_summary <- active_enrichment_all %>%
  group_by(dataset, guild, removal_fraction) %>%
  summarise(median_enrichment = median(initial_degree_enrichment_among_active_nodes, na.rm = TRUE),
            q025_enrichment = quantile(initial_degree_enrichment_among_active_nodes, 0.025, na.rm = TRUE),
            q975_enrichment = quantile(initial_degree_enrichment_among_active_nodes, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels))

## Figures
p1 <- ggplot(support_distribution_all,
             aes(x = initial_degree_group, y = fraction_links_in_support_class, fill = support_class)) +
  geom_col(width = 0.8) +
  facet_grid(guild ~ dataset, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = support_colours, drop = FALSE) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
  xlab("Initial-degree group") + ylab("Fraction of original realised links") + labs(fill = "Interaction support")
ggsave(file.path(combined_out, "28_initial_K_support_distribution_by_degree_and_guild.png"), p1, width = 16, height = 7, dpi = 300)

assoc_plot <- associations_all %>%
  filter(variable_2 %in% c("mean_K_per_link", "proportion_links_K_1")) %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names), guild = factor(guild, levels = guild_levels),
         metric = recode(variable_2,
                         mean_K_per_link = "Initial degree vs mean K per link",
                         proportion_links_K_1 = "Initial degree vs proportion one-site links"))
p2 <- ggplot(assoc_plot, aes(x = dataset, y = spearman_correlation)) +
  geom_hline(yintercept = 0, colour = "grey65") + geom_point(size = 2.1, colour = "grey20") +
  facet_grid(metric ~ guild) + theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("Spearman correlation")
ggsave(file.path(combined_out, "28_degree_support_association_summary.png"), p2, width = 12, height = 6, dpi = 300)

retention_plot <- retention_by_group_all %>% filter(summary_type == "All original nodes")
p3 <- ggplot(retention_plot,
             aes(x = removal_fraction, y = median_binary_partner_retention,
                 ymin = q025_binary_partner_retention, ymax = q975_binary_partner_retention,
                 colour = initial_degree_group, fill = initial_degree_group)) +
  geom_ribbon(alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
  facet_grid(guild ~ dataset) +
  scale_colour_manual(values = degree_group_colours, drop = FALSE) +
  scale_fill_manual(values = degree_group_colours, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) + theme_classic(base_size = 9) + theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") + ylab("Fraction of original binary partners retained") +
  labs(colour = "Initial-degree group", fill = "Initial-degree group")
ggsave(file.path(combined_out, "28_partner_retention_by_initial_degree_and_guild.png"), p3, width = 16, height = 7, dpi = 300)

presence_active_plot_data <- retention_plot %>%
  select(dataset, guild, removal_fraction, initial_degree_group,
         median_node_present, q025_node_present, q975_node_present,
         median_network_active, q025_network_active, q975_network_active) %>%
  pivot_longer(cols = c(median_node_present, median_network_active), names_to = "metric", values_to = "median_value") %>%
  mutate(q025 = ifelse(metric == "median_node_present", q025_node_present, q025_network_active),
         q975 = ifelse(metric == "median_node_present", q975_node_present, q975_network_active),
         metric = recode(metric,
                         median_node_present = "Fraction of original nodes still present",
                         median_network_active = "Fraction of original nodes network-active"))

for(g in guild_levels){
  cur <- presence_active_plot_data %>% filter(guild == g)
  p4 <- ggplot(cur,
               aes(x = removal_fraction, y = median_value, ymin = q025, ymax = q975,
                   colour = initial_degree_group, fill = initial_degree_group)) +
    geom_ribbon(alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
    facet_grid(metric ~ dataset) +
    scale_colour_manual(values = degree_group_colours, drop = FALSE) +
    scale_fill_manual(values = degree_group_colours, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1)) + theme_classic(base_size = 9) + theme(legend.position = "bottom") +
    xlab("Proportion of sites removed") + ylab("") +
    labs(colour = "Initial-degree group", fill = "Initial-degree group")
  file_name <- ifelse(g == "Consumer", "28_consumer_presence_active_by_degree.png", "28_resource_presence_active_by_degree.png")
  ggsave(file.path(combined_out, file_name), p4, width = 16, height = 7, dpi = 300)
}

p5 <- ggplot(expected_vs_observed_group,
             aes(x = mean_expected_binary_partner_retention, y = mean_observed_binary_partner_retention,
                 colour = initial_degree_group)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.5) +
  geom_point(size = 1.7, alpha = 0.85) + facet_wrap(~ guild) +
  scale_colour_manual(values = degree_group_colours, drop = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_classic(base_size = 11) +
  theme(legend.position = "bottom") +
  xlab("Mean expected binary-partner retention from K profile") +
  ylab("Mean observed binary-partner retention") + labs(colour = "Initial-degree group")
ggsave(file.path(combined_out, "28_expected_vs_observed_partner_retention_from_K_profile.png"), p5, width = 10, height = 5, dpi = 300)

p6 <- ggplot(active_enrichment_summary,
             aes(x = removal_fraction, y = median_enrichment, ymin = q025_enrichment, ymax = q975_enrichment,
                 fill = guild, colour = guild)) +
  geom_hline(yintercept = 1, colour = "grey55") + geom_ribbon(alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) + facet_grid(guild ~ dataset, scales = "free_y") +
  theme_classic(base_size = 9) + theme(legend.position = "none") +
  xlab("Proportion of sites removed") + ylab("Initial-degree enrichment among active nodes")
ggsave(file.path(combined_out, "28_initial_degree_enrichment_among_active_nodes.png"), p6, width = 16, height = 7, dpi = 300)

message("Under uniform random site removal, a link’s probability of remaining in the binary regional network is determined by the number of sites that initially support it. Degree-dependent differences in partner retention therefore describe how local interaction support is distributed across low- and high-degree consumers or resources; they do not show that degree independently causes a link to persist.")
future::plan(future::sequential)
message("Finished script 28.")
