## ------------------------------------------------------------
## Script: All/scripts/29_revisiting_degree_structure_under_site_removal.R
##
## Purpose:
## Empirical re-analysis of binary regional interaction degree under
## random site removal, separately for consumers and resources.
##
## This script:
##   - uses only original realised interaction links (full_K >= 1);
##   - does not fit models or degree distributions;
##   - compares cumulative retained-degree distributions among all
##     original nodes and among active nodes only;
##   - tracks node-level degree transitions and the opposing processes
##     of active-node enrichment and degree compression.
##
## Run from the parent repository folder.
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

## ---------------------------
## Settings
## ---------------------------

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500
min_sites_retained <- 1

result_type <- "script29_revisiting_degree_structure_under_site_removal"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

selected_removal_levels_for_cdf <- c(0, 0.4, 0.8)

## Parallelise across datasets by default.
n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")

## ---------------------------
## Small helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

make_node_id <- function(guild, node){
  paste(guild, node, sep = "___")
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

## Tied tertile-style grouping.
## This keeps identical degree values together whenever possible.
assign_tied_tertiles <- function(df,
                                 value_col = "initial_degree",
                                 group_col = "initial_degree_group"){
  values <- sort(unique(df[[value_col]]))
  if(length(values) == 0){
    df[[group_col]] <- character()
    return(df)
  }
  if(length(values) == 1){
    map <- data.frame(
      value = values,
      group = "Middle initial degree"
    )
  } else {
    counts_by_value <- df %>%
      count(.data[[value_col]], name = "n_nodes") %>%
      arrange(.data[[value_col]]) %>%
      mutate(
        cum_before = lag(cumsum(n_nodes), default = 0),
        cum_mid = cum_before + n_nodes / 2,
        frac_mid = cum_mid / sum(n_nodes),
        group = case_when(
          frac_mid <= 1/3 ~ "Lower initial degree",
          frac_mid <= 2/3 ~ "Middle initial degree",
          TRUE ~ "Higher initial degree"
        )
      )
    names(counts_by_value)[names(counts_by_value) == value_col] <- "value"
    map <- counts_by_value %>% select(value, group)
  }
  out <- df %>%
    left_join(map, by = setNames("value", value_col))
  names(out)[names(out) == "group"] <- group_col
  out
}

make_site_subsets_local <- function(all_sites){
  all_sites <- sort(unique(as.character(all_sites)))
  n_sites <- length(all_sites)
  subset_list <- list()
  subset_index <- list()
  counter <- 1L

  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_sites * (1 - removal)))
    n_keep <- min(n_keep, n_sites)
    reps_here <- ifelse(removal == 0, 1, n_site_reps)

    for(r in seq_len(reps_here)){
      sites_keep <- if(removal == 0) all_sites else sample(all_sites, size = n_keep, replace = FALSE)
      subset_list[[counter]] <- sort(as.character(sites_keep))
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        replicate = r,
        n_sites_total = n_sites,
        n_sites_kept = length(sites_keep),
        n_sites_removed = n_sites - length(sites_keep),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }

  list(index = bind_rows(subset_index), subsets = subset_list)
}

make_original_links <- function(dataset, cooc_triples, interactions){
  cooc_counts <- cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource)) %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  int_counts <- interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource)) %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  original_links <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource", "pair_id")) %>%
    mutate(
      full_K = replace_na(full_K, 0L),
      dataset = dataset
    ) %>%
    filter(full_K >= 1) %>%
    select(dataset, consumer, resource, pair_id, full_n, full_K)

  if(any(original_links$full_K > original_links$full_n)){
    stop("Found original links with full_K > full_n in ", dataset)
  }

  original_links
}

make_node_tables <- function(dataset, original_links){
  consumer_nodes <- original_links %>%
    group_by(node = consumer) %>%
    summarise(initial_degree = n_distinct(resource), .groups = "drop") %>%
    mutate(dataset = dataset, guild = "Consumer")

  resource_nodes <- original_links %>%
    group_by(node = resource) %>%
    summarise(initial_degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(dataset = dataset, guild = "Resource")

  nodes <- bind_rows(consumer_nodes, resource_nodes) %>%
    select(dataset, guild, node, initial_degree) %>%
    group_by(dataset, guild) %>%
    group_modify(~ assign_tied_tertiles(.x, value_col = "initial_degree",
                                        group_col = "initial_degree_group")) %>%
    ungroup()

  degree_groups <- nodes %>%
    group_by(dataset, guild, initial_degree_group) %>%
    summarise(
      number_of_nodes = n(),
      minimum_initial_degree = min(initial_degree, na.rm = TRUE),
      median_initial_degree = median(initial_degree, na.rm = TRUE),
      maximum_initial_degree = max(initial_degree, na.rm = TRUE),
      .groups = "drop"
    )

  list(nodes = nodes, degree_groups = degree_groups)
}

## Convert a retained set of original links into retained degree for both guilds.
retained_degrees_for_subset <- function(original_links,
                                        retained_interactions,
                                        nodes,
                                        removal_fraction,
                                        replicate){

  retained_pairs <- retained_interactions %>%
    distinct(consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource))

  link_state <- original_links %>%
    mutate(retained_link = pair_id %in% retained_pairs$pair_id)

  consumer_retained <- link_state %>%
    group_by(node = consumer) %>%
    summarise(retained_degree = n_distinct(resource[retained_link]), .groups = "drop") %>%
    mutate(guild = "Consumer")

  resource_retained <- link_state %>%
    group_by(node = resource) %>%
    summarise(retained_degree = n_distinct(consumer[retained_link]), .groups = "drop") %>%
    mutate(guild = "Resource")

  bind_rows(consumer_retained, resource_retained) %>%
    right_join(nodes, by = c("guild", "node")) %>%
    mutate(
      retained_degree = replace_na(retained_degree, 0L),
      removal_fraction = removal_fraction,
      replicate = replicate,
      network_active = retained_degree > 0,
      retained_degree_fraction = retained_degree / initial_degree
    ) %>%
    select(dataset, guild, node, initial_degree, initial_degree_group,
           removal_fraction, replicate, retained_degree, network_active,
           retained_degree_fraction)
}

cdf_from_degrees <- function(deg_vec, denominator_type){
  deg_vec <- as.integer(deg_vec)

  if(denominator_type == "All original nodes"){
    denom <- length(deg_vec)
    max_deg <- max(deg_vec, na.rm = TRUE)
  } else {
    deg_vec <- deg_vec[deg_vec > 0]
    denom <- length(deg_vec)
    max_deg <- ifelse(denom == 0, 0L, max(deg_vec, na.rm = TRUE))
  }

  if(is.na(max_deg) || max_deg < 1 || denom == 0){
    return(data.frame(
      degree_threshold = integer(),
      cumulative_probability = numeric()
    ))
  }

  x <- seq_len(max_deg)
  data.frame(
    degree_threshold = x,
    cumulative_probability = vapply(x, function(xx){
      mean(deg_vec >= xx)
    }, numeric(1))
  )
}

make_cdf_rows <- function(node_transitions_subset){
  node_counts <- node_transitions_subset %>%
    group_by(dataset, guild, removal_fraction, replicate) %>%
    summarise(
      number_of_original_nodes = n(),
      number_of_active_nodes = sum(network_active, na.rm = TRUE),
      .groups = "drop"
    )

  all_cdf <- node_transitions_subset %>%
    group_by(dataset, guild, removal_fraction, replicate) %>%
    group_modify(~ cdf_from_degrees(.x$retained_degree, "All original nodes")) %>%
    ungroup() %>%
    mutate(distribution_type = "All original nodes")

  active_cdf <- node_transitions_subset %>%
    group_by(dataset, guild, removal_fraction, replicate) %>%
    group_modify(~ cdf_from_degrees(.x$retained_degree, "Active nodes only")) %>%
    ungroup() %>%
    mutate(distribution_type = "Active nodes only")

  bind_rows(all_cdf, active_cdf) %>%
    left_join(node_counts, by = c("dataset", "guild", "removal_fraction", "replicate")) %>%
    select(dataset, guild, removal_fraction, replicate, distribution_type,
           degree_threshold, cumulative_probability,
           number_of_original_nodes, number_of_active_nodes)
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  message("Running script 29 degree structure under site removal: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource)) %>%
    distinct(site, consumer, resource)

  interactions <- site_tables$empirical_site_interactions %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource)) %>%
    distinct(site, consumer, resource)

  ## Be defensive: an observed interaction necessarily implies an operational
  ## co-occurrence record at that site.
  cooc_triples <- bind_rows(
    cooc_triples,
    interactions %>% select(site, consumer, resource)
  ) %>%
    distinct(site, consumer, resource)

  original_links <- make_original_links(dataset, cooc_triples, interactions)
  node_obj <- make_node_tables(dataset, original_links)
  nodes <- node_obj$nodes
  degree_groups <- node_obj$degree_groups

  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets_local(all_sites)

  node_transition_list <- vector("list", nrow(subset_object$index))
  composition_list <- vector("list", nrow(subset_object$index))

  mean_initial_all <- nodes %>%
    group_by(dataset, guild) %>%
    summarise(mean_initial_degree_all_nodes = mean(initial_degree), .groups = "drop")

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    removal <- subset_object$index$removal_fraction[i]
    rep_id <- subset_object$index$replicate[i]

    retained_interactions <- interactions %>%
      filter(site %in% sites_keep) %>%
      distinct(site, consumer, resource)

    trans_i <- retained_degrees_for_subset(
      original_links = original_links,
      retained_interactions = retained_interactions,
      nodes = nodes,
      removal_fraction = removal,
      replicate = rep_id
    )

    node_transition_list[[i]] <- trans_i

    composition_list[[i]] <- trans_i %>%
      group_by(dataset, guild, removal_fraction, replicate) %>%
      summarise(
        fraction_active = mean(network_active, na.rm = TRUE),
        mean_retained_degree_among_active_nodes =
          ifelse(any(network_active), mean(retained_degree[network_active], na.rm = TRUE), NA_real_),
        median_retained_degree_among_active_nodes =
          ifelse(any(network_active), median(retained_degree[network_active], na.rm = TRUE), NA_real_),
        mean_initial_degree_among_active_nodes =
          ifelse(any(network_active), mean(initial_degree[network_active], na.rm = TRUE), NA_real_),
        median_initial_degree_among_active_nodes =
          ifelse(any(network_active), median(initial_degree[network_active], na.rm = TRUE), NA_real_),
        .groups = "drop"
      ) %>%
      left_join(mean_initial_all, by = c("dataset", "guild")) %>%
      mutate(
        initial_degree_enrichment_among_active_nodes =
          mean_initial_degree_among_active_nodes / mean_initial_degree_all_nodes,
        degree_compression_among_active_nodes =
          mean_retained_degree_among_active_nodes / mean_initial_degree_among_active_nodes
      )
  }

  node_transitions <- bind_rows(node_transition_list)
  composition_compression <- bind_rows(composition_list)

  cumulative_degree_distributions <- make_cdf_rows(node_transitions)

  cumulative_degree_summary <- cumulative_degree_distributions %>%
    group_by(dataset, guild, removal_fraction, distribution_type, degree_threshold) %>%
    summarise(
      median_cumulative_probability = median(cumulative_probability, na.rm = TRUE),
      q025 = quantile(cumulative_probability, 0.025, na.rm = TRUE),
      q975 = quantile(cumulative_probability, 0.975, na.rm = TRUE),
      number_of_original_nodes = first(number_of_original_nodes),
      median_number_of_active_nodes = median(number_of_active_nodes, na.rm = TRUE),
      .groups = "drop"
    )

  transition_rep_group <- node_transitions %>%
    group_by(dataset, guild, removal_fraction, replicate, initial_degree_group) %>%
    summarise(
      number_of_original_nodes = n(),
      number_active = sum(network_active, na.rm = TRUE),
      fraction_retained_degree_zero = mean(retained_degree == 0, na.rm = TRUE),
      median_retained_degree_fraction = median(retained_degree_fraction, na.rm = TRUE),
      median_retained_degree_among_active_nodes =
        ifelse(any(network_active), median(retained_degree[network_active], na.rm = TRUE), NA_real_),
      median_initial_degree_among_active_nodes =
        ifelse(any(network_active), median(initial_degree[network_active], na.rm = TRUE), NA_real_),
      .groups = "drop"
    )

  transition_group_summary <- transition_rep_group %>%
    group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
    summarise(
      number_of_original_nodes = first(number_of_original_nodes),
      median_number_active = median(number_active, na.rm = TRUE),
      median_fraction_retained_degree_zero = median(fraction_retained_degree_zero, na.rm = TRUE),
      q025_fraction_retained_degree_zero = quantile(fraction_retained_degree_zero, 0.025, na.rm = TRUE),
      q975_fraction_retained_degree_zero = quantile(fraction_retained_degree_zero, 0.975, na.rm = TRUE),
      median_retained_degree_fraction = median(median_retained_degree_fraction, na.rm = TRUE),
      q025_retained_degree_fraction = quantile(median_retained_degree_fraction, 0.025, na.rm = TRUE),
      q975_retained_degree_fraction = quantile(median_retained_degree_fraction, 0.975, na.rm = TRUE),
      median_retained_degree_among_active_nodes =
        median(median_retained_degree_among_active_nodes, na.rm = TRUE),
      median_initial_degree_among_active_nodes =
        median(median_initial_degree_among_active_nodes, na.rm = TRUE),
      .groups = "drop"
    )

  ## Validation checks
  check_1_sep_guilds <- all(c("Consumer", "Resource") %in% unique(nodes$guild))
  check_2_initial_degree <- all(nodes$initial_degree >= 1)
  check_3_retained_not_exceed_initial <- all(node_transitions$retained_degree <= node_transitions$initial_degree)

  cdf_check1 <- cumulative_degree_distributions %>%
    filter(distribution_type == "All original nodes", degree_threshold == 1) %>%
    left_join(
      node_transitions %>%
        group_by(dataset, guild, removal_fraction, replicate) %>%
        summarise(fraction_active = mean(network_active), .groups = "drop"),
      by = c("dataset", "guild", "removal_fraction", "replicate")
    ) %>%
    mutate(ok = abs(cumulative_probability - fraction_active) < 1e-12)

  check_4_all_node_cdf_x1_equals_fraction_active <- all(cdf_check1$ok)

  cdf_check2 <- cumulative_degree_distributions %>%
    filter(distribution_type == "Active nodes only", degree_threshold == 1) %>%
    mutate(ok = abs(cumulative_probability - 1) < 1e-12)

  check_5_active_cdf_x1_equals_1 <- all(cdf_check2$ok)

  cdf_monotonic <- cumulative_degree_distributions %>%
    arrange(dataset, guild, removal_fraction, replicate, distribution_type, degree_threshold) %>%
    group_by(dataset, guild, removal_fraction, replicate, distribution_type) %>%
    summarise(
      ok = all(diff(cumulative_probability) <= 1e-12),
      .groups = "drop"
    )
  check_6_cdf_monotonic <- all(cdf_monotonic$ok)

  zero_comp <- composition_compression %>%
    filter(removal_fraction == 0)
  check_7_enrichment_zero <- all(abs(zero_comp$initial_degree_enrichment_among_active_nodes - 1) < 1e-12)
  check_8_compression_zero <- all(abs(zero_comp$degree_compression_among_active_nodes - 1) < 1e-12)

  check_rows <- data.frame(
    dataset = dataset,
    consumers_and_resources_analysed_separately = check_1_sep_guilds,
    all_original_nodes_initial_degree_at_least_1 = check_2_initial_degree,
    retained_degree_never_exceeds_initial_degree = check_3_retained_not_exceed_initial,
    all_node_cdf_threshold_1_equals_fraction_active =
      check_4_all_node_cdf_x1_equals_fraction_active,
    active_only_cdf_threshold_1_equals_1_when_active_nodes_exist =
      check_5_active_cdf_x1_equals_1,
    cumulative_probabilities_nonincreasing_with_degree_threshold =
      check_6_cdf_monotonic,
    initial_degree_enrichment_equals_1_at_zero_removal = check_7_enrichment_zero,
    degree_compression_equals_1_at_zero_removal = check_8_compression_zero,
    zero_degree_nodes_only_in_all_node_denominator = TRUE,
    no_parametric_degree_distribution_fitting_or_model_comparison = TRUE,
    stringsAsFactors = FALSE
  )

  if(!all(unlist(check_rows[,-1]))){
    print(check_rows)
    stop("Validation failed for dataset: ", dataset)
  }

  write.csv2(cumulative_degree_distributions,
             file.path(out_dir, paste0(dataset, "_29_cumulative_degree_distributions.csv")),
             row.names = FALSE)
  write.csv2(node_transitions,
             file.path(out_dir, paste0(dataset, "_29_node_degree_transitions.csv")),
             row.names = FALSE)

  list(
    nodes = nodes,
    degree_groups = degree_groups,
    cumulative_degree_distributions = cumulative_degree_distributions,
    cumulative_degree_summary = cumulative_degree_summary,
    node_degree_transitions = node_transitions,
    transition_rep_group = transition_rep_group,
    transition_group_summary = transition_group_summary,
    composition_compression = composition_compression,
    checks = check_rows
  )
}

## ---------------------------
## Run
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)

names(all_outputs) <- all_dataset_names
future::plan(future::sequential)

nodes_all <- bind_rows(lapply(all_outputs, `[[`, "nodes"))
degree_groups_all <- bind_rows(lapply(all_outputs, `[[`, "degree_groups"))
cumulative_degree_distributions_all <- bind_rows(lapply(all_outputs, `[[`, "cumulative_degree_distributions"))
cumulative_degree_summary_all <- bind_rows(lapply(all_outputs, `[[`, "cumulative_degree_summary"))
node_degree_transitions_all <- bind_rows(lapply(all_outputs, `[[`, "node_degree_transitions"))
transition_rep_group_all <- bind_rows(lapply(all_outputs, `[[`, "transition_rep_group"))
transition_group_summary_all <- bind_rows(lapply(all_outputs, `[[`, "transition_group_summary"))
composition_compression_all <- bind_rows(lapply(all_outputs, `[[`, "composition_compression"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

dataset_levels <- all_dataset_names
guild_levels <- c("Consumer", "Resource")
degree_group_levels <- c("Lower initial degree", "Middle initial degree", "Higher initial degree")

## ---------------------------
## Save required tables
## ---------------------------

write.csv2(
  cumulative_degree_distributions_all,
  file.path(combined_out, "29_cumulative_degree_distributions.csv"),
  row.names = FALSE
)

write.csv2(
  cumulative_degree_summary_all,
  file.path(combined_out, "29_cumulative_degree_distribution_summary.csv"),
  row.names = FALSE
)

write.csv2(
  node_degree_transitions_all,
  file.path(combined_out, "29_node_degree_transitions.csv"),
  row.names = FALSE
)

write.csv2(
  transition_group_summary_all,
  file.path(combined_out, "29_degree_transition_group_summary.csv"),
  row.names = FALSE
)

write.csv2(
  composition_compression_all,
  file.path(combined_out, "29_degree_composition_compression_summary.csv"),
  row.names = FALSE
)

write.csv2(
  checks_all,
  file.path(combined_out, "29_degree_structure_checks.csv"),
  row.names = FALSE
)

## Extra useful table with group definitions.
write.csv2(
  degree_groups_all,
  file.path(combined_out, "29_initial_degree_groups_by_guild.csv"),
  row.names = FALSE
)

message("")
message("Validation summary for script 29:")
print(checks_all)

## ---------------------------
## Plot helpers and themes
## ---------------------------

removal_colours <- c(
  "0" = "#1b9e77",
  "0.4" = "#d95f02",
  "0.8" = "#7570b3"
)

degree_group_colours <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#d95f02",
  "Higher initial degree" = "#7570b3"
)

base_theme_29 <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom"
  )

## ---------------------------
## Figure 1: CDF all nodes
## ---------------------------

cdf_all_plot <- cumulative_degree_summary_all %>%
  filter(
    distribution_type == "All original nodes",
    removal_fraction %in% selected_removal_levels_for_cdf,
    median_cumulative_probability > 0
  ) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    guild = factor(guild, levels = guild_levels),
    removal_fraction_label = factor(as.character(removal_fraction),
                                    levels = as.character(selected_removal_levels_for_cdf))
  )

p1 <- ggplot(
  cdf_all_plot,
  aes(x = degree_threshold,
      y = median_cumulative_probability,
      colour = removal_fraction_label,
      fill = removal_fraction_label,
      group = removal_fraction_label)
) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 0.9, na.rm = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = removal_colours, name = "Sites removed") +
  scale_fill_manual(values = removal_colours, name = "Sites removed") +
  facet_grid(guild ~ dataset, scales = "free") +
  base_theme_29 +
  xlab("Retained binary degree") +
  ylab("P(retained degree >= x) among all original nodes")

ggsave(
  file.path(combined_out, "29_cumulative_degree_distributions_all_nodes.png"),
  p1,
  width = 16,
  height = 6.5,
  dpi = 300
)

## ---------------------------
## Figure 2: CDF active nodes
## ---------------------------

cdf_active_plot <- cumulative_degree_summary_all %>%
  filter(
    distribution_type == "Active nodes only",
    removal_fraction %in% selected_removal_levels_for_cdf,
    median_cumulative_probability > 0
  ) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    guild = factor(guild, levels = guild_levels),
    removal_fraction_label = factor(as.character(removal_fraction),
                                    levels = as.character(selected_removal_levels_for_cdf))
  )

p2 <- ggplot(
  cdf_active_plot,
  aes(x = degree_threshold,
      y = median_cumulative_probability,
      colour = removal_fraction_label,
      fill = removal_fraction_label,
      group = removal_fraction_label)
) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 0.9, na.rm = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = removal_colours, name = "Sites removed") +
  scale_fill_manual(values = removal_colours, name = "Sites removed") +
  facet_grid(guild ~ dataset, scales = "free") +
  base_theme_29 +
  xlab("Retained binary degree") +
  ylab("P(retained degree >= x | retained degree > 0)")

ggsave(
  file.path(combined_out, "29_cumulative_degree_distributions_active_nodes.png"),
  p2,
  width = 16,
  height = 6.5,
  dpi = 300
)

## ---------------------------
## Figure 3: degree transitions by initial degree
## ---------------------------

transition_plot <- transition_group_summary_all %>%
  select(dataset, guild, removal_fraction, initial_degree_group,
         median_fraction_retained_degree_zero,
         q025_fraction_retained_degree_zero,
         q975_fraction_retained_degree_zero,
         median_retained_degree_fraction,
         q025_retained_degree_fraction,
         q975_retained_degree_fraction) %>%
  pivot_longer(
    cols = c(median_fraction_retained_degree_zero,
             median_retained_degree_fraction),
    names_to = "metric",
    values_to = "median_value"
  ) %>%
  mutate(
    q025 = ifelse(metric == "median_fraction_retained_degree_zero",
                  q025_fraction_retained_degree_zero,
                  q025_retained_degree_fraction),
    q975 = ifelse(metric == "median_fraction_retained_degree_zero",
                  q975_fraction_retained_degree_zero,
                  q975_retained_degree_fraction),
    metric = recode(
      metric,
      median_fraction_retained_degree_zero = "Fraction with retained degree = 0",
      median_retained_degree_fraction = "Median retained-degree fraction"
    ),
    dataset = factor(dataset, levels = dataset_levels),
    guild = factor(guild, levels = guild_levels),
    initial_degree_group = factor(initial_degree_group, levels = degree_group_levels)
  )

p3 <- ggplot(
  transition_plot,
  aes(x = removal_fraction,
      y = median_value,
      colour = initial_degree_group,
      fill = initial_degree_group)
) +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.2, na.rm = TRUE) +
  scale_colour_manual(values = degree_group_colours, name = "Initial degree group") +
  scale_fill_manual(values = degree_group_colours, name = "Initial degree group") +
  facet_grid(guild + metric ~ dataset) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme_29 +
  xlab("Proportion of sites removed") +
  ylab("Fraction")

ggsave(
  file.path(combined_out, "29_degree_transition_by_initial_degree.png"),
  p3,
  width = 16,
  height = 9,
  dpi = 300
)

## ---------------------------
## Figure 4: composition and compression
## ---------------------------

comp_summary <- composition_compression_all %>%
  select(dataset, guild, removal_fraction, replicate,
         initial_degree_enrichment_among_active_nodes,
         degree_compression_among_active_nodes) %>%
  pivot_longer(
    cols = c(initial_degree_enrichment_among_active_nodes,
             degree_compression_among_active_nodes),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(dataset, guild, removal_fraction, metric) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    metric = recode(
      metric,
      initial_degree_enrichment_among_active_nodes = "Initial-degree enrichment among active nodes",
      degree_compression_among_active_nodes = "Degree compression among active nodes"
    ),
    dataset = factor(dataset, levels = dataset_levels),
    guild = factor(guild, levels = guild_levels)
  )

p4 <- ggplot(
  comp_summary,
  aes(x = removal_fraction,
      y = median_value,
      fill = guild,
      colour = guild)
) +
  geom_hline(yintercept = 1, colour = "grey60") +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.2, na.rm = TRUE) +
  facet_grid(metric + guild ~ dataset, scales = "free_y") +
  base_theme_29 +
  xlab("Proportion of sites removed") +
  ylab("Value") +
  guides(fill = "none", colour = "none")

ggsave(
  file.path(combined_out, "29_degree_composition_and_compression.png"),
  p4,
  width = 16,
  height = 8.5,
  dpi = 300
)

## ---------------------------
## Interpretation reminder
## ---------------------------

message("")
message("Changes in the cumulative degree distribution must be interpreted together with node loss and degree compression. A broadly similar distribution among active nodes does not imply that the original network structure is unchanged, because low-degree nodes may disappear while initially high-degree nodes remain active but lose partners.")

message("")
message("Finished script 29. Outputs saved in: ", combined_out)
