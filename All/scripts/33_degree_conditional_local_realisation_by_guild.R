## ------------------------------------------------------------
## Script: All/scripts/33_degree_conditional_local_realisation_by_guild.R
##
## Purpose:
## Degree-conditional local realisation by guild.
##
## Question:
## At the same exact number of sites where a consumer-resource pair is
## recorded together, do high-degree consumers/resources realise a smaller
## fraction of those co-occurrence sites as observed interactions?
##
## Notes:
## - Uses all co-occurring pairs (full_n >= 1), including full_K = 0.
## - Analyses consumers and resources separately.
## - No models, null models, regressions, GAMs, site removal, or p-values.
## - Bootstrap unit is the focal node.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr", "tidyr", "ggplot2", "tibble", "purrr",
  "scales", "patchwork", "future", "future.apply"
)

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---------------------------
## User controls
## ---------------------------

out_dir <- "All/outputs/33_degree_conditional_local_realisation_by_guild"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

n_boot <- 1000
min_pairs_per_group <- 10
min_nodes_per_group <- 3

## Use parallel across datasets. Set to FALSE for debugging.
use_parallel <- TRUE
n_workers <- max(1, min(7, future::availableCores() - 1))

if(use_parallel){
  future::plan(future::multisession, workers = n_workers)
} else {
  future::plan(future::sequential)
}

## ---------------------------
## Plot style
## ---------------------------

degree_levels <- c(
  "Lower initial degree",
  "Middle initial degree",
  "Higher initial degree"
)

degree_cols <- c(
  "Lower initial degree"  = "#0072B2",
  "Middle initial degree" = "#E69F00",
  "Higher initial degree" = "#D55E00"
)

dataset_order <- all_dataset_names

save_both <- function(p, name, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height)
}

theme_33 <- function(base_size = 10){
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

## ---------------------------
## Helpers
## ---------------------------

assign_tied_tertiles <- function(df, value_col, group_col = "degree_group"){
  stopifnot(value_col %in% names(df))

  tmp <- df %>%
    distinct(.data[[value_col]]) %>%
    arrange(.data[[value_col]]) %>%
    mutate(
      value_rank = row_number(),
      n_values = n(),
      tercile_raw = ceiling(3 * value_rank / n_values),
      tercile_raw = pmin(pmax(tercile_raw, 1), 3),
      "{group_col}" := factor(
        case_when(
          tercile_raw == 1 ~ "Lower initial degree",
          tercile_raw == 2 ~ "Middle initial degree",
          tercile_raw == 3 ~ "Higher initial degree"
        ),
        levels = degree_levels
      )
    ) %>%
    select(all_of(value_col), all_of(group_col))

  df %>% left_join(tmp, by = value_col)
}

make_pair_table <- function(dataset){
  message("Loading all co-occurring pairs: ", dataset)

  st <- get_dataset_site_tables(dataset)

  cooc <- st$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )

  ints <- st$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )

  ## Defensive consistency: an observed interaction is also a co-occurrence.
  ## This only adds missing operational co-occurrence rows if a source table
  ## failed to record them explicitly.
  cooc <- bind_rows(cooc, ints) %>%
    distinct(site, consumer, resource)

  cooc_counts <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  int_counts <- ints %>%
    group_by(consumer, resource) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  pair_table <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      dataset = dataset,
      full_K = tidyr::replace_na(full_K, 0L),
      ever_interacted = full_K >= 1,
      local_realisation = full_K / full_n
    ) %>%
    select(dataset, consumer, resource, full_n, full_K,
           ever_interacted, local_realisation)

  if(any(pair_table$full_n < 1)){
    stop("full_n < 1 detected in ", dataset)
  }
  if(any(pair_table$full_K < 0 | pair_table$full_K > pair_table$full_n)){
    bad <- pair_table %>% filter(full_K < 0 | full_K > full_n)
    print(head(bad, 20))
    stop("Invalid full_K values in ", dataset)
  }

  pair_table
}

make_guild_long_table <- function(pair_table){
  consumer_degree <- pair_table %>%
    filter(ever_interacted) %>%
    group_by(consumer) %>%
    summarise(initial_degree = n_distinct(resource), .groups = "drop") %>%
    mutate(guild = "Consumer", focal_node = consumer) %>%
    select(guild, focal_node, initial_degree)

  resource_degree <- pair_table %>%
    filter(ever_interacted) %>%
    group_by(resource) %>%
    summarise(initial_degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(guild = "Resource", focal_node = resource) %>%
    select(guild, focal_node, initial_degree)

  degree_table <- bind_rows(consumer_degree, resource_degree) %>%
    group_by(guild) %>%
    group_modify(~ assign_tied_tertiles(.x, "initial_degree", "degree_group")) %>%
    ungroup()

  consumer_long <- pair_table %>%
    transmute(
      dataset, guild = "Consumer",
      focal_node = consumer,
      partner_node = resource,
      full_n, full_K, ever_interacted, local_realisation
    )

  resource_long <- pair_table %>%
    transmute(
      dataset, guild = "Resource",
      focal_node = resource,
      partner_node = consumer,
      full_n, full_K, ever_interacted, local_realisation
    )

  bind_rows(consumer_long, resource_long) %>%
    inner_join(degree_table, by = c("guild", "focal_node")) %>%
    mutate(
      degree_group = factor(degree_group, levels = degree_levels),
      outcome_A_eligible = TRUE,
      outcome_B_eligible = TRUE,
      outcome_C_eligible = ever_interacted & full_n >= 2
    ) %>%
    select(dataset, guild, focal_node, partner_node,
           initial_degree, degree_group,
           full_n, full_K, ever_interacted, local_realisation,
           outcome_A_eligible, outcome_B_eligible, outcome_C_eligible)
}

summarise_exact_n <- function(long_pairs){
  outcomes <- list(
    "Fraction of co-occurrence sites with interaction" = long_pairs %>%
      filter(outcome_A_eligible) %>%
      group_by(dataset, guild, full_n, degree_group) %>%
      summarise(
        number_pair_records = n(),
        number_focal_nodes = n_distinct(focal_node),
        estimate = mean(local_realisation, na.rm = TRUE),
        mean_full_K = mean(full_K, na.rm = TRUE),
        median_full_K = median(full_K, na.rm = TRUE),
        .groups = "drop"
      ),
    "Pairs ever observed interacting" = long_pairs %>%
      filter(outcome_B_eligible) %>%
      group_by(dataset, guild, full_n, degree_group) %>%
      summarise(
        number_pair_records = n(),
        number_focal_nodes = n_distinct(focal_node),
        estimate = mean(ever_interacted, na.rm = TRUE),
        mean_full_K = mean(full_K, na.rm = TRUE),
        median_full_K = median(full_K, na.rm = TRUE),
        .groups = "drop"
      ),
    "Interaction repetition after a link exists" = long_pairs %>%
      filter(outcome_C_eligible) %>%
      group_by(dataset, guild, full_n, degree_group) %>%
      summarise(
        number_pair_records = n(),
        number_focal_nodes = n_distinct(focal_node),
        estimate = mean(local_realisation, na.rm = TRUE),
        mean_full_K = mean(full_K, na.rm = TRUE),
        median_full_K = median(full_K, na.rm = TRUE),
        .groups = "drop"
      )
  )

  bind_rows(lapply(names(outcomes), function(nm){
    outcomes[[nm]] %>% mutate(outcome = nm)
  })) %>%
    group_by(dataset, guild, outcome, full_n) %>%
    mutate(
      n_groups_present = n_distinct(degree_group),
      eligible_exact_n_stratum =
        first(n_groups_present) == 3 &&
        all(number_pair_records >= min_pairs_per_group) &&
        all(number_focal_nodes >= min_nodes_per_group)
    ) %>%
    ungroup() %>%
    select(dataset, guild, full_n, degree_group, outcome,
           number_pair_records, number_focal_nodes, estimate,
           mean_full_K, median_full_K, eligible_exact_n_stratum)
}

reference_weights <- function(long_pairs, exact_summary){
  eligible_n <- exact_summary %>%
    filter(eligible_exact_n_stratum) %>%
    distinct(dataset, guild, outcome, full_n)

  out_A <- long_pairs %>%
    filter(outcome_A_eligible) %>%
    mutate(outcome = "Fraction of co-occurrence sites with interaction")

  out_B <- long_pairs %>%
    filter(outcome_B_eligible) %>%
    mutate(outcome = "Pairs ever observed interacting")

  out_C <- long_pairs %>%
    filter(outcome_C_eligible) %>%
    mutate(outcome = "Interaction repetition after a link exists")

  bind_rows(out_A, out_B, out_C) %>%
    semi_join(eligible_n, by = c("dataset", "guild", "outcome", "full_n")) %>%
    count(dataset, guild, outcome, full_n, name = "n_pairs_ref") %>%
    group_by(dataset, guild, outcome) %>%
    mutate(weight = n_pairs_ref / sum(n_pairs_ref)) %>%
    ungroup()
}

compute_standardized_once <- function(long_pairs, weights, weighting_type = "Pair-weighted primary"){
  if(nrow(weights) == 0){
    return(tibble(
      dataset = character(),
      guild = character(),
      outcome = character(),
      weighting_type = character(),
      result_type = character(),
      degree_group = character(),
      estimate = numeric(),
      number_eligible_exact_n_strata = integer(),
      minimum_eligible_n = integer(),
      maximum_eligible_n = integer(),
      total_eligible_pair_records = integer(),
      total_eligible_focal_nodes = integer()
    ))
  }

  make_outcome_data <- function(outcome_name){
    if(outcome_name == "Fraction of co-occurrence sites with interaction"){
      long_pairs %>%
        filter(outcome_A_eligible) %>%
        mutate(outcome = outcome_name, outcome_value = local_realisation)
    } else if(outcome_name == "Pairs ever observed interacting"){
      long_pairs %>%
        filter(outcome_B_eligible) %>%
        mutate(outcome = outcome_name, outcome_value = as.numeric(ever_interacted))
    } else {
      long_pairs %>%
        filter(outcome_C_eligible) %>%
        mutate(outcome = outcome_name, outcome_value = local_realisation)
    }
  }

  outcome_data <- bind_rows(lapply(unique(weights$outcome), make_outcome_data)) %>%
    semi_join(weights, by = c("dataset", "guild", "outcome", "full_n"))

  if(nrow(outcome_data) == 0){
    return(tibble())
  }

  if(weighting_type == "Pair-weighted primary"){
    exact_est <- outcome_data %>%
      group_by(dataset, guild, outcome, full_n, degree_group) %>%
      summarise(exact_estimate = mean(outcome_value, na.rm = TRUE), .groups = "drop")
  } else {
    exact_est <- outcome_data %>%
      group_by(dataset, guild, outcome, full_n, degree_group, focal_node) %>%
      summarise(node_estimate = mean(outcome_value, na.rm = TRUE), .groups = "drop") %>%
      group_by(dataset, guild, outcome, full_n, degree_group) %>%
      summarise(exact_estimate = mean(node_estimate, na.rm = TRUE), .groups = "drop")
  }

  group_est <- exact_est %>%
    left_join(weights, by = c("dataset", "guild", "outcome", "full_n")) %>%
    group_by(dataset, guild, outcome, degree_group) %>%
    summarise(
      estimate = sum(weight * exact_estimate, na.rm = TRUE),
      number_eligible_exact_n_strata = n_distinct(full_n),
      minimum_eligible_n = min(full_n, na.rm = TRUE),
      maximum_eligible_n = max(full_n, na.rm = TRUE),
      .groups = "drop"
    )

  metadata <- outcome_data %>%
    group_by(dataset, guild, outcome) %>%
    summarise(
      total_eligible_pair_records = n(),
      total_eligible_focal_nodes = n_distinct(focal_node),
      .groups = "drop"
    )

  group_rows <- group_est %>%
    left_join(metadata, by = c("dataset", "guild", "outcome")) %>%
    mutate(
      weighting_type = weighting_type,
      result_type = "Degree group estimate"
    )

  contrasts <- group_rows %>%
    select(dataset, guild, outcome, weighting_type,
           degree_group, estimate,
           number_eligible_exact_n_strata,
           minimum_eligible_n, maximum_eligible_n,
           total_eligible_pair_records, total_eligible_focal_nodes) %>%
    tidyr::pivot_wider(names_from = degree_group, values_from = estimate) %>%
    mutate(
      `Higher minus lower` = `Higher initial degree` - `Lower initial degree`,
      `Middle minus lower` = `Middle initial degree` - `Lower initial degree`
    ) %>%
    select(dataset, guild, outcome, weighting_type,
           number_eligible_exact_n_strata,
           minimum_eligible_n, maximum_eligible_n,
           total_eligible_pair_records, total_eligible_focal_nodes,
           `Higher minus lower`, `Middle minus lower`) %>%
    pivot_longer(
      cols = c(`Higher minus lower`, `Middle minus lower`),
      names_to = "result_type",
      values_to = "estimate"
    ) %>%
    mutate(degree_group = NA_character_)

  bind_rows(
    group_rows %>%
      select(dataset, guild, outcome, weighting_type, result_type,
             degree_group, estimate,
             number_eligible_exact_n_strata,
             minimum_eligible_n, maximum_eligible_n,
             total_eligible_pair_records, total_eligible_focal_nodes),
    contrasts
  )
}

bootstrap_standardized <- function(long_pairs, weights, n_boot, weighting_type){
  if(nrow(weights) == 0){
    return(tibble())
  }

  groups <- long_pairs %>%
    semi_join(weights, by = c("dataset", "guild", "full_n")) %>%
    distinct(dataset, guild, degree_group, focal_node)

  boot_one <- function(b){
    sampled_nodes <- groups %>%
      group_by(dataset, guild, degree_group) %>%
      summarise(
        focal_node = list(sample(focal_node, size = n(), replace = TRUE)),
        .groups = "drop"
      ) %>%
      tidyr::unnest(focal_node) %>%
      group_by(dataset, guild, degree_group, focal_node) %>%
      mutate(.copy_id = row_number()) %>%
      ungroup()

    boot_pairs <- sampled_nodes %>%
      left_join(long_pairs,
                by = c("dataset", "guild", "degree_group", "focal_node"),
                relationship = "many-to-many")

    compute_standardized_once(
      boot_pairs,
      weights,
      weighting_type = weighting_type
    ) %>%
      mutate(bootstrap_replicate = b)
  }

  bind_rows(lapply(seq_len(n_boot), boot_one))
}

make_exact_boot_intervals <- function(long_pairs, exact_summary, guild_filter){
  focal_data <- long_pairs %>%
    filter(guild == guild_filter)

  eligible <- exact_summary %>%
    filter(
      guild == guild_filter,
      outcome == "Fraction of co-occurrence sites with interaction",
      eligible_exact_n_stratum
    ) %>%
    distinct(dataset, guild, full_n)

  if(nrow(eligible) == 0){
    return(tibble(
      dataset = factor(character(), levels = dataset_order),
      guild = character(),
      full_n = integer(),
      degree_group = factor(character(), levels = degree_levels),
      q025 = numeric(),
      q975 = numeric()
    ))
  }

  focal_data <- focal_data %>%
    filter(outcome_A_eligible) %>%
    semi_join(eligible, by = c("dataset", "guild", "full_n"))

  node_groups <- focal_data %>%
    distinct(dataset, guild, degree_group, focal_node)

  boot_one <- function(b){
    sampled <- node_groups %>%
      group_by(dataset, guild, degree_group) %>%
      summarise(
        focal_node = list(sample(focal_node, size = n(), replace = TRUE)),
        .groups = "drop"
      ) %>%
      tidyr::unnest(focal_node) %>%
      group_by(dataset, guild, degree_group, focal_node) %>%
      mutate(.copy_id = row_number()) %>%
      ungroup()

    sampled %>%
      left_join(focal_data,
                by = c("dataset", "guild", "degree_group", "focal_node"),
                relationship = "many-to-many") %>%
      group_by(dataset, guild, full_n, degree_group) %>%
      summarise(estimate = mean(local_realisation, na.rm = TRUE), .groups = "drop") %>%
      mutate(bootstrap_replicate = b)
  }

  bind_rows(lapply(seq_len(n_boot), boot_one)) %>%
    group_by(dataset, guild, full_n, degree_group) %>%
    summarise(
      q025 = quantile(estimate, 0.025, na.rm = TRUE),
      q975 = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

run_one_dataset <- function(dataset){
  message("Running script 33: ", dataset)

  pair_table <- make_pair_table(dataset)
  long_pairs <- make_guild_long_table(pair_table)

  degree_groups <- long_pairs %>%
    distinct(dataset, guild, focal_node, initial_degree, degree_group) %>%
    group_by(dataset, guild, degree_group) %>%
    summarise(
      number_focal_nodes = n_distinct(focal_node),
      minimum_initial_degree = min(initial_degree),
      median_initial_degree = median(initial_degree),
      maximum_initial_degree = max(initial_degree),
      .groups = "drop"
    )

  exact_summary <- summarise_exact_n(long_pairs)
  weights <- reference_weights(long_pairs, exact_summary)

  std_pair <- compute_standardized_once(
    long_pairs,
    weights,
    weighting_type = "Pair-weighted primary"
  )

  std_node <- compute_standardized_once(
    long_pairs,
    weights,
    weighting_type = "Focal-node-weighted sensitivity"
  )

  boot_pair <- bootstrap_standardized(
    long_pairs,
    weights,
    n_boot = n_boot,
    weighting_type = "Pair-weighted primary"
  )

  boot_node <- bootstrap_standardized(
    long_pairs,
    weights,
    n_boot = n_boot,
    weighting_type = "Focal-node-weighted sensitivity"
  )

  boot_all <- bind_rows(boot_pair, boot_node)

  boot_ci <- boot_all %>%
    group_by(dataset, guild, outcome, weighting_type, result_type, degree_group) %>%
    summarise(
      bootstrap_q025 = quantile(estimate, 0.025, na.rm = TRUE),
      bootstrap_q975 = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  std_all <- bind_rows(std_pair, std_node) %>%
    left_join(
      boot_ci,
      by = c("dataset", "guild", "outcome", "weighting_type",
             "result_type", "degree_group")
    )

  checks <- tibble(
    dataset = dataset,
    includes_full_K_zero_pairs = any(long_pairs$full_K == 0),
    all_full_n_ge_1 = all(long_pairs$full_n >= 1),
    all_K_between_0_and_n = all(long_pairs$full_K >= 0 & long_pairs$full_K <= long_pairs$full_n),
    consumers_and_resources_present = all(c("Consumer", "Resource") %in% unique(long_pairs$guild)),
    outcome_A_includes_n1 = any(long_pairs$outcome_A_eligible & long_pairs$full_n == 1),
    outcome_B_includes_n1 = any(long_pairs$outcome_B_eligible & long_pairs$full_n == 1),
    outcome_C_excludes_n1 = !any(long_pairs$outcome_C_eligible & long_pairs$full_n == 1),
    outcome_C_excludes_K0 = !any(long_pairs$outcome_C_eligible & long_pairs$full_K == 0),
    eligible_exact_n_min_pairs_ok =
      all(exact_summary$number_pair_records[exact_summary$eligible_exact_n_stratum] >= min_pairs_per_group),
    eligible_exact_n_min_nodes_ok =
      all(exact_summary$number_focal_nodes[exact_summary$eligible_exact_n_stratum] >= min_nodes_per_group),
    weights_sum_to_1 = {
      ws <- weights %>%
        group_by(dataset, guild, outcome) %>%
        summarise(s = sum(weight), .groups = "drop")
      if(nrow(ws) == 0) TRUE else all(abs(ws$s - 1) < 1e-8)
    },
    bootstrap_at_focal_node_level = TRUE,
    no_model_null_site_removal_or_regression_used = TRUE
  )

  list(
    pair_data = long_pairs,
    exact_summary = exact_summary,
    standardized = std_all,
    degree_groups = degree_groups,
    checks = checks
  )
}

## ---------------------------
## Run
## ---------------------------

if(use_parallel){
  outputs <- future.apply::future_lapply(
    all_dataset_names,
    run_one_dataset,
    future.seed = TRUE
  )
} else {
  outputs <- lapply(all_dataset_names, run_one_dataset)
}

pair_data <- bind_rows(lapply(outputs, `[[`, "pair_data")) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    degree_group = factor(degree_group, levels = degree_levels)
  )

exact_summary <- bind_rows(lapply(outputs, `[[`, "exact_summary")) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    degree_group = factor(degree_group, levels = degree_levels)
  )

standardized <- bind_rows(lapply(outputs, `[[`, "standardized")) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    degree_group = factor(degree_group, levels = degree_levels)
  )

degree_groups <- bind_rows(lapply(outputs, `[[`, "degree_groups")) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_order),
    degree_group = factor(degree_group, levels = degree_levels)
  )

checks <- bind_rows(lapply(outputs, `[[`, "checks"))

## ---------------------------
## Write tables
## ---------------------------

write.csv2(
  pair_data,
  file.path(out_dir, "33_all_cooccurrence_pairs_by_degree_group.csv"),
  row.names = FALSE
)

write.csv2(
  exact_summary,
  file.path(out_dir, "33_exact_n_degree_conditional_outcomes.csv"),
  row.names = FALSE
)

write.csv2(
  standardized,
  file.path(out_dir, "33_standardized_degree_conditional_outcomes.csv"),
  row.names = FALSE
)

write.csv2(
  degree_groups,
  file.path(out_dir, "33_initial_degree_groups_by_guild.csv"),
  row.names = FALSE
)

write.csv2(
  checks,
  file.path(out_dir, "33_degree_conditional_realisation_checks.csv"),
  row.names = FALSE
)

## ---------------------------
## Figure 1: exact-n local realisation, consumers/resources separately
## ---------------------------

plot_exact_guild <- function(guild_name, out_name){
  plot_dat <- exact_summary %>%
    filter(
      guild == guild_name,
      outcome == "Fraction of co-occurrence sites with interaction",
      eligible_exact_n_stratum
    )

  if(nrow(plot_dat) == 0){
    warning("No eligible exact-n strata for ", guild_name, ". Skipping figure.")
    return(invisible(NULL))
  }

  boot_int <- bind_rows(lapply(all_dataset_names, function(ds){
    message("Exact-n bootstrap intervals for ", guild_name, ": ", ds)
    make_exact_boot_intervals(
      pair_data %>% filter(dataset == ds),
      exact_summary %>% filter(dataset == ds),
      guild_name
    )
  })) %>%
    mutate(dataset = factor(dataset, levels = dataset_order),
           degree_group = factor(degree_group, levels = degree_levels))

  plot_dat <- plot_dat %>%
    left_join(
      boot_int,
      by = c("dataset", "guild", "full_n", "degree_group")
    )

  p <- ggplot(
    plot_dat,
    aes(x = full_n, y = estimate, colour = degree_group, fill = degree_group)
  ) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      alpha = 0.12,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
    scale_colour_manual(values = degree_cols, drop = FALSE) +
    scale_fill_manual(values = degree_cols, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_33(base_size = 9) +
    xlab("Sites where pair is recorded together") +
    ylab("Fraction of co-occurrence sites with observed interaction") +
    labs(
      colour = "Initial-degree group",
      fill = "Initial-degree group",
      title = paste0(guild_name, ": local realisation at the same co-occurrence support"),
      subtitle = "Pairs are compared only when they are recorded together in the same number of sites."
    )

  save_both(p, out_name, width = 13, height = 8)
  invisible(p)
}

plot_exact_guild(
  "Consumer",
  "33_consumers_local_realisation_given_same_cooccurrence_support"
)

plot_exact_guild(
  "Resource",
  "33_resources_local_realisation_given_same_cooccurrence_support"
)

## ---------------------------
## Figure 2: higher-minus-lower summary
## ---------------------------

summary_plot_data <- standardized %>%
  filter(
    weighting_type == "Pair-weighted primary",
    result_type == "Higher minus lower"
  ) %>%
  mutate(
    guild = factor(guild, levels = c("Consumer", "Resource"),
                   labels = c("Consumers", "Resources")),
    outcome = factor(
      outcome,
      levels = c(
        "Fraction of co-occurrence sites with interaction",
        "Pairs ever observed interacting",
        "Interaction repetition after a link exists"
      ),
      labels = c(
        "Fraction of co-occurrence sites with interaction",
        "Pairs ever observed interacting",
        "Interaction repetition after a link exists"
      )
    )
  )

p_summary <- ggplot(
  summary_plot_data,
  aes(x = dataset, y = estimate)
) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
    width = 0,
    colour = "grey45",
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_point(size = 2, colour = "#222222", na.rm = TRUE) +
  facet_grid(outcome ~ guild, scales = "free_y") +
  theme_33(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab("Higher initial degree minus lower initial degree") +
  ggtitle(
    "Degree-conditional local realisation",
    subtitle = "Negative values mean that high-degree nodes realise less of the same co-occurrence opportunity."
  )

save_both(
  p_summary,
  "33_degree_conditional_realisation_summary",
  width = 13,
  height = 9
)

## ---------------------------
## Figure 3: decomposition standardized values
## ---------------------------

decomp_data <- standardized %>%
  filter(
    weighting_type == "Pair-weighted primary",
    result_type == "Degree group estimate"
  ) %>%
  mutate(
    guild = factor(guild, levels = c("Consumer", "Resource"),
                   labels = c("Consumers", "Resources")),
    outcome = factor(
      outcome,
      levels = c(
        "Fraction of co-occurrence sites with interaction",
        "Pairs ever observed interacting",
        "Interaction repetition after a link exists"
      ),
      labels = c(
        "Local realisation across all co-occurring pairs",
        "Entry into regional interaction network",
        "Repetition after a link exists"
      )
    )
  )

p_decomp <- ggplot(
  decomp_data,
  aes(x = degree_group, y = estimate, colour = degree_group)
) +
  geom_errorbar(
    aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
    width = 0.12,
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_point(size = 1.8, na.rm = TRUE) +
  facet_grid(outcome ~ guild) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_33(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  xlab(NULL) +
  ylab("Standardized value") +
  labs(
    colour = "Initial-degree group",
    title = "Decomposition of degree-conditional realisation",
    subtitle = "Values are matched to a common distribution of co-occurrence support within each dataset and guild."
  )

save_both(
  p_decomp,
  "33_degree_realisation_decomposition",
  width = 12,
  height = 9
)

## ---------------------------
## Console summary
## ---------------------------

message("\nValidation summary:")
print(checks)

message("\nEligible exact-n strata by dataset/guild/outcome:")
eligible_summary <- exact_summary %>%
  filter(eligible_exact_n_stratum) %>%
  distinct(dataset, guild, outcome, full_n) %>%
  count(dataset, guild, outcome, name = "n_eligible_exact_n")
print(eligible_summary, n = Inf)

message("\nInterpretation boundary:")
message(
  "Initial degree is calculated from observed regional interactions, so this analysis describes how local interaction realisation is associated with low- and high-degree nodes. Matching co-occurrence support removes the mechanical effect that pairs recorded together in more sites have more scope for a low realised fraction. The results do not establish that degree itself causes interactions to be more or less locally realised."
)

message("\nSaved script 33 outputs in: ", out_dir)
