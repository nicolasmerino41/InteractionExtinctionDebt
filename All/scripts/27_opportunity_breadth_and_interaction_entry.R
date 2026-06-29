## ------------------------------------------------------------
## Script: All/scripts/27_opportunity_breadth_and_interaction_entry.R
##
## Purpose:
## Purely empirical full-network analysis of whether consumers with
## broader pre-interaction opportunity have co-occurring pairs that are
## more likely to enter the realised regional interaction network.
##
## Uses all full-network co-occurring consumer-resource pairs, including
## pairs with full_K = 0. No models, null models, site-removal analyses,
## or refitting are used.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

## ---------------------------
## Settings
## ---------------------------

n_boot <- 1000
min_pairs_per_group_exact_n <- 10

## Parallelisation across datasets.
## This keeps consumer-level bootstrap resampling inside each dataset,
## avoiding nested parallelism while still using multiple cores.
n_workers <- max(1, parallelly::availableCores() - 1)

result_type <- "script27_opportunity_breadth_and_interaction_entry"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

## ---------------------------
## Extra packages
## ---------------------------

extra_packages <- c("purrr", "future", "future.apply", "parallelly")
for(pkg in extra_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---------------------------
## Labels and colours
## ---------------------------

cooc_group_levels <- c(
  "Lower co-occurrence degree",
  "Middle co-occurrence degree",
  "Higher co-occurrence degree"
)

occupancy_group_levels <- c(
  "Lower occupancy",
  "Middle occupancy",
  "Higher occupancy"
)

predictor_types <- c("Co-occurrence degree", "Consumer occupancy")

outcomes <- c(
  "Ever interacted",
  "Occupancy across all co-occurring pairs",
  "Occupancy conditional on ever interacting"
)

cooc_colours <- c(
  "Lower co-occurrence degree" = "#1b9e77",
  "Middle co-occurrence degree" = "#7570b3",
  "Higher co-occurrence degree" = "#d95f02"
)

occupancy_colours <- c(
  "Lower occupancy" = "#1b9e77",
  "Middle occupancy" = "#7570b3",
  "Higher occupancy" = "#d95f02"
)

## ---------------------------
## Helper functions
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2){
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = "spearman"))
}

## Tertile-style grouping based on predictor values.
## This keeps identical predictor values in the same group by cutting the
## distribution of unique predictor values, rather than ntile()-splitting rows.
assign_three_groups <- function(values, labels){
  values <- as.numeric(values)
  out <- rep(NA_character_, length(values))
  ok <- is.finite(values)
  if(!any(ok)) return(out)

  unique_vals <- sort(unique(values[ok]))

  if(length(unique_vals) == 1){
    out[ok] <- labels[2]
    return(out)
  }

  if(length(unique_vals) == 2){
    out[ok & values == unique_vals[1]] <- labels[1]
    out[ok & values == unique_vals[2]] <- labels[3]
    return(out)
  }

  q <- as.numeric(stats::quantile(unique_vals, probs = c(1/3, 2/3), type = 1, na.rm = TRUE))
  out[ok & values <= q[1]] <- labels[1]
  out[ok & values > q[1] & values <= q[2]] <- labels[2]
  out[ok & values > q[2]] <- labels[3]
  out
}

read_if_exists <- function(file){
  if(!file.exists(file)) return(NULL)
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1) return(x)
  tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
}

empty_standardized_table <- function(){
  data.frame(
    dataset = character(),
    predictor_type = character(),
    outcome = character(),
    result_type = character(),
    group = character(),
    estimate = numeric(),
    bootstrap_q025 = numeric(),
    bootstrap_q975 = numeric(),
    n_eligible_exact_n_strata = integer(),
    minimum_eligible_n = integer(),
    maximum_eligible_n = integer(),
    total_eligible_cooccurring_pairs = integer(),
    total_eligible_consumers = integer(),
    stringsAsFactors = FALSE
  )
}

empty_exact_boot_table <- function(){
  data.frame(
    dataset = character(),
    predictor_type = character(),
    full_n = integer(),
    group = character(),
    q025_proportion_ever_interacted = numeric(),
    q975_proportion_ever_interacted = numeric(),
    stringsAsFactors = FALSE
  )
}

## ---------------------------
## Data construction
## ---------------------------

build_full_pair_data <- function(dataset){

  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    ) %>%
    distinct(site, consumer, resource, pair_id)

  interactions <- site_tables$empirical_site_interactions %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    ) %>%
    distinct(site, consumer, resource, pair_id)

  ## Interactions are treated as observed co-occurrence cells if a loader has
  ## omitted them from cooc_triples. This is a defensive consistency step: an
  ## observed interaction at a site implies the pair was recorded together there.
  cooc <- bind_rows(
    cooc,
    interactions %>% select(site, consumer, resource, pair_id)
  ) %>%
    distinct(site, consumer, resource, pair_id)

  pair_n <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  pair_k <- interactions %>%
    group_by(pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  pair_data <- pair_n %>%
    left_join(pair_k, by = "pair_id") %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      ever_interacted = full_K >= 1,
      interaction_occupancy = full_K / full_n,
      n_is_one = full_n == 1
    ) %>%
    select(dataset, consumer, resource, pair_id,
           full_n, full_K, ever_interacted,
           interaction_occupancy, n_is_one)

  if(any(pair_data$full_n < 1, na.rm = TRUE)) stop("full_n < 1 in ", dataset)
  if(any(pair_data$full_K < 0 | pair_data$full_K > pair_data$full_n, na.rm = TRUE)){
    stop("Invalid full_K/full_n relationship in ", dataset)
  }

  ## Consumer occupancy from helper occupancy table. If a consumer is absent
  ## there for any reason, fall back to the number of sites in cooc_triples.
  consumer_occ_from_occupancy <- site_tables$occupancy %>%
    mutate(site = as.character(site), species = as.character(species)) %>%
    filter(trophic_level == "consumer") %>%
    group_by(consumer = species) %>%
    summarise(consumer_occupancy = n_distinct(site), .groups = "drop")

  consumer_occ_from_cooc <- cooc %>%
    group_by(consumer) %>%
    summarise(consumer_occupancy_from_cooc = n_distinct(site), .groups = "drop")

  consumer_table <- pair_data %>%
    group_by(consumer) %>%
    summarise(
      cooccurrence_degree = n_distinct(resource),
      realised_interaction_degree = n_distinct(resource[ever_interacted]),
      fraction_cooccurring_partners_ever_interacted =
        realised_interaction_degree / cooccurrence_degree,
      mean_interaction_occupancy_all_pairs = mean(interaction_occupancy, na.rm = TRUE),
      mean_interaction_occupancy_ever_interacted =
        ifelse(any(ever_interacted),
               mean(interaction_occupancy[ever_interacted], na.rm = TRUE),
               NA_real_),
      .groups = "drop"
    ) %>%
    left_join(consumer_occ_from_occupancy, by = "consumer") %>%
    left_join(consumer_occ_from_cooc, by = "consumer") %>%
    mutate(
      consumer_occupancy = ifelse(is.na(consumer_occupancy),
                                  consumer_occupancy_from_cooc,
                                  consumer_occupancy),
      cooccurrence_degree_group = assign_three_groups(cooccurrence_degree, cooc_group_levels),
      occupancy_group = assign_three_groups(consumer_occupancy, occupancy_group_levels)
    ) %>%
    select(-consumer_occupancy_from_cooc)

  pair_data <- pair_data %>%
    left_join(
      consumer_table %>%
        select(consumer, cooccurrence_degree, cooccurrence_degree_group,
               consumer_occupancy, occupancy_group,
               realised_interaction_degree),
      by = "consumer"
    )

  list(
    pair_data = pair_data,
    consumer_table = consumer_table
  )
}

make_group_definitions <- function(dataset, consumer_table){
  bind_rows(
    consumer_table %>%
      group_by(group = cooccurrence_degree_group) %>%
      summarise(
        n_consumers = n_distinct(consumer),
        minimum_predictor_value = min(cooccurrence_degree, na.rm = TRUE),
        median_predictor_value = median(cooccurrence_degree, na.rm = TRUE),
        maximum_predictor_value = max(cooccurrence_degree, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(dataset = dataset, predictor_type = "Co-occurrence degree"),
    consumer_table %>%
      group_by(group = occupancy_group) %>%
      summarise(
        n_consumers = n_distinct(consumer),
        minimum_predictor_value = min(consumer_occupancy, na.rm = TRUE),
        median_predictor_value = median(consumer_occupancy, na.rm = TRUE),
        maximum_predictor_value = max(consumer_occupancy, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(dataset = dataset, predictor_type = "Consumer occupancy")
  ) %>%
    filter(!is.na(group)) %>%
    select(dataset, predictor_type, group, n_consumers,
           minimum_predictor_value, median_predictor_value, maximum_predictor_value)
}

## ---------------------------
## Exact-n summaries
## ---------------------------

make_exact_summary_one_predictor <- function(pair_data, predictor_type){

  if(predictor_type == "Co-occurrence degree"){
    dat <- pair_data %>%
      mutate(group = cooccurrence_degree_group) %>%
      filter(!is.na(group))
    group_levels <- cooc_group_levels
  } else {
    dat <- pair_data %>%
      mutate(group = occupancy_group) %>%
      filter(!is.na(group))
    group_levels <- occupancy_group_levels
  }

  dat$group <- factor(dat$group, levels = group_levels)

  out <- dat %>%
    group_by(dataset, full_n, group) %>%
    summarise(
      n_cooccurring_pairs = n(),
      n_consumers = n_distinct(consumer),
      n_ever_interacted = sum(ever_interacted, na.rm = TRUE),
      proportion_ever_interacted = mean(ever_interacted, na.rm = TRUE),
      mean_interaction_occupancy_all_pairs = mean(interaction_occupancy, na.rm = TRUE),
      mean_interaction_occupancy_conditional =
        ifelse(any(ever_interacted),
               mean(interaction_occupancy[ever_interacted], na.rm = TRUE),
               NA_real_),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    mutate(
      predictor_type = predictor_type,
      group = as.character(group)
    )

  eligibility <- out %>%
    group_by(dataset, full_n) %>%
    summarise(
      eligible_exact_n_stratum =
        all(group_levels %in% group) &&
        all(n_cooccurring_pairs[match(group_levels, group)] >= min_pairs_per_group_exact_n),
      .groups = "drop"
    )

  out %>%
    left_join(eligibility, by = c("dataset", "full_n")) %>%
    select(dataset, predictor_type, full_n, group,
           n_cooccurring_pairs, n_consumers, n_ever_interacted,
           proportion_ever_interacted,
           mean_interaction_occupancy_all_pairs,
           mean_interaction_occupancy_conditional,
           eligible_exact_n_stratum)
}

make_exact_summary <- function(pair_data){
  bind_rows(
    make_exact_summary_one_predictor(pair_data, "Co-occurrence degree"),
    make_exact_summary_one_predictor(pair_data, "Consumer occupancy")
  )
}

## ---------------------------
## Standardized summaries and bootstrap
## ---------------------------

estimate_standardized <- function(dat, predictor_type, eligible_ns, weights){

  if(length(eligible_ns) == 0){
    return(empty_standardized_table())
  }

  if(predictor_type == "Co-occurrence degree"){
    group_var <- "cooccurrence_degree_group"
    group_levels <- cooc_group_levels
  } else {
    group_var <- "occupancy_group"
    group_levels <- occupancy_group_levels
  }

  d <- dat %>%
    filter(full_n %in% eligible_ns, !is.na(.data[[group_var]])) %>%
    mutate(group = .data[[group_var]])

  if(nrow(d) == 0){
    return(empty_standardized_table())
  }

  n_stats <- d %>%
    group_by(dataset, group, full_n) %>%
    summarise(
      prop_ever = mean(ever_interacted, na.rm = TRUE),
      occ_all = mean(interaction_occupancy, na.rm = TRUE),
      occ_cond = ifelse(any(ever_interacted),
                        mean(interaction_occupancy[ever_interacted], na.rm = TRUE),
                        NA_real_),
      .groups = "drop"
    ) %>%
    left_join(weights, by = "full_n")

  group_values <- bind_rows(lapply(group_levels, function(g){
    cur <- n_stats %>% filter(group == g)

    value_one <- function(col){
      if(nrow(cur) == 0) return(NA_real_)
      if(any(is.na(cur[[col]])) || any(!eligible_ns %in% cur$full_n)) return(NA_real_)
      sum(cur[[col]] * cur$reference_weight, na.rm = TRUE)
    }

    data.frame(
      dataset = unique(dat$dataset)[1],
      predictor_type = predictor_type,
      outcome = outcomes,
      result_type = "Group standardized value",
      group = g,
      estimate = c(value_one("prop_ever"), value_one("occ_all"), value_one("occ_cond")),
      stringsAsFactors = FALSE
    )
  }))

  get_group_est <- function(g, outcome_name){
    group_values$estimate[group_values$group == g & group_values$outcome == outcome_name][1]
  }

  contrasts <- bind_rows(lapply(outcomes, function(outcome_name){
    lower <- if(predictor_type == "Co-occurrence degree") cooc_group_levels[1] else occupancy_group_levels[1]
    middle <- if(predictor_type == "Co-occurrence degree") cooc_group_levels[2] else occupancy_group_levels[2]
    higher <- if(predictor_type == "Co-occurrence degree") cooc_group_levels[3] else occupancy_group_levels[3]

    data.frame(
      dataset = unique(dat$dataset)[1],
      predictor_type = predictor_type,
      outcome = outcome_name,
      result_type = c("Higher minus lower", "Middle minus lower"),
      group = NA_character_,
      estimate = c(
        get_group_est(higher, outcome_name) - get_group_est(lower, outcome_name),
        get_group_est(middle, outcome_name) - get_group_est(lower, outcome_name)
      ),
      stringsAsFactors = FALSE
    )
  }))

  bind_rows(group_values, contrasts)
}

bootstrap_standardized <- function(dat, predictor_type, eligible_ns, weights, n_boot){

  if(length(eligible_ns) == 0){
    return(empty_standardized_table())
  }

  if(predictor_type == "Co-occurrence degree"){
    group_var <- "cooccurrence_degree_group"
    group_levels <- cooc_group_levels
  } else {
    group_var <- "occupancy_group"
    group_levels <- occupancy_group_levels
  }

  d <- dat %>%
    filter(full_n %in% eligible_ns, !is.na(.data[[group_var]])) %>%
    mutate(group = .data[[group_var]])

  if(nrow(d) == 0){
    return(empty_standardized_table())
  }

  consumers_by_group <- d %>%
    distinct(group, consumer) %>%
    split(.$group)

  boot_list <- vector("list", n_boot)

  for(b in seq_len(n_boot)){
    sampled <- lapply(group_levels, function(g){
      cur_consumers <- consumers_by_group[[g]]$consumer
      if(is.null(cur_consumers) || length(cur_consumers) == 0){
        return(NULL)
      }
      sampled_consumers <- sample(cur_consumers, size = length(cur_consumers), replace = TRUE)
      bind_rows(lapply(seq_along(sampled_consumers), function(i){
        d %>%
          filter(group == g, consumer == sampled_consumers[i]) %>%
          mutate(bootstrap_consumer_copy = i)
      }))
    })

    sampled_d <- bind_rows(sampled)

    boot_list[[b]] <- estimate_standardized(
      dat = sampled_d,
      predictor_type = predictor_type,
      eligible_ns = eligible_ns,
      weights = weights
    ) %>%
      mutate(bootstrap_replicate = b)
  }

  boot_values <- bind_rows(boot_list)

  if(nrow(boot_values) == 0){
    return(empty_standardized_table())
  }

  boot_values %>%
    group_by(dataset, predictor_type, outcome, result_type, group) %>%
    summarise(
      bootstrap_q025 = quantile(estimate, 0.025, na.rm = TRUE, names = FALSE),
      bootstrap_q975 = quantile(estimate, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

bootstrap_exact_proportions <- function(dat, predictor_type, eligible_ns, n_boot){

  if(length(eligible_ns) == 0){
    return(empty_exact_boot_table())
  }

  if(predictor_type == "Co-occurrence degree"){
    group_var <- "cooccurrence_degree_group"
    group_levels <- cooc_group_levels
  } else {
    group_var <- "occupancy_group"
    group_levels <- occupancy_group_levels
  }

  d <- dat %>%
    filter(full_n %in% eligible_ns, !is.na(.data[[group_var]])) %>%
    mutate(group = .data[[group_var]])

  if(nrow(d) == 0){
    return(empty_exact_boot_table())
  }

  consumers_by_group <- d %>%
    distinct(group, consumer) %>%
    split(.$group)

  boot_list <- vector("list", n_boot)

  for(b in seq_len(n_boot)){
    sampled <- lapply(group_levels, function(g){
      cur_consumers <- consumers_by_group[[g]]$consumer
      if(is.null(cur_consumers) || length(cur_consumers) == 0){
        return(NULL)
      }
      sampled_consumers <- sample(cur_consumers, size = length(cur_consumers), replace = TRUE)
      bind_rows(lapply(seq_along(sampled_consumers), function(i){
        d %>%
          filter(group == g, consumer == sampled_consumers[i]) %>%
          mutate(bootstrap_consumer_copy = i)
      }))
    })

    sampled_d <- bind_rows(sampled)

    boot_list[[b]] <- sampled_d %>%
      group_by(dataset, predictor_type = predictor_type, full_n, group) %>%
      summarise(
        proportion_ever_interacted = mean(ever_interacted, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(bootstrap_replicate = b)
  }

  bind_rows(boot_list) %>%
    group_by(dataset, predictor_type, full_n, group) %>%
    summarise(
      q025_proportion_ever_interacted = quantile(proportion_ever_interacted, 0.025, na.rm = TRUE, names = FALSE),
      q975_proportion_ever_interacted = quantile(proportion_ever_interacted, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

make_standardized_tables <- function(pair_data, exact_summary){

  out_all <- list()
  exact_boot_all <- list()

  for(pred_type in predictor_types){

    eligible_ns <- exact_summary %>%
      filter(predictor_type == pred_type, eligible_exact_n_stratum) %>%
      distinct(full_n) %>%
      arrange(full_n) %>%
      pull(full_n)

    if(length(eligible_ns) == 0){
      next
    }

    weights <- pair_data %>%
      filter(full_n %in% eligible_ns) %>%
      count(full_n, name = "n_pairs") %>%
      mutate(reference_weight = n_pairs / sum(n_pairs)) %>%
      select(full_n, reference_weight)

    estimates <- estimate_standardized(pair_data, pred_type, eligible_ns, weights)
    boot_intervals <- bootstrap_standardized(pair_data, pred_type, eligible_ns, weights, n_boot)

    if(nrow(boot_intervals) > 0){
      estimates <- estimates %>%
        left_join(
          boot_intervals,
          by = c("dataset", "predictor_type", "outcome", "result_type", "group")
        )
    } else {
      estimates$bootstrap_q025 <- NA_real_
      estimates$bootstrap_q975 <- NA_real_
    }

    estimates <- estimates %>%
      mutate(
        n_eligible_exact_n_strata = length(eligible_ns),
        minimum_eligible_n = min(eligible_ns),
        maximum_eligible_n = max(eligible_ns),
        total_eligible_cooccurring_pairs = sum(weights$n_pairs),
        total_eligible_consumers = pair_data %>%
          filter(full_n %in% eligible_ns) %>%
          summarise(n = n_distinct(consumer)) %>%
          pull(n)
      )

    out_all[[pred_type]] <- estimates
    exact_boot_all[[pred_type]] <- bootstrap_exact_proportions(pair_data, pred_type, eligible_ns, n_boot)
  }

  list(
    standardized = bind_rows(out_all),
    exact_boot = bind_rows(exact_boot_all)
  )
}

## ---------------------------
## Consumer-level correlations
## ---------------------------

make_consumer_correlations <- function(dataset, consumer_table){

  corr_vars <- c(
    "cooccurrence_degree",
    "consumer_occupancy",
    "realised_interaction_degree",
    "fraction_cooccurring_partners_ever_interacted",
    "mean_interaction_occupancy_all_pairs",
    "mean_interaction_occupancy_ever_interacted"
  )

  pairs <- combn(corr_vars, 2, simplify = FALSE)

  bind_rows(lapply(pairs, function(v){
    data.frame(
      dataset = dataset,
      variable_1 = v[1],
      variable_2 = v[2],
      spearman_correlation = safe_spearman(consumer_table[[v[1]]], consumer_table[[v[2]]]),
      n_consumers = sum(is.finite(consumer_table[[v[1]]]) & is.finite(consumer_table[[v[2]]])),
      stringsAsFactors = FALSE
    )
  }))
}

## ---------------------------
## Run one dataset
## ---------------------------

run_one_dataset <- function(dataset){

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tibble)
    library(purrr)
  })

  message("Running script 27 opportunity breadth and interaction entry: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  built <- build_full_pair_data(dataset)
  pair_data <- built$pair_data
  consumer_table <- built$consumer_table %>% mutate(dataset = dataset, .before = 1)

  group_defs <- make_group_definitions(dataset, consumer_table)

  exact_summary <- make_exact_summary(pair_data)

  ## Add exact-n eligibility flags to pair-level data.
  eligible_cooc_ns <- exact_summary %>%
    filter(predictor_type == "Co-occurrence degree", eligible_exact_n_stratum) %>%
    distinct(full_n) %>%
    pull(full_n)

  eligible_occ_ns <- exact_summary %>%
    filter(predictor_type == "Consumer occupancy", eligible_exact_n_stratum) %>%
    distinct(full_n) %>%
    pull(full_n)

  pair_data <- pair_data %>%
    mutate(
      eligible_exact_n_cooccurrence_degree = full_n %in% eligible_cooc_ns,
      eligible_exact_n_occupancy = full_n %in% eligible_occ_ns
    )

  std <- make_standardized_tables(pair_data, exact_summary)
  standardized <- std$standardized
  exact_boot <- std$exact_boot

  exact_summary_plot <- exact_summary %>%
    left_join(
      exact_boot,
      by = c("dataset", "predictor_type", "full_n", "group")
    )

  consumer_correlations <- make_consumer_correlations(dataset, consumer_table)

  n1_structure <- bind_rows(
    pair_data %>%
      filter(full_n == 1) %>%
      count(dataset, predictor_type = "Co-occurrence degree",
            group = cooccurrence_degree_group, full_K, name = "n_pairs") %>%
      group_by(dataset, predictor_type, group) %>%
      mutate(fraction = n_pairs / sum(n_pairs)) %>%
      ungroup(),
    pair_data %>%
      filter(full_n == 1) %>%
      count(dataset, predictor_type = "Consumer occupancy",
            group = occupancy_group, full_K, name = "n_pairs") %>%
      group_by(dataset, predictor_type, group) %>%
      mutate(fraction = n_pairs / sum(n_pairs)) %>%
      ungroup()
  ) %>%
    mutate(outcome = ifelse(full_K == 0, "full_K = 0", "full_K = 1"))

  checks <- data.frame(
    dataset = dataset,
    includes_full_K_zero = any(pair_data$full_K == 0),
    all_full_n_ge_1 = all(pair_data$full_n >= 1),
    all_K_between_0_and_n = all(pair_data$full_K >= 0 & pair_data$full_K <= pair_data$full_n),
    n1_has_K0 = any(pair_data$full_n == 1 & pair_data$full_K == 0),
    n1_has_K1 = any(pair_data$full_n == 1 & pair_data$full_K == 1),
    conditional_occupancy_excludes_K0 = TRUE,
    exact_n_min10_all_groups_cooccurrence_degree =
      ifelse(length(eligible_cooc_ns) == 0, TRUE,
             exact_summary %>%
               filter(predictor_type == "Co-occurrence degree",
                      eligible_exact_n_stratum) %>%
               summarise(ok = all(n_cooccurring_pairs >= min_pairs_per_group_exact_n)) %>%
               pull(ok)),
    exact_n_min10_all_groups_occupancy =
      ifelse(length(eligible_occ_ns) == 0, TRUE,
             exact_summary %>%
               filter(predictor_type == "Consumer occupancy",
                      eligible_exact_n_stratum) %>%
               summarise(ok = all(n_cooccurring_pairs >= min_pairs_per_group_exact_n)) %>%
               pull(ok)),
    weights_sum_to_1_cooccurrence_degree =
      ifelse(length(eligible_cooc_ns) == 0, TRUE,
             abs(pair_data %>%
                   filter(full_n %in% eligible_cooc_ns) %>%
                   count(full_n) %>%
                   mutate(w = n / sum(n)) %>%
                   summarise(s = sum(w)) %>%
                   pull(s) - 1) < 1e-8),
    weights_sum_to_1_occupancy =
      ifelse(length(eligible_occ_ns) == 0, TRUE,
             abs(pair_data %>%
                   filter(full_n %in% eligible_occ_ns) %>%
                   count(full_n) %>%
                   mutate(w = n / sum(n)) %>%
                   summarise(s = sum(w)) %>%
                   pull(s) - 1) < 1e-8),
    bootstrap_is_consumer_level = TRUE,
    grouping_avoids_realised_interaction_degree = TRUE,
    no_models_or_site_removal_used = TRUE
  )

  write.csv2(pair_data,
             file.path(out_dir, paste0(dataset, "_27_all_cooccurrence_pairs_with_opportunity_predictors.csv")),
             row.names = FALSE)
  write.csv2(consumer_table,
             file.path(out_dir, paste0(dataset, "_27_consumer_opportunity_breadth.csv")),
             row.names = FALSE)
  write.csv2(exact_summary,
             file.path(out_dir, paste0(dataset, "_27_ever_interacted_by_exact_n_and_group.csv")),
             row.names = FALSE)
  write.csv2(standardized,
             file.path(out_dir, paste0(dataset, "_27_standardized_opportunity_breadth_comparisons.csv")),
             row.names = FALSE)

  list(
    pair_data = pair_data,
    consumer_table = consumer_table,
    exact_summary = exact_summary,
    exact_summary_plot = exact_summary_plot,
    standardized = standardized,
    consumer_correlations = consumer_correlations,
    group_defs = group_defs,
    n1_structure = n1_structure,
    checks = checks
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

message("Using ", n_workers, " parallel workers across datasets.")
future::plan(future::multisession, workers = n_workers)

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)
names(all_outputs) <- all_dataset_names

future::plan(future::sequential)

pair_data_all <- bind_rows(lapply(all_outputs, `[[`, "pair_data"))
consumer_table_all <- bind_rows(lapply(all_outputs, `[[`, "consumer_table"))
exact_summary_all <- bind_rows(lapply(all_outputs, `[[`, "exact_summary"))
exact_summary_plot_all <- bind_rows(lapply(all_outputs, `[[`, "exact_summary_plot"))
standardized_all <- bind_rows(lapply(all_outputs, `[[`, "standardized"))
consumer_correlations_all <- bind_rows(lapply(all_outputs, `[[`, "consumer_correlations"))
group_defs_all <- bind_rows(lapply(all_outputs, `[[`, "group_defs"))
n1_structure_all <- bind_rows(lapply(all_outputs, `[[`, "n1_structure"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

pair_data_all$dataset <- factor(pair_data_all$dataset, levels = all_dataset_names)
exact_summary_all$dataset <- factor(exact_summary_all$dataset, levels = all_dataset_names)
exact_summary_plot_all$dataset <- factor(exact_summary_plot_all$dataset, levels = all_dataset_names)
standardized_all$dataset <- factor(standardized_all$dataset, levels = all_dataset_names)
consumer_table_all$dataset <- factor(consumer_table_all$dataset, levels = all_dataset_names)
group_defs_all$dataset <- factor(group_defs_all$dataset, levels = all_dataset_names)
n1_structure_all$dataset <- factor(n1_structure_all$dataset, levels = all_dataset_names)

## ---------------------------
## Save combined tables
## ---------------------------

pair_data_out <- pair_data_all %>%
  select(dataset, consumer, resource, full_n, full_K,
         ever_interacted, interaction_occupancy,
         cooccurrence_degree,
         `cooccurrence-degree group` = cooccurrence_degree_group,
         consumer_occupancy,
         `occupancy group` = occupancy_group,
         `full realised interaction degree` = realised_interaction_degree,
         eligible_exact_n_cooccurrence_degree,
         eligible_exact_n_occupancy)

consumer_out <- consumer_table_all %>%
  select(dataset, consumer,
         cooccurrence_degree,
         `cooccurrence-degree group` = cooccurrence_degree_group,
         consumer_occupancy,
         `occupancy group` = occupancy_group,
         `realised interaction degree` = realised_interaction_degree,
         `fraction of co-occurring partners ever observed interacting` = fraction_cooccurring_partners_ever_interacted,
         `mean interaction occupancy across all co-occurring pairs` = mean_interaction_occupancy_all_pairs,
         `mean interaction occupancy among ever-interacted pairs` = mean_interaction_occupancy_ever_interacted)

write.csv2(pair_data_out,
           file.path(combined_out, "27_all_cooccurrence_pairs_with_opportunity_predictors.csv"),
           row.names = FALSE)

write.csv2(consumer_out,
           file.path(combined_out, "27_consumer_opportunity_breadth.csv"),
           row.names = FALSE)

write.csv2(exact_summary_all,
           file.path(combined_out, "27_ever_interacted_by_exact_n_and_group.csv"),
           row.names = FALSE)

write.csv2(standardized_all,
           file.path(combined_out, "27_standardized_opportunity_breadth_comparisons.csv"),
           row.names = FALSE)

write.csv2(consumer_correlations_all,
           file.path(combined_out, "27_consumer_summary_correlations.csv"),
           row.names = FALSE)

write.csv2(group_defs_all,
           file.path(combined_out, "27_opportunity_breadth_group_definitions.csv"),
           row.names = FALSE)

write.csv2(checks_all,
           file.path(combined_out, "27_opportunity_breadth_interaction_entry_checks.csv"),
           row.names = FALSE)

## Save n = 1 table in requested wide format.
n1_out <- n1_structure_all %>%
  select(dataset, predictor_type, group, outcome, n_pairs, fraction) %>%
  pivot_wider(
    names_from = outcome,
    values_from = c(n_pairs, fraction),
    values_fill = 0
  )

write.csv2(n1_out,
           file.path(combined_out, "27_n1_outcomes_by_degree.csv"),
           row.names = FALSE)

## ---------------------------
## Figures
## ---------------------------

## Figure 1: exact-n ever-interacted by co-occurrence degree.
plot_cooc <- exact_summary_plot_all %>%
  filter(predictor_type == "Co-occurrence degree", eligible_exact_n_stratum) %>%
  mutate(group = factor(group, levels = cooc_group_levels))

if(nrow(plot_cooc) > 0){
  p1 <- ggplot(plot_cooc,
               aes(x = full_n,
                   y = proportion_ever_interacted,
                   colour = group,
                   fill = group)) +
    geom_ribbon(aes(ymin = q025_proportion_ever_interacted,
                    ymax = q975_proportion_ever_interacted),
                alpha = 0.15, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 1.7, na.rm = TRUE) +
    scale_colour_manual(values = cooc_colours, drop = FALSE) +
    scale_fill_manual(values = cooc_colours, drop = FALSE) +
    facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 10) +
    xlab("Number of sites where pair was recorded together") +
    ylab("Fraction of co-occurring pairs ever observed interacting") +
    labs(
      colour = "Co-occurrence degree group",
      fill = "Co-occurrence degree group",
      subtitle = "All co-occurring pairs included; groups defined from potential partner breadth only."
    )

  ggsave(file.path(combined_out, "27_ever_interacted_by_cooccurrence_degree.png"),
         p1, width = 14, height = 7, dpi = 300)
}

## Figure 2: exact-n ever-interacted by consumer occupancy.
plot_occ <- exact_summary_plot_all %>%
  filter(predictor_type == "Consumer occupancy", eligible_exact_n_stratum) %>%
  mutate(group = factor(group, levels = occupancy_group_levels))

if(nrow(plot_occ) > 0){
  p2 <- ggplot(plot_occ,
               aes(x = full_n,
                   y = proportion_ever_interacted,
                   colour = group,
                   fill = group)) +
    geom_ribbon(aes(ymin = q025_proportion_ever_interacted,
                    ymax = q975_proportion_ever_interacted),
                alpha = 0.15, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 1.7, na.rm = TRUE) +
    scale_colour_manual(values = occupancy_colours, drop = FALSE) +
    scale_fill_manual(values = occupancy_colours, drop = FALSE) +
    facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 10) +
    xlab("Number of sites where pair was recorded together") +
    ylab("Fraction of co-occurring pairs ever observed interacting") +
    labs(
      colour = "Occupancy group",
      fill = "Occupancy group",
      subtitle = "All co-occurring pairs included; groups defined from consumer occurrence across sampled sites only."
    )

  ggsave(file.path(combined_out, "27_ever_interacted_by_consumer_occupancy.png"),
         p2, width = 14, height = 7, dpi = 300)
}

## Figure 3: standardized interaction entry by opportunity breadth.
std_entry_groups <- standardized_all %>%
  filter(outcome == "Ever interacted",
         result_type == "Group standardized value") %>%
  mutate(panel = "Standardized probability of ever interacting")

std_entry_contrast <- standardized_all %>%
  filter(outcome == "Ever interacted",
         result_type == "Higher minus lower") %>%
  mutate(group = "Higher minus lower",
         panel = "Higher minus lower")

std_entry_plot <- bind_rows(std_entry_groups, std_entry_contrast) %>%
  mutate(
    panel = factor(panel, levels = c("Standardized probability of ever interacting", "Higher minus lower")),
    dataset = factor(dataset, levels = all_dataset_names)
  )

if(nrow(std_entry_plot) > 0){
  p3 <- ggplot(std_entry_plot,
               aes(x = dataset,
                   y = estimate,
                   ymin = bootstrap_q025,
                   ymax = bootstrap_q975,
                   colour = group)) +
    geom_hline(data = data.frame(panel = "Higher minus lower"),
               aes(yintercept = 0), colour = "grey50", inherit.aes = FALSE) +
    geom_pointrange(position = position_dodge(width = 0.55), na.rm = TRUE) +
    facet_grid(panel ~ predictor_type, scales = "free_y") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("Estimate") +
    labs(colour = "Group or contrast")

  ggsave(file.path(combined_out, "27_standardized_interaction_entry_by_opportunity_breadth.png"),
         p3, width = 14, height = 8, dpi = 300)
}

## Figure 4: standardized occupancy after entry.
std_occ_plot <- standardized_all %>%
  filter(outcome %in% c("Occupancy across all co-occurring pairs",
                        "Occupancy conditional on ever interacting"),
         result_type == "Group standardized value") %>%
  mutate(
    outcome = recode(outcome,
                     "Occupancy across all co-occurring pairs" = "Across all co-occurring pairs",
                     "Occupancy conditional on ever interacting" = "Among pairs ever observed interacting"),
    dataset = factor(dataset, levels = all_dataset_names)
  )

if(nrow(std_occ_plot) > 0){
  p4 <- ggplot(std_occ_plot,
               aes(x = dataset,
                   y = estimate,
                   ymin = bootstrap_q025,
                   ymax = bootstrap_q975,
                   colour = group)) +
    geom_pointrange(position = position_dodge(width = 0.55), na.rm = TRUE) +
    facet_grid(outcome ~ predictor_type, scales = "free_y") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("Standardized mean interaction occupancy") +
    labs(colour = "Group")

  ggsave(file.path(combined_out, "27_standardized_occupancy_after_entry.png"),
         p4, width = 14, height = 8, dpi = 300)
}

## Figure 5: n = 1 outcomes.
n1_plot <- n1_structure_all %>%
  filter(predictor_type %in% predictor_types) %>%
  mutate(
    group = factor(group, levels = c(cooc_group_levels, occupancy_group_levels)),
    outcome = factor(outcome, levels = c("full_K = 0", "full_K = 1"))
  )

if(nrow(n1_plot) > 0){
  p5 <- ggplot(n1_plot,
               aes(x = group, y = fraction, fill = outcome)) +
    geom_col(position = "fill") +
    facet_grid(predictor_type ~ dataset, scales = "free_x", space = "free_x") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("Fraction of pairs with full_n = 1") +
    labs(fill = "Outcome")

  ggsave(file.path(combined_out, "27_n1_outcomes_by_initial_degree.png"),
         p5, width = 16, height = 7, dpi = 300)
}

## ---------------------------
## Console summary and reminder
## ---------------------------

message("\nValidation summary for script 27:")
print(checks_all)

message("\nDegree is defined from the observed regional interaction network. Differences among degree groups therefore describe the structure of co-occurrence opportunities associated with low- and high-degree consumers; they do not establish that degree causes an interaction to occur.")
message("\nGroups are defined from co-occurrence opportunity or consumer occupancy, not from realised interactions. Results describe whether consumers with broader opportunity are associated with a higher probability that co-occurring pairs are ever observed interacting; they do not establish a causal effect of breadth or occupancy.")

message("Finished script 27 opportunity breadth and interaction entry.")
