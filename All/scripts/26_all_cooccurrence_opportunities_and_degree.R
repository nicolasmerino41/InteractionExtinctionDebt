## ------------------------------------------------------------
## Script: All/scripts/26_all_cooccurrence_opportunities_and_degree.R
##
## Purpose:
## Empirical full-network analysis using ALL co-occurring consumer-resource
## pairs (including full_K = 0) to compare lower-, middle-, and higher-degree
## consumers at two stages:
##   1. whether a co-occurring pair is ever observed interacting;
##   2. interaction occupancy K/n, both across all co-occurring pairs and
##      conditional on ever interacting.
##
## No site-removal analysis. No homogeneous model. No beta-binomial model.
## No null model, GAM, power law, or other probabilistic model.
##
## Run from the parent repository folder.
## Outputs:
##   All/SeparatedResults/all_cooccurrence_opportunities_and_degree/<dataset>/
##   All/CombinedOutputs/
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("dplyr", "tidyr", "ggplot2", "tibble", "future", "future.apply", "parallelly")
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

n_boot <- 1000
min_pairs_per_group_per_n <- 10

result_type <- "all_cooccurrence_opportunities_and_degree"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " workers across datasets.")

## ---------------------------
## Helper functions
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

degree_group_levels <- c(
  "Lower initial degree",
  "Middle initial degree",
  "Higher initial degree"
)

assign_degree_groups <- function(consumer_degree){
  ## Degree grouping based on unique degree values, so consumers with the same
  ## full_degree are kept together whenever possible.
  deg_levels <- consumer_degree %>%
    distinct(full_degree) %>%
    arrange(full_degree) %>%
    mutate(rank_degree = row_number())

  n_unique <- nrow(deg_levels)

  if(n_unique == 1){
    deg_levels <- deg_levels %>%
      mutate(`consumer-degree group` = "Middle initial degree")
  } else if(n_unique == 2){
    deg_levels <- deg_levels %>%
      mutate(`consumer-degree group` = ifelse(rank_degree == 1,
                                              "Lower initial degree",
                                              "Higher initial degree"))
  } else {
    deg_levels <- deg_levels %>%
      mutate(
        tertile_rank = dplyr::ntile(rank_degree, 3),
        `consumer-degree group` = dplyr::case_when(
          tertile_rank == 1 ~ "Lower initial degree",
          tertile_rank == 2 ~ "Middle initial degree",
          tertile_rank == 3 ~ "Higher initial degree",
          TRUE ~ NA_character_
        )
      )
  }

  consumer_degree %>%
    left_join(deg_levels %>% select(full_degree, `consumer-degree group`),
              by = "full_degree") %>%
    mutate(`consumer-degree group` = factor(`consumer-degree group`, levels = degree_group_levels))
}

safe_weighted_sum <- function(x, w){
  ok <- is.finite(x) & is.finite(w)
  if(!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok], na.rm = TRUE)
}

## Build all co-occurring pairs, including full_K = 0.
build_all_pair_table <- function(dataset){
  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  ## Operational consistency: an observed interaction must also count as a
  ## recorded-together cell. This avoids full_K > full_n due to loader quirks.
  missing_int_cells <- anti_join(ints, cooc, by = c("site", "consumer", "resource", "pair_id"))
  if(nrow(missing_int_cells) > 0){
    warning(dataset, ": adding ", nrow(missing_int_cells),
            " observed interaction cells to co-occurrence table for consistency.")
    cooc <- bind_rows(cooc, missing_int_cells) %>%
      distinct(site, consumer, resource, pair_id)
  }

  n_table <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  k_table <- ints %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  pair_table <- n_table %>%
    left_join(k_table, by = c("consumer", "resource", "pair_id")) %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      ever_interacted = full_K >= 1,
      interaction_occupancy = full_K / full_n,
      n_is_one = full_n == 1
    )

  consumer_degree <- pair_table %>%
    filter(ever_interacted) %>%
    group_by(dataset, consumer) %>%
    summarise(full_degree = n_distinct(resource), .groups = "drop") %>%
    assign_degree_groups()

  ## Keep the analysis focused on consumers represented in the observed regional
  ## interaction network, because realised degree is undefined for consumers with
  ## no realised partner.
  pair_table <- pair_table %>%
    inner_join(consumer_degree, by = c("dataset", "consumer")) %>%
    select(dataset, consumer, resource, full_degree, `consumer-degree group`,
           full_n, full_K, ever_interacted, interaction_occupancy, n_is_one)

  degree_groups <- consumer_degree %>%
    group_by(dataset, `consumer-degree group`) %>%
    summarise(
      number_of_consumers = n_distinct(consumer),
      minimum_full_degree = min(full_degree, na.rm = TRUE),
      median_full_degree = median(full_degree, na.rm = TRUE),
      maximum_full_degree = max(full_degree, na.rm = TRUE),
      .groups = "drop"
    )

  list(pair_table = pair_table, degree_groups = degree_groups)
}

summarise_by_exact_n <- function(pair_table){
  group_counts <- pair_table %>%
    group_by(dataset, full_n, `consumer-degree group`) %>%
    summarise(number_of_cooccurring_pairs = n(), .groups = "drop")

  eligible_n <- group_counts %>%
    group_by(dataset, full_n) %>%
    summarise(
      n_groups_present = n_distinct(`consumer-degree group`),
      min_pairs_any_group = min(number_of_cooccurring_pairs),
      has_lower = any(`consumer-degree group` == "Lower initial degree"),
      has_middle = any(`consumer-degree group` == "Middle initial degree"),
      has_higher = any(`consumer-degree group` == "Higher initial degree"),
      eligible_exact_n = n_groups_present == 3 &&
        has_lower && has_middle && has_higher &&
        min_pairs_any_group >= min_pairs_per_group_per_n,
      .groups = "drop"
    )

  exact_summary <- pair_table %>%
    group_by(dataset, full_n, `consumer-degree group`) %>%
    summarise(
      number_of_cooccurring_pairs = n(),
      number_of_consumers = n_distinct(consumer),
      number_ever_interacted = sum(ever_interacted, na.rm = TRUE),
      proportion_ever_interacted = mean(ever_interacted, na.rm = TRUE),
      mean_full_K = mean(full_K, na.rm = TRUE),
      median_full_K = median(full_K, na.rm = TRUE),
      mean_occupancy_across_all_cooccurring_pairs = mean(interaction_occupancy, na.rm = TRUE),
      median_occupancy_across_all_cooccurring_pairs = median(interaction_occupancy, na.rm = TRUE),
      mean_occupancy_among_ever_interacted_pairs = ifelse(any(ever_interacted),
                                                          mean(interaction_occupancy[ever_interacted], na.rm = TRUE),
                                                          NA_real_),
      median_occupancy_among_ever_interacted_pairs = ifelse(any(ever_interacted),
                                                            median(interaction_occupancy[ever_interacted], na.rm = TRUE),
                                                            NA_real_),
      .groups = "drop"
    ) %>%
    left_join(eligible_n %>% select(dataset, full_n, eligible_exact_n),
              by = c("dataset", "full_n")) %>%
    mutate(eligible_exact_n = replace_na(eligible_exact_n, FALSE))

  list(exact_summary = exact_summary, eligible_n = eligible_n)
}

make_reference_weights <- function(pair_table, eligible_n){
  eligible_values <- eligible_n %>%
    filter(eligible_exact_n) %>%
    pull(full_n)

  pair_table %>%
    filter(full_n %in% eligible_values) %>%
    count(dataset, full_n, name = "n_pairs_reference") %>%
    group_by(dataset) %>%
    mutate(reference_weight = n_pairs_reference / sum(n_pairs_reference)) %>%
    ungroup()
}

standardise_from_pair_table <- function(pair_table, reference_weights){
  if(nrow(reference_weights) == 0){
    return(tibble(
      dataset = character(),
      outcome = character(),
      result_type = character(),
      `consumer-degree group` = character(),
      estimate = numeric(),
      bootstrap_replicate = integer()
    ))
  }

  by_n_group <- pair_table %>%
    semi_join(reference_weights, by = c("dataset", "full_n")) %>%
    group_by(dataset, full_n, `consumer-degree group`) %>%
    summarise(
      value_ever = mean(ever_interacted, na.rm = TRUE),
      value_occupancy_all = mean(interaction_occupancy, na.rm = TRUE),
      value_occupancy_conditional = ifelse(any(ever_interacted),
                                           mean(interaction_occupancy[ever_interacted], na.rm = TRUE),
                                           NA_real_),
      .groups = "drop"
    ) %>%
    tidyr::complete(dataset, full_n, `consumer-degree group` = factor(degree_group_levels, levels = degree_group_levels)) %>%
    left_join(reference_weights %>% select(dataset, full_n, reference_weight),
              by = c("dataset", "full_n"))

  bind_rows(
    by_n_group %>%
      group_by(dataset, `consumer-degree group`) %>%
      summarise(estimate = safe_weighted_sum(value_ever, reference_weight), .groups = "drop") %>%
      mutate(outcome = "Ever interacted"),
    by_n_group %>%
      group_by(dataset, `consumer-degree group`) %>%
      summarise(estimate = safe_weighted_sum(value_occupancy_all, reference_weight), .groups = "drop") %>%
      mutate(outcome = "Occupancy across all co-occurring pairs"),
    by_n_group %>%
      group_by(dataset, `consumer-degree group`) %>%
      summarise(estimate = safe_weighted_sum(value_occupancy_conditional, reference_weight), .groups = "drop") %>%
      mutate(outcome = "Occupancy conditional on ever interacting")
  ) %>%
    mutate(result_type = "Degree group standardized value") %>%
    select(dataset, outcome, result_type, `consumer-degree group`, estimate)
}

make_contrasts <- function(std_values){
  if(nrow(std_values) == 0) return(tibble())

  wide <- std_values %>%
    filter(result_type == "Degree group standardized value") %>%
    select(dataset, outcome, `consumer-degree group`, estimate) %>%
    tidyr::pivot_wider(names_from = `consumer-degree group`, values_from = estimate)

  bind_rows(
    wide %>%
      transmute(dataset, outcome,
                result_type = "Higher minus lower",
                `consumer-degree group` = NA_character_,
                estimate = `Higher initial degree` - `Lower initial degree`),
    wide %>%
      transmute(dataset, outcome,
                result_type = "Middle minus lower",
                `consumer-degree group` = NA_character_,
                estimate = `Middle initial degree` - `Lower initial degree`)
  )
}

bootstrap_one_dataset <- function(pair_table, reference_weights, n_boot){
  if(nrow(reference_weights) == 0){
    return(tibble(
      dataset = character(),
      outcome = character(),
      result_type = character(),
      `consumer-degree group` = character(),
      estimate = numeric(),
      bootstrap_replicate = integer()
    ))
  }

  consumers_by_group <- pair_table %>%
    distinct(`consumer-degree group`, consumer) %>%
    filter(!is.na(`consumer-degree group`)) %>%
    group_by(`consumer-degree group`) %>%
    summarise(consumers = list(consumer), n_consumers = n(), .groups = "drop")

  if(nrow(consumers_by_group) == 0){
    return(tibble(
      dataset = character(),
      outcome = character(),
      result_type = character(),
      `consumer-degree group` = character(),
      estimate = numeric(),
      bootstrap_replicate = integer()
    ))
  }

  boot_rows <- vector("list", n_boot)

  for(b in seq_len(n_boot)){
    sampled <- bind_rows(lapply(seq_len(nrow(consumers_by_group)), function(i){
      g <- consumers_by_group$`consumer-degree group`[i]
      cur_consumers <- consumers_by_group$consumers[[i]]
      sampled_consumers <- sample(cur_consumers, size = length(cur_consumers), replace = TRUE)

      bind_rows(lapply(seq_along(sampled_consumers), function(j){
        pair_table %>%
          filter(`consumer-degree group` == g, consumer == sampled_consumers[j]) %>%
          mutate(bootstrap_consumer_id = paste0(as.character(g), "__", j, "__", sampled_consumers[j]))
      }))
    }))

    std <- standardise_from_pair_table(sampled, reference_weights)
    boot_rows[[b]] <- bind_rows(std, make_contrasts(std)) %>%
      mutate(bootstrap_replicate = b)
  }

  bind_rows(boot_rows)
}

bootstrap_exact_n_intervals <- function(pair_table, eligible_n, n_boot){
  eligible_values <- eligible_n %>%
    filter(eligible_exact_n) %>%
    pull(full_n)

  if(length(eligible_values) == 0){
    return(tibble(
      dataset = character(),
      full_n = integer(),
      `consumer-degree group` = character(),
      boot_q025 = numeric(),
      boot_q975 = numeric()
    ))
  }

  pair_table_eligible <- pair_table %>%
    filter(full_n %in% eligible_values)

  consumers_by_group <- pair_table_eligible %>%
    distinct(`consumer-degree group`, consumer) %>%
    filter(!is.na(`consumer-degree group`)) %>%
    group_by(`consumer-degree group`) %>%
    summarise(consumers = list(consumer), .groups = "drop")

  if(nrow(consumers_by_group) == 0){
    return(tibble(
      dataset = character(),
      full_n = integer(),
      `consumer-degree group` = character(),
      boot_q025 = numeric(),
      boot_q975 = numeric()
    ))
  }

  boot_rows <- vector("list", n_boot)

  for(b in seq_len(n_boot)){
    sampled <- bind_rows(lapply(seq_len(nrow(consumers_by_group)), function(i){
      g <- consumers_by_group$`consumer-degree group`[i]
      cur_consumers <- consumers_by_group$consumers[[i]]
      sampled_consumers <- sample(cur_consumers, size = length(cur_consumers), replace = TRUE)

      bind_rows(lapply(seq_along(sampled_consumers), function(j){
        pair_table_eligible %>%
          filter(`consumer-degree group` == g, consumer == sampled_consumers[j])
      }))
    }))

    boot_rows[[b]] <- sampled %>%
      group_by(dataset, full_n, `consumer-degree group`) %>%
      summarise(proportion_ever_interacted = mean(ever_interacted, na.rm = TRUE), .groups = "drop") %>%
      mutate(bootstrap_replicate = b)
  }

  bind_rows(boot_rows) %>%
    group_by(dataset, full_n, `consumer-degree group`) %>%
    summarise(
      boot_q025 = quantile(proportion_ever_interacted, 0.025, na.rm = TRUE),
      boot_q975 = quantile(proportion_ever_interacted, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_bootstrap <- function(observed_std, boot_std, reference_weights, pair_table){
  if(nrow(observed_std) == 0){
    return(tibble())
  }

  observed_all <- bind_rows(observed_std, make_contrasts(observed_std))

  boot_ci <- boot_std %>%
    group_by(dataset, outcome, result_type, `consumer-degree group`) %>%
    summarise(
      bootstrap_q025 = quantile(estimate, 0.025, na.rm = TRUE),
      bootstrap_q975 = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  ref_info <- reference_weights %>%
    group_by(dataset) %>%
    summarise(
      number_of_eligible_exact_n_strata = n_distinct(full_n),
      minimum_eligible_n = min(full_n, na.rm = TRUE),
      maximum_eligible_n = max(full_n, na.rm = TRUE),
      .groups = "drop"
    )

  eligible_pairs <- pair_table %>%
    semi_join(reference_weights, by = c("dataset", "full_n")) %>%
    group_by(dataset) %>%
    summarise(
      total_eligible_cooccurring_pairs = n(),
      total_eligible_consumers = n_distinct(consumer),
      .groups = "drop"
    )

  observed_all %>%
    left_join(boot_ci, by = c("dataset", "outcome", "result_type", "consumer-degree group")) %>%
    left_join(ref_info, by = "dataset") %>%
    left_join(eligible_pairs, by = "dataset")
}

make_n1_table <- function(pair_table){
  pair_table %>%
    filter(full_n == 1) %>%
    group_by(dataset, `consumer-degree group`) %>%
    summarise(
      number_of_cooccurring_pairs_with_full_n_1 = n(),
      number_with_full_K_0 = sum(full_K == 0, na.rm = TRUE),
      fraction_with_full_K_0 = number_with_full_K_0 / number_of_cooccurring_pairs_with_full_n_1,
      number_with_full_K_1 = sum(full_K == 1, na.rm = TRUE),
      fraction_with_full_K_1 = number_with_full_K_1 / number_of_cooccurring_pairs_with_full_n_1,
      .groups = "drop"
    )
}

validate_dataset <- function(dataset, pair_table, exact_summary, reference_weights, std_summary){
  includes_K0 <- any(pair_table$full_K == 0)
  full_n_ok <- all(pair_table$full_n >= 1)
  K_bounds_ok <- all(pair_table$full_K >= 0 & pair_table$full_K <= pair_table$full_n)

  n1 <- pair_table %>% filter(full_n == 1)
  n1_has_K0 <- ifelse(nrow(n1) == 0, NA, any(n1$full_K == 0))
  n1_has_K1 <- ifelse(nrow(n1) == 0, NA, any(n1$full_K == 1))

  conditional_has_K0 <- pair_table %>%
    filter(ever_interacted) %>%
    summarise(any_K0 = any(full_K == 0)) %>%
    pull(any_K0)

  eligible_counts_ok <- exact_summary %>%
    filter(eligible_exact_n) %>%
    group_by(full_n) %>%
    summarise(min_pairs = min(number_of_cooccurring_pairs), n_groups = n_distinct(`consumer-degree group`), .groups = "drop") %>%
    summarise(ok = ifelse(n() == 0, TRUE, all(min_pairs >= min_pairs_per_group_per_n & n_groups == 3))) %>%
    pull(ok)

  weights_sum_ok <- if(nrow(reference_weights) == 0){
    TRUE
  } else {
    reference_weights %>%
      group_by(dataset) %>%
      summarise(weight_sum = sum(reference_weight), .groups = "drop") %>%
      summarise(ok = all(abs(weight_sum - 1) < 1e-8)) %>%
      pull(ok)
  }

  tibble(
    dataset = dataset,
    includes_pairs_with_full_K_0 = includes_K0,
    all_pairs_have_full_n_at_least_1 = full_n_ok,
    all_pairs_satisfy_0_le_full_K_le_full_n = K_bounds_ok,
    n1_pairs_have_full_K_0_where_present = n1_has_K0,
    n1_pairs_have_full_K_1_where_present = n1_has_K1,
    no_full_K_0_in_conditional_occupancy_calculations = !conditional_has_K0,
    all_eligible_exact_n_have_at_least_10_pairs_in_every_degree_group = eligible_counts_ok,
    standardized_reference_weights_sum_to_1 = weights_sum_ok,
    bootstrap_resampling_consumer_level_retains_all_pairs = TRUE,
    no_site_removal_model_null_or_refitting_used = TRUE,
    all_checks_pass = includes_K0 && full_n_ok && K_bounds_ok && isTRUE(!conditional_has_K0) &&
      isTRUE(eligible_counts_ok) && isTRUE(weights_sum_ok)
  )
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

  message("Running script 26 all co-occurrence opportunities and degree: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  built <- build_all_pair_table(dataset)
  pair_table <- built$pair_table
  degree_groups <- built$degree_groups

  exact_objects <- summarise_by_exact_n(pair_table)
  exact_summary <- exact_objects$exact_summary
  eligible_n <- exact_objects$eligible_n

  reference_weights <- make_reference_weights(pair_table, eligible_n)
  observed_std <- standardise_from_pair_table(pair_table, reference_weights)
  boot_std <- bootstrap_one_dataset(pair_table, reference_weights, n_boot)
  standardized_summary <- summarise_bootstrap(observed_std, boot_std, reference_weights, pair_table)

  exact_boot_intervals <- bootstrap_exact_n_intervals(pair_table, eligible_n, n_boot)
  exact_summary_plot <- exact_summary %>%
    left_join(exact_boot_intervals, by = c("dataset", "full_n", "consumer-degree group"))

  n1_table <- make_n1_table(pair_table)

  pair_table_out <- pair_table %>%
    left_join(eligible_n %>% select(dataset, full_n, eligible_exact_n), by = c("dataset", "full_n")) %>%
    mutate(eligible_exact_n_main_analysis = replace_na(eligible_exact_n, FALSE)) %>%
    select(dataset, consumer, resource, full_degree, `consumer-degree group`,
           full_n, full_K, ever_interacted, interaction_occupancy,
           n_is_one, eligible_exact_n_main_analysis)

  checks <- validate_dataset(dataset, pair_table, exact_summary, reference_weights, standardized_summary)

  write.csv2(pair_table_out, file.path(out_dir, paste0(dataset, "_26_all_cooccurrence_pair_data.csv")), row.names = FALSE)
  write.csv2(exact_summary, file.path(out_dir, paste0(dataset, "_26_ever_interacted_and_occupancy_by_exact_n.csv")), row.names = FALSE)
  write.csv2(standardized_summary, file.path(out_dir, paste0(dataset, "_26_standardized_degree_comparisons.csv")), row.names = FALSE)
  write.csv2(n1_table, file.path(out_dir, paste0(dataset, "_26_n1_outcomes_by_degree.csv")), row.names = FALSE)
  write.csv2(degree_groups, file.path(out_dir, paste0(dataset, "_26_initial_consumer_degree_groups.csv")), row.names = FALSE)
  write.csv2(checks, file.path(out_dir, paste0(dataset, "_26_all_cooccurrence_degree_checks.csv")), row.names = FALSE)

  list(
    pair_table = pair_table_out,
    exact_summary = exact_summary,
    exact_summary_plot = exact_summary_plot,
    standardized_summary = standardized_summary,
    n1_table = n1_table,
    degree_groups = degree_groups,
    checks = checks
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)
names(all_outputs) <- all_dataset_names
future::plan(future::sequential)

pair_data_all <- bind_rows(lapply(all_outputs, `[[`, "pair_table"))
exact_summary_all <- bind_rows(lapply(all_outputs, `[[`, "exact_summary"))
exact_summary_plot_all <- bind_rows(lapply(all_outputs, `[[`, "exact_summary_plot"))
standardized_summary_all <- bind_rows(lapply(all_outputs, `[[`, "standardized_summary"))
n1_table_all <- bind_rows(lapply(all_outputs, `[[`, "n1_table"))
degree_groups_all <- bind_rows(lapply(all_outputs, `[[`, "degree_groups"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(pair_data_all, file.path(combined_out, "26_all_cooccurrence_pair_data.csv"), row.names = FALSE)
write.csv2(exact_summary_all, file.path(combined_out, "26_ever_interacted_and_occupancy_by_exact_n.csv"), row.names = FALSE)
write.csv2(standardized_summary_all, file.path(combined_out, "26_standardized_degree_comparisons.csv"), row.names = FALSE)
write.csv2(n1_table_all, file.path(combined_out, "26_n1_outcomes_by_degree.csv"), row.names = FALSE)
write.csv2(degree_groups_all, file.path(combined_out, "26_initial_consumer_degree_groups.csv"), row.names = FALSE)
write.csv2(checks_all, file.path(combined_out, "26_all_cooccurrence_degree_checks.csv"), row.names = FALSE)

## ---------------------------
## Figures
## ---------------------------

degree_colours <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#7570b3",
  "Higher initial degree" = "#d95f02"
)

## Figure 1: exact-n ever interacted
plot_exact <- exact_summary_plot_all %>%
  filter(eligible_exact_n) %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = degree_group_levels)
  )

p1 <- ggplot(plot_exact,
             aes(x = full_n,
                 y = proportion_ever_interacted,
                 colour = `consumer-degree group`,
                 group = `consumer-degree group`)) +
  geom_ribbon(aes(ymin = boot_q025, ymax = boot_q975, fill = `consumer-degree group`),
              alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = degree_colours, drop = FALSE) +
  scale_fill_manual(values = degree_colours, drop = FALSE) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 10) +
  xlab("Number of sites where pair was recorded together") +
  ylab("Fraction of co-occurring pairs ever observed interacting") +
  labs(colour = "Consumer-degree group", fill = "Consumer-degree group",
       subtitle = "All co-occurring pairs included; comparisons shown only where every degree group has at least 10 pairs.")

ggsave(file.path(combined_out, "26_ever_interacted_given_exact_cooccurrence_support.png"),
       p1, width = 14, height = 7, dpi = 300)

## Figure 2: standardized ever interacted and higher-minus-lower
std_ever <- standardized_summary_all %>%
  filter(outcome == "Ever interacted") %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = degree_group_levels)
  )

p2_top <- std_ever %>%
  filter(result_type == "Degree group standardized value") %>%
  ggplot(aes(x = dataset, y = estimate, colour = `consumer-degree group`)) +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
                width = 0.15, position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = degree_colours, drop = FALSE) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Standardized probability of ever interacting") +
  labs(colour = "Consumer-degree group")

p2_bottom <- std_ever %>%
  filter(result_type == "Higher minus lower") %>%
  ggplot(aes(x = dataset, y = estimate)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
                width = 0.15, colour = "grey25") +
  geom_point(size = 2, colour = "grey10") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Higher minus lower")

if(requireNamespace("patchwork", quietly = TRUE)){
  p2 <- p2_top / p2_bottom
  ggsave(file.path(combined_out, "26_standardized_ever_interacted_by_initial_degree.png"),
         p2, width = 12, height = 8, dpi = 300)
} else {
  ggsave(file.path(combined_out, "26_standardized_ever_interacted_by_initial_degree.png"),
         p2_top, width = 12, height = 5, dpi = 300)
}

## Figure 3: standardized occupancy across all and conditional
occ_plot <- standardized_summary_all %>%
  filter(result_type == "Degree group standardized value",
         outcome %in% c("Occupancy across all co-occurring pairs",
                        "Occupancy conditional on ever interacting")) %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = degree_group_levels),
    outcome_label = recode(
      outcome,
      "Occupancy across all co-occurring pairs" = "Across all co-occurring pairs",
      "Occupancy conditional on ever interacting" = "Among pairs ever observed interacting"
    )
  )

p3 <- ggplot(occ_plot,
             aes(x = dataset, y = estimate, colour = `consumer-degree group`)) +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
                width = 0.15, position = position_dodge(width = 0.65)) +
  geom_point(size = 2, position = position_dodge(width = 0.65)) +
  scale_colour_manual(values = degree_colours, drop = FALSE) +
  facet_grid(outcome_label ~ ., scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Standardized mean interaction occupancy") +
  labs(colour = "Consumer-degree group")

ggsave(file.path(combined_out, "26_standardized_occupancy_all_and_conditional.png"),
       p3, width = 12, height = 7, dpi = 300)

## Figure 4: n = 1 outcomes
n1_plot <- n1_table_all %>%
  select(dataset, `consumer-degree group`, fraction_with_full_K_0, fraction_with_full_K_1) %>%
  pivot_longer(cols = c(fraction_with_full_K_0, fraction_with_full_K_1),
               names_to = "outcome", values_to = "fraction") %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = degree_group_levels),
    outcome = recode(outcome,
                     fraction_with_full_K_0 = "full_K = 0",
                     fraction_with_full_K_1 = "full_K = 1")
  )

p4 <- ggplot(n1_plot,
             aes(x = `consumer-degree group`, y = fraction, fill = outcome)) +
  geom_col(position = "fill") +
  facet_wrap(~ dataset, ncol = 5) +
  scale_fill_manual(values = c("full_K = 0" = "grey70", "full_K = 1" = "grey20")) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Consumer-degree group") +
  ylab("Fraction of full_n = 1 co-occurring pairs") +
  labs(fill = "Outcome")

ggsave(file.path(combined_out, "26_n1_outcomes_by_initial_degree.png"),
       p4, width = 14, height = 7, dpi = 300)

## ---------------------------
## Console output
## ---------------------------

message("\nValidation summary:")
print(checks_all %>% select(dataset, all_checks_pass, includes_pairs_with_full_K_0,
                            all_pairs_satisfy_0_le_full_K_le_full_n,
                            standardized_reference_weights_sum_to_1), n = Inf)

message("\nDegree is defined from the observed regional interaction network. Differences among degree groups therefore describe the structure of co-occurrence opportunities associated with low- and high-degree consumers; they do not establish that degree causes an interaction to occur.")

message("\nFinished script 26. Outputs saved in:")
message("  ", combined_out)
