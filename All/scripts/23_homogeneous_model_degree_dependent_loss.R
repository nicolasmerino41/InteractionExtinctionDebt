## ------------------------------------------------------------
## Script: All/scripts/23_homogeneous_model_degree_dependent_loss.R
##
## Purpose:
## Compare empirical degree-dependent consumer loss patterns from script 22
## with the original homogeneous Galiana co-occurrence-based baseline.
##
## Uses the 10 original datasets only.
## Does not modify earlier scripts.
## Uses the shared helper script for dataset loading.
##
## Outputs:
##   All/CombinedOutputs/23_model_vs_empirical_*.png/csv
##   All/CombinedOutputs/23_homogeneous_model_degree_loss_checks.csv
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
n_model_sims <- 500
min_sites_retained <- 1

result_type <- "script23_homogeneous_model_degree_dependent_loss"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " workers.")

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script is intended for the original 10 datasets.")
}

## ---------------------------
## General helpers
## ---------------------------

read_semicolon_or_comma <- function(file){
  if(!file.exists(file)) return(NULL)
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1) return(x)
  tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
}

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

safe_cor <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

summ_quant <- function(x){
  x <- x[is.finite(x)]
  if(length(x) == 0){
    return(data.frame(median = NA_real_, q025 = NA_real_, q975 = NA_real_))
  }
  data.frame(
    median = median(x, na.rm = TRUE),
    q025 = as.numeric(quantile(x, 0.025, na.rm = TRUE)),
    q975 = as.numeric(quantile(x, 0.975, na.rm = TRUE))
  )
}

conditioned_link_prob <- function(n, N_alpha, p){
  denom <- 1 - (1 - p)^N_alpha
  out <- (1 - (1 - p)^n) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

expected_conditioned_links <- function(pair_table, p){
  sum(conditioned_link_prob(pair_table$n, pair_table$N_alpha, p), na.rm = TRUE)
}

fit_conditioned_p <- function(pair_table, observed_links){
  if(nrow(pair_table) == 0 || is.na(observed_links)) return(NA_real_)
  if(observed_links <= 0) return(0)
  max_expected <- expected_conditioned_links(pair_table, 1 - 1e-12)
  if(observed_links >= max_expected) return(1 - 1e-12)
  f <- function(p){ expected_conditioned_links(pair_table, p) - observed_links }
  tryCatch(uniroot(f, interval = c(1e-12, 1 - 1e-12), tol = 1e-10)$root,
           error = function(e) NA_real_)
}

get_fitted_p_from_script16 <- function(dataset, pair_table, observed_links){
  f16 <- file.path(combined_out, "16_consumer_generalism_repeatability_summary.csv")
  ## script 16 is usually in All/CombinedOutputs, not this result folder
  f16_global <- file.path("All", "CombinedOutputs", "16_consumer_generalism_repeatability_summary.csv")
  tab <- read_semicolon_or_comma(f16_global)
  if(!is.null(tab) && "dataset" %in% names(tab)){
    p_cols <- c("fitted_conditioned_p", "conditioned_p", "p_conditioned", "fitted_p")
    p_col <- p_cols[p_cols %in% names(tab)][1]
    if(!is.na(p_col)){
      val <- tab[tab$dataset == dataset, p_col][1]
      if(length(val) == 1 && is.finite(as.numeric(val))) return(as.numeric(val))
    }
  }
  fit_conditioned_p(pair_table, observed_links)
}

assign_degree_groups <- function(degree_table){
  degree_table <- degree_table %>%
    arrange(full_degree, consumer) %>%
    mutate(rank_index = row_number(), n_consumers = n())

  if(nrow(degree_table) == 0){
    degree_table$initial_degree_group <- character()
    return(degree_table)
  }

  unique_degrees <- sort(unique(degree_table$full_degree))

  if(length(unique_degrees) == 1){
    degree_table$initial_degree_group <- "Middle initial degree"
    return(degree_table %>% select(-rank_index, -n_consumers))
  }

  degree_counts <- degree_table %>%
    count(full_degree, name = "n_degree") %>%
    arrange(full_degree) %>%
    mutate(cum_n = cumsum(n_degree), midpoint = cum_n - n_degree / 2,
           rel_mid = midpoint / sum(n_degree))

  degree_counts <- degree_counts %>%
    mutate(initial_degree_group = case_when(
      rel_mid <= 1/3 ~ "Lower initial degree",
      rel_mid <= 2/3 ~ "Middle initial degree",
      TRUE ~ "Higher initial degree"
    ))

  degree_table %>%
    left_join(degree_counts %>% select(full_degree, initial_degree_group), by = "full_degree") %>%
    select(-rank_index, -n_consumers)
}

make_group_definition <- function(grouped_degree, dataset, replicate = NA_integer_, source = "Empirical data"){
  grouped_degree %>%
    group_by(initial_degree_group) %>%
    summarise(
      dataset = dataset,
      source = source,
      replicate = replicate,
      number_of_consumers = n_distinct(consumer),
      minimum_full_degree = min(full_degree, na.rm = TRUE),
      maximum_full_degree = max(full_degree, na.rm = TRUE),
      median_full_degree = median(full_degree, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(dataset, source, replicate, initial_degree_group,
           number_of_consumers, minimum_full_degree, maximum_full_degree,
           median_full_degree)
}

make_pair_table <- function(dataset, site_tables){
  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource))

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource))

  ## Operational safeguard: an observed interaction implies recorded together.
  cooc <- bind_rows(cooc, ints %>% select(site, consumer, resource, pair_id)) %>%
    distinct(site, consumer, resource, pair_id)

  pair_table <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(n = n_distinct(site), .groups = "drop") %>%
    left_join(
      ints %>%
        group_by(pair_id) %>%
        summarise(K = n_distinct(site), .groups = "drop"),
      by = "pair_id"
    ) %>%
    mutate(K = replace_na(K, 0L),
           is_realised_link = K > 0) %>%
    group_by(consumer) %>%
    mutate(N_alpha = sum(n, na.rm = TRUE)) %>%
    ungroup()

  list(pair_table = pair_table, cooc = cooc, interactions = ints)
}

make_or_reuse_site_subsets <- function(dataset, all_sites){
  ## Reuse script 22 subsets only if a saved subset file exists. Most previous scripts
  ## did not save subsets, so the fallback is the original reproducible helper.
  candidate_files <- c(
    file.path("All", "SeparatedResults", "script22_degree_dependent_opportunity_and_interaction_loss", dataset,
              paste0(dataset, "_22_site_removal_subsets.rds")),
    file.path("All", "SeparatedResults", "script19_describe_layers_of_interaction_loss", dataset,
              paste0(dataset, "_19_site_removal_subsets.rds"))
  )
  existing <- candidate_files[file.exists(candidate_files)][1]
  if(!is.na(existing)){
    return(readRDS(existing))
  }
  make_site_subsets(all_sites, removal_levels, n_site_reps)
}

## ---------------------------
## Simulation helpers
## ---------------------------

simulate_conditioned_full_network <- function(cooc, consumers, p){
  out_list <- vector("list", length(consumers))
  names(out_list) <- consumers

  for(con in consumers){
    cells <- cooc %>% filter(consumer == con)
    if(nrow(cells) == 0) next

    tries <- 0L
    repeat{
      tries <- tries + 1L
      draw <- rbinom(nrow(cells), size = 1, prob = p)
      if(any(draw == 1) || tries > 10000L) break
    }

    if(any(draw == 1)){
      out_list[[con]] <- cells[draw == 1, c("site", "consumer", "resource", "pair_id")]
    } else {
      ## Extremely unlikely fallback for numerical edge cases: force one observed cell.
      idx <- sample(seq_len(nrow(cells)), 1)
      out_list[[con]] <- cells[idx, c("site", "consumer", "resource", "pair_id")]
    }
  }

  bind_rows(out_list) %>% distinct(site, consumer, resource, pair_id)
}

consumer_metrics_for_subsets <- function(dataset, source, replicate, subset_object,
                                         cooc, interactions, full_links_grouped){
  consumers <- full_links_grouped$consumer
  full_degree_lookup <- full_links_grouped %>%
    select(consumer, full_degree, initial_degree_group)

  out <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    rem <- subset_object$index$removal_fraction[i]
    site_rep <- subset_object$index$site_rep[i]

    present_consumers <- cooc %>%
      filter(site %in% sites_keep, consumer %in% consumers) %>%
      distinct(consumer) %>%
      mutate(consumer_present = 1L)

    retained_n <- cooc %>%
      filter(site %in% sites_keep) %>%
      semi_join(full_links_grouped %>% select(pair_id), by = "pair_id") %>%
      group_by(pair_id) %>%
      summarise(retained_n = n_distinct(site), .groups = "drop")

    retained_K <- interactions %>%
      filter(site %in% sites_keep) %>%
      semi_join(full_links_grouped %>% select(pair_id), by = "pair_id") %>%
      group_by(pair_id) %>%
      summarise(retained_K = n_distinct(site), .groups = "drop")

    link_states <- full_links_grouped %>%
      select(consumer, resource, pair_id, full_degree, initial_degree_group) %>%
      left_join(retained_n, by = "pair_id") %>%
      left_join(retained_K, by = "pair_id") %>%
      mutate(
        retained_n = replace_na(retained_n, 0L),
        retained_K = replace_na(retained_K, 0L),
        opportunity_retained = retained_n > 0,
        interaction_retained = retained_K > 0,
        missing_despite_opportunity = retained_n > 0 & retained_K == 0
      )

    consumer_rows <- link_states %>%
      group_by(consumer, full_degree, initial_degree_group) %>%
      summarise(
        opportunity_retention = mean(opportunity_retained),
        interaction_retention = mean(interaction_retained),
        interaction_missing_despite_opportunity = mean(missing_despite_opportunity),
        .groups = "drop"
      ) %>%
      left_join(present_consumers, by = "consumer") %>%
      mutate(
        consumer_present = replace_na(consumer_present, 0L),
        dataset = dataset,
        source = source,
        replicate = replicate,
        removal_fraction = rem,
        site_rep = site_rep
      ) %>%
      select(dataset, source, replicate, removal_fraction, site_rep,
             consumer, initial_degree_group, full_degree, consumer_present,
             opportunity_retention, interaction_retention,
             interaction_missing_despite_opportunity)

    out[[i]] <- consumer_rows
  }

  bind_rows(out)
}

summarise_group_by_subset <- function(consumer_subset_data){
  consumer_subset_data %>%
    group_by(dataset, source, replicate, removal_fraction, site_rep, initial_degree_group) %>%
    summarise(
      number_of_consumers_contributing = n_distinct(consumer),
      `consumer presence` = mean(consumer_present, na.rm = TRUE),
      `opportunity retention` = mean(opportunity_retention, na.rm = TRUE),
      `interaction retention` = mean(interaction_retention, na.rm = TRUE),
      `interaction missing despite opportunity` = mean(interaction_missing_despite_opportunity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(`consumer presence`, `opportunity retention`,
               `interaction retention`, `interaction missing despite opportunity`),
      names_to = "metric",
      values_to = "group_mean"
    )
}

summarise_across_replicates <- function(group_by_subset){
  group_by_subset %>%
    group_by(dataset, source, removal_fraction, initial_degree_group, metric) %>%
    summarise(
      median = median(group_mean, na.rm = TRUE),
      q025 = as.numeric(quantile(group_mean, 0.025, na.rm = TRUE)),
      q975 = as.numeric(quantile(group_mean, 0.975, na.rm = TRUE)),
      number_of_consumers_contributing = median(number_of_consumers_contributing, na.rm = TRUE),
      number_of_simulation_networks = ifelse(first(source) == "Homogeneous model",
                                             n_distinct(replicate), NA_integer_),
      .groups = "drop"
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
  })

  message("Running script 23 homogeneous degree-dependent loss: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  built <- make_pair_table(dataset, site_tables)
  pair_table <- built$pair_table
  cooc <- built$cooc
  interactions <- built$interactions

  observed_links <- sum(pair_table$is_realised_link)
  p_fit <- get_fitted_p_from_script16(dataset, pair_table, observed_links)
  expected_links <- expected_conditioned_links(pair_table, p_fit)

  empirical_full_links <- pair_table %>%
    filter(is_realised_link) %>%
    select(consumer, resource, pair_id) %>%
    group_by(consumer) %>%
    mutate(full_degree = n_distinct(resource)) %>%
    ungroup() %>%
    left_join(
      pair_table %>% select(pair_id, n, K),
      by = "pair_id"
    )

  empirical_degrees <- empirical_full_links %>%
    distinct(consumer, full_degree) %>%
    assign_degree_groups()

  empirical_full_links_grouped <- empirical_full_links %>%
    left_join(empirical_degrees %>% select(consumer, initial_degree_group), by = "consumer")

  empirical_group_def <- make_group_definition(empirical_degrees, dataset,
                                               replicate = NA_integer_,
                                               source = "Empirical data")

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_or_reuse_site_subsets(dataset, all_sites)
  ## Save actual subsets used so future scripts can reuse exactly these.
  saveRDS(subset_object, file.path(out_dir, paste0(dataset, "_23_site_removal_subsets.rds")))

  empirical_consumer_data <- consumer_metrics_for_subsets(
    dataset = dataset,
    source = "Empirical data",
    replicate = NA_integer_,
    subset_object = subset_object,
    cooc = cooc,
    interactions = interactions,
    full_links_grouped = empirical_full_links_grouped
  )

  empirical_group_subset <- summarise_group_by_subset(empirical_consumer_data)
  empirical_summary <- summarise_across_replicates(empirical_group_subset)

  consumers <- sort(unique(cooc$consumer))

  sim_outputs <- future.apply::future_lapply(
    seq_len(n_model_sims),
    function(sim_id){
      suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
      })

      sim_interactions <- simulate_conditioned_full_network(cooc, consumers, p_fit)

      sim_full_links <- sim_interactions %>%
        distinct(consumer, resource, pair_id) %>%
        group_by(consumer) %>%
        mutate(full_degree = n_distinct(resource)) %>%
        ungroup()

      sim_degrees <- sim_full_links %>%
        distinct(consumer, full_degree) %>%
        assign_degree_groups()

      sim_full_links_grouped <- sim_full_links %>%
        left_join(sim_degrees %>% select(consumer, initial_degree_group), by = "consumer")

      sim_consumer_data <- consumer_metrics_for_subsets(
        dataset = dataset,
        source = "Homogeneous model",
        replicate = sim_id,
        subset_object = subset_object,
        cooc = cooc,
        interactions = sim_interactions,
        full_links_grouped = sim_full_links_grouped
      )

      sim_group_subset <- summarise_group_by_subset(sim_consumer_data)

      list(
        group_subset = sim_group_subset,
        group_definition = make_group_definition(sim_degrees, dataset,
                                                 replicate = sim_id,
                                                 source = "Homogeneous model"),
        sim_check = data.frame(
          dataset = dataset,
          replicate = sim_id,
          all_consumers_have_partner_full_network =
            all(sim_degrees$full_degree >= 1) && nrow(sim_degrees) == length(consumers),
          n_consumers_with_zero_partners_after_removal =
            sum(sim_consumer_data$interaction_retention == 0 & sim_consumer_data$removal_fraction > 0, na.rm = TRUE),
          duplicate_site_cell_detected = any(duplicated(sim_interactions[, c("site", "consumer", "resource")])),
          interactions_outside_cooccurrence_detected =
            nrow(anti_join(sim_interactions, cooc, by = c("site", "consumer", "resource", "pair_id"))) > 0,
          interaction_retention_exceeds_opportunity =
            any(sim_consumer_data$interaction_retention > sim_consumer_data$opportunity_retention + 1e-12, na.rm = TRUE),
          decomposition_ok = all(abs((1 - sim_consumer_data$opportunity_retention) +
                                       sim_consumer_data$interaction_retention +
                                       sim_consumer_data$interaction_missing_despite_opportunity - 1) < 1e-8,
                                 na.rm = TRUE)
        )
      )
    },
    future.seed = TRUE
  )

  model_group_subset <- bind_rows(lapply(sim_outputs, `[[`, "group_subset"))
  model_group_def <- bind_rows(lapply(sim_outputs, `[[`, "group_definition"))
  sim_checks <- bind_rows(lapply(sim_outputs, `[[`, "sim_check"))
  model_summary <- summarise_across_replicates(model_group_subset)

  summary_all <- bind_rows(empirical_summary, model_summary)

  difference <- empirical_summary %>%
    select(dataset, removal_fraction, initial_degree_group, metric,
           empirical_median = median) %>%
    left_join(
      model_summary %>%
        select(dataset, removal_fraction, initial_degree_group, metric,
               model_median = median, model_q025 = q025, model_q975 = q975),
      by = c("dataset", "removal_fraction", "initial_degree_group", "metric")
    ) %>%
    mutate(
      empirical_minus_homogeneous_model_median = empirical_median - model_median,
      empirical_position_relative_to_model_range = case_when(
        is.na(empirical_median) | is.na(model_q025) | is.na(model_q975) ~ NA_character_,
        empirical_median < model_q025 ~ "below model 95% range",
        empirical_median > model_q975 ~ "above model 95% range",
        TRUE ~ "inside model 95% range"
      )
    )

  empirical_check_data <- empirical_consumer_data

  checks <- data.frame(
    dataset = dataset,
    fitted_p = p_fit,
    observed_full_network_link_number = observed_links,
    expected_full_network_link_number = expected_links,
    link_number_calibration_ok = abs(observed_links - expected_links) < 1e-6,
    every_consumer_has_partner_in_every_simulated_full_network =
      all(sim_checks$all_consumers_have_partner_full_network),
    consumers_can_have_zero_partners_after_removal =
      any(sim_checks$n_consumers_with_zero_partners_after_removal > 0),
    no_interaction_outside_cooccurrence =
      !any(sim_checks$interactions_outside_cooccurrence_detected),
    no_duplicate_site_cells =
      !any(sim_checks$duplicate_site_cell_detected),
    interaction_retention_never_exceeds_opportunity =
      !any(sim_checks$interaction_retention_exceeds_opportunity),
    decomposition_ok_every_consumer = all(sim_checks$decomposition_ok),
    p_not_refitted_after_removal = TRUE,
    empirical_and_model_use_identical_site_removal_subsets = TRUE,
    no_beta_binomial_or_null_or_degree_distribution_model_used = TRUE
  )

  write.csv2(empirical_group_subset,
             file.path(out_dir, paste0(dataset, "_23_empirical_group_by_subset.csv")),
             row.names = FALSE)
  write.csv2(model_group_subset,
             file.path(out_dir, paste0(dataset, "_23_model_group_by_subset.csv")),
             row.names = FALSE)
  write.csv2(difference,
             file.path(out_dir, paste0(dataset, "_23_curve_difference.csv")),
             row.names = FALSE)

  list(
    summary = summary_all,
    difference = difference,
    checks = checks,
    group_definitions = bind_rows(empirical_group_def, model_group_def)
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

summary_all <- bind_rows(lapply(all_outputs, `[[`, "summary"))
difference_all <- bind_rows(lapply(all_outputs, `[[`, "difference"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))
group_defs_all <- bind_rows(lapply(all_outputs, `[[`, "group_definitions"))

summary_all <- summary_all %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    initial_degree_group = factor(initial_degree_group,
                                  levels = c("Lower initial degree",
                                             "Middle initial degree",
                                             "Higher initial degree")),
    source = factor(source, levels = c("Homogeneous model", "Empirical data")),
    metric = factor(metric,
                    levels = c("consumer presence",
                               "opportunity retention",
                               "interaction retention",
                               "interaction missing despite opportunity"))
  )

difference_all <- difference_all %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    initial_degree_group = factor(initial_degree_group,
                                  levels = c("Lower initial degree",
                                             "Middle initial degree",
                                             "Higher initial degree")),
    metric = factor(metric,
                    levels = c("consumer presence",
                               "opportunity retention",
                               "interaction retention",
                               "interaction missing despite opportunity"))
  )

write.csv2(summary_all,
           file.path(combined_out, "23_model_vs_empirical_degree_loss_summary.csv"),
           row.names = FALSE)
write.csv2(difference_all,
           file.path(combined_out, "23_model_vs_empirical_curve_difference.csv"),
           row.names = FALSE)
write.csv2(checks_all,
           file.path(combined_out, "23_homogeneous_model_degree_loss_checks.csv"),
           row.names = FALSE)
write.csv2(group_defs_all,
           file.path(combined_out, "23_initial_degree_group_definitions_empirical_and_model.csv"),
           row.names = FALSE)

message("Validation summary for script 23:")
print(checks_all)

## ---------------------------
## Plots
## ---------------------------

group_colours <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#7570b3",
  "Higher initial degree" = "#d95f02"
)

plot_model_empirical <- function(data, metrics_keep, ylab_text, file_name,
                                 facet_formula = NULL, width = 14, height = 7){
  model_data <- data %>% filter(source == "Homogeneous model", metric %in% metrics_keep)
  emp_data <- data %>% filter(source == "Empirical data", metric %in% metrics_keep)

  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_ribbon(
      data = model_data,
      aes(x = removal_fraction, ymin = q025, ymax = q975,
          fill = initial_degree_group),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(
      data = emp_data,
      aes(x = removal_fraction, y = median,
          colour = initial_degree_group),
      linewidth = 0.9
    ) +
    scale_colour_manual(values = group_colours, name = "Initial degree group") +
    scale_fill_manual(values = group_colours, name = "Homogeneous-model range") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom") +
    xlab("Proportion of sites removed") +
    ylab(ylab_text)

  if(is.null(facet_formula)){
    p <- p + facet_wrap(~ dataset, ncol = 5)
  } else {
    p <- p + facet_grid(facet_formula)
  }

  ggsave(file.path(combined_out, file_name), p, width = width, height = height, dpi = 300)
  p
}

p_presence <- plot_model_empirical(
  data = summary_all,
  metrics_keep = "consumer presence",
  ylab_text = "Fraction of consumers still present",
  file_name = "23_model_vs_empirical_consumer_presence.png",
  width = 14,
  height = 7
)

retention_data <- summary_all %>%
  filter(metric %in% c("opportunity retention", "interaction retention")) %>%
  mutate(metric_label = recode(as.character(metric),
                               `opportunity retention` = "Recorded together",
                               `interaction retention` = "Observed interacting"))

model_ret <- retention_data %>% filter(source == "Homogeneous model")
emp_ret <- retention_data %>% filter(source == "Empirical data")

p_ret <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_ribbon(
    data = model_ret,
    aes(x = removal_fraction, ymin = q025, ymax = q975,
        fill = initial_degree_group),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(
    data = emp_ret,
    aes(x = removal_fraction, y = median,
        colour = initial_degree_group),
    linewidth = 0.85
  ) +
  scale_colour_manual(values = group_colours, name = "Initial degree group") +
  scale_fill_manual(values = group_colours, name = "Homogeneous-model range") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(metric_label ~ dataset) +
  theme_classic(base_size = 8.5) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Proportion of sites removed") +
  ylab("Fraction of original interaction partners")

ggsave(file.path(combined_out, "23_model_vs_empirical_opportunity_and_interaction_retention.png"),
       p_ret, width = 16, height = 7.5, dpi = 300)

p_missing <- plot_model_empirical(
  data = summary_all,
  metrics_keep = "interaction missing despite opportunity",
  ylab_text = "Fraction recorded together but no longer interacting",
  file_name = "23_model_vs_empirical_interaction_missing_despite_opportunity.png",
  width = 14,
  height = 7
)

future::plan(future::sequential)

message("Finished script 23 homogeneous model degree-dependent loss comparison.")
