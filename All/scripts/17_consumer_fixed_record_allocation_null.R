## ------------------------------------------------------------
## Script: All/scripts/17_consumer_fixed_record_allocation_null.R
##
## Purpose:
## Consumer-specific fixed-record allocation null model following script 16.
##
## Main question:
## Could the relationship between consumer partner breadth and interaction
## repeatability arise simply because each consumer has a fixed number of
## observed local interaction records that must be distributed among its
## possible consumer-resource-site co-occurrence cells?
##
## Important limits:
##   - original 10 datasets only
##   - consumers only
##   - same conditioned original Galiana IR model as script 16
##   - no beta-binomial models
##   - no site-removal outcome
##   - no ecological/causal interpretation
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

packages_extra <- c("dplyr", "tidyr", "tibble", "ggplot2", "purrr")
for(pkg in packages_extra){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---------------------------
## Settings
## ---------------------------

n_allocation_sims <- 2000

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script is intended for the original 10 Galiana datasets only.")
}

dirs <- make_output_dirs("script17_consumer_fixed_record_allocation_null")
sep_out <- dirs$separated
combined_out <- dirs$combined

script16_summary_file <- file.path(
  "All", "outputs", "script16_consumer_generalism_and_repeatability",
  "combined", "16_consumer_generalism_repeatability_summary.csv"
)

## ---------------------------
## Helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

make_cell_id <- function(consumer, resource, site){
  paste(consumer, resource, site, sep = "___")
}

read_semicolon_or_comma <- function(file){
  if(!file.exists(file)) return(NULL)
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1) return(x)
  tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
}

plink_unconditioned <- function(n, p){
  1 - (1 - p)^n
}

p_consumer_has_link <- function(N_alpha, p){
  1 - (1 - p)^N_alpha
}

plink_conditioned_consumer <- function(n_ai, N_alpha, p){
  denom <- p_consumer_has_link(N_alpha, p)
  ifelse(denom > 0,
         plink_unconditioned(n_ai, p) / denom,
         NA_real_)
}

expected_repeatability_given_realised <- function(n, p){
  denom <- plink_unconditioned(n, p)
  ifelse(denom > 0, p / denom, NA_real_)
}

spearman_safe <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

mc_p_one_sided_lower <- function(empirical, null_values){
  null_values <- null_values[is.finite(null_values)]
  if(!is.finite(empirical) || length(null_values) == 0) return(NA_real_)
  (sum(null_values <= empirical) + 1) / (length(null_values) + 1)
}

consumer_N_table <- function(pair_table){
  pair_table %>%
    group_by(consumer) %>%
    summarise(
      N_alpha = sum(n_cooccurring_sites),
      number_cooccurring_resources = n_distinct(resource),
      .groups = "drop"
    )
}

expected_total_conditioned_links <- function(pair_table, p){
  N_tab <- consumer_N_table(pair_table)
  pair_table %>%
    left_join(N_tab, by = "consumer") %>%
    mutate(expected_link_probability = plink_conditioned_consumer(
      n_ai = n_cooccurring_sites,
      N_alpha = N_alpha,
      p = p
    )) %>%
    summarise(expected_total = sum(expected_link_probability, na.rm = TRUE)) %>%
    pull(expected_total)
}

fit_conditioned_p <- function(pair_table, observed_total_links){

  if(observed_total_links <= 0) return(0)
  if(observed_total_links >= nrow(pair_table)) return(1)

  f <- function(p){
    expected_total_conditioned_links(pair_table, p) - observed_total_links
  }

  f_low <- f(1e-12)
  f_high <- f(1 - 1e-12)

  if(!is.finite(f_low) || !is.finite(f_high)){
    return(NA_real_)
  }

  if(f_low * f_high > 0){
    opt <- optimize(function(p) abs(f(p)), interval = c(1e-12, 1 - 1e-12))
    return(opt$minimum)
  }

  uniroot(f, interval = c(1e-12, 1 - 1e-12), tol = 1e-12)$root
}

make_dataset_tables <- function(dataset){

  site_tables <- get_dataset_site_tables(dataset)

  cooc_cells <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      pair_id = make_pair_id(consumer, resource),
      cell_id = make_cell_id(consumer, resource, site)
    )

  empirical_cells <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      pair_id = make_pair_id(consumer, resource),
      cell_id = make_cell_id(consumer, resource, site)
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
      n_interacting_sites = replace_na(n_interacting_sites, 0L),
      is_realised_link = n_interacting_sites > 0L,
      observed_repeatability = ifelse(is_realised_link,
                                      n_interacting_sites / n_cooccurring_sites,
                                      NA_real_)
    ) %>%
    select(dataset, consumer, resource, pair_id,
           n_cooccurring_sites, n_interacting_sites,
           is_realised_link, observed_repeatability)

  list(
    pair_table = pair_table,
    cooc_cells = cooc_cells,
    empirical_cells = empirical_cells
  )
}

add_fixed_expectations <- function(pair_table, p){
  N_tab <- consumer_N_table(pair_table)

  pair_table %>%
    left_join(N_tab, by = "consumer") %>%
    mutate(
      expected_link_probability = plink_conditioned_consumer(
        n_ai = n_cooccurring_sites,
        N_alpha = N_alpha,
        p = p
      ),
      expected_repeatability = expected_repeatability_given_realised(
        n = n_cooccurring_sites,
        p = p
      )
    )
}

build_consumer_metrics_from_expected_table <- function(pair_expected){

  pair_with_excess <- pair_expected %>%
    mutate(
      repeatability_excess = observed_repeatability - expected_repeatability
    )

  partner_metrics <- pair_with_excess %>%
    group_by(dataset, consumer) %>%
    summarise(
      number_cooccurring_resources = first(number_cooccurring_resources),
      observed_realised_partners = sum(is_realised_link, na.rm = TRUE),
      expected_realised_partners = sum(expected_link_probability, na.rm = TRUE),
      `Realised partners above model expectation` =
        observed_realised_partners - expected_realised_partners,
      .groups = "drop"
    )

  repeatability_metrics <- pair_with_excess %>%
    filter(is_realised_link, n_cooccurring_sites >= 2) %>%
    group_by(dataset, consumer) %>%
    summarise(
      number_realised_links_n_ge_2 = n(),
      mean_observed_repeatability = mean(observed_repeatability, na.rm = TRUE),
      mean_expected_repeatability = mean(expected_repeatability, na.rm = TRUE),
      `Extra repeatability of realised links` =
        mean(repeatability_excess, na.rm = TRUE),
      .groups = "drop"
    )

  partner_metrics %>%
    left_join(repeatability_metrics, by = c("dataset", "consumer")) %>%
    filter(!is.na(`Extra repeatability of realised links`))
}

compute_consumer_correlation <- function(species_metrics){
  spearman_safe(
    species_metrics$`Realised partners above model expectation`,
    species_metrics$`Extra repeatability of realised links`
  )
}

simulate_fixed_record_allocation <- function(pair_expected, cooc_cells, empirical_cells){

  ## M_alpha is fixed from the observed empirical interaction records.
  consumers <- unique(pair_expected$consumer)

  M_by_consumer <- data.frame(consumer = consumers) %>%
    left_join(
      empirical_cells %>%
        group_by(consumer) %>%
        summarise(M_alpha = n_distinct(cell_id), .groups = "drop"),
      by = "consumer"
    ) %>%
    mutate(M_alpha = replace_na(M_alpha, 0L))

  selected_cells <- bind_rows(lapply(consumers, function(cons){
    possible <- cooc_cells %>% filter(consumer == cons)
    M_alpha <- M_by_consumer %>%
      filter(consumer == cons) %>%
      pull(M_alpha)

    if(length(M_alpha) == 0) M_alpha <- 0L

    if(M_alpha < 0 || M_alpha > nrow(possible)){
      stop("Invalid M_alpha for consumer ", cons,
           ": M_alpha = ", M_alpha,
           ", possible cells = ", nrow(possible))
    }

    if(M_alpha == 0){
      return(possible[0, ])
    }

    possible[sample(seq_len(nrow(possible)), size = M_alpha, replace = FALSE), ]
  }))

  K_sim <- selected_cells %>%
    group_by(pair_id) %>%
    summarise(n_interacting_sites_sim = n_distinct(cell_id), .groups = "drop")

  sim_pair <- pair_expected %>%
    select(-n_interacting_sites, -is_realised_link,
           -observed_repeatability) %>%
    left_join(K_sim, by = "pair_id") %>%
    mutate(
      n_interacting_sites = replace_na(n_interacting_sites_sim, 0L),
      is_realised_link = n_interacting_sites > 0L,
      observed_repeatability = ifelse(is_realised_link,
                                      n_interacting_sites / n_cooccurring_sites,
                                      NA_real_)
    ) %>%
    select(-n_interacting_sites_sim)

  list(
    sim_pair = sim_pair,
    selected_cells = selected_cells,
    M_by_consumer = M_by_consumer
  )
}

validate_simulation <- function(sim_object, pair_expected, cooc_cells, expected_partners_reference){

  sim_pair <- sim_object$sim_pair
  selected_cells <- sim_object$selected_cells
  M_by_consumer <- sim_object$M_by_consumer

  ## 1. Each consumer's simulated total equals observed M_alpha.
  sim_M <- selected_cells %>%
    group_by(consumer) %>%
    summarise(M_simulated = n_distinct(cell_id), .groups = "drop")

  M_check <- M_by_consumer %>%
    left_join(sim_M, by = "consumer") %>%
    mutate(
      M_simulated = replace_na(M_simulated, 0L),
      total_records_match = M_alpha == M_simulated
    )

  total_records_match_all <- all(M_check$total_records_match)

  ## 2. Every simulated K_ai is between 0 and n_ai.
  K_bounds_ok <- all(sim_pair$n_interacting_sites >= 0 &
                       sim_pair$n_interacting_sites <= sim_pair$n_cooccurring_sites)

  ## 3. No assignment outside allowed co-occurrence cells.
  outside_allowed_count <- selected_cells %>%
    anti_join(cooc_cells %>% select(cell_id), by = "cell_id") %>%
    nrow()
  outside_allowed_ok <- outside_allowed_count == 0L

  ## 4. No selected cell more than once.
  duplicate_selected_cell_count <- selected_cells %>%
    count(cell_id, name = "n_selected") %>%
    filter(n_selected > 1L) %>%
    nrow()
  no_duplicate_cells_ok <- duplicate_selected_cell_count == 0L

  ## 5. Expected-partner values remain unchanged.
  current_expected <- sim_pair %>%
    group_by(consumer) %>%
    summarise(expected_realised_partners = sum(expected_link_probability, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(consumer)

  reference <- expected_partners_reference %>% arrange(consumer)

  expected_partners_unchanged <- isTRUE(all.equal(
    current_expected$consumer,
    reference$consumer,
    check.attributes = FALSE
  )) && isTRUE(all.equal(
    current_expected$expected_realised_partners,
    reference$expected_realised_partners,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))

  list(
    ok = total_records_match_all && K_bounds_ok && outside_allowed_ok &&
      no_duplicate_cells_ok && expected_partners_unchanged,
    M_check = M_check,
    total_records_match_all = total_records_match_all,
    K_bounds_ok = K_bounds_ok,
    outside_allowed_ok = outside_allowed_ok,
    no_duplicate_cells_ok = no_duplicate_cells_ok,
    expected_partners_unchanged = expected_partners_unchanged,
    outside_allowed_count = outside_allowed_count,
    duplicate_selected_cell_count = duplicate_selected_cell_count
  )
}

get_script16_p_if_available <- function(dataset){
  s16 <- read_semicolon_or_comma(script16_summary_file)
  if(is.null(s16)) return(NA_real_)
  if(!all(c("dataset", "fitted_conditioned_p") %in% names(s16))) return(NA_real_)
  val <- s16 %>%
    filter(dataset == !!dataset) %>%
    slice(1) %>%
    pull(fitted_conditioned_p)
  if(length(val) == 0) NA_real_ else as.numeric(val)
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){

  message("Running script 17 fixed-record allocation null: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tables <- make_dataset_tables(dataset)
  pair_table <- tables$pair_table
  cooc_cells <- tables$cooc_cells
  empirical_cells <- tables$empirical_cells

  observed_total_links <- sum(pair_table$is_realised_link, na.rm = TRUE)

  p_from_script16 <- get_script16_p_if_available(dataset)
  p_hat <- if(is.finite(p_from_script16)) {
    p_from_script16
  } else {
    fit_conditioned_p(pair_table, observed_total_links)
  }

  expected_total_links <- expected_total_conditioned_links(pair_table, p_hat)
  pair_expected <- add_fixed_expectations(pair_table, p_hat)

  expected_partners_reference <- pair_expected %>%
    group_by(consumer) %>%
    summarise(expected_realised_partners = sum(expected_link_probability, na.rm = TRUE),
              .groups = "drop")

  empirical_species <- build_consumer_metrics_from_expected_table(pair_expected)
  empirical_rho <- compute_consumer_correlation(empirical_species)
  empirical_n_included <- nrow(empirical_species)

  consumer_check_base <- cooc_cells %>%
    group_by(consumer) %>%
    summarise(number_possible_interaction_cells = n_distinct(cell_id),
              .groups = "drop") %>%
    left_join(
      empirical_cells %>%
        group_by(consumer) %>%
        summarise(observed_total_interaction_records_M_alpha = n_distinct(cell_id),
                  .groups = "drop"),
      by = "consumer"
    ) %>%
    mutate(observed_total_interaction_records_M_alpha =
             replace_na(observed_total_interaction_records_M_alpha, 0L)) %>%
    left_join(
      pair_table %>%
        group_by(consumer) %>%
        summarise(maximum_allowed_n_ai = max(n_cooccurring_sites, na.rm = TRUE),
                  .groups = "drop"),
      by = "consumer"
    )

  correlations <- vector("numeric", n_allocation_sims)
  n_included <- vector("integer", n_allocation_sims)

  total_match_all_reps_by_consumer <- setNames(rep(TRUE, nrow(consumer_check_base)),
                                               consumer_check_base$consumer)
  max_simulated_K_by_consumer <- setNames(rep(0, nrow(consumer_check_base)),
                                          consumer_check_base$consumer)

  validation_rows <- vector("list", n_allocation_sims)

  for(rep_i in seq_len(n_allocation_sims)){

    sim_object <- simulate_fixed_record_allocation(
      pair_expected = pair_expected,
      cooc_cells = cooc_cells,
      empirical_cells = empirical_cells
    )

    validation <- validate_simulation(
      sim_object = sim_object,
      pair_expected = pair_expected,
      cooc_cells = cooc_cells,
      expected_partners_reference = expected_partners_reference
    )

    validation_rows[[rep_i]] <- data.frame(
      dataset = dataset,
      replicate = rep_i,
      total_records_match_all_consumers = validation$total_records_match_all,
      K_bounds_ok = validation$K_bounds_ok,
      no_assignments_outside_allowed_cells = validation$outside_allowed_ok,
      no_duplicate_selected_cells = validation$no_duplicate_cells_ok,
      expected_partners_unchanged = validation$expected_partners_unchanged,
      outside_allowed_count = validation$outside_allowed_count,
      duplicate_selected_cell_count = validation$duplicate_selected_cell_count
    )

    if(!validation$ok){
      print(validation_rows[[rep_i]])
      stop("Validation failed for dataset ", dataset, " replicate ", rep_i)
    }

    M_check <- validation$M_check
    total_match_all_reps_by_consumer[M_check$consumer] <-
      total_match_all_reps_by_consumer[M_check$consumer] & M_check$total_records_match

    max_K_here <- sim_object$sim_pair %>%
      group_by(consumer) %>%
      summarise(max_K = max(n_interacting_sites, na.rm = TRUE), .groups = "drop")
    max_simulated_K_by_consumer[max_K_here$consumer] <- pmax(
      max_simulated_K_by_consumer[max_K_here$consumer],
      max_K_here$max_K
    )

    sim_species <- build_consumer_metrics_from_expected_table(sim_object$sim_pair)
    correlations[rep_i] <- compute_consumer_correlation(sim_species)
    n_included[rep_i] <- nrow(sim_species)
  }

  validation_table <- bind_rows(validation_rows)

  finite_corr <- correlations[is.finite(correlations)]

  summary_row <- data.frame(
    dataset = dataset,
    fitted_conditioned_p = p_hat,
    observed_total_regional_links = observed_total_links,
    expected_total_regional_links = expected_total_links,
    calibration_error = observed_total_links - expected_total_links,
    empirical_spearman_correlation = empirical_rho,
    null_mean = mean(finite_corr, na.rm = TRUE),
    null_q025 = as.numeric(quantile(finite_corr, 0.025, na.rm = TRUE)),
    null_q975 = as.numeric(quantile(finite_corr, 0.975, na.rm = TRUE)),
    monte_carlo_p_one_sided_lower = mc_p_one_sided_lower(empirical_rho, finite_corr),
    empirical_number_consumers_included = empirical_n_included,
    median_simulated_number_consumers_included = median(n_included, na.rm = TRUE),
    minimum_simulated_number_consumers_included = min(n_included, na.rm = TRUE),
    maximum_simulated_number_consumers_included = max(n_included, na.rm = TRUE),
    n_allocation_sims_requested = n_allocation_sims,
    n_allocation_sims_finite_correlations = length(finite_corr)
  )

  correlation_table <- data.frame(
    dataset = dataset,
    replicate = seq_len(n_allocation_sims),
    simulated_spearman_correlation = correlations,
    number_consumers_included = n_included
  )

  consumer_checks <- consumer_check_base %>%
    mutate(
      simulated_total_interaction_records_matched_M_alpha_every_replicate =
        as.logical(total_match_all_reps_by_consumer[consumer]),
      maximum_simulated_K_ai = as.numeric(max_simulated_K_by_consumer[consumer])
    ) %>%
    select(
      consumer,
      number_possible_interaction_cells,
      observed_total_interaction_records_M_alpha,
      simulated_total_interaction_records_matched_M_alpha_every_replicate,
      maximum_simulated_K_ai,
      maximum_allowed_n_ai
    ) %>%
    mutate(dataset = dataset, .before = 1)

  write.csv2(correlation_table,
             file.path(out_dir, paste0(dataset, "_17_fixed_record_allocation_correlations.csv")),
             row.names = FALSE)
  write.csv2(consumer_checks,
             file.path(out_dir, paste0(dataset, "_17_fixed_record_allocation_checks.csv")),
             row.names = FALSE)
  write.csv2(validation_table,
             file.path(out_dir, paste0(dataset, "_17_validation_by_replicate.csv")),
             row.names = FALSE)

  message("Validation summary for ", dataset, ": ",
          "records match = ", all(validation_table$total_records_match_all_consumers), "; ",
          "K bounds = ", all(validation_table$K_bounds_ok), "; ",
          "inside allowed cells = ", all(validation_table$no_assignments_outside_allowed_cells), "; ",
          "no duplicate cells = ", all(validation_table$no_duplicate_selected_cells), "; ",
          "expected partners unchanged = ", all(validation_table$expected_partners_unchanged))

  list(
    summary_row = summary_row,
    correlation_table = correlation_table,
    consumer_checks = consumer_checks,
    validation_table = validation_table
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- lapply(all_dataset_names, run_one_dataset)
names(all_outputs) <- all_dataset_names

summary_all <- bind_rows(lapply(all_outputs, `[[`, "summary_row"))
correlations_all <- bind_rows(lapply(all_outputs, `[[`, "correlation_table"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "consumer_checks"))
validation_all <- bind_rows(lapply(all_outputs, `[[`, "validation_table"))

## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(
  summary_all,
  file.path(combined_out, "17_consumer_fixed_record_allocation_summary.csv"),
  row.names = FALSE
)

write.csv2(
  correlations_all,
  file.path(combined_out, "17_consumer_fixed_record_allocation_correlations.csv"),
  row.names = FALSE
)

write.csv2(
  checks_all,
  file.path(combined_out, "17_consumer_fixed_record_allocation_checks.csv"),
  row.names = FALSE
)

write.csv2(
  validation_all,
  file.path(combined_out, "17_consumer_fixed_record_allocation_validation_by_replicate.csv"),
  row.names = FALSE
)

## ---------------------------
## Console summaries
## ---------------------------

message("\nScript 17 validation summary:\n")
print(
  validation_all %>%
    group_by(dataset) %>%
    summarise(
      all_consumer_totals_match_M_alpha = all(total_records_match_all_consumers),
      all_K_values_within_bounds = all(K_bounds_ok),
      all_assignments_inside_allowed_cells = all(no_assignments_outside_allowed_cells),
      no_duplicate_cells_in_any_replicate = all(no_duplicate_selected_cells),
      expected_partners_unchanged_in_all_replicates = all(expected_partners_unchanged),
      .groups = "drop"
    ),
  row.names = FALSE
)

message("\nScript 17 summary table:\n")
print(
  summary_all %>%
    select(dataset,
           fitted_conditioned_p,
           empirical_spearman_correlation,
           null_mean,
           null_q025,
           null_q975,
           monte_carlo_p_one_sided_lower,
           empirical_number_consumers_included,
           median_simulated_number_consumers_included,
           minimum_simulated_number_consumers_included,
           maximum_simulated_number_consumers_included),
  row.names = FALSE
)

## ---------------------------
## Plots
## ---------------------------

plot_summary <- summary_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p1 <- ggplot(plot_summary, aes(x = dataset)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_linerange(
    aes(ymin = null_q025,
        ymax = null_q975,
        colour = "Consumer-specific random allocation"),
    linewidth = 2.2,
    alpha = 0.85
  ) +
  geom_point(
    aes(y = empirical_spearman_correlation,
        colour = "Empirical data"),
    size = 2.7
  ) +
  scale_colour_manual(
    name = "",
    values = c(
      "Consumer-specific random allocation" = "#56B4E9",
      "Empirical data" = "#D55E00"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  xlab("") +
  ylab("Spearman correlation")

ggsave(
  file.path(combined_out, "17_consumer_fixed_record_allocation_null_test.png"),
  p1,
  width = 10,
  height = 5.5,
  dpi = 300
)

plot_distributions <- correlations_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names)) %>%
  left_join(
    summary_all %>% select(dataset, empirical_spearman_correlation),
    by = "dataset"
  )

p2 <- ggplot(plot_distributions,
             aes(x = simulated_spearman_correlation)) +
  geom_histogram(aes(fill = "Consumer-specific random allocation"),
                 bins = 30,
                 alpha = 0.75,
                 colour = "white",
                 na.rm = TRUE) +
  geom_vline(aes(xintercept = empirical_spearman_correlation,
                 colour = "Empirical data"),
             linewidth = 0.8) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_fill_manual(
    name = "",
    values = c("Consumer-specific random allocation" = "#56B4E9")
  ) +
  scale_colour_manual(
    name = "",
    values = c("Empirical data" = "#D55E00")
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank()
  ) +
  xlab("Spearman correlation") +
  ylab("Simulation count")

ggsave(
  file.path(combined_out, "17_consumer_fixed_record_allocation_null_distributions.png"),
  p2,
  width = 14,
  height = 7,
  dpi = 300
)

## ---------------------------
## Interpretation notes
## ---------------------------

notes <- c(
  "Script 17: consumer fixed-record allocation null",
  "",
  "This script follows script 16 and uses the same conditioned original Galiana IR model.",
  "For each consumer, the observed total number of local interaction records M_alpha is fixed.",
  "Each simulation randomly allocates exactly M_alpha records without replacement among that consumer's actual co-occurring consumer-resource-site cells.",
  "The script does not preserve resource totals, site totals, pair totals, or network-wide totals except as implied by fixed consumer M_alpha.",
  "The fitted conditioned p is not refitted in simulations.",
  "The one-sided Monte Carlo p-value is P(null correlation <= empirical correlation), testing whether the empirical relationship is more negative than this null expectation.",
  "No ecological mechanism or causal interpretation is inferred here."
)

writeLines(notes, con = file.path(combined_out, "17_interpretation_notes.txt"))

message("Finished script 17 consumer fixed-record allocation null.")
