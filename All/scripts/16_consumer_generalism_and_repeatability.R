## ------------------------------------------------------------
## Script: All/scripts/16_consumer_generalism_and_repeatability.R
##
## Purpose:
## Link the original conditioned consumer-level IR model idea
## to pair-level repeatability among realised links.
##
## Question:
## Do consumers with more realised resource partners than expected
## by the conditioned homogeneous co-occurrence-frequency model also
## have realised links that occur more repeatedly across sites than
## expected by the same homogeneous model?
##
## Important limits:
##   - original 10 datasets only
##   - consumers only
##   - no beta-binomial models
##   - no site-removal outcome
##   - descriptive association, not causal interpretation
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

n_null_sims <- 1000
max_conditioning_attempts <- 10000

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script is intended for the original 10 Galiana datasets only.")
}

dirs <- make_output_dirs("script16_consumer_generalism_and_repeatability")
sep_out <- dirs$separated
combined_out <- dirs$combined

## ---------------------------
## Helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

readable_bool <- function(x){
  ifelse(isTRUE(x), "TRUE", "FALSE")
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

mc_p_two_sided <- function(empirical, null_values){
  null_values <- null_values[is.finite(null_values)]
  if(!is.finite(empirical) || length(null_values) == 0) return(NA_real_)
  ## Two-sided randomisation-style p-value using absolute deviation from null mean.
  centre <- mean(null_values)
  (sum(abs(null_values - centre) >= abs(empirical - centre)) + 1) /
    (length(null_values) + 1)
}

make_pair_table <- function(dataset, site_tables){

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource))

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource))

  cooc_counts <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

  int_counts <- ints %>%
    group_by(pair_id) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop")

  cooc_counts %>%
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

  ## In rare edge cases the conditioned expectation may not bracket exactly.
  ## Use optimize fallback rather than silently failing.
  if(f_low * f_high > 0){
    opt <- optimize(function(p) abs(f(p)), interval = c(1e-12, 1 - 1e-12))
    return(opt$minimum)
  }

  uniroot(f, interval = c(1e-12, 1 - 1e-12), tol = 1e-12)$root
}

build_consumer_metrics <- function(pair_table, p){

  N_tab <- consumer_N_table(pair_table)

  pair_expected <- pair_table %>%
    left_join(N_tab, by = "consumer") %>%
    mutate(
      expected_link_probability = plink_conditioned_consumer(
        n_ai = n_cooccurring_sites,
        N_alpha = N_alpha,
        p = p
      ),
      expected_repeatability = expected_repeatability_given_realised(n_cooccurring_sites, p),
      repeatability_excess = observed_repeatability - expected_repeatability
    )

  partner_metrics <- pair_expected %>%
    group_by(dataset, consumer) %>%
    summarise(
      number_cooccurring_resources = first(number_cooccurring_resources),
      observed_realised_partners = sum(is_realised_link, na.rm = TRUE),
      expected_realised_partners = sum(expected_link_probability, na.rm = TRUE),
      `Realised partners above model expectation` =
        observed_realised_partners - expected_realised_partners,
      .groups = "drop"
    )

  repeatability_metrics <- pair_expected %>%
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

simulate_conditioned_K_for_consumer <- function(n_vec, p){
  if(length(n_vec) == 0) return(integer(0))

  for(attempt in seq_len(max_conditioning_attempts)){
    K <- rbinom(length(n_vec), size = n_vec, prob = p)
    if(any(K > 0L)) return(K)
  }

  ## Extremely unlikely fallback: force one realised link, still respecting n.
  K <- rbinom(length(n_vec), size = n_vec, prob = p)
  possible <- which(n_vec > 0)
  chosen <- sample(possible, 1)
  K[chosen] <- max(K[chosen], 1L)
  K
}

simulate_pair_table_conditioned <- function(pair_table, p){

  pair_table %>%
    group_by(consumer) %>%
    group_modify(function(.x, .y){
      K_sim <- simulate_conditioned_K_for_consumer(.x$n_cooccurring_sites, p)
      .x$n_interacting_sites <- K_sim
      .x$is_realised_link <- K_sim > 0L
      .x$observed_repeatability <- ifelse(.x$is_realised_link,
                                          K_sim / .x$n_cooccurring_sites,
                                          NA_real_)
      .x
    }) %>%
    ungroup()
}

run_null_sims <- function(pair_table, p, n_sims){
  replicate(n_sims, {
    sim_pair_table <- simulate_pair_table_conditioned(pair_table, p)
    sim_species <- build_consumer_metrics(sim_pair_table, p)
    spearman_safe(
      sim_species$`Realised partners above model expectation`,
      sim_species$`Extra repeatability of realised links`
    )
  })
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){

  message("Running script 16 consumer generalism-repeatability analysis: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  pair_table <- make_pair_table(dataset, site_tables)

  observed_total_links <- sum(pair_table$is_realised_link, na.rm = TRUE)
  p_hat <- fit_conditioned_p(pair_table, observed_total_links)
  expected_total_links <- expected_total_conditioned_links(pair_table, p_hat)

  species_metrics <- build_consumer_metrics(pair_table, p_hat)

  empirical_rho <- spearman_safe(
    species_metrics$`Realised partners above model expectation`,
    species_metrics$`Extra repeatability of realised links`
  )

  null_rhos <- run_null_sims(pair_table, p_hat, n_null_sims)
  null_rhos_finite <- null_rhos[is.finite(null_rhos)]

  summary_row <- data.frame(
    dataset = dataset,
    fitted_conditioned_p = p_hat,
    observed_total_regional_links = observed_total_links,
    expected_total_regional_links = expected_total_links,
    calibration_error = observed_total_links - expected_total_links,
    number_consumers_total = n_distinct(pair_table$consumer),
    number_consumers_included_repeatability_analysis = nrow(species_metrics),
    empirical_spearman_correlation = empirical_rho,
    null_mean = mean(null_rhos_finite, na.rm = TRUE),
    null_q025 = as.numeric(quantile(null_rhos_finite, 0.025, na.rm = TRUE)),
    null_q975 = as.numeric(quantile(null_rhos_finite, 0.975, na.rm = TRUE)),
    monte_carlo_p_two_sided = mc_p_two_sided(empirical_rho, null_rhos_finite),
    n_null_sims_requested = n_null_sims,
    n_null_sims_finite = length(null_rhos_finite)
  )

  null_table <- data.frame(
    dataset = dataset,
    simulation = seq_along(null_rhos),
    spearman_correlation = null_rhos
  )

  write.csv2(species_metrics,
             file.path(out_dir, paste0(dataset, "_16_consumer_metrics.csv")),
             row.names = FALSE)
  write.csv2(null_table,
             file.path(out_dir, paste0(dataset, "_16_null_spearman_correlations.csv")),
             row.names = FALSE)

  list(
    species_metrics = species_metrics,
    summary_row = summary_row,
    null_table = null_table
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- lapply(all_dataset_names, run_one_dataset)
names(all_outputs) <- all_dataset_names

species_all <- bind_rows(lapply(all_outputs, `[[`, "species_metrics"))
summary_all <- bind_rows(lapply(all_outputs, `[[`, "summary_row"))
null_all <- bind_rows(lapply(all_outputs, `[[`, "null_table"))

## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(
  species_all,
  file.path(combined_out, "16_consumer_generalism_repeatability_species.csv"),
  row.names = FALSE
)

write.csv2(
  summary_all,
  file.path(combined_out, "16_consumer_generalism_repeatability_summary.csv"),
  row.names = FALSE
)

write.csv2(
  null_all,
  file.path(combined_out, "16_consumer_generalism_repeatability_null_correlations.csv"),
  row.names = FALSE
)

## ---------------------------
## Console summary
## ---------------------------

message("\nScript 16 summary table:\n")
print(
  summary_all %>%
    select(dataset,
           fitted_conditioned_p,
           observed_total_regional_links,
           expected_total_regional_links,
           calibration_error,
           number_consumers_total,
           number_consumers_included_repeatability_analysis,
           empirical_spearman_correlation,
           null_mean,
           null_q025,
           null_q975,
           monte_carlo_p_two_sided),
  row.names = FALSE
)

## ---------------------------
## Plots
## ---------------------------

plot_species <- species_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p1 <- ggplot(
  plot_species,
  aes(x = `Realised partners above model expectation`,
      y = `Extra repeatability of realised links`)
) +
  geom_hline(aes(colour = "Zero reference"), yintercept = 0,
             linewidth = 0.35, alpha = 0.7) +
  geom_vline(aes(colour = "Zero reference"), xintercept = 0,
             linewidth = 0.35, alpha = 0.7) +
  geom_point(aes(colour = "Consumers"), alpha = 0.75, size = 1.8) +
  geom_smooth(aes(colour = "Trend"), method = "lm", formula = y ~ x,
              se = TRUE, linewidth = 0.8) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  scale_colour_manual(
    name = "",
    values = c(
      "Consumers" = "grey20",
      "Trend" = "#0072B2",
      "Zero reference" = "grey75"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank()
  ) +
  xlab("Realised partners above model expectation") +
  ylab("Extra repeatability of realised links")

ggsave(
  file.path(combined_out, "16_consumer_generalism_repeatability.png"),
  p1,
  width = 14,
  height = 7,
  dpi = 300
)

plot_summary <- summary_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p2 <- ggplot(plot_summary, aes(x = dataset)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_linerange(
    aes(ymin = null_q025,
        ymax = null_q975,
        colour = "Homogeneous-model simulations"),
    linewidth = 2.2,
    alpha = 0.8
  ) +
  geom_point(
    aes(y = empirical_spearman_correlation,
        colour = "Empirical data"),
    size = 2.6
  ) +
  scale_colour_manual(
    name = "",
    values = c(
      "Homogeneous-model simulations" = "#56B4E9",
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
  file.path(combined_out, "16_consumer_generalism_repeatability_null_test.png"),
  p2,
  width = 10,
  height = 5.5,
  dpi = 300
)

## ---------------------------
## Interpretation notes
## ---------------------------
notes <- c(
  "Script 16: consumer generalism and repeatability",
  "",
  "This script analyses consumers only, matching the direction closest to the original Figure 4 idea.",
  "The fitted model is the original conditioned homogeneous IR formulation.",
  "The x-variable is observed realised partners minus expected realised partners under the conditioned model.",
  "The y-variable is mean observed repeatability minus homogeneous expected repeatability among realised links with n >= 2.",
  "Links with n = 1 are excluded from repeatability because K/n is necessarily 1 among realised links.",
  "The null comparison simulates Binomial(n_ai, p) interaction counts while conditioning each consumer to have at least one realised link.",
  "The fitted p is not refitted in simulations.",
  "No causal interpretation is made here; this is a focused descriptive diagnostic."
)

writeLines(notes, con = file.path(combined_out, "16_interpretation_notes.txt"))

message("Finished script 16 consumer generalism-repeatability analysis.")
