## ------------------------------------------------------------
## Script: All/scripts/14a_validate_beta_binomial_repeatability.R
##
## Purpose:
## Validation/debugging script for the beta-binomial repeatability
## calculation used in:
##
##   All/scripts/14_beta_binomial_occupancy_model_all.R
##
## Main validation target:
##
##   E[K / n | K > 0, n]
##
## Compare:
##   A. analytical shortcut: mu / P_link(n)
##   B. exact beta-binomial PMF calculation
##   C. Monte Carlo simulation from fitted Beta-Binomial model
##
## Uses write.csv2().
## Does not modify script 14.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

## ---------------------------
## Settings
## ---------------------------

candidate_n <- c(1, 2, 3, 5, 10, 20, 30, 40)
n_pair_sims <- 200000
tolerance <- 0.005

dirs <- make_output_dirs("script14a_validate_beta_binomial_repeatability")
sep_out <- dirs$separated
combined_out <- dirs$combined

fit_file <- file.path(
  combined_out,
  "14_model_parameter_and_fit_summary.csv"
)

repeatability_file <- file.path(
  combined_out,
  "14_repeatability_by_ncooccurrence.csv"
)


## ---------------------------
## Helpers
## ---------------------------

read_semicolon_or_comma <- function(file){
  if(!file.exists(file)){
    return(NULL)
  }
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1){
    return(x)
  }
  read.csv(file, stringsAsFactors = FALSE)
}

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

log_betabinom_pmf <- function(k, n, alpha, beta){
  lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta)
}

plink_MBB <- function(n, alpha, beta){
  1 - exp(lbeta(alpha, beta + n) - lbeta(alpha, beta))
}

repeatability_shortcut <- function(n, mu, alpha, beta){
  pl <- plink_MBB(n, alpha, beta)
  ifelse(pl > 0, mu / pl, NA_real_)
}

repeatability_exact_pmf <- function(n, alpha, beta){

  k <- 0:n
  log_p <- log_betabinom_pmf(k, n, alpha, beta)
  p <- exp(log_p)

  p_link_exact <- 1 - p[k == 0]

  mean_repeatability_exact <- sum((k[k > 0] / n) * p[k > 0]) / p_link_exact

  data.frame(
    p_link_exact = p_link_exact,
    mean_repeatability_exact = mean_repeatability_exact
  )
}

repeatability_simulation <- function(n, alpha, beta, n_pair_sims){

  p_ij <- rbeta(n_pair_sims, shape1 = alpha, shape2 = beta)
  k <- rbinom(n_pair_sims, size = n, prob = p_ij)

  realised <- k > 0
  repeatability <- k[realised] / n

  data.frame(
    p_link_simulated = mean(realised),
    mean_repeatability_simulated = mean(repeatability, na.rm = TRUE),
    median_repeatability_simulated = median(repeatability, na.rm = TRUE),
    q025_repeatability_simulated = quantile(repeatability, 0.025, na.rm = TRUE),
    q975_repeatability_simulated = quantile(repeatability, 0.975, na.rm = TRUE),
    fraction_K0_simulated = mean(k == 0),
    fraction_K_equals_n_simulated = mean(k[realised] == n, na.rm = TRUE)
  )
}

beta_mode <- function(alpha, beta){
  if(is.na(alpha) || is.na(beta)) return(NA_real_)
  if(alpha > 1 && beta > 1){
    return((alpha - 1) / (alpha + beta - 2))
  }
  NA_real_
}

make_pair_table <- function(dataset){

  site_tables <- get_dataset_site_tables(dataset)

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
      is_realised_link = n_interacting_sites > 0,
      empirical_repeatability = n_interacting_sites / n_cooccurring_sites
    ) %>%
    select(dataset, consumer, resource, pair_id,
           n_cooccurring_sites, n_interacting_sites,
           is_realised_link, empirical_repeatability)
}

## Minimal refit helpers, used only if script 14 output is unavailable.
log_kappa_bounds <- c(log(1e-4), log(1e8))
mu_bounds <- c(1e-8, 1 - 1e-8)

expected_links_MBB <- function(n_vec, alpha, beta){
  sum(plink_MBB(n_vec, alpha, beta), na.rm = TRUE)
}

solve_mu_for_kappa <- function(n_vec, observed_links, kappa){
  if(length(n_vec) == 0 || is.na(kappa)) return(NA_real_)
  if(observed_links <= 0) return(mu_bounds[1])
  if(observed_links >= length(n_vec)) return(mu_bounds[2])

  f <- function(mu){
    alpha <- mu * kappa
    beta <- (1 - mu) * kappa
    expected_links_MBB(n_vec, alpha, beta) - observed_links
  }

  f_low <- f(mu_bounds[1])
  f_high <- f(mu_bounds[2])

  if(!is.finite(f_low) || !is.finite(f_high) || f_low > 0 || f_high < 0){
    return(NA_real_)
  }

  tryCatch(
    uniroot(f, interval = mu_bounds, tol = 1e-10)$root,
    error = function(e) NA_real_
  )
}

fit_MBB_minimal <- function(n_vec, k_vec, observed_links){

  objective <- function(log_kappa){
    kappa <- exp(log_kappa)
    mu <- solve_mu_for_kappa(n_vec, observed_links, kappa)

    if(is.na(mu)){
      return(1e100)
    }

    alpha <- mu * kappa
    beta <- (1 - mu) * kappa

    ll <- sum(log_betabinom_pmf(k_vec, n_vec, alpha, beta), na.rm = TRUE)

    if(!is.finite(ll)){
      return(1e100)
    }

    -ll
  }

  fit <- tryCatch(
    optimize(objective, interval = log_kappa_bounds),
    error = function(e) e
  )

  if(inherits(fit, "error")){
    return(data.frame(
      alpha_MBB = NA_real_,
      beta_MBB = NA_real_,
      mu_MBB = NA_real_,
      kappa_MBB = NA_real_,
      rho_MBB = NA_real_,
      MBB_success = FALSE,
      MBB_message = fit$message
    ))
  }

  kappa <- exp(fit$minimum)
  mu <- solve_mu_for_kappa(n_vec, observed_links, kappa)

  if(is.na(mu)){
    return(data.frame(
      alpha_MBB = NA_real_,
      beta_MBB = NA_real_,
      mu_MBB = NA_real_,
      kappa_MBB = kappa,
      rho_MBB = NA_real_,
      MBB_success = FALSE,
      MBB_message = "mu calibration failed"
    ))
  }

  data.frame(
    alpha_MBB = mu * kappa,
    beta_MBB = (1 - mu) * kappa,
    mu_MBB = mu,
    kappa_MBB = kappa,
    rho_MBB = 1 / (kappa + 1),
    MBB_success = TRUE,
    MBB_message = ""
  )
}


## ---------------------------
## Load fitted parameters if possible
## ---------------------------

fit_summary <- read_semicolon_or_comma(fit_file)

if(is.null(fit_summary)){
  message("Script 14 fit summary not found. Reproducing only MBB parameter fitting.")

  fit_summary <- bind_rows(lapply(all_dataset_names, function(dataset){
    message("Refitting minimal MBB for ", dataset)
    pair_table <- make_pair_table(dataset)
    n_vec <- pair_table$n_cooccurring_sites
    k_vec <- pair_table$n_interacting_sites
    observed_links <- sum(pair_table$is_realised_link)

    fit <- fit_MBB_minimal(n_vec, k_vec, observed_links)
    data.frame(dataset = dataset, fit)
  }))
} else {
  required_cols <- c("dataset", "alpha_MBB", "beta_MBB", "mu_MBB",
                     "kappa_MBB", "rho_MBB", "MBB_success")
  missing_cols <- setdiff(required_cols, names(fit_summary))

  if(length(missing_cols) > 0){
    stop("Fit summary exists but is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
}

## Only original 10 datasets.
fit_summary <- fit_summary %>%
  filter(dataset %in% all_dataset_names)

if(length(unique(fit_summary$dataset)) != length(all_dataset_names)){
  warning("Not all datasets in all_dataset_names were found in fit summary.")
}


## ---------------------------
## Build selected n values and validation table
## ---------------------------

validation_rows <- list()
shape_rows <- list()

for(dataset in all_dataset_names){

  message("Validating beta-binomial repeatability: ", dataset)

  params <- fit_summary %>%
    filter(dataset == !!dataset) %>%
    slice(1)

  if(nrow(params) == 0 || !isTRUE(as.logical(params$MBB_success))){
    validation_rows[[dataset]] <- data.frame(
      dataset = dataset,
      alpha_MBB = NA_real_,
      beta_MBB = NA_real_,
      mu_MBB = NA_real_,
      kappa_MBB = NA_real_,
      rho_MBB = NA_real_,
      n_cooccurring_sites = NA_integer_,
      p_link_shortcut = NA_real_,
      p_link_exact = NA_real_,
      p_link_simulated = NA_real_,
      mean_repeatability_shortcut = NA_real_,
      mean_repeatability_exact = NA_real_,
      mean_repeatability_simulated = NA_real_,
      median_repeatability_simulated = NA_real_,
      q025_repeatability_simulated = NA_real_,
      q975_repeatability_simulated = NA_real_,
      fraction_K0_simulated = NA_real_,
      fraction_K_equals_n_simulated = NA_real_,
      absolute_difference_shortcut_vs_exact = NA_real_,
      absolute_difference_exact_vs_simulated = NA_real_,
      validation_pass = FALSE
    )
    next
  }

  pair_table <- make_pair_table(dataset)
  max_n <- max(pair_table$n_cooccurring_sites, na.rm = TRUE)
  median_n <- as.integer(round(median(pair_table$n_cooccurring_sites, na.rm = TRUE)))

  n_values <- sort(unique(c(candidate_n, median_n, max_n)))
  n_values <- n_values[n_values >= 1 & n_values <= max_n]

  alpha <- as.numeric(params$alpha_MBB)
  beta <- as.numeric(params$beta_MBB)
  mu <- as.numeric(params$mu_MBB)
  kappa <- as.numeric(params$kappa_MBB)
  rho <- as.numeric(params$rho_MBB)

  ## Parameter-shape table
  shape_rows[[dataset]] <- data.frame(
    dataset = dataset,
    alpha_MBB = alpha,
    beta_MBB = beta,
    mu_MBB = mu,
    kappa_MBB = kappa,
    rho_MBB = rho,
    beta_mode = beta_mode(alpha, beta),
    probability_pij_less_0_01 = pbeta(0.01, shape1 = alpha, shape2 = beta),
    probability_pij_greater_0_5 = 1 - pbeta(0.5, shape1 = alpha, shape2 = beta),
    probability_pij_greater_0_9 = 1 - pbeta(0.9, shape1 = alpha, shape2 = beta),
    q01_pij = qbeta(0.01, shape1 = alpha, shape2 = beta),
    q05_pij = qbeta(0.05, shape1 = alpha, shape2 = beta),
    q50_pij = qbeta(0.50, shape1 = alpha, shape2 = beta),
    q95_pij = qbeta(0.95, shape1 = alpha, shape2 = beta),
    q99_pij = qbeta(0.99, shape1 = alpha, shape2 = beta)
  )

  validation_rows[[dataset]] <- bind_rows(lapply(n_values, function(n){

    p_link_short <- plink_MBB(n, alpha, beta)
    mean_short <- repeatability_shortcut(n, mu, alpha, beta)

    exact <- repeatability_exact_pmf(n, alpha, beta)
    sim <- repeatability_simulation(n, alpha, beta, n_pair_sims)

    diff_short_exact <- abs(mean_short - exact$mean_repeatability_exact)
    diff_exact_sim <- abs(exact$mean_repeatability_exact - sim$mean_repeatability_simulated)

    data.frame(
      dataset = dataset,
      alpha_MBB = alpha,
      beta_MBB = beta,
      mu_MBB = mu,
      kappa_MBB = kappa,
      rho_MBB = rho,
      n_cooccurring_sites = n,
      p_link_shortcut = p_link_short,
      p_link_exact = exact$p_link_exact,
      p_link_simulated = sim$p_link_simulated,
      mean_repeatability_shortcut = mean_short,
      mean_repeatability_exact = exact$mean_repeatability_exact,
      mean_repeatability_simulated = sim$mean_repeatability_simulated,
      median_repeatability_simulated = sim$median_repeatability_simulated,
      q025_repeatability_simulated = sim$q025_repeatability_simulated,
      q975_repeatability_simulated = sim$q975_repeatability_simulated,
      fraction_K0_simulated = sim$fraction_K0_simulated,
      fraction_K_equals_n_simulated = sim$fraction_K_equals_n_simulated,
      absolute_difference_shortcut_vs_exact = diff_short_exact,
      absolute_difference_exact_vs_simulated = diff_exact_sim,
      validation_pass = diff_short_exact < tolerance & diff_exact_sim < tolerance
    )
  }))
}

validation_table <- bind_rows(validation_rows)
shape_table <- bind_rows(shape_rows)


## ---------------------------
## Save outputs
## ---------------------------

write.csv2(
  validation_table,
  file.path(combined_out, "14a_beta_binomial_repeatability_validation.csv"),
  row.names = FALSE
)

write.csv2(
  shape_table,
  file.path(combined_out, "14a_fitted_beta_distribution_summary.csv"),
  row.names = FALSE
)


## ---------------------------
## Compact validation figure
## ---------------------------

plot_data <- validation_table %>%
  filter(!is.na(n_cooccurring_sites)) %>%
  select(dataset,
         n_cooccurring_sites,
         mean_repeatability_shortcut,
         mean_repeatability_exact,
         mean_repeatability_simulated) %>%
  pivot_longer(
    cols = c(mean_repeatability_shortcut,
             mean_repeatability_exact,
             mean_repeatability_simulated),
    names_to = "calculation",
    values_to = "conditional_mean_repeatability"
  ) %>%
  mutate(
    calculation = recode(
      calculation,
      mean_repeatability_shortcut = "Analytical shortcut",
      mean_repeatability_exact = "Exact PMF",
      mean_repeatability_simulated = "Simulation"
    ),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p <- ggplot(plot_data,
            aes(x = n_cooccurring_sites,
                y = conditional_mean_repeatability,
                colour = calculation)) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  facet_wrap(~ dataset, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Number of co-occurring sites") +
  ylab("E[K / n | K > 0, n]") +
  ggtitle("Validation of beta-binomial conditional repeatability",
          subtitle = "Shortcut, exact PMF calculation, and direct simulation should overlap")

ggsave(
  file.path(combined_out, "14a_beta_binomial_repeatability_validation.png"),
  p,
  width = 14,
  height = 7,
  dpi = 300
)


## ---------------------------
## Interpretation note
## ---------------------------

all_pass <- all(validation_table$validation_pass, na.rm = TRUE)
max_short_exact <- max(validation_table$absolute_difference_shortcut_vs_exact, na.rm = TRUE)
max_exact_sim <- max(validation_table$absolute_difference_exact_vs_simulated, na.rm = TRUE)

if(all_pass){
  interpretation <- c(
    "Beta-binomial repeatability validation",
    "",
    "Result: shortcut, exact PMF, and Monte Carlo simulation agree within the requested tolerance.",
    "",
    paste0("Maximum absolute difference between shortcut and exact PMF: ", signif(max_short_exact, 4)),
    paste0("Maximum absolute difference between exact PMF and simulation: ", signif(max_exact_sim, 4)),
    "",
    "Diagnosis:",
    "The near-perfect MBB conditional repeatability is a genuine implication of the fitted beta-binomial model, not a plotting or conditioning error.",
    "",
    "Interpretation:",
    "High fitted overdispersion can produce a latent propensity distribution with many pairs near zero interaction probability and a smaller set of pairs with high interaction probability.",
    "Because repeatability is conditioned on K > 0, realised links are disproportionately drawn from the high-propensity part of the fitted beta distribution.",
    "Those realised links are therefore predicted to occupy many or most of their co-occurring sites."
  )
} else {
  failed <- validation_table %>%
    filter(!validation_pass) %>%
    select(dataset, n_cooccurring_sites,
           absolute_difference_shortcut_vs_exact,
           absolute_difference_exact_vs_simulated)

  fail_file <- file.path(combined_out, "14a_failed_validation_cases.csv")
  write.csv2(failed, fail_file, row.names = FALSE)

  interpretation <- c(
    "Beta-binomial repeatability validation",
    "",
    "Result: at least one dataset/n combination failed the requested tolerance.",
    "",
    paste0("Maximum absolute difference between shortcut and exact PMF: ", signif(max_short_exact, 4)),
    paste0("Maximum absolute difference between exact PMF and simulation: ", signif(max_exact_sim, 4)),
    "",
    "Diagnosis:",
    "The validation failure means the near-perfect MBB repeatability may involve a formula, conditioning, plotting, parameter extraction, or numerical issue.",
    "",
    "Next step:",
    paste0("Inspect failed cases saved in: ", fail_file),
    "If shortcut differs from exact PMF, the analytical formula in script 14 is wrong.",
    "If exact PMF differs from simulation, the simulation or numerical PMF implementation should be checked.",
    "If both match but the plotted script-14 curve differs, the problem is likely plotting or parameter extraction."
  )
}

writeLines(
  interpretation,
  con = file.path(combined_out, "14a_validation_interpretation.txt")
)

message("Finished script 14a beta-binomial repeatability validation.")
