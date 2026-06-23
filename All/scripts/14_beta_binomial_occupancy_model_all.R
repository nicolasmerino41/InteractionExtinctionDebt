## Script: All/scripts/14_beta_binomial_occupancy_model_all.R
##
## Purpose:
## Test whether an exchangeable beta-binomial occupancy model improves
## over the calibrated homogeneous co-occurrence-to-interaction model.
##
## This is the SHORT version:
##   - no relative-error extra figure
##   - no model-generated envelope figure
##   - no posterior-predictive simulations
##
## Important fixes in this version:
##   - Avoid vector-collapsing `ifelse(mbb$success, vector, NA)` bugs.
##   - The MBB repeatability curve is now built from exact
##   beta-binomial PMF conditional means:
##
##     E[K / n | K > 0, n]
##
##   not from any previously misconstructed plotting variable.
##
## Uses original 10 Galiana datasets only.
## Source:
##   All/scripts/00_dataset_loaders_and_helpers_all.R
##
## Uses write.csv2().
## Parallelised across datasets.
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
n_site_reps <- 100
min_sites_retained <- 2

log_kappa_bounds <- c(log(1e-4), log(1e8))
mu_bounds <- c(1e-8, 1 - 1e-8)

dirs <- make_output_dirs("script14_beta_binomial_occupancy_model")
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script assumes the original 10 Galiana datasets only.")
}


## ---------------------------
## Helpers
## ---------------------------

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0) return(NA_real_)
  (observed - expected) / expected
}

auc_trapezoid <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 2) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

log_betabinom_pmf <- function(k, n, alpha, beta){
  lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta)
}

plink_M0 <- function(n, p){
  1 - (1 - p)^n
}

plink_MBB <- function(n, alpha, beta){
  1 - exp(lbeta(alpha, beta + n) - lbeta(alpha, beta))
}

expected_links_M0 <- function(n_vec, p){
  sum(plink_M0(n_vec, p), na.rm = TRUE)
}

expected_links_MBB <- function(n_vec, alpha, beta){
  sum(plink_MBB(n_vec, alpha, beta), na.rm = TRUE)
}

fit_p_M0 <- function(n_vec, observed_links){
  if(length(n_vec) == 0 || is.na(observed_links)) return(NA_real_)
  if(observed_links <= 0) return(0)
  if(observed_links >= length(n_vec)) return(1)
  
  f <- function(p){
    expected_links_M0(n_vec, p) - observed_links
  }
  
  tryCatch(
    uniroot(f, interval = c(0, 1), tol = 1e-10)$root,
    error = function(e) NA_real_
  )
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

fit_MBB <- function(n_vec, k_vec, observed_links){
  
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
    return(list(success = FALSE, message = fit$message,
                mu = NA_real_, kappa = NA_real_, alpha = NA_real_, beta = NA_real_,
                rho = NA_real_, logLik = NA_real_, boundary_warning = NA_character_))
  }
  
  log_kappa <- fit$minimum
  kappa <- exp(log_kappa)
  mu <- solve_mu_for_kappa(n_vec, observed_links, kappa)
  
  if(is.na(mu)){
    return(list(success = FALSE, message = "mu calibration failed at optimum kappa",
                mu = NA_real_, kappa = kappa, alpha = NA_real_, beta = NA_real_,
                rho = NA_real_, logLik = NA_real_, boundary_warning = NA_character_))
  }
  
  alpha <- mu * kappa
  beta <- (1 - mu) * kappa
  ll <- sum(log_betabinom_pmf(k_vec, n_vec, alpha, beta), na.rm = TRUE)
  rho <- 1 / (kappa + 1)
  
  boundary_warning <- NA_character_
  if(abs(log_kappa - log_kappa_bounds[1]) < 1e-3){
    boundary_warning <- "kappa_at_lower_bound_high_overdispersion"
  }
  if(abs(log_kappa - log_kappa_bounds[2]) < 1e-3){
    boundary_warning <- "kappa_at_upper_bound_near_binomial"
  }
  
  list(success = TRUE, message = "", mu = mu, kappa = kappa,
       alpha = alpha, beta = beta, rho = rho, logLik = ll,
       boundary_warning = boundary_warning)
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
  
  pair_table <- cooc_counts %>%
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
  
  list(pair_table = pair_table, cooc = cooc, interactions = ints)
}

make_site_subsets_local <- function(sites){
  sites <- sort(unique(as.character(sites)))
  n_sites <- length(sites)
  subset_list <- list()
  subset_index <- list()
  counter <- 1
  
  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_sites * (1 - removal)))
    n_keep <- min(n_keep, n_sites)
    reps_here <- ifelse(removal == 0, 1, n_site_reps)
    
    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) sites else sample(sites, n_keep, replace = FALSE)
      
      subset_list[[counter]] <- sort(keep)
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        actual_removal_fraction = 1 - length(keep) / n_sites,
        site_rep = r,
        n_sites_kept = length(keep)
      )
      counter <- counter + 1
    }
  }
  
  list(subsets = subset_list, index = bind_rows(subset_index))
}

## Validated exact conditional repeatability functions.
expected_conditional_repeatability_M0_exact <- function(n, p){
  k <- 0:n
  pk <- dbinom(k, size = n, prob = p)
  p_link <- 1 - pk[k == 0]
  if(is.na(p_link) || p_link <= 0) return(NA_real_)
  sum((k[k > 0] / n) * pk[k > 0]) / p_link
}

expected_conditional_repeatability_MBB_exact <- function(n, alpha, beta){
  k <- 0:n
  log_pk <- log_betabinom_pmf(k, n, alpha, beta)
  pk <- exp(log_pk)
  p_link <- 1 - pk[k == 0]
  if(is.na(p_link) || p_link <= 0) return(NA_real_)
  sum((k[k > 0] / n) * pk[k > 0]) / p_link
}

prob_one_site_link_M0 <- function(n, p){
  dbinom(1, size = n, prob = p)
}

prob_one_site_link_MBB <- function(n, alpha, beta){
  exp(log_betabinom_pmf(1, n, alpha, beta))
}

occupancy_distribution_rows <- function(dataset, pair_table, params){
  n_values <- sort(unique(pair_table$n_cooccurring_sites))
  
  bind_rows(lapply(n_values, function(n){
    obs_tab <- pair_table %>%
      filter(n_cooccurring_sites == n) %>%
      count(k = n_interacting_sites, name = "observed_frequency")
    
    all_k <- data.frame(k = 0:n)
    
    all_k %>%
      left_join(obs_tab, by = "k") %>%
      mutate(
        observed_frequency = replace_na(observed_frequency, 0L),
        dataset = dataset,
        n_cooccurring_sites = n,
        n_pairs_with_this_n = sum(observed_frequency),
        expected_frequency_M0 =
          n_pairs_with_this_n * dbinom(k, size = n, prob = params$p_M0),
        expected_frequency_MBB =
          if (isTRUE(params$MBB_success)) {
            n_pairs_with_this_n * exp(log_betabinom_pmf(k, n, params$alpha_MBB, params$beta_MBB))
          } else {
            rep(NA_real_, length(k))
          }
      )
  }))
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
  
  message("Running beta-binomial occupancy model: ", dataset)
  
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  site_tables <- get_dataset_site_tables(dataset)
  built <- make_pair_table(dataset, site_tables)
  
  pair_table <- built$pair_table
  cooc <- built$cooc
  interactions <- built$interactions
  
  n_vec <- pair_table$n_cooccurring_sites
  k_vec <- pair_table$n_interacting_sites
  observed_full_links <- sum(pair_table$is_realised_link)
  
  ## Fit M0
  p_M0 <- fit_p_M0(n_vec, observed_full_links)
  expected_full_M0 <- expected_links_M0(n_vec, p_M0)
  ll_M0 <- sum(dbinom(k_vec, size = n_vec, prob = p_M0, log = TRUE), na.rm = TRUE)
  
  ## Fit MBB
  mbb <- fit_MBB(n_vec, k_vec, observed_full_links)
  
  expected_full_MBB <- if(mbb$success){
    expected_links_MBB(n_vec, mbb$alpha, mbb$beta)
  } else {
    NA_real_
  }
  
  ## Full-network repeatability checks.
  ## These use exact PMF conditional means, not the old plotting shortcut.
  realised <- pair_table %>% filter(is_realised_link)
  obs_mean_rep <- mean(realised$empirical_repeatability, na.rm = TRUE)
  obs_prop_one <- mean(realised$n_interacting_sites == 1, na.rm = TRUE)
  
  exp_mean_rep_M0_by_pair <- sapply(
    n_vec,
    expected_conditional_repeatability_M0_exact,
    p = p_M0
  )
  
  exp_mean_rep_M0 <- weighted.mean(
    exp_mean_rep_M0_by_pair,
    w = plink_M0(n_vec, p_M0),
    na.rm = TRUE
  )
  
  exp_prop_one_M0 <- sum(prob_one_site_link_M0(n_vec, p_M0), na.rm = TRUE) /
    sum(plink_M0(n_vec, p_M0), na.rm = TRUE)
  
  exp_mean_rep_MBB <- NA_real_
  exp_prop_one_MBB <- NA_real_
  
  if(mbb$success){
    exp_mean_rep_MBB_by_pair <- sapply(
      n_vec,
      expected_conditional_repeatability_MBB_exact,
      alpha = mbb$alpha,
      beta = mbb$beta
    )
    
    exp_mean_rep_MBB <- weighted.mean(
      exp_mean_rep_MBB_by_pair,
      w = plink_MBB(n_vec, mbb$alpha, mbb$beta),
      na.rm = TRUE
    )
    
    exp_prop_one_MBB <- sum(prob_one_site_link_MBB(n_vec, mbb$alpha, mbb$beta), na.rm = TRUE) /
      sum(plink_MBB(n_vec, mbb$alpha, mbb$beta), na.rm = TRUE)
  }
  
  fit_summary <- data.frame(
    dataset = dataset,
    observed_full_links = observed_full_links,
    expected_full_links_M0 = expected_full_M0,
    calibration_error_M0 = observed_full_links - expected_full_M0,
    expected_full_links_MBB = expected_full_MBB,
    calibration_error_MBB = observed_full_links - expected_full_MBB,
    logLik_M0 = ll_M0,
    logLik_MBB = mbb$logLik,
    delta_logLik_MBB_minus_M0 = mbb$logLik - ll_M0,
    p_M0 = p_M0,
    mu_MBB = mbb$mu,
    kappa_MBB = mbb$kappa,
    rho_MBB = mbb$rho,
    alpha_MBB = mbb$alpha,
    beta_MBB = mbb$beta,
    MBB_success = mbb$success,
    MBB_message = mbb$message,
    MBB_boundary_warning = mbb$boundary_warning,
    observed_mean_repeatability_realised_links = obs_mean_rep,
    expected_mean_repeatability_realised_links_M0 = exp_mean_rep_M0,
    expected_mean_repeatability_realised_links_MBB = exp_mean_rep_MBB,
    observed_proportion_one_site_realised_links = obs_prop_one,
    expected_proportion_one_site_realised_links_M0 = exp_prop_one_M0,
    expected_proportion_one_site_realised_links_MBB = exp_prop_one_MBB,
    n_cooccurrence_pairs = nrow(pair_table),
    n_site_level_cooccurrences = nrow(cooc),
    n_site_level_interactions = nrow(interactions),
    repeatability_calculation_note =
      "MBB repeatability uses exact PMF E[K/n | K > 0, n], validated against script 14a."
  )
  
  ## Repeatability by exact n.
  repeatability_by_n <- pair_table %>%
    filter(is_realised_link) %>%
    group_by(dataset, n_cooccurring_sites) %>%
    summarise(
      n_empirical_realised_links = n(),
      empirical_mean_repeatability = mean(empirical_repeatability, na.rm = TRUE),
      empirical_median_repeatability = median(empirical_repeatability, na.rm = TRUE),
      .groups = "drop"
    )
  
  all_n <- data.frame(
    dataset = dataset,
    n_cooccurring_sites = sort(unique(n_vec))
  )
  
  repeatability_by_n <- all_n %>%
    left_join(repeatability_by_n,
              by = c("dataset", "n_cooccurring_sites")) %>%
    mutate(
      n_empirical_realised_links = replace_na(n_empirical_realised_links, 0L),
      M0_expected_conditional_mean_repeatability =
        sapply(n_cooccurring_sites,
               expected_conditional_repeatability_M0_exact,
               p = p_M0),
      MBB_expected_conditional_mean_repeatability =
        if (isTRUE(mbb$success)) {
          sapply(n_cooccurring_sites,
                 expected_conditional_repeatability_MBB_exact,
                 alpha = mbb$alpha,
                 beta = mbb$beta)
        } else {
          rep(NA_real_, length(n_cooccurring_sites))
        },
      MBB_expected_conditional_mean_repeatability_source =
        ifelse(mbb$success, "exact_beta_binomial_pmf_validated_14a", NA_character_),
      M0_expected_probability_realised_link = plink_M0(n_cooccurring_sites, p_M0),
      MBB_expected_probability_realised_link =
        if (isTRUE(mbb$success)) {
          plink_MBB(n_cooccurring_sites, mbb$alpha, mbb$beta)
        } else {
          rep(NA_real_, length(n_cooccurring_sites))
        }
    )
  
  occupancy_distribution <- occupancy_distribution_rows(
    dataset = dataset,
    pair_table = pair_table,
    params = list(
      p_M0 = p_M0,
      alpha_MBB = mbb$alpha,
      beta_MBB = mbb$beta,
      MBB_success = mbb$success
    )
  )
  
  ## Site-removal comparison
  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets_local(all_sites)
  
  removal_rows <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    sites_keep <- subset_object$subsets[[i]]
    subset_row <- subset_object$index[i, ]
    
    cooc_retained_counts <- cooc %>%
      filter(site %in% sites_keep) %>%
      group_by(pair_id) %>%
      summarise(n_retained = n_distinct(site), .groups = "drop")
    
    empirical_retained_links <- interactions %>%
      filter(site %in% sites_keep) %>%
      distinct(pair_id) %>%
      nrow()
    
    n_ret <- cooc_retained_counts$n_retained
    
    expected_M0 <- expected_links_M0(n_ret, p_M0)
    expected_MBB <- if(mbb$success){
      expected_links_MBB(n_ret, mbb$alpha, mbb$beta)
    } else NA_real_
    
    bind_rows(
      data.frame(
        dataset = dataset,
        model_type = "M0_homogeneous_binomial",
        subset_id = subset_row$subset_id,
        removal_fraction = subset_row$removal_fraction,
        actual_removal_fraction = subset_row$actual_removal_fraction,
        site_rep = subset_row$site_rep,
        n_sites_kept = subset_row$n_sites_kept,
        empirical_retained_links = empirical_retained_links,
        expected_retained_links = expected_M0,
        error_observed_minus_expected = empirical_retained_links - expected_M0,
        relative_error = safe_relative_difference(empirical_retained_links, expected_M0),
        empirical_retained_link_proportion = empirical_retained_links / observed_full_links,
        expected_retained_link_proportion = expected_M0 / expected_full_M0,
        absolute_error_retained_link_proportion =
          abs(empirical_retained_links / observed_full_links - expected_M0 / expected_full_M0)
      ),
      data.frame(
        dataset = dataset,
        model_type = "MBB_beta_binomial",
        subset_id = subset_row$subset_id,
        removal_fraction = subset_row$removal_fraction,
        actual_removal_fraction = subset_row$actual_removal_fraction,
        site_rep = subset_row$site_rep,
        n_sites_kept = subset_row$n_sites_kept,
        empirical_retained_links = empirical_retained_links,
        expected_retained_links = expected_MBB,
        error_observed_minus_expected = empirical_retained_links - expected_MBB,
        relative_error = safe_relative_difference(empirical_retained_links, expected_MBB),
        empirical_retained_link_proportion = empirical_retained_links / observed_full_links,
        expected_retained_link_proportion = expected_MBB / expected_full_MBB,
        absolute_error_retained_link_proportion =
          abs(empirical_retained_links / observed_full_links - expected_MBB / expected_full_MBB)
      )
    )
  }))
  
  removal_summary <- removal_rows %>%
    group_by(dataset, model_type, removal_fraction, actual_removal_fraction) %>%
    summarise(
      mean_empirical_retained_links = mean(empirical_retained_links, na.rm = TRUE),
      mean_expected_retained_links = mean(expected_retained_links, na.rm = TRUE),
      mean_error = mean(error_observed_minus_expected, na.rm = TRUE),
      mean_relative_error = mean(relative_error, na.rm = TRUE),
      RMSE = sqrt(mean(error_observed_minus_expected^2, na.rm = TRUE)),
      MAE = mean(abs(error_observed_minus_expected), na.rm = TRUE),
      mean_empirical_retained_link_proportion = mean(empirical_retained_link_proportion, na.rm = TRUE),
      mean_expected_retained_link_proportion = mean(expected_retained_link_proportion, na.rm = TRUE),
      mean_absolute_error_retained_link_proportion =
        mean(absolute_error_retained_link_proportion, na.rm = TRUE),
      n_subsets = n(),
      .groups = "drop"
    )
  
  performance <- removal_rows %>%
    group_by(dataset, model_type) %>%
    summarise(
      RMSE_retained_link_proportion =
        sqrt(mean((empirical_retained_link_proportion - expected_retained_link_proportion)^2, na.rm = TRUE)),
      MAE_retained_link_proportion =
        mean(abs(empirical_retained_link_proportion - expected_retained_link_proportion), na.rm = TRUE),
      signed_bias_retained_link_proportion =
        mean(empirical_retained_link_proportion - expected_retained_link_proportion, na.rm = TRUE),
      absolute_error_AUC =
        auc_trapezoid(actual_removal_fraction,
                      abs(empirical_retained_link_proportion - expected_retained_link_proportion)),
      .groups = "drop"
    ) %>%
    group_by(dataset) %>%
    mutate(
      RMSE_M0 = RMSE_retained_link_proportion[model_type == "M0_homogeneous_binomial"][1],
      MAE_M0 = MAE_retained_link_proportion[model_type == "M0_homogeneous_binomial"][1],
      AUC_M0 = absolute_error_AUC[model_type == "M0_homogeneous_binomial"][1],
      RMSE_relative_to_M0 = RMSE_retained_link_proportion / RMSE_M0,
      MAE_relative_to_M0 = MAE_retained_link_proportion / MAE_M0,
      AUC_relative_to_M0 = absolute_error_AUC / AUC_M0,
      RMSE_improvement_relative_to_M0 = 1 - RMSE_relative_to_M0,
      MAE_improvement_relative_to_M0 = 1 - MAE_relative_to_M0,
      AUC_improvement_relative_to_M0 = 1 - AUC_relative_to_M0
    ) %>%
    ungroup()
  
  write.csv2(pair_table, file.path(out_dir, paste0(dataset, "_14_pair_occupancy_table.csv")), row.names = FALSE)
  write.csv2(removal_summary, file.path(out_dir, paste0(dataset, "_14_site_removal_summary.csv")), row.names = FALSE)
  
  list(
    fit_summary = fit_summary,
    pair_table = pair_table,
    repeatability_by_n = repeatability_by_n,
    occupancy_distribution = occupancy_distribution,
    removal_rows = removal_rows,
    removal_summary = removal_summary,
    performance = performance
  )
}


## ---------------------------
## Run in parallel
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)

names(all_outputs) <- all_dataset_names

fit_summary_all <- bind_rows(lapply(all_outputs, `[[`, "fit_summary"))
pair_table_all <- bind_rows(lapply(all_outputs, `[[`, "pair_table"))
repeatability_by_n_all <- bind_rows(lapply(all_outputs, `[[`, "repeatability_by_n"))
occupancy_distribution_all <- bind_rows(lapply(all_outputs, `[[`, "occupancy_distribution"))
removal_rows_all <- bind_rows(lapply(all_outputs, `[[`, "removal_rows"))
removal_summary_all <- bind_rows(lapply(all_outputs, `[[`, "removal_summary"))
performance_all <- bind_rows(lapply(all_outputs, `[[`, "performance"))


## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(
  fit_summary_all,
  file.path(combined_out, "14_model_parameter_and_fit_summary.csv"),
  row.names = FALSE
)

write.csv2(
  pair_table_all,
  file.path(combined_out, "14_pair_occupancy_table_combined.csv"),
  row.names = FALSE
)

write.csv2(
  removal_rows_all,
  file.path(combined_out, "14_site_removal_model_comparison_by_subset.csv"),
  row.names = FALSE
)

write.csv2(
  removal_summary_all,
  file.path(combined_out, "14_site_removal_model_comparison_summary.csv"),
  row.names = FALSE
)

write.csv2(
  performance_all,
  file.path(combined_out, "14_site_removal_performance_summary.csv"),
  row.names = FALSE
)

write.csv2(
  repeatability_by_n_all,
  file.path(combined_out, "14_repeatability_by_ncooccurrence.csv"),
  row.names = FALSE
)

write.csv2(
  occupancy_distribution_all,
  file.path(combined_out, "14_occupancy_count_distribution.csv"),
  row.names = FALSE
)


## ---------------------------
## Combined plots
## ---------------------------

model_colours <- c(
  "empirical" = "black",
  "M0_homogeneous_binomial" = "orange",
  "MBB_beta_binomial" = "blue"
)

## Figure 1: repeatability by co-occurrence frequency
repeat_plot <- repeatability_by_n_all %>%
  select(dataset, n_cooccurring_sites,
         n_empirical_realised_links,
         empirical_mean_repeatability,
         M0_expected_conditional_mean_repeatability,
         MBB_expected_conditional_mean_repeatability) %>%
  pivot_longer(
    cols = c(empirical_mean_repeatability,
             M0_expected_conditional_mean_repeatability,
             MBB_expected_conditional_mean_repeatability),
    names_to = "curve",
    values_to = "mean_repeatability"
  ) %>%
  mutate(
    model_type = case_when(
      curve == "empirical_mean_repeatability" ~ "empirical",
      curve == "M0_expected_conditional_mean_repeatability" ~ "M0_homogeneous_binomial",
      curve == "MBB_expected_conditional_mean_repeatability" ~ "MBB_beta_binomial"
    ),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p1 <- ggplot(repeat_plot,
             aes(x = n_cooccurring_sites,
                 y = mean_repeatability,
                 colour = model_type)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  geom_point(aes(size = ifelse(model_type == "empirical", n_empirical_realised_links, NA_real_)),
             alpha = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = model_colours) +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  facet_wrap(~ dataset, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Number of co-occurring sites") +
  ylab("Mean repeatability among realised links") +
  ggtitle("Repeatability by co-occurrence frequency",
          subtitle = "MBB curve uses exact beta-binomial PMF conditional means validated in script 14a")

ggsave(
  file.path(combined_out, "14_repeatability_by_ncooccurrence_models.png"),
  p1,
  width = 14,
  height = 7,
  dpi = 300
)

## Figure 2: site-removal prediction error
error_plot <- removal_summary_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p2 <- ggplot(error_plot,
             aes(x = actual_removal_fraction,
                 y = mean_error,
                 colour = model_type)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = model_colours[c("M0_homogeneous_binomial", "MBB_beta_binomial")]) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Actual fraction of sites removed") +
  ylab("Observed - expected retained links") +
  ggtitle("Site-removal prediction error",
          subtitle = "Both models are calibrated to match full-network link number at zero removal")

ggsave(
  file.path(combined_out, "14_site_removal_prediction_error_models.png"),
  p2,
  width = 14,
  height = 7,
  dpi = 300
)

## Figure 3: cross-dataset model improvement
improvement_plot <- performance_all %>%
  select(dataset, model_type, RMSE_relative_to_M0) %>%
  left_join(fit_summary_all %>% select(dataset, rho_MBB), by = "dataset") %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p3a <- ggplot(improvement_plot,
              aes(x = dataset,
                  y = RMSE_relative_to_M0,
                  colour = model_type)) +
  geom_hline(yintercept = 1, colour = "grey40") +
  geom_point(size = 2.5) +
  scale_colour_manual(values = model_colours[c("M0_homogeneous_binomial", "MBB_beta_binomial")]) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Retention-proportion RMSE relative to M0") +
  ggtitle("A. Retention prediction performance")

p3b <- fit_summary_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names)) %>%
  ggplot(aes(x = dataset, y = rho_MBB)) +
  geom_point(size = 2.5, colour = "blue") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Fitted beta-binomial overdispersion, rho") +
  ggtitle("B. Fitted exchangeable occupancy overdispersion")

if(requireNamespace("patchwork", quietly = TRUE)){
  p3 <- p3a / p3b
  ggsave(
    file.path(combined_out, "14_beta_binomial_model_improvement_summary.png"),
    p3,
    width = 12,
    height = 8,
    dpi = 300
  )
} else {
  ggsave(
    file.path(combined_out, "14_beta_binomial_model_improvement_summary.png"),
    p3a,
    width = 12,
    height = 5,
    dpi = 300
  )
}

## ---------------------------
## Interpretation notes
## ---------------------------

notes <- c(
  "Script 14: beta-binomial occupancy model",
  "",
  "1. MBB is a parsimonious exchangeable heterogeneity model, not a model that identifies which pair has high interaction propensity.",
  "2. Both M0 and MBB are calibrated to match the observed full metaweb link count. Differences after removal therefore reflect predicted occupancy structure, not baseline link-number mismatch.",
  "3. No parameters are refitted after site removal.",
  "4. A better beta-binomial fit means that a single global p is insufficient to describe the distribution of interaction occurrence across operational co-occurrence opportunities.",
  "5. A reduction in removal error means that pair-to-pair heterogeneity in local interaction propensity explains part of the empirical link-retention excess without needing to predict specific pair identities.",
  "6. Site removal remains a structural/sampling diagnostic, not a habitat-loss forecast.",
  "7. The model uses the raw/unconditioned formulation only; no consumer-level conditioned formula is used.",
  "8. The MBB repeatability curve in this rebuilt version uses exact beta-binomial PMF conditional means E[K/n | K > 0, n], matching the validation logic from script 14a."
)

writeLines(
  notes,
  con = file.path(combined_out, "14_interpretation_notes.txt")
)

future::plan(future::sequential)

message("Finished rebuilt short script 14 beta-binomial occupancy model.")