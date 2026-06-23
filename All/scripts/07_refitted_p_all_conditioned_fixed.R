## ------------------------------------------------------------
## Script: All/scripts/07_refitted_p_all_original_style.R
##
## Refit p under site removal using the original Galiana-style
## conditioned expectation and original-style model object.
##
## Key points:
## 1. The original scripts compute:
##      prob_expected = [1 - (1 - p)^n_cooc] /
##                      [1 - (1 - p)^N_alpha]
## 2. They check fit with:
##      sum(expected_int_model$prob_expected)
##      sum(number_interactions_pred$interactions)
## 3. They build expected_int_model after:
##      interaction_proportion <- merge(number_interactions,
##                                      proportion_coocurrence,
##                                      by = "species")
##
## This script mimics that logic for each retained-site subset.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script7_refitted_p")
sep_out <- dirs$separated
combined_out <- dirs$combined

## ------------------------------------------------------------
## Expected links
## ------------------------------------------------------------

expected_raw_from_p <- function(cooc_pairs, p){
  
  if(nrow(cooc_pairs) == 0 || is.na(p)){
    return(NA_real_)
  }
  
  sum(1 - (1 - p)^cooc_pairs$n_cooc, na.rm = TRUE)
}

expected_conditioned_from_p <- function(cooc_pairs, p){
  
  if(nrow(cooc_pairs) == 0 || is.na(p)){
    return(NA_real_)
  }
  
  ## Original N_alpha: sum of all co-occurrence frequencies
  ## for a given consumer/species.
  tmp <- cooc_pairs %>%
    group_by(consumer) %>%
    mutate(N_alpha = sum(n_cooc, na.rm = TRUE)) %>%
    ungroup()
  
  ## The conditioned expression is undefined exactly at p = 0.
  ## Its limit as p -> 0 is n_cooc / N_alpha.
  if(p <= 0){
    tmp <- tmp %>%
      mutate(
        prob_expected = ifelse(N_alpha > 0, n_cooc / N_alpha, NA_real_)
      )
    
    return(sum(tmp$prob_expected, na.rm = TRUE))
  }
  
  tmp <- tmp %>%
    mutate(
      numerator = 1 - (1 - p)^n_cooc,
      denominator = 1 - (1 - p)^N_alpha,
      prob_expected = ifelse(denominator > 0,
                             numerator / denominator,
                             NA_real_)
    )
  
  sum(tmp$prob_expected, na.rm = TRUE)
}

fit_p_conditioned_to_links <- function(cooc_pairs, L_emp){
  
  if(nrow(cooc_pairs) == 0 || is.na(L_emp)){
    return(NA_real_)
  }
  
  if(L_emp <= 0){
    return(0)
  }
  
  lower_expected <- expected_conditioned_from_p(cooc_pairs, 0)
  upper_expected <- expected_conditioned_from_p(cooc_pairs, 1)
  
  ## With the conditioned model, the lower bound is approximately
  ## one expected link per included consumer. If the observed number
  ## of links is at or below this, the best possible p is 0.
  if(L_emp <= lower_expected){
    return(0)
  }
  
  if(L_emp >= upper_expected){
    return(1)
  }
  
  f <- function(p){
    expected_conditioned_from_p(cooc_pairs, p) - L_emp
  }
  
  uniroot(f, interval = c(0, 1), tol = 1e-8)$root
}

## ------------------------------------------------------------
## Build original-style model object for one site subset
## ------------------------------------------------------------

build_original_style_model_object <- function(site_tables, sites_keep){
  
  ## Equivalent to the empirical metaweb for the retained sites.
  ## This corresponds to the original number_interactions_pred total:
  ## one unique consumer-resource interaction link counted from
  ## the consumer side.
  realised_pairs <- site_tables$empirical_site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource)
  
  L_emp <- nrow(realised_pairs)
  
  ## Co-occurrence table before the original-style merge.
  ## This corresponds to proportion_coocurrence before merging with
  ## number_interactions.
  cooc_pairs_all <- site_tables$cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")
  
  ## Original-style analogue of:
  ## number_interactions <- rbind(number_interactions_pred,
  ##                              number_interactions_prey)
  ##
  ## Then:
  ## interaction_proportion <- merge(number_interactions,
  ##                                 proportion_coocurrence,
  ##                                 by = "species")
  ##
  ## Since cooc_pairs_all$consumer is the original "species" column
  ## in proportion_coocurrence, this keeps co-occurrence rows whose
  ## consumer/species is present in the realised interaction network.
  realised_species_for_merge <- unique(c(realised_pairs$consumer,
                                         realised_pairs$resource))
  
  cooc_pairs_model <- cooc_pairs_all %>%
    filter(consumer %in% realised_species_for_merge)
  
  ## Diagnostics: interactions should normally be a subset of
  ## co-occurrences. If not, there is a loader/data issue.
  missing_realised_from_cooc <- realised_pairs %>%
    anti_join(cooc_pairs_all, by = c("consumer", "resource"))
  
  diagnostics <- data.frame(
    n_realised_links = L_emp,
    n_realised_consumers = length(unique(realised_pairs$consumer)),
    n_realised_resources = length(unique(realised_pairs$resource)),
    n_realised_species_for_merge = length(realised_species_for_merge),
    n_cooc_pairs_all = nrow(cooc_pairs_all),
    n_cooc_pairs_model = nrow(cooc_pairs_model),
    n_model_consumers = length(unique(cooc_pairs_model$consumer)),
    n_realised_links_missing_from_cooc = nrow(missing_realised_from_cooc)
  )
  
  list(
    realised_pairs = realised_pairs,
    cooc_pairs_all = cooc_pairs_all,
    cooc_pairs_model = cooc_pairs_model,
    L_emp = L_emp,
    diagnostics = diagnostics
  )
}

## ------------------------------------------------------------
## Refit one subset
## ------------------------------------------------------------

refit_for_subset <- function(dataset, site_tables, sites_keep, subset_row, p_fixed){
  
  obj <- build_original_style_model_object(site_tables, sites_keep)
  
  cooc_pairs <- obj$cooc_pairs_model
  L_emp <- obj$L_emp
  diag <- obj$diagnostics
  
  if(nrow(cooc_pairs) == 0){
    return(data.frame(
      dataset = dataset,
      subset_id = subset_row$subset_id,
      removal_fraction = subset_row$removal_fraction,
      site_rep = subset_row$site_rep,
      n_sites_kept = subset_row$n_sites_kept,
      p_fixed_original = p_fixed,
      p_refit_conditioned = NA_real_,
      n_realised_links = L_emp,
      expected_links_fixed_p_conditioned = NA_real_,
      expected_links_refitted_p_conditioned = NA_real_,
      expected_links_fixed_p_raw = NA_real_,
      fixed_p_divergence_conditioned = NA_real_,
      relative_fixed_p_divergence_conditioned = NA_real_,
      fixed_p_divergence_raw = NA_real_,
      relative_fixed_p_divergence_raw = NA_real_,
      n_realised_consumers = diag$n_realised_consumers,
      n_realised_resources = diag$n_realised_resources,
      n_realised_species_for_merge = diag$n_realised_species_for_merge,
      n_cooc_pairs_all = diag$n_cooc_pairs_all,
      n_cooc_pairs_model = diag$n_cooc_pairs_model,
      n_model_consumers = diag$n_model_consumers,
      n_realised_links_missing_from_cooc = diag$n_realised_links_missing_from_cooc
    ))
  }
  
  p_refit <- fit_p_conditioned_to_links(cooc_pairs, L_emp)
  
  exp_fixed_conditioned <- expected_conditioned_from_p(cooc_pairs, p_fixed)
  exp_refit_conditioned <- expected_conditioned_from_p(cooc_pairs, p_refit)
  
  ## Keep raw values only as diagnostics/comparison with older raw outputs.
  exp_fixed_raw <- expected_raw_from_p(cooc_pairs, p_fixed)
  
  data.frame(
    dataset = dataset,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_kept = subset_row$n_sites_kept,
    
    p_fixed_original = p_fixed,
    p_refit_conditioned = p_refit,
    
    n_realised_links = L_emp,
    
    expected_links_fixed_p_conditioned = exp_fixed_conditioned,
    expected_links_refitted_p_conditioned = exp_refit_conditioned,
    expected_links_fixed_p_raw = exp_fixed_raw,
    
    fixed_p_divergence_conditioned = L_emp - exp_fixed_conditioned,
    relative_fixed_p_divergence_conditioned =
      safe_relative_difference(L_emp, exp_fixed_conditioned),
    
    fixed_p_divergence_raw = L_emp - exp_fixed_raw,
    relative_fixed_p_divergence_raw =
      safe_relative_difference(L_emp, exp_fixed_raw),
    
    n_realised_consumers = diag$n_realised_consumers,
    n_realised_resources = diag$n_realised_resources,
    n_realised_species_for_merge = diag$n_realised_species_for_merge,
    n_cooc_pairs_all = diag$n_cooc_pairs_all,
    n_cooc_pairs_model = diag$n_cooc_pairs_model,
    n_model_consumers = diag$n_model_consumers,
    n_realised_links_missing_from_cooc = diag$n_realised_links_missing_from_cooc
  )
}

## ------------------------------------------------------------
## Run one dataset
## ------------------------------------------------------------

run_one_dataset <- function(dataset){
  
  message("Running script 7 original-style conditioned refit: ", dataset)
  
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_fixed <- get_p_fixed(dataset)
  site_tables <- get_dataset_site_tables(dataset)
  
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  
  by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    refit_for_subset(
      dataset = dataset,
      site_tables = site_tables,
      sites_keep = subset_object$subsets[[i]],
      subset_row = subset_object$index[i, ],
      p_fixed = p_fixed
    )
  }))
  
  ## Full-network baseline from the reconstructed/original-style object.
  ## This is the safest baseline for site-removal change.
  baseline <- by_subset %>%
    filter(removal_fraction == 0) %>%
    summarise(
      p_full_conditioned = mean(p_refit_conditioned, na.rm = TRUE),
      expected_full_fixed_p_conditioned =
        mean(expected_links_fixed_p_conditioned, na.rm = TRUE),
      realised_links_full = mean(n_realised_links, na.rm = TRUE),
      fixed_p_divergence_full_conditioned =
        mean(fixed_p_divergence_conditioned, na.rm = TRUE),
      relative_fixed_p_divergence_full_conditioned =
        mean(relative_fixed_p_divergence_conditioned, na.rm = TRUE),
      .groups = "drop"
    )
  
  ## Convert NaN to NA if needed.
  baseline[is.nan(as.matrix(baseline))] <- NA
  
  by_subset <- by_subset %>%
    mutate(
      p_full_conditioned = baseline$p_full_conditioned,
      
      ## Change relative to original hardcoded p.
      p_abs_change_vs_original_fixed =
        p_refit_conditioned - p_fixed_original,
      p_relative_change_vs_original_fixed =
        ifelse(p_fixed_original > 0,
               p_refit_conditioned / p_fixed_original,
               NA_real_),
      
      ## Change relative to the full reconstructed network.
      ## This should equal 1 at removal_fraction == 0.
      p_abs_change_vs_full =
        p_refit_conditioned - p_full_conditioned,
      p_relative_change_vs_full =
        ifelse(p_full_conditioned > 0,
               p_refit_conditioned / p_full_conditioned,
               NA_real_)
    )
  
  summary <- by_subset %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      p_fixed_original = first(p_fixed_original),
      p_full_conditioned = first(p_full_conditioned),
      
      mean_p_refit_conditioned = mean(p_refit_conditioned, na.rm = TRUE),
      sd_p_refit_conditioned = sd(p_refit_conditioned, na.rm = TRUE),
      
      mean_p_abs_change_vs_original_fixed =
        mean(p_abs_change_vs_original_fixed, na.rm = TRUE),
      mean_p_relative_change_vs_original_fixed =
        mean(p_relative_change_vs_original_fixed, na.rm = TRUE),
      
      mean_p_abs_change_vs_full =
        mean(p_abs_change_vs_full, na.rm = TRUE),
      mean_p_relative_change_vs_full =
        mean(p_relative_change_vs_full, na.rm = TRUE),
      
      mean_realised_links = mean(n_realised_links, na.rm = TRUE),
      
      mean_expected_links_fixed_p_conditioned =
        mean(expected_links_fixed_p_conditioned, na.rm = TRUE),
      mean_expected_links_refitted_p_conditioned =
        mean(expected_links_refitted_p_conditioned, na.rm = TRUE),
      mean_expected_links_fixed_p_raw =
        mean(expected_links_fixed_p_raw, na.rm = TRUE),
      
      mean_fixed_p_divergence_conditioned =
        mean(fixed_p_divergence_conditioned, na.rm = TRUE),
      mean_relative_fixed_p_divergence_conditioned =
        mean(relative_fixed_p_divergence_conditioned, na.rm = TRUE),
      
      mean_fixed_p_divergence_raw =
        mean(fixed_p_divergence_raw, na.rm = TRUE),
      mean_relative_fixed_p_divergence_raw =
        mean(relative_fixed_p_divergence_raw, na.rm = TRUE),
      
      mean_n_realised_consumers =
        mean(n_realised_consumers, na.rm = TRUE),
      mean_n_realised_resources =
        mean(n_realised_resources, na.rm = TRUE),
      mean_n_cooc_pairs_all =
        mean(n_cooc_pairs_all, na.rm = TRUE),
      mean_n_cooc_pairs_model =
        mean(n_cooc_pairs_model, na.rm = TRUE),
      mean_n_model_consumers =
        mean(n_model_consumers, na.rm = TRUE),
      mean_n_realised_links_missing_from_cooc =
        mean(n_realised_links_missing_from_cooc, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  baseline_check <- summary %>%
    filter(removal_fraction == 0) %>%
    mutate(
      full_refit_vs_original_fixed =
        p_full_conditioned / p_fixed_original,
      full_expected_fixed_vs_realised =
        mean_expected_links_fixed_p_conditioned / mean_realised_links
    )
  
  write.csv2(
    by_subset,
    file.path(out_dir, paste0(dataset, "_script7_refitted_p_by_subset.csv")),
    row.names = FALSE
  )
  
  write.csv2(
    summary,
    file.path(out_dir, paste0(dataset, "_script7_refitted_p_summary.csv")),
    row.names = FALSE
  )
  
  write.csv2(
    baseline_check,
    file.path(out_dir, paste0(dataset, "_script7_full_network_baseline_check.csv")),
    row.names = FALSE
  )
  
  list(
    by_subset = by_subset,
    summary = summary,
    baseline_check = baseline_check
  )
}

## ------------------------------------------------------------
## Run all datasets
## ------------------------------------------------------------

outs <- lapply(all_dataset_names, run_one_dataset)
names(outs) <- all_dataset_names

combined <- bind_rows(lapply(outs, `[[`, "by_subset"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary"))
combined_baseline_check <- bind_rows(lapply(outs, `[[`, "baseline_check"))

write.csv2(
  combined,
  file.path(combined_out, "script7_refitted_p_combined.csv"),
  row.names = FALSE
)

write.csv2(
  combined_summary,
  file.path(combined_out, "script7_refitted_p_summary_combined.csv"),
  row.names = FALSE
)

write.csv2(
  combined_baseline_check,
  file.path(combined_out, "script7_full_network_baseline_check_combined.csv"),
  row.names = FALSE
)

combined_summary$dataset <- factor(combined_summary$dataset,
                                   levels = all_dataset_names)

## ------------------------------------------------------------
## Plots
## ------------------------------------------------------------

p_refit_plot <- ggplot(
  combined_summary,
  aes(x = removal_fraction,
      y = mean_p_refit_conditioned)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_hline(aes(yintercept = p_fixed_original), linetype = 2) +
  geom_hline(aes(yintercept = p_full_conditioned), linetype = 3) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") +
  ylab("Conditioned refitted p") +
  ggtitle(
    "Conditioned refitted p under site removal",
    subtitle = "Dashed = original fixed p; dotted = full-network refit from current reconstruction"
  )

ggsave(
  file.path(combined_out, "07_refitted_p_vs_site_removal.png"),
  p_refit_plot,
  width = 14,
  height = 7,
  dpi = 300
)

p_rel_full <- ggplot(
  combined_summary,
  aes(x = removal_fraction,
      y = mean_p_relative_change_vs_full)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") +
  ylab("Conditioned refitted p / full-network refitted p") +
  ggtitle(
    "Relative change in conditioned p under site removal",
    subtitle = "Baseline is the full reconstructed network, so removal 0 should equal 1"
  )

ggsave(
  file.path(combined_out, "07_refitted_p_relative_to_full_vs_site_removal.png"),
  p_rel_full,
  width = 14,
  height = 7,
  dpi = 300
)

p_rel_original <- ggplot(
  combined_summary,
  aes(x = removal_fraction,
      y = mean_p_relative_change_vs_original_fixed)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") +
  ylab("Conditioned refitted p / original fixed p") +
  ggtitle(
    "Relative change in conditioned p compared with original fixed p",
    subtitle = "Use with baseline-check table; mismatch at 0 indicates object or rounding differences"
  )

ggsave(
  file.path(combined_out, "07_refitted_p_relative_to_original_fixed_vs_site_removal.png"),
  p_rel_original,
  width = 14,
  height = 7,
  dpi = 300
)

p_scatter <- ggplot(
  combined_summary,
  aes(x = mean_relative_fixed_p_divergence_conditioned,
      y = mean_p_relative_change_vs_full)
) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  facet_wrap(~ dataset, ncol = 5, scales = "free") +
  theme_classic(base_size = 11) +
  xlab("Conditioned fixed-p relative divergence") +
  ylab("Conditioned refitted p / full-network refitted p") +
  ggtitle("Conditioned fixed-p divergence vs refitted-p change")

ggsave(
  file.path(combined_out, "07_divergence_vs_refitted_p_change.png"),
  p_scatter,
  width = 14,
  height = 7,
  dpi = 300
)

message("Finished script 7 original-style conditioned refit.")