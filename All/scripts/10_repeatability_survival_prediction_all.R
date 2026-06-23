## ------------------------------------------------------------
## Script: All/scripts/10_repeatability_survival_prediction_all.R
## ------------------------------------------------------------
source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script10_repeatability_survival_prediction")
sep_out <- dirs$separated
combined_out <- dirs$combined

read_semicolon_or_comma <- function(file){
  if(!file.exists(file)) return(NULL)
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1) return(x)
  read.csv(file, stringsAsFactors = FALSE)
}

retention_for_dataset <- function(dataset){
  message("Running script 10: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p_fixed <- get_p_fixed(dataset)
  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples
  ints <- site_tables$empirical_site_interactions
  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  full_links <- ints %>% distinct(consumer, resource)
  L_full <- nrow(full_links)
  full_expected <- expected_links_from_subset(cooc_triples, all_sites, p_fixed)$expected_links_raw
  link_k <- ints %>% group_by(consumer, resource) %>% summarise(k_interacting_sites_full = n_distinct(site), .groups = "drop")
  by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    sites_keep <- subset_object$subsets[[i]]
    rem <- subset_object$index$removal_fraction[i]
    observed_links <- ints %>% filter(site %in% sites_keep) %>% distinct(consumer, resource) %>% nrow()
    expected_fixed <- expected_links_from_subset(cooc_triples, sites_keep, p_fixed)$expected_links_raw
    repeatability_predicted_links <- sum(1 - rem^link_k$k_interacting_sites_full, na.rm = TRUE)
    data.frame(
      dataset = dataset,
      subset_id = subset_object$index$subset_id[i],
      removal_fraction = rem,
      site_rep = subset_object$index$site_rep[i],
      n_sites_kept = subset_object$index$n_sites_kept[i],
      observed_retained_links = observed_links,
      observed_retained_proportion = observed_links / L_full,
      repeatability_predicted_links = repeatability_predicted_links,
      repeatability_predicted_proportion = repeatability_predicted_links / L_full,
      fixed_p_model_expected_links = expected_fixed,
      fixed_p_model_retained_proportion = expected_fixed / full_expected,
      observed_minus_repeatability_prediction = observed_links - repeatability_predicted_links,
      observed_minus_fixed_p_model = observed_links - expected_fixed
    )
  }))
  summary <- by_subset %>% group_by(dataset, removal_fraction) %>% summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  write.csv2(by_subset, file.path(out_dir, paste0(dataset, "_script10_repeatability_survival_prediction.csv")), row.names = FALSE)
  write.csv2(summary, file.path(out_dir, paste0(dataset, "_script10_repeatability_survival_prediction_summary.csv")), row.names = FALSE)
  list(by_subset = by_subset, summary = summary)
}

outs <- lapply(all_dataset_names, retention_for_dataset)
names(outs) <- all_dataset_names
combined <- bind_rows(lapply(outs, `[[`, "by_subset"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary"))
write.csv2(combined, file.path(combined_out, "script10_repeatability_survival_prediction_combined.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "script10_repeatability_survival_prediction_summary_combined.csv"), row.names = FALSE)
combined_summary$dataset <- factor(combined_summary$dataset, levels = all_dataset_names)

long <- combined_summary %>%
  select(dataset, removal_fraction, observed_retained_proportion, repeatability_predicted_proportion, fixed_p_model_retained_proportion) %>%
  pivot_longer(cols = -c(dataset, removal_fraction), names_to = "curve", values_to = "retention")

p_ret <- ggplot(long, aes(removal_fraction, retention, linetype = curve)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Retained-link proportion") +
  ggtitle("Retained-link proportion: empirical, repeatability-only, and fixed-p model")

ggsave(file.path(combined_out, "10_repeatability_survival_prediction.png"), p_ret, width = 14, height = 7, dpi = 300)

p_diff_rep <- ggplot(combined_summary, aes(removal_fraction, observed_minus_repeatability_prediction)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") + theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Observed - repeatability prediction") +
  ggtitle("Difference between observed links and repeatability-only prediction")

ggsave(file.path(combined_out, "10_observed_minus_repeatability_prediction.png"), p_diff_rep, width = 14, height = 7, dpi = 300)

p_diff_model <- ggplot(combined_summary, aes(removal_fraction, observed_minus_fixed_p_model)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") + theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Observed - fixed-p model") +
  ggtitle("Difference between observed links and fixed-p model")

ggsave(file.path(combined_out, "10_observed_minus_fixed_p_model.png"), p_diff_model, width = 14, height = 7, dpi = 300)

## Combined synthesis table and plots

div <- read_semicolon_or_comma(file.path(combined_out, "site_removal_summary_all_datasets.csv"))
conv <- read_semicolon_or_comma(file.path(combined_out, "script6_conversion_ratios_summary_combined.csv"))
refp <- read_semicolon_or_comma(file.path(combined_out, "script7_refitted_p_summary_combined.csv"))
rephet <- read_semicolon_or_comma(file.path(combined_out, "script8_repeatability_heterogeneity_summary_combined.csv"))
beta <- read_semicolon_or_comma(file.path(combined_out, "script9_interaction_beta_diversity_summary_combined.csv"))
surv <- combined_summary

synth <- surv %>% select(dataset, removal_fraction, repeatability_prediction_error = observed_minus_repeatability_prediction)
if(!is.null(div) && "mean_relative_divergence_raw" %in% names(div)) synth <- synth %>% left_join(div %>% select(dataset, removal_fraction, fixed_p_relative_divergence = mean_relative_divergence_raw), by = c("dataset", "removal_fraction"))
if(!is.null(conv)) synth <- synth %>% left_join(conv %>% select(dataset, removal_fraction, metaweb_conversion = mean_metaweb_conversion, site_level_conversion = mean_site_level_conversion), by = c("dataset", "removal_fraction"))
if(!is.null(refp)) synth <- synth %>% left_join(refp %>% select(dataset, removal_fraction, refitted_p = mean_p_refit, relative_change_p = mean_p_relative_change), by = c("dataset", "removal_fraction"))
if(!is.null(rephet)) synth <- synth %>% left_join(rephet %>% select(dataset, removal_fraction, mean_repeatability, proportion_one_site_realised_links, pair_level_interaction_rate_variance = variance_pair_interaction_rate), by = c("dataset", "removal_fraction"))
if(!is.null(beta)) synth <- synth %>% left_join(beta %>% select(dataset, removal_fraction, interaction_beta_diversity = interaction_mean_jaccard), by = c("dataset", "removal_fraction"))
write.csv2(synth, file.path(combined_out, "10_conversion_spatial_arrangement_summary.csv"), row.names = FALSE)

plot_scatter <- function(x, y, xlab, ylab, filename){
  if(!(x %in% names(synth)) || !(y %in% names(synth))) return(NULL)
  p <- ggplot(synth, aes(.data[[x]], .data[[y]])) + geom_point() +
    facet_wrap(~dataset, ncol = 5, scales = "free") + theme_classic(base_size = 11) +
    xlab(xlab) + ylab(ylab) + ggtitle(paste(ylab, "vs", xlab))
  ggsave(file.path(combined_out, filename), p, width = 14, height = 7, dpi = 300)
}

plot_scatter("fixed_p_relative_divergence", "relative_change_p", "Fixed-p relative divergence", "Refitted p / fixed p", "10_divergence_vs_refitted_p_change_synthesis.png")
plot_scatter("fixed_p_relative_divergence", "mean_repeatability", "Fixed-p relative divergence", "Mean repeatability", "10_divergence_vs_repeatability.png")
plot_scatter("fixed_p_relative_divergence", "proportion_one_site_realised_links", "Fixed-p relative divergence", "Proportion one-site links", "10_divergence_vs_one_site_links.png")
plot_scatter("fixed_p_relative_divergence", "interaction_beta_diversity", "Fixed-p relative divergence", "Interaction beta diversity", "10_divergence_vs_interaction_beta_diversity.png")
plot_scatter("relative_change_p", "mean_repeatability", "Refitted p / fixed p", "Mean repeatability", "10_refitted_p_change_vs_repeatability.png")
plot_scatter("metaweb_conversion", "site_level_conversion", "Metaweb conversion", "Site-level conversion", "10_metaweb_conversion_vs_site_level_conversion.png")

message("Finished script 10 and combined synthesis.")
