
## ------------------------------------------------------------
## Script: All/scripts/11.2_summarise_community_predictors_all.R
##
## Purpose:
## Dataset-level summary of script 11 outputs.
##
## Main question:
## Across datasets, are communities with stronger fixed-p divergence also
## those with higher repeatability, stronger interaction heterogeneity,
## or stronger hub/generalist dependence?
##
## Outputs:
##   All/SeparatedResults/script11.2_dataset_level_predictors/
##   All/CombinedOutputs/
##
## Rules:
##   - reuse existing outputs
##   - no site-removal simulations
##   - combined plots only
##   - write.csv2() for all tables
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("ggrepel")
for(pkg in packages_extra){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
if(!require(patchwork)){
  install.packages("patchwork")
}
library(patchwork)
dirs <- make_output_dirs("script11.2_dataset_level_predictors")
sep_out <- dirs$separated
combined_out <- dirs$combined

## ---------------------------
## 1. File reading helpers
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

require_file <- function(file, label){
  x <- read_semicolon_or_comma(file)
  if(is.null(x)){
    stop("Missing required file for ", label, ": ", file)
  }
  x
}

get_value_at_removal <- function(df, value_col, removal_value){
  df %>%
    filter(removal_fraction == removal_value) %>%
    select(dataset, all_of(value_col)) %>%
    rename(!!paste0(value_col, "_", gsub("\\.", "_", as.character(removal_value))) := all_of(value_col))
}

slope_by_dataset <- function(df, value_col, out_col){
  df %>%
    filter(!is.na(.data[[value_col]]), !is.na(removal_fraction)) %>%
    group_by(dataset) %>%
    summarise(
      !!out_col := ifelse(
        n_distinct(removal_fraction) >= 2,
        coef(lm(.data[[value_col]] ~ removal_fraction))[2],
        NA_real_
      ),
      .groups = "drop"
    )
}

safe_cor <- function(df, response, predictor){
  cur <- df %>%
    select(dataset, response_value = all_of(response), predictor_value = all_of(predictor)) %>%
    filter(is.finite(response_value), is.finite(predictor_value))

  if(nrow(cur) < 4 || n_distinct(cur$response_value) < 2 || n_distinct(cur$predictor_value) < 2){
    return(data.frame(
      response_variable = response,
      predictor_variable = predictor,
      spearman_rho = NA_real_,
      p_value = NA_real_,
      n_datasets_used = nrow(cur)
    ))
  }

  test <- suppressWarnings(cor.test(cur$response_value, cur$predictor_value, method = "spearman", exact = FALSE))

  data.frame(
    response_variable = response,
    predictor_variable = predictor,
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    n_datasets_used = nrow(cur)
  )
}

## ---------------------------
## 2. Load script 11 outputs
## ---------------------------

community_file <- file.path(combined_out, "script11_community_predictors_summary_combined.csv")
hub_file <- file.path(combined_out, "script11_specialist_generalist_link_retention_combined.csv")
concentration_file <- file.path(combined_out, "script11_retained_link_concentration_combined.csv")

community <- require_file(community_file, "script 11 community predictors")
hub <- require_file(hub_file, "script 11 specialist/hub retention")
concentration <- require_file(concentration_file, "script 11 retained-link concentration")

## Ensure numeric removal fractions after csv2 import.
community$removal_fraction <- as.numeric(as.character(community$removal_fraction))
hub$removal_fraction <- as.numeric(as.character(hub$removal_fraction))
concentration$removal_fraction <- as.numeric(as.character(concentration$removal_fraction))

## ---------------------------
## 3. Choose divergence and refitted-p columns
## ---------------------------

divergence_candidates <- c(
  "conditioned_fixed_p_relative_divergence",
  "raw_fixed_p_relative_divergence",
  "relative_divergence_conditioned",
  "relative_divergence_raw",
  "mean_relative_divergence_conditioned",
  "mean_relative_divergence_raw"
)

divergence_col <- divergence_candidates[divergence_candidates %in% names(community)][1]

if(is.na(divergence_col)){
  stop("No usable divergence column found in script11_community_predictors_summary_combined.csv")
}

community <- community %>%
  mutate(
    divergence_metric_used = divergence_col,
    preferred_relative_divergence = .data[[divergence_col]]
  )

refit_candidates <- c(
  "refitted_p_relative_to_full_network",
  "refitted_p_relative_to_fixed_p",
  "mean_p_relative_change_conditioned",
  "mean_p_relative_change"
)

refit_col <- refit_candidates[refit_candidates %in% names(community)][1]

if(is.na(refit_col)){
  community$preferred_refitted_p_ratio <- NA_real_
  refit_col <- NA_character_
} else {
  community <- community %>%
    mutate(preferred_refitted_p_ratio = .data[[refit_col]])
}

## ---------------------------
## 4. Dataset-level community summary
## ---------------------------

base_dataset <- data.frame(dataset = sort(unique(community$dataset)))

div_summary <- community %>%
  select(dataset, removal_fraction, preferred_relative_divergence) %>%
  group_by(dataset) %>%
  summarise(
    divergence_at_0 = preferred_relative_divergence[match(0, removal_fraction)],
    divergence_at_0_4 = preferred_relative_divergence[match(0.4, removal_fraction)],
    divergence_at_0_8 = preferred_relative_divergence[match(0.8, removal_fraction)],
    divergence_change_0_to_0_8 = divergence_at_0_8 - divergence_at_0,
    divergence_slope = ifelse(
      n_distinct(removal_fraction[!is.na(preferred_relative_divergence)]) >= 2,
      coef(lm(preferred_relative_divergence ~ removal_fraction))[2],
      NA_real_
    ),
    maximum_divergence = max(preferred_relative_divergence, na.rm = TRUE),
    divergence_metric_used = first(community$divergence_metric_used[community$dataset == first(dataset)]),
    .groups = "drop"
  ) %>%
  mutate(maximum_divergence = ifelse(is.infinite(maximum_divergence), NA_real_, maximum_divergence))

refit_summary <- community %>%
  select(dataset, removal_fraction, preferred_refitted_p_ratio) %>%
  group_by(dataset) %>%
  summarise(
    refitted_p_ratio_at_0 = preferred_refitted_p_ratio[match(0, removal_fraction)],
    refitted_p_ratio_at_0_8 = preferred_refitted_p_ratio[match(0.8, removal_fraction)],
    refitted_p_ratio_change_0_to_0_8 = refitted_p_ratio_at_0_8 - refitted_p_ratio_at_0,
    refitted_p_ratio_slope = ifelse(
      n_distinct(removal_fraction[!is.na(preferred_refitted_p_ratio)]) >= 2,
      coef(lm(preferred_refitted_p_ratio ~ removal_fraction))[2],
      NA_real_
    ),
    refitted_p_ratio_metric_used = ifelse(is.na(refit_col), NA_character_, refit_col),
    .groups = "drop"
  )

repeatability_summary <- community %>%
  group_by(dataset) %>%
  summarise(
    mean_repeatability_at_0 = mean_repeatability[match(0, removal_fraction)],
    mean_repeatability_at_0_8 = mean_repeatability[match(0.8, removal_fraction)],
    mean_repeatability_change_0_to_0_8 = mean_repeatability_at_0_8 - mean_repeatability_at_0,
    mean_repeatability_slope = ifelse(
      "mean_repeatability" %in% names(cur_data()) &&
        n_distinct(removal_fraction[!is.na(mean_repeatability)]) >= 2,
      coef(lm(mean_repeatability ~ removal_fraction))[2],
      NA_real_
    ),
    one_site_links_at_0 = proportion_one_site_realised_links[match(0, removal_fraction)],
    one_site_links_at_0_8 = proportion_one_site_realised_links[match(0.8, removal_fraction)],
    one_site_links_change_0_to_0_8 = one_site_links_at_0_8 - one_site_links_at_0,
    one_site_links_slope = ifelse(
      "proportion_one_site_realised_links" %in% names(cur_data()) &&
        n_distinct(removal_fraction[!is.na(proportion_one_site_realised_links)]) >= 2,
      coef(lm(proportion_one_site_realised_links ~ removal_fraction))[2],
      NA_real_
    ),
    .groups = "drop"
  )

heterogeneity_summary <- community %>%
  group_by(dataset) %>%
  summarise(
    interaction_rate_variance_at_0 = variance_pair_interaction_rate[match(0, removal_fraction)],
    interaction_rate_variance_at_0_8 = variance_pair_interaction_rate[match(0.8, removal_fraction)],
    interaction_rate_variance_change_0_to_0_8 = interaction_rate_variance_at_0_8 - interaction_rate_variance_at_0,
    interaction_rate_variance_slope = ifelse(
      "variance_pair_interaction_rate" %in% names(cur_data()) &&
        n_distinct(removal_fraction[!is.na(variance_pair_interaction_rate)]) >= 2,
      coef(lm(variance_pair_interaction_rate ~ removal_fraction))[2],
      NA_real_
    ),
    zero_interaction_pairs_at_0 = proportion_zero_interaction_pairs[match(0, removal_fraction)],
    zero_interaction_pairs_at_0_8 = proportion_zero_interaction_pairs[match(0.8, removal_fraction)],
    zero_interaction_pairs_change_0_to_0_8 = zero_interaction_pairs_at_0_8 - zero_interaction_pairs_at_0,
    zero_interaction_pairs_slope = ifelse(
      "proportion_zero_interaction_pairs" %in% names(cur_data()) &&
        n_distinct(removal_fraction[!is.na(proportion_zero_interaction_pairs)]) >= 2,
      coef(lm(proportion_zero_interaction_pairs ~ removal_fraction))[2],
      NA_real_
    ),
    .groups = "drop"
  )

## ---------------------------
## 5. Hub/generalist dataset summaries
## ---------------------------

hub_summary <- hub %>%
  group_by(dataset, trophic_level, removal_fraction) %>%
  summarise(
    proportion_retained_links_attached_to_hubs = mean(proportion_retained_links_attached_to_hubs, na.rm = TRUE),
    proportion_retained_links_attached_to_specialists = mean(proportion_retained_links_attached_to_specialists, na.rm = TRUE),
    survival_hub_attached_links = mean(survival_hub_attached_links, na.rm = TRUE),
    survival_specialist_attached_links = mean(survival_specialist_attached_links, na.rm = TRUE),
    .groups = "drop"
  )

hub_dataset_level <- hub_summary %>%
  group_by(dataset, trophic_level) %>%
  summarise(
    hub_attached_link_proportion_at_0 = proportion_retained_links_attached_to_hubs[match(0, removal_fraction)],
    hub_attached_link_proportion_at_0_8 = proportion_retained_links_attached_to_hubs[match(0.8, removal_fraction)],
    hub_attached_link_proportion_change_0_to_0_8 =
      hub_attached_link_proportion_at_0_8 - hub_attached_link_proportion_at_0,
    hub_attached_link_proportion_slope = ifelse(
      n_distinct(removal_fraction[!is.na(proportion_retained_links_attached_to_hubs)]) >= 2,
      coef(lm(proportion_retained_links_attached_to_hubs ~ removal_fraction))[2],
      NA_real_
    ),
    hub_attached_link_survival_at_0_8 = survival_hub_attached_links[match(0.8, removal_fraction)],
    specialist_attached_link_survival_at_0_8 = survival_specialist_attached_links[match(0.8, removal_fraction)],
    hub_minus_specialist_survival_at_0_8 =
      hub_attached_link_survival_at_0_8 - specialist_attached_link_survival_at_0_8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = dataset,
    names_from = trophic_level,
    values_from = c(
      hub_attached_link_proportion_at_0,
      hub_attached_link_proportion_at_0_8,
      hub_attached_link_proportion_change_0_to_0_8,
      hub_attached_link_proportion_slope,
      hub_attached_link_survival_at_0_8,
      specialist_attached_link_survival_at_0_8,
      hub_minus_specialist_survival_at_0_8
    ),
    names_glue = "{.value}_{trophic_level}"
  )

## ---------------------------
## 6. Link concentration summaries
## ---------------------------

concentration_summary <- concentration %>%
  group_by(dataset, trophic_level, removal_fraction) %>%
  summarise(
    retained_link_gini = mean(retained_link_gini, na.rm = TRUE),
    top10_percent_species_link_share = mean(top10_percent_species_link_share, na.rm = TRUE),
    .groups = "drop"
  )

concentration_dataset_level <- concentration_summary %>%
  group_by(dataset, trophic_level) %>%
  summarise(
    retained_link_gini_at_0 = retained_link_gini[match(0, removal_fraction)],
    retained_link_gini_at_0_8 = retained_link_gini[match(0.8, removal_fraction)],
    retained_link_gini_change_0_to_0_8 = retained_link_gini_at_0_8 - retained_link_gini_at_0,
    top10_link_share_at_0 = top10_percent_species_link_share[match(0, removal_fraction)],
    top10_link_share_at_0_8 = top10_percent_species_link_share[match(0.8, removal_fraction)],
    top10_link_share_change_0_to_0_8 = top10_link_share_at_0_8 - top10_link_share_at_0,
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = dataset,
    names_from = trophic_level,
    values_from = c(
      retained_link_gini_at_0,
      retained_link_gini_at_0_8,
      retained_link_gini_change_0_to_0_8,
      top10_link_share_at_0,
      top10_link_share_at_0_8,
      top10_link_share_change_0_to_0_8
    ),
    names_glue = "{.value}_{trophic_level}"
  )

## ---------------------------
## 7. Final dataset-level table
## ---------------------------

dataset_summary <- base_dataset %>%
  left_join(div_summary, by = "dataset") %>%
  left_join(refit_summary, by = "dataset") %>%
  left_join(repeatability_summary, by = "dataset") %>%
  left_join(heterogeneity_summary, by = "dataset") %>%
  left_join(hub_dataset_level, by = "dataset") %>%
  left_join(concentration_dataset_level, by = "dataset")

write.csv2(
  dataset_summary,
  file.path(combined_out, "script11.2_dataset_level_predictor_summary.csv"),
  row.names = FALSE
)

write.csv2(
  dataset_summary,
  file.path(combined_out, "script11.2_dataset_level_predictor_matrix.csv"),
  row.names = FALSE
)

## ---------------------------
## 8. Spearman correlation table
## ---------------------------

responses <- c(
  "divergence_at_0_8",
  "divergence_slope",
  "maximum_divergence",
  "refitted_p_ratio_at_0_8"
)

predictors <- c(
  "mean_repeatability_at_0",
  "mean_repeatability_slope",
  "interaction_rate_variance_at_0",
  "interaction_rate_variance_slope",
  "one_site_links_at_0",
  "one_site_links_slope",
  "hub_attached_link_proportion_at_0_consumer",
  "hub_attached_link_proportion_at_0_resource",
  "hub_minus_specialist_survival_at_0_8_consumer",
  "hub_minus_specialist_survival_at_0_8_resource",
  "retained_link_gini_at_0_consumer",
  "retained_link_gini_at_0_resource",
  "top10_link_share_at_0_consumer",
  "top10_link_share_at_0_resource"
)

correlation_table <- bind_rows(lapply(responses, function(resp){
  bind_rows(lapply(predictors, function(pred){
    if(!(resp %in% names(dataset_summary)) || !(pred %in% names(dataset_summary))){
      return(data.frame(
        response_variable = resp,
        predictor_variable = pred,
        spearman_rho = NA_real_,
        p_value = NA_real_,
        n_datasets_used = 0
      ))
    }
    safe_cor(dataset_summary, resp, pred)
  }))
})) %>%
  arrange(desc(abs(spearman_rho)))

write.csv2(
  correlation_table,
  file.path(combined_out, "script11.2_dataset_level_spearman_correlations.csv"),
  row.names = FALSE
)

## ---------------------------
## 9. Ranked summary table
## ---------------------------

ranking_metrics <- c(
  "divergence_at_0_8",
  "divergence_slope",
  "mean_repeatability_at_0",
  "mean_repeatability_slope",
  "interaction_rate_variance_at_0",
  "interaction_rate_variance_slope",
  "hub_attached_link_proportion_at_0_consumer",
  "hub_attached_link_proportion_at_0_resource"
)

rankings <- dataset_summary %>%
  select(dataset, all_of(intersect(ranking_metrics, names(dataset_summary))))

for(metric in intersect(ranking_metrics, names(rankings))){
  rankings[[paste0(metric, "_rank")]] <- rank(-rankings[[metric]], ties.method = "min", na.last = "keep")
}

write.csv2(
  rankings,
  file.path(combined_out, "script11.2_dataset_rankings.csv"),
  row.names = FALSE
)

## ---------------------------
## 10. Combined dataset-level plots
## ---------------------------

plot_labeled_scatter <- function(df, x, y, xlab, ylab, title){

  if(!(x %in% names(df)) || !(y %in% names(df))){
    return(ggplot() + theme_void() + ggtitle(paste("Missing:", x, "or", y)))
  }

  plot_df <- df %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]]))

  ggplot(plot_df, aes(x = .data[[x]], y = .data[[y]], label = dataset)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, linetype = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    theme_classic(base_size = 11) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)
}

## Figure 1: divergence vs repeatability
p1a <- plot_labeled_scatter(
  dataset_summary,
  "mean_repeatability_at_0",
  "divergence_at_0_8",
  "Baseline mean repeatability",
  "Divergence at 0.8",
  "Divergence vs baseline repeatability"
)

p1b <- plot_labeled_scatter(
  dataset_summary,
  "mean_repeatability_slope",
  "divergence_at_0_8",
  "Repeatability slope",
  "Divergence at 0.8",
  "Divergence vs repeatability slope"
)

p1 <- p1a + p1b
ggsave(
  file.path(combined_out, "script11.2_divergence_vs_repeatability_dataset_level.png"),
  p1,
  width = 12,
  height = 5.5,
  dpi = 300
)

## Figure 2: divergence vs interaction heterogeneity
p2a <- plot_labeled_scatter(
  dataset_summary,
  "interaction_rate_variance_at_0",
  "divergence_at_0_8",
  "Baseline pair-level interaction-rate variance",
  "Divergence at 0.8",
  "Divergence vs baseline interaction heterogeneity"
)

p2b <- plot_labeled_scatter(
  dataset_summary,
  "interaction_rate_variance_slope",
  "divergence_at_0_8",
  "Interaction-rate variance slope",
  "Divergence at 0.8",
  "Divergence vs interaction heterogeneity slope"
)

p2 <- p2a + p2b
ggsave(
  file.path(combined_out, "script11.2_divergence_vs_interaction_heterogeneity_dataset_level.png"),
  p2,
  width = 12,
  height = 5.5,
  dpi = 300
)

## Figure 3: divergence vs hub dependence
p3a <- plot_labeled_scatter(
  dataset_summary,
  "hub_attached_link_proportion_at_0_consumer",
  "divergence_at_0_8",
  "Consumer hub-attached link proportion at 0",
  "Divergence at 0.8",
  "Divergence vs consumer hub dependence"
)

p3b <- plot_labeled_scatter(
  dataset_summary,
  "hub_attached_link_proportion_at_0_resource",
  "divergence_at_0_8",
  "Resource hub-attached link proportion at 0",
  "Divergence at 0.8",
  "Divergence vs resource hub dependence"
)

p3c <- plot_labeled_scatter(
  dataset_summary,
  "hub_minus_specialist_survival_at_0_8_consumer",
  "divergence_at_0_8",
  "Consumer hub survival minus specialist survival at 0.8",
  "Divergence at 0.8",
  "Divergence vs consumer hub survival advantage"
)

p3d <- plot_labeled_scatter(
  dataset_summary,
  "hub_minus_specialist_survival_at_0_8_resource",
  "divergence_at_0_8",
  "Resource hub survival minus specialist survival at 0.8",
  "Divergence at 0.8",
  "Divergence vs resource hub survival advantage"
)

p3 <- (p3a + p3b) / (p3c + p3d)
ggsave(
  file.path(combined_out, "script11.2_divergence_vs_hub_dependence_dataset_level.png"),
  p3,
  width = 12,
  height = 10,
  dpi = 300
)

## Figure 4: divergence vs link concentration
p4a <- plot_labeled_scatter(
  dataset_summary,
  "retained_link_gini_at_0_consumer",
  "divergence_at_0_8",
  "Consumer retained-link Gini at 0",
  "Divergence at 0.8",
  "Divergence vs consumer link concentration"
)

p4b <- plot_labeled_scatter(
  dataset_summary,
  "retained_link_gini_at_0_resource",
  "divergence_at_0_8",
  "Resource retained-link Gini at 0",
  "Divergence at 0.8",
  "Divergence vs resource link concentration"
)

p4 <- p4a + p4b
ggsave(
  file.path(combined_out, "script11.2_divergence_vs_link_concentration_dataset_level.png"),
  p4,
  width = 12,
  height = 5.5,
  dpi = 300
)

## Figure 5: refitted p ratio vs repeatability/heterogeneity
p5a <- plot_labeled_scatter(
  dataset_summary,
  "mean_repeatability_at_0",
  "refitted_p_ratio_at_0_8",
  "Baseline mean repeatability",
  "Refitted p ratio at 0.8",
  "Refitted p ratio vs baseline repeatability"
)

p5b <- plot_labeled_scatter(
  dataset_summary,
  "interaction_rate_variance_at_0",
  "refitted_p_ratio_at_0_8",
  "Baseline pair-level interaction-rate variance",
  "Refitted p ratio at 0.8",
  "Refitted p ratio vs interaction heterogeneity"
)

p5 <- p5a + p5b
ggsave(
  file.path(combined_out, "script11.2_refitted_p_vs_repeatability_heterogeneity_dataset_level.png"),
  p5,
  width = 12,
  height = 5.5,
  dpi = 300
)

message("Finished script 11.2 dataset-level predictor summary.")
