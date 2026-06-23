## ------------------------------------------------------------
## Script: All/scripts/11_community_predictors_and_link_concentration_all.R
## ------------------------------------------------------------
source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script11_community_predictors")
sep_out <- dirs$separated
combined_out <- dirs$combined

read_semicolon_or_comma <- function(file){
  if(!file.exists(file)) return(NULL)
  x <- tryCatch(read.csv2(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if(!is.null(x) && ncol(x) > 1) return(x)
  read.csv(file, stringsAsFactors = FALSE)
}

coalesce_col <- function(df, candidates, default = NA_real_){
  hit <- candidates[candidates %in% names(df)]
  if(length(hit) == 0) return(rep(default, nrow(df)))
  df[[hit[1]]]
}

## ------------------------------------------------------------
## Part A. Combined predictor summary from previous scripts
## ------------------------------------------------------------
script7 <- read_semicolon_or_comma(file.path(combined_out, "script7_refitted_p_summary_combined.csv"))
script8 <- read_semicolon_or_comma(file.path(combined_out, "script8_repeatability_heterogeneity_summary_combined.csv"))
degree <- read_semicolon_or_comma(file.path(combined_out, "degree_metrics_empirical_vs_model_envelope_all_datasets.csv"))
script10 <- read_semicolon_or_comma(file.path(combined_out, "script10_repeatability_survival_prediction_summary_combined.csv"))
div1 <- read_semicolon_or_comma(file.path(combined_out, "site_removal_summary_all_datasets.csv"))

if(!is.null(script7)){
  divergence_summary <- script7 %>%
    mutate(
      observed_realised_links = coalesce_col(., c("mean_realised_links")),
      expected_fixed_p_links_conditioned = coalesce_col(., c("mean_expected_links_fixed_p_conditioned", "mean_expected_links_fixed_p")),
      conditioned_fixed_p_absolute_divergence = coalesce_col(., c("mean_fixed_p_divergence_conditioned", "mean_fixed_p_divergence")),
      conditioned_fixed_p_relative_divergence = coalesce_col(., c("mean_relative_fixed_p_divergence_conditioned", "mean_relative_fixed_p_divergence")),
      raw_fixed_p_relative_divergence = coalesce_col(., c("mean_relative_fixed_p_divergence_raw", "mean_relative_fixed_p_divergence")),
      refitted_p = coalesce_col(., c("mean_p_refit_conditioned", "mean_p_refit")),
      refitted_p_relative_to_fixed_p = coalesce_col(., c("mean_p_relative_change_conditioned", "mean_p_relative_change"))
    ) %>%
    select(dataset, removal_fraction, observed_realised_links,
           expected_fixed_p_links_conditioned,
           conditioned_fixed_p_absolute_divergence,
           conditioned_fixed_p_relative_divergence,
           raw_fixed_p_relative_divergence,
           refitted_p,
           refitted_p_relative_to_fixed_p)
  
  full_p <- divergence_summary %>%
    filter(removal_fraction == 0) %>%
    select(dataset, full_network_refitted_p = refitted_p)
  
  divergence_summary <- divergence_summary %>%
    left_join(full_p, by = "dataset") %>%
    mutate(refitted_p_relative_to_full_network = refitted_p / full_network_refitted_p)
  
} else if(!is.null(div1)){
  divergence_summary <- div1 %>%
    transmute(
      dataset,
      removal_fraction,
      observed_realised_links = mean_realised_links,
      expected_fixed_p_links_conditioned = mean_expected_links_conditioned,
      conditioned_fixed_p_absolute_divergence = mean_divergence_conditioned,
      conditioned_fixed_p_relative_divergence = mean_relative_divergence_conditioned,
      raw_fixed_p_relative_divergence = mean_relative_divergence_raw,
      refitted_p = NA_real_,
      refitted_p_relative_to_fixed_p = NA_real_,
      full_network_refitted_p = NA_real_,
      refitted_p_relative_to_full_network = NA_real_
    )
} else {
  stop("No divergence file found. Run script 1 or script 7 first.")
}

community_summary <- divergence_summary

if(!is.null(script8)){
  rep_cols <- intersect(
    c("dataset", "removal_fraction", "mean_repeatability", "median_repeatability",
      "proportion_one_site_realised_links", "mean_interacting_sites_realised",
      "gini_interacting_sites_realised", "variance_pair_interaction_rate",
      "proportion_zero_interaction_pairs"),
    names(script8)
  )
  community_summary <- community_summary %>%
    left_join(script8[, rep_cols, drop = FALSE], by = c("dataset", "removal_fraction"))
}

if(!is.null(script10)){
  keep10 <- intersect(c("dataset", "removal_fraction", "observed_retained_proportion",
                        "repeatability_predicted_proportion", "fixed_p_model_retained_proportion",
                        "observed_minus_repeatability_prediction", "observed_minus_fixed_p_model"), names(script10))
  community_summary <- community_summary %>%
    left_join(script10[, keep10, drop = FALSE], by = c("dataset", "removal_fraction"))
}

if(!is.null(degree)){
  degree_summary_long <- degree %>%
    filter(metric %in% c("mean_degree", "maximum_degree", "gini_degree",
                         "proportion_degree_1", "n_species_degree_gt0")) %>%
    select(dataset, removal_fraction, trophic_level, metric, empirical_mean) %>%
    pivot_wider(names_from = metric, values_from = empirical_mean) %>%
    rename(
      degree_mean = mean_degree,
      degree_maximum = maximum_degree,
      degree_gini = gini_degree,
      degree1_proportion = proportion_degree_1,
      n_species_degree_gt0 = n_species_degree_gt0
    )
  
  full_deg <- degree_summary_long %>%
    filter(removal_fraction == 0) %>%
    select(dataset, trophic_level,
           full_degree_mean = degree_mean,
           full_degree_maximum = degree_maximum,
           full_degree_gini = degree_gini,
           full_degree1_proportion = degree1_proportion,
           full_n_species_degree_gt0 = n_species_degree_gt0)
  
  degree_summary_long <- degree_summary_long %>%
    left_join(full_deg, by = c("dataset", "trophic_level")) %>%
    mutate(
      relative_degree_mean = degree_mean / full_degree_mean,
      relative_degree_maximum = degree_maximum / full_degree_maximum,
      relative_degree_gini = degree_gini / full_degree_gini,
      relative_degree1_proportion = degree1_proportion / full_degree1_proportion,
      relative_n_species_degree_gt0 = n_species_degree_gt0 / full_n_species_degree_gt0
    )
  
  degree_wide <- degree_summary_long %>%
    pivot_wider(
      id_cols = c(dataset, removal_fraction),
      names_from = trophic_level,
      values_from = c(degree_mean, degree_maximum, degree_gini,
                      degree1_proportion, n_species_degree_gt0,
                      relative_degree_mean, relative_degree_maximum,
                      relative_degree_gini, relative_degree1_proportion,
                      relative_n_species_degree_gt0),
      names_glue = "{.value}_{trophic_level}"
    )
  
  community_summary <- community_summary %>%
    left_join(degree_wide, by = c("dataset", "removal_fraction"))
}

write.csv2(community_summary,
           file.path(combined_out, "script11_community_predictors_summary_combined.csv"),
           row.names = FALSE)

## ------------------------------------------------------------
## Parts C-D. Specialist/generalist/hub retention and concentration
## ------------------------------------------------------------
classify_full_metaweb_species <- function(site_tables){
  full_pairs <- site_tables$empirical_site_interactions %>% distinct(consumer, resource)
  consumer_deg <- full_pairs %>% count(species = consumer, name = "degree") %>% mutate(trophic_level = "consumer")
  resource_deg <- full_pairs %>% count(species = resource, name = "degree") %>% mutate(trophic_level = "resource")
  bind_rows(consumer_deg, resource_deg) %>%
    group_by(trophic_level) %>%
    mutate(
      top10_cutoff = quantile(degree, 0.90, na.rm = TRUE),
      top25_cutoff = quantile(degree, 0.75, na.rm = TRUE),
      n_top10 = sum(degree >= top10_cutoff, na.rm = TRUE),
      hub_cutoff = ifelse(n_top10 >= 2, top10_cutoff, top25_cutoff),
      species_class = case_when(
        degree == 1 ~ "specialist",
        degree >= hub_cutoff ~ "hub",
        TRUE ~ "intermediate"
      )
    ) %>%
    ungroup() %>%
    select(trophic_level, species, full_degree = degree, species_class)
}

link_class_summary_one_subset <- function(dataset, site_tables, species_classes, sites_keep, subset_row){
  retained_pairs <- site_tables$empirical_site_interactions %>% filter(site %in% sites_keep) %>% distinct(consumer, resource)
  full_pairs <- site_tables$empirical_site_interactions %>% distinct(consumer, resource)
  summarise_level <- function(level){
    if(level == "consumer"){
      retained_species_links <- retained_pairs %>% count(species = consumer, name = "retained_links")
      full_species_links <- full_pairs %>% count(species = consumer, name = "full_links")
    } else {
      retained_species_links <- retained_pairs %>% count(species = resource, name = "retained_links")
      full_species_links <- full_pairs %>% count(species = resource, name = "full_links")
    }
    classes <- species_classes %>% filter(trophic_level == level)
    joined <- full_species_links %>%
      full_join(retained_species_links, by = "species") %>%
      mutate(full_links = replace_na(full_links, 0), retained_links = replace_na(retained_links, 0)) %>%
      left_join(classes, by = "species")
    total_retained <- sum(joined$retained_links, na.rm = TRUE)
    data.frame(
      dataset = dataset,
      subset_id = subset_row$subset_id,
      removal_fraction = subset_row$removal_fraction,
      site_rep = subset_row$site_rep,
      n_sites_kept = subset_row$n_sites_kept,
      trophic_level = level,
      retained_links_total = total_retained,
      proportion_retained_links_attached_to_specialists = sum(joined$retained_links[joined$species_class == "specialist"], na.rm = TRUE) / total_retained,
      proportion_retained_links_attached_to_hubs = sum(joined$retained_links[joined$species_class == "hub"], na.rm = TRUE) / total_retained,
      survival_specialist_attached_links = sum(joined$retained_links[joined$species_class == "specialist"], na.rm = TRUE) / sum(joined$full_links[joined$species_class == "specialist"], na.rm = TRUE),
      survival_hub_attached_links = sum(joined$retained_links[joined$species_class == "hub"], na.rm = TRUE) / sum(joined$full_links[joined$species_class == "hub"], na.rm = TRUE)
    )
  }
  bind_rows(summarise_level("consumer"), summarise_level("resource")) %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.x), NA_real_, .x)))
}

link_concentration_one_subset <- function(dataset, site_tables, sites_keep, subset_row){
  retained_pairs <- site_tables$empirical_site_interactions %>% filter(site %in% sites_keep) %>% distinct(consumer, resource)
  summarise_level <- function(level){
    counts <- if(level == "consumer") retained_pairs %>% count(species = consumer, name = "link_count") else retained_pairs %>% count(species = resource, name = "link_count")
    if(nrow(counts) == 0){
      return(data.frame(dataset = dataset, subset_id = subset_row$subset_id,
                        removal_fraction = subset_row$removal_fraction, site_rep = subset_row$site_rep,
                        n_sites_kept = subset_row$n_sites_kept, trophic_level = level,
                        retained_link_gini = NA_real_, top1_species_link_share = NA_real_,
                        top5_percent_species_link_share = NA_real_, top10_percent_species_link_share = NA_real_,
                        effective_number_species = NA_real_))
    }
    x <- sort(counts$link_count, decreasing = TRUE)
    total <- sum(x)
    p <- x / total
    n_top5 <- max(1, ceiling(0.05 * length(x)))
    n_top10 <- max(1, ceiling(0.10 * length(x)))
    data.frame(
      dataset = dataset,
      subset_id = subset_row$subset_id,
      removal_fraction = subset_row$removal_fraction,
      site_rep = subset_row$site_rep,
      n_sites_kept = subset_row$n_sites_kept,
      trophic_level = level,
      retained_link_gini = gini_coefficient(x),
      top1_species_link_share = x[1] / total,
      top5_percent_species_link_share = sum(x[seq_len(n_top5)]) / total,
      top10_percent_species_link_share = sum(x[seq_len(n_top10)]) / total,
      effective_number_species = 1 / sum(p^2)
    )
  }
  bind_rows(summarise_level("consumer"), summarise_level("resource"))
}

run_link_class_and_concentration <- function(dataset){
  message("Script 11 new calculations: ", dataset)
  site_tables <- get_dataset_site_tables(dataset)
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  species_classes <- classify_full_metaweb_species(site_tables)
  link_class <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    link_class_summary_one_subset(dataset, site_tables, species_classes, subset_object$subsets[[i]], subset_object$index[i,])
  }))
  concentration <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    link_concentration_one_subset(dataset, site_tables, subset_object$subsets[[i]], subset_object$index[i,])
  }))
  list(link_class = link_class, concentration = concentration)
}

new_outputs <- lapply(all_dataset_names, run_link_class_and_concentration)
names(new_outputs) <- all_dataset_names
specialist_hub_retention <- bind_rows(lapply(new_outputs, `[[`, "link_class"))
retained_link_concentration <- bind_rows(lapply(new_outputs, `[[`, "concentration"))

write.csv2(specialist_hub_retention, file.path(combined_out, "script11_specialist_generalist_link_retention_combined.csv"), row.names = FALSE)
write.csv2(retained_link_concentration, file.path(combined_out, "script11_retained_link_concentration_combined.csv"), row.names = FALSE)

## ------------------------------------------------------------
## Combined plots only
## ------------------------------------------------------------
community_summary$dataset <- factor(community_summary$dataset, levels = all_dataset_names)

plot_scatter_dataset <- function(df, x, y, filename, xlab, ylab){
  if(!(x %in% names(df)) || !(y %in% names(df))) return(NULL)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet_wrap(~dataset, ncol = 5, scales = "free") +
    theme_classic(base_size = 11) + xlab(xlab) + ylab(ylab)
  ggsave(file.path(combined_out, filename), p, width = 14, height = 7, dpi = 300)
}

plot_scatter_dataset(community_summary, "mean_repeatability", "conditioned_fixed_p_relative_divergence", "script11_divergence_vs_repeatability.png", "Mean repeatability", "Conditioned fixed-p relative divergence")
plot_scatter_dataset(community_summary, "proportion_one_site_realised_links", "conditioned_fixed_p_relative_divergence", "script11_divergence_vs_one_site_links.png", "Proportion one-site realised links", "Conditioned fixed-p relative divergence")
plot_scatter_dataset(community_summary, "variance_pair_interaction_rate", "conditioned_fixed_p_relative_divergence", "script11_divergence_vs_interaction_rate_variance.png", "Pair-level interaction-rate variance", "Conditioned fixed-p relative divergence")
plot_scatter_dataset(community_summary, "mean_repeatability", "refitted_p_relative_to_full_network", "script11_refitted_p_change_vs_repeatability.png", "Mean repeatability", "Refitted p / full-network refitted p")
plot_scatter_dataset(community_summary, "variance_pair_interaction_rate", "refitted_p_relative_to_full_network", "script11_refitted_p_change_vs_interaction_rate_variance.png", "Pair-level interaction-rate variance", "Refitted p / full-network refitted p")

if(exists("degree_summary_long")){
  degree_plot_data <- degree_summary_long %>%
    left_join(divergence_summary %>% select(dataset, removal_fraction, conditioned_fixed_p_relative_divergence), by = c("dataset", "removal_fraction")) %>%
    mutate(dataset = factor(dataset, levels = all_dataset_names))
  
  plot_degree_scatter <- function(x, filename, xlab){
    p <- ggplot(degree_plot_data, aes(x = .data[[x]], y = conditioned_fixed_p_relative_divergence)) +
      geom_point(alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
      facet_grid(trophic_level ~ dataset, scales = "free") +
      theme_classic(base_size = 10) + xlab(xlab) + ylab("Conditioned fixed-p relative divergence")
    ggsave(file.path(combined_out, filename), p, width = 16, height = 7, dpi = 300)
  }
  plot_degree_scatter("degree_gini", "script11_divergence_vs_degree_gini.png", "Degree Gini")
  plot_degree_scatter("degree_maximum", "script11_divergence_vs_max_degree.png", "Maximum degree")
  plot_degree_scatter("degree1_proportion", "script11_divergence_vs_degree1_species.png", "Proportion degree-1 species")
}

specialist_hub_summary <- specialist_hub_retention %>%
  group_by(dataset, removal_fraction, trophic_level) %>%
  summarise(
    proportion_retained_links_attached_to_specialists = mean(proportion_retained_links_attached_to_specialists, na.rm = TRUE),
    proportion_retained_links_attached_to_hubs = mean(proportion_retained_links_attached_to_hubs, na.rm = TRUE),
    survival_specialist_attached_links = mean(survival_specialist_attached_links, na.rm = TRUE),
    survival_hub_attached_links = mean(survival_hub_attached_links, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(dataset = factor(dataset, levels = all_dataset_names))

p_hub <- ggplot(specialist_hub_summary, aes(removal_fraction, proportion_retained_links_attached_to_hubs)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_grid(trophic_level ~ dataset, scales = "free_y") + theme_classic(base_size = 10) + xlab("Fraction of sites removed") + ylab("Proportion retained links attached to hubs")
ggsave(file.path(combined_out, "script11_hub_attached_links_vs_site_removal.png"), p_hub, width = 16, height = 7, dpi = 300)

p_spec <- ggplot(specialist_hub_summary, aes(removal_fraction, proportion_retained_links_attached_to_specialists)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_grid(trophic_level ~ dataset, scales = "free_y") + theme_classic(base_size = 10) + xlab("Fraction of sites removed") + ylab("Proportion retained links attached to specialists")
ggsave(file.path(combined_out, "script11_specialist_attached_links_vs_site_removal.png"), p_spec, width = 16, height = 7, dpi = 300)

p_hub_spec <- ggplot(specialist_hub_summary, aes(survival_specialist_attached_links, survival_hub_attached_links)) + geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(alpha = 0.8) + facet_grid(trophic_level ~ dataset, scales = "free") + theme_classic(base_size = 10) + xlab("Specialist-attached link survival") + ylab("Hub-attached link survival")
ggsave(file.path(combined_out, "script11_hub_vs_specialist_link_survival.png"), p_hub_spec, width = 16, height = 7, dpi = 300)

hub_div <- specialist_hub_summary %>% left_join(divergence_summary %>% select(dataset, removal_fraction, conditioned_fixed_p_relative_divergence), by = c("dataset", "removal_fraction"))
p_hub_div <- ggplot(hub_div, aes(proportion_retained_links_attached_to_hubs, conditioned_fixed_p_relative_divergence)) + geom_point(alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) + facet_grid(trophic_level ~ dataset, scales = "free") + theme_classic(base_size = 10) + xlab("Proportion retained links attached to hubs") + ylab("Conditioned fixed-p relative divergence")
ggsave(file.path(combined_out, "script11_divergence_vs_hub_concentration.png"), p_hub_div, width = 16, height = 7, dpi = 300)

concentration_summary <- retained_link_concentration %>%
  group_by(dataset, removal_fraction, trophic_level) %>%
  summarise(
    retained_link_gini = mean(retained_link_gini, na.rm = TRUE),
    top10_percent_species_link_share = mean(top10_percent_species_link_share, na.rm = TRUE),
    effective_number_species = mean(effective_number_species, na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(dataset = factor(dataset, levels = all_dataset_names))

p_conc_gini <- ggplot(concentration_summary, aes(removal_fraction, retained_link_gini)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_grid(trophic_level ~ dataset, scales = "free_y") + theme_classic(base_size = 10) + xlab("Fraction of sites removed") + ylab("Gini of retained-link counts")
ggsave(file.path(combined_out, "script11_retained_link_gini_vs_site_removal.png"), p_conc_gini, width = 16, height = 7, dpi = 300)

p_top10 <- ggplot(concentration_summary, aes(removal_fraction, top10_percent_species_link_share)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_grid(trophic_level ~ dataset, scales = "free_y") + theme_classic(base_size = 10) + xlab("Fraction of sites removed") + ylab("Share links attached to top 10% species")
ggsave(file.path(combined_out, "script11_top10_link_share_vs_site_removal.png"), p_top10, width = 16, height = 7, dpi = 300)

conc_div <- concentration_summary %>% left_join(divergence_summary %>% select(dataset, removal_fraction, conditioned_fixed_p_relative_divergence), by = c("dataset", "removal_fraction"))
p_conc_div <- ggplot(conc_div, aes(retained_link_gini, conditioned_fixed_p_relative_divergence)) + geom_point(alpha = 0.8) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) + facet_grid(trophic_level ~ dataset, scales = "free") + theme_classic(base_size = 10) + xlab("Gini of retained-link counts") + ylab("Conditioned fixed-p relative divergence")
ggsave(file.path(combined_out, "script11_retained_link_concentration_vs_divergence.png"), p_conc_div, width = 16, height = 7, dpi = 300)

message("Finished script 11.")
