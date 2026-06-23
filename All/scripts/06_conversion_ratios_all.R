## ------------------------------------------------------------
## Script: All/scripts/06_conversion_ratios_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script6_conversion_ratios")
sep_out <- dirs$separated
combined_out <- dirs$combined

conversion_for_subset <- function(dataset, site_tables, sites_keep, subset_row){
  cooc <- site_tables$cooc_triples %>% filter(site %in% sites_keep)
  ints <- site_tables$empirical_site_interactions %>% filter(site %in% sites_keep)
  
  n_cooc_pairs <- cooc %>% distinct(consumer, resource) %>% nrow()
  n_realised_links <- ints %>% distinct(consumer, resource) %>% nrow()
  n_cooc_site_opportunities <- cooc %>% distinct(site, consumer, resource) %>% nrow()
  n_interaction_site_presences <- ints %>% distinct(site, consumer, resource) %>% nrow()
  
  data.frame(
    dataset = dataset,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_kept = subset_row$n_sites_kept,
    n_realised_links = n_realised_links,
    n_cooc_pairs = n_cooc_pairs,
    metaweb_conversion = ifelse(n_cooc_pairs > 0, n_realised_links / n_cooc_pairs, NA_real_),
    n_interaction_site_presences = n_interaction_site_presences,
    n_cooc_site_opportunities = n_cooc_site_opportunities,
    site_level_conversion = ifelse(n_cooc_site_opportunities > 0,
                                   n_interaction_site_presences / n_cooc_site_opportunities,
                                   NA_real_)
  )
}

run_one_dataset <- function(dataset){
  message("Running script 6: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  site_tables <- get_dataset_site_tables(dataset)
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  
  by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    conversion_for_subset(dataset, site_tables, subset_object$subsets[[i]], subset_object$index[i,])
  }))
  
  summary <- by_subset %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      mean_metaweb_conversion = mean(metaweb_conversion, na.rm = TRUE),
      sd_metaweb_conversion = sd(metaweb_conversion, na.rm = TRUE),
      mean_site_level_conversion = mean(site_level_conversion, na.rm = TRUE),
      sd_site_level_conversion = sd(site_level_conversion, na.rm = TRUE),
      mean_realised_links = mean(n_realised_links, na.rm = TRUE),
      mean_cooc_pairs = mean(n_cooc_pairs, na.rm = TRUE),
      mean_interaction_site_presences = mean(n_interaction_site_presences, na.rm = TRUE),
      mean_cooc_site_opportunities = mean(n_cooc_site_opportunities, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv2(by_subset, file.path(out_dir, paste0(dataset, "_script6_conversion_ratios_by_subset.csv")), row.names = FALSE)
  write.csv2(summary, file.path(out_dir, paste0(dataset, "_script6_conversion_ratios_summary.csv")), row.names = FALSE)
  list(by_subset = by_subset, summary = summary)
}

outs <- lapply(all_dataset_names, run_one_dataset)
names(outs) <- all_dataset_names
combined <- bind_rows(lapply(outs, `[[`, "by_subset"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary"))

write.csv2(combined, file.path(combined_out, "script6_conversion_ratios_combined.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "script6_conversion_ratios_summary_combined.csv"), row.names = FALSE)

combined_summary$dataset <- factor(combined_summary$dataset, levels = all_dataset_names)

p_meta <- ggplot(combined_summary, aes(removal_fraction, mean_metaweb_conversion)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Metaweb conversion") +
  ggtitle("Metaweb co-occurrence-to-interaction conversion")

ggsave(file.path(combined_out, "06_metaweb_conversion_vs_site_removal.png"),
       p_meta, width = 14, height = 7, dpi = 300)

p_site <- ggplot(combined_summary, aes(removal_fraction, mean_site_level_conversion)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Site-level conversion") +
  ggtitle("Site-level co-occurrence-to-interaction conversion")

ggsave(file.path(combined_out, "06_site_level_conversion_vs_site_removal.png"),
       p_site, width = 14, height = 7, dpi = 300)

long <- combined_summary %>%
  select(dataset, removal_fraction, mean_metaweb_conversion, mean_site_level_conversion) %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "conversion_type", values_to = "conversion")

p_both <- ggplot(long, aes(removal_fraction, conversion, linetype = conversion_type)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Conversion ratio") +
  ggtitle("Metaweb and site-level conversion under site removal")

ggsave(file.path(combined_out, "06_conversion_ratios_vs_site_removal.png"),
       p_both, width = 14, height = 7, dpi = 300)

message("Finished script 6.")
