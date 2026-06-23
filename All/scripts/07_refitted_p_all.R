## ------------------------------------------------------------
## Script: All/scripts/07_refitted_p_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script7_refitted_p")
sep_out <- dirs$separated
combined_out <- dirs$combined

expected_raw_from_p <- function(n_cooc_vec, p){
  sum(1 - (1 - p)^n_cooc_vec, na.rm = TRUE)
}

expected_conditioned_from_p <- function(cooc_pairs, p){
  
  tmp <- cooc_pairs %>%
    group_by(consumer) %>%
    mutate(
      N_alpha = sum(n_cooc),
      prob_expected = (1 - (1 - p)^n_cooc) /
        (1 - (1 - p)^N_alpha)
    ) %>%
    ungroup()
  
  sum(tmp$prob_expected, na.rm = TRUE)
}

fit_p_to_links <- function(n_cooc_vec, L_emp){
  if(length(n_cooc_vec) == 0 || is.na(L_emp)) return(NA_real_)
  max_links <- length(n_cooc_vec)
  if(L_emp <= 0) return(0)
  if(L_emp >= max_links) return(1)
  f <- function(p) expected_conditioned_from_p(n_cooc_vec, p) - L_emp
  if(f(0) > 0 || f(1) < 0) return(NA_real_)
  uniroot(f, interval = c(0, 1), tol = 1e-8)$root
}

refit_for_subset <- function(dataset, site_tables, sites_keep, subset_row, p_fixed){
  cooc_pairs <- site_tables$cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")
  
  L_emp <- site_tables$empirical_site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource) %>%
    nrow()
  
  if(nrow(cooc_pairs) == 0){
    return(data.frame(dataset = dataset, subset_id = subset_row$subset_id,
                      removal_fraction = subset_row$removal_fraction, site_rep = subset_row$site_rep,
                      n_sites_kept = subset_row$n_sites_kept, p_fixed = p_fixed,
                      p_refit = NA_real_, p_abs_change = NA_real_, p_relative_change = NA_real_,
                      n_realised_links = L_emp, expected_links_fixed_p = NA_real_,
                      expected_links_refitted_p = NA_real_, fixed_p_divergence = NA_real_,
                      relative_fixed_p_divergence = NA_real_))
  }
  
  p_refit <- fit_p_to_links(cooc_pairs$n_cooc, L_emp)
  exp_fixed <- expected_raw_from_p(cooc_pairs$n_cooc, p_fixed)
  exp_refit <- expected_raw_from_p(cooc_pairs$n_cooc, p_refit)
  
  data.frame(
    dataset = dataset,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_kept = subset_row$n_sites_kept,
    p_fixed = p_fixed,
    p_refit = p_refit,
    p_abs_change = p_refit - p_fixed,
    p_relative_change = p_refit / p_fixed,
    n_realised_links = L_emp,
    expected_links_fixed_p = exp_fixed,
    expected_links_refitted_p = exp_refit,
    fixed_p_divergence = L_emp - exp_fixed,
    relative_fixed_p_divergence = safe_relative_difference(L_emp, exp_fixed)
  )
}

run_one_dataset <- function(dataset){
  message("Running script 7: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  p_fixed <- get_p_fixed(dataset)
  site_tables <- get_dataset_site_tables(dataset)
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  
  by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    refit_for_subset(dataset, site_tables, subset_object$subsets[[i]], subset_object$index[i,], p_fixed)
  }))
  
  summary <- by_subset %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      p_fixed = first(p_fixed),
      mean_p_refit = mean(p_refit, na.rm = TRUE),
      sd_p_refit = sd(p_refit, na.rm = TRUE),
      mean_p_abs_change = mean(p_abs_change, na.rm = TRUE),
      mean_p_relative_change = mean(p_relative_change, na.rm = TRUE),
      mean_realised_links = mean(n_realised_links, na.rm = TRUE),
      mean_expected_links_fixed_p = mean(expected_links_fixed_p, na.rm = TRUE),
      mean_fixed_p_divergence = mean(fixed_p_divergence, na.rm = TRUE),
      mean_relative_fixed_p_divergence = mean(relative_fixed_p_divergence, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv2(by_subset, file.path(out_dir, paste0(dataset, "_script7_refitted_p_by_subset.csv")), row.names = FALSE)
  write.csv2(summary, file.path(out_dir, paste0(dataset, "_script7_refitted_p_summary.csv")), row.names = FALSE)
  list(by_subset = by_subset, summary = summary)
}

outs <- lapply(all_dataset_names, run_one_dataset)
names(outs) <- all_dataset_names
combined <- bind_rows(lapply(outs, `[[`, "by_subset"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary"))

write.csv2(combined, file.path(combined_out, "script7_refitted_p_combined.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "script7_refitted_p_summary_combined.csv"), row.names = FALSE)

combined_summary$dataset <- factor(combined_summary$dataset, levels = all_dataset_names)

p_refit_plot <- ggplot(combined_summary, aes(removal_fraction, mean_p_refit)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  geom_hline(aes(yintercept = p_fixed), linetype = 2) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Refitted p") +
  ggtitle("Refitted p under site removal")

ggsave(file.path(combined_out, "07_refitted_p_vs_site_removal.png"),
       p_refit_plot, width = 14, height = 7, dpi = 300)

p_rel <- ggplot(combined_summary, aes(removal_fraction, mean_p_relative_change)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 11) +
  xlab("Fraction of sites removed") + ylab("Refitted p / fixed p") +
  ggtitle("Relative change in p under site removal")

ggsave(file.path(combined_out, "07_refitted_p_relative_change_vs_site_removal.png"),
       p_rel, width = 14, height = 7, dpi = 300)

p_scatter <- ggplot(combined_summary, aes(mean_relative_fixed_p_divergence, mean_p_relative_change)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point() +
  facet_wrap(~dataset, ncol = 5, scales = "free") +
  theme_classic(base_size = 11) +
  xlab("Fixed-p relative divergence") + ylab("Refitted p / fixed p") +
  ggtitle("Fixed-p divergence vs refitted-p change")

ggsave(file.path(combined_out, "07_divergence_vs_refitted_p_change.png"),
       p_scatter, width = 14, height = 7, dpi = 300)

message("Finished script 7.")
