## ------------------------------------------------------------
## Script: All/scripts/09_interaction_beta_diversity_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100

dirs <- make_output_dirs("script9_interaction_beta_diversity")
sep_out <- dirs$separated
combined_out <- dirs$combined

pairwise_set_dissimilarity <- function(df, sites_keep, max_site_pairs = 5000){
  
  site_sets <- lapply(sites_keep, function(s){
    cur <- df %>% filter(site == s)
    unique(paste(cur$consumer, cur$resource, sep = "___"))
  })
  names(site_sets) <- sites_keep
  
  n_sites <- length(site_sets)
  
  if(n_sites < 2){
    return(data.frame(
      mean_jaccard = NA_real_,
      median_jaccard = NA_real_,
      mean_sorensen = NA_real_,
      median_sorensen = NA_real_
    ))
  }
  
  site_pairs <- combn(seq_len(n_sites), 2)
  
  if(ncol(site_pairs) > max_site_pairs){
    keep_cols <- sample(seq_len(ncol(site_pairs)), max_site_pairs)
    site_pairs <- site_pairs[, keep_cols, drop = FALSE]
  }
  
  vals_j <- numeric(ncol(site_pairs))
  vals_s <- numeric(ncol(site_pairs))
  
  for(k in seq_len(ncol(site_pairs))){
    
    a <- site_sets[[site_pairs[1, k]]]
    b <- site_sets[[site_pairs[2, k]]]
    
    union_n <- length(union(a, b))
    inter_n <- length(intersect(a, b))
    sum_n <- length(a) + length(b)
    
    vals_j[k] <- ifelse(union_n > 0, 1 - inter_n / union_n, NA_real_)
    vals_s[k] <- ifelse(sum_n > 0, 1 - (2 * inter_n) / sum_n, NA_real_)
  }
  
  data.frame(
    mean_jaccard = mean(vals_j, na.rm = TRUE),
    median_jaccard = median(vals_j, na.rm = TRUE),
    mean_sorensen = mean(vals_s, na.rm = TRUE),
    median_sorensen = median(vals_s, na.rm = TRUE)
  )
}

beta_for_subset <- function(dataset, site_tables, sites_keep, subset_row){
  message(dataset, " subset ", subset_row$subset_id,
          " removal ", subset_row$removal_fraction,
          " sites ", length(sites_keep))
  ints <- site_tables$empirical_site_interactions %>% filter(site %in% sites_keep)
  cooc <- site_tables$cooc_triples %>% filter(site %in% sites_keep)
  int_beta <- pairwise_set_dissimilarity(ints, sites_keep)
  cooc_beta <- pairwise_set_dissimilarity(cooc, sites_keep)
  data.frame(
    dataset = dataset,
    subset_id = subset_row$subset_id,
    removal_fraction = subset_row$removal_fraction,
    site_rep = subset_row$site_rep,
    n_sites_retained = length(sites_keep),
    n_sites_with_interactions = length(unique(ints$site)),
    interaction_mean_jaccard = int_beta$mean_jaccard,
    interaction_median_jaccard = int_beta$median_jaccard,
    interaction_mean_sorensen = int_beta$mean_sorensen,
    interaction_median_sorensen = int_beta$median_sorensen,
    cooccurrence_mean_jaccard = cooc_beta$mean_jaccard,
    cooccurrence_median_jaccard = cooc_beta$median_jaccard,
    cooccurrence_mean_sorensen = cooc_beta$mean_sorensen,
    cooccurrence_median_sorensen = cooc_beta$median_sorensen
  )
}

run_one_dataset <- function(dataset){
  message("Running script 9: ", dataset)
  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  site_tables <- get_dataset_site_tables(dataset)
  all_sites <- sort(unique(site_tables$cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    beta_for_subset(dataset, site_tables, subset_object$subsets[[i]], subset_object$index[i,])
  }))
  summary <- by_subset %>%
    group_by(dataset, removal_fraction) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  write.csv2(by_subset, file.path(out_dir, paste0(dataset, "_script9_interaction_beta_diversity_by_subset.csv")), row.names = FALSE)
  write.csv2(summary, file.path(out_dir, paste0(dataset, "_script9_interaction_beta_diversity_summary.csv")), row.names = FALSE)
  list(by_subset = by_subset, summary = summary)
}

outs <- lapply(all_dataset_names, run_one_dataset)
names(outs) <- all_dataset_names
combined <- bind_rows(lapply(outs, `[[`, "by_subset"))
combined_summary <- bind_rows(lapply(outs, `[[`, "summary"))
write.csv2(combined, file.path(combined_out, "script9_interaction_beta_diversity_combined.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "script9_interaction_beta_diversity_summary_combined.csv"), row.names = FALSE)
combined_summary$dataset <- factor(combined_summary$dataset, levels = all_dataset_names)

plot_beta <- function(metric, ylab, filename){
  p <- ggplot(combined_summary, aes(removal_fraction, .data[[metric]])) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    facet_wrap(~dataset, ncol = 5, scales = "free_y") +
    theme_classic(base_size = 11) +
    xlab("Fraction of sites removed") + ylab(ylab) + ggtitle(ylab)
  ggsave(file.path(combined_out, filename), p, width = 14, height = 7, dpi = 300)
}
plot_beta("interaction_mean_jaccard", "Mean interaction Jaccard dissimilarity", "combined_interaction_beta_diversity_vs_site_removal.png")
plot_beta("cooccurrence_mean_jaccard", "Mean co-occurrence Jaccard dissimilarity", "combined_cooccurrence_beta_diversity_vs_site_removal.png")

message("Finished script 9.")
