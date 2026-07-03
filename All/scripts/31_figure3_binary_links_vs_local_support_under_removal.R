## ------------------------------------------------------------
## Script: All/scripts/31_figure3_binary_links_vs_local_support_under_removal.R
## Binary links versus local interaction support under random site removal.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "patchwork", "scales", "future", "future.apply", "parallelly")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 50
out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, name, width = 11, height = 7){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}
make_pair_id <- function(consumer, resource) paste(consumer, resource, sep = "___")

run_one_dataset <- function(dataset){
  message("Running Figure 3 removal: ", dataset)
  site_tables <- get_dataset_site_tables(dataset)
  cooc <- site_tables$cooc_triples %>% distinct(site, consumer, resource) %>% mutate(pair_id = make_pair_id(consumer, resource))
  ints <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource) %>% mutate(pair_id = make_pair_id(consumer, resource))
  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)

  full_pairs <- cooc %>%
    group_by(pair_id, consumer, resource) %>% summarise(full_n = n_distinct(site), .groups = "drop") %>%
    left_join(ints %>% group_by(pair_id) %>% summarise(full_K = n_distinct(site), .groups = "drop"), by = "pair_id") %>%
    mutate(full_K = replace_na(full_K, 0L), full_B = full_K > 0)
  realised <- full_pairs %>% filter(full_K >= 1)
  full_sum_K <- sum(realised$full_K)
  full_links <- nrow(realised)

  rows <- lapply(seq_len(nrow(subset_object$index)), function(i){
    idx <- subset_object$index[i,]
    sites_keep <- subset_object$subsets[[idx$subset_id]]
    ret_n <- cooc %>% filter(site %in% sites_keep, pair_id %in% realised$pair_id) %>% group_by(pair_id) %>% summarise(retained_n = n_distinct(site), .groups = "drop")
    ret_K <- ints %>% filter(site %in% sites_keep, pair_id %in% realised$pair_id) %>% group_by(pair_id) %>% summarise(retained_K = n_distinct(site), .groups = "drop")
    link_level <- realised %>%
      select(pair_id, full_n, full_K) %>% left_join(ret_n, by = "pair_id") %>% left_join(ret_K, by = "pair_id") %>%
      mutate(retained_n = replace_na(retained_n, 0L), retained_K = replace_na(retained_K, 0L),
             state = case_when(retained_K > 0 ~ "Interaction still present", retained_n > 0 ~ "Species still recorded together, interaction not observed", TRUE ~ "Species no longer recorded together"))
    retention <- data.frame(dataset = dataset, removal_fraction = idx$removal_fraction, replicate = idx$site_rep, retained_site_count = idx$n_sites_kept,
                            local_interaction_support_remaining = sum(link_level$retained_K) / full_sum_K,
                            regional_interaction_links_still_present = sum(link_level$retained_K > 0) / full_links)
    states <- link_level %>% count(state, name = "n_links") %>% mutate(dataset = dataset, removal_fraction = idx$removal_fraction, replicate = idx$site_rep, fraction = n_links / full_links)
    list(retention = retention, states = states)
  })
  list(retention = bind_rows(lapply(rows, `[[`, "retention")), states = bind_rows(lapply(rows, `[[`, "states")))
}

workers <- max(1, min(length(all_dataset_names), parallelly::availableCores() - 1))
future::plan(future::multisession, workers = workers)
outputs <- future.apply::future_lapply(all_dataset_names, run_one_dataset, future.seed = TRUE)
future::plan(future::sequential)

retention_reps <- bind_rows(lapply(outputs, `[[`, "retention"))
state_reps <- bind_rows(lapply(outputs, `[[`, "states"))
write.csv2(retention_reps, file.path(out_dir, "31_site_removal_link_and_support_retention.csv"), row.names = FALSE)
write.csv2(state_reps, file.path(out_dir, "31_site_removal_original_link_states.csv"), row.names = FALSE)

retention_long <- retention_reps %>%
  pivot_longer(cols = c(local_interaction_support_remaining, regional_interaction_links_still_present), names_to = "quantity", values_to = "fraction_remaining") %>%
  mutate(quantity = recode(quantity, local_interaction_support_remaining = "Sites still supporting interactions", regional_interaction_links_still_present = "Regional interaction links still present"))

dataset_curves <- retention_long %>% group_by(dataset, removal_fraction, quantity) %>% summarise(median_value = median(fraction_remaining, na.rm = TRUE), .groups = "drop")
cross_curves <- dataset_curves %>% group_by(removal_fraction, quantity) %>% summarise(median_value = median(median_value, na.rm = TRUE), q25 = quantile(median_value, 0.25, na.rm = TRUE), q75 = quantile(median_value, 0.75, na.rm = TRUE), .groups = "drop")
cols <- c("Regional interaction links still present" = "#E69F00", "Sites still supporting interactions" = "#0072B2")

panel_a <- ggplot() +
  geom_ribbon(data = cross_curves, aes(x = removal_fraction, ymin = q25, ymax = q75, fill = quantity), alpha = 0.13) +
  geom_line(data = dataset_curves, aes(x = removal_fraction, y = median_value, group = interaction(dataset, quantity), colour = quantity), alpha = 0.20, linewidth = 0.45) +
  geom_line(data = cross_curves, aes(x = removal_fraction, y = median_value, colour = quantity), linewidth = 1.25) +
  geom_point(data = cross_curves, aes(x = removal_fraction, y = median_value, colour = quantity), size = 1.8) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) + coord_cartesian(ylim = c(0,1)) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") + ylab("Fraction remaining") +
  ggtitle("A. Regional links remain while local support is eroded")

state_levels <- c("Interaction still present", "Species still recorded together, interaction not observed", "Species no longer recorded together")
state_cols <- c("Interaction still present" = "#009E73", "Species still recorded together, interaction not observed" = "#CC79A7", "Species no longer recorded together" = "#999999")
state_cross <- state_reps %>%
  group_by(dataset, removal_fraction, replicate, state) %>% summarise(fraction = sum(fraction), .groups = "drop") %>%
  complete(dataset, removal_fraction, replicate, state = state_levels, fill = list(fraction = 0)) %>%
  group_by(dataset, removal_fraction, state) %>% summarise(dataset_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  group_by(removal_fraction, state) %>% summarise(fraction = mean(dataset_fraction, na.rm = TRUE), .groups = "drop") %>%
  mutate(state = factor(state, levels = state_levels)) %>% filter(removal_fraction > 0)

panel_b <- ggplot(state_cross, aes(x = removal_fraction, y = fraction, fill = state)) +
  geom_area(position = "fill", colour = "white", linewidth = 0.15, alpha = 0.95) +
  scale_fill_manual(values = state_cols, drop = FALSE) + scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom") +
  xlab("Proportion of sites removed") + ylab("Share of original interactions") +
  ggtitle("B. What happens to original interactions?", subtitle = "Species can still meet in the retained landscape after their observed interaction has disappeared.")

fig <- panel_a / panel_b + plot_annotation(title = "Binary links outlast local interaction support under site removal")
save_both(fig, "31_Figure3_binary_links_vs_local_support_under_removal", width = 11, height = 8)

checks <- data.frame(check = c("zero_removal_link_retention_is_1", "zero_removal_support_retention_is_1", "state_fractions_sum_to_1"),
                     pass = c(all(abs(retention_reps$regional_interaction_links_still_present[retention_reps$removal_fraction == 0] - 1) < 1e-10), all(abs(retention_reps$local_interaction_support_remaining[retention_reps$removal_fraction == 0] - 1) < 1e-10), state_reps %>% group_by(dataset, removal_fraction, replicate) %>% summarise(s = sum(fraction), .groups = "drop") %>% summarise(ok = all(abs(s - 1) < 1e-8)) %>% pull(ok)))
write.csv2(checks, file.path(out_dir, "31_Figure3_checks.csv"), row.names = FALSE)
message("Saved Figure 3 outputs in: ", out_dir)
