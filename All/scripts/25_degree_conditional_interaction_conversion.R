## ------------------------------------------------------------
## Script: All/scripts/25_degree_conditional_interaction_conversion.R
##
## Purpose:
## Empirical full-network analysis asking whether, for realised links
## with the same exact co-occurrence support n, local interaction
## conversion K/n differs among links attached to lower-, middle-, and
## higher-degree consumers.
##
## Uses original 10 Galiana datasets only.
## No site-removal simulation. No probabilistic model. No null model.
## No regressions/GAMs/power laws/beta-binomial models.
##
## Run from the parent repository folder.
## Outputs:
##   All/SeparatedResults/degree_conditional_interaction_conversion/<dataset>/
##   All/CombinedOutputs/

    )
  }

  bind_rows(boot_rows) %>%
    group_by(dataset, result_type, `consumer-degree group`, weighting) %>%
    summarise(
      bootstrap_q025 = quantile(value, 0.025, na.rm = TRUE),
      bootstrap_q975 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_n1_table <- function(link_table){
  link_table %>%
    group_by(dataset, `consumer-degree group`) %>%
    summarise(
      number_of_observed_links = n(),
      number_with_full_n_1 = sum(full_n == 1, na.rm = TRUE),
      proportion_with_full_n_1 = number_with_full_n_1 / number_of_observed_links,
      number_with_full_n_ge_2 = sum(full_n >= 2, na.rm = TRUE),
      proportion_with_full_n_ge_2 = number_with_full_n_ge_2 / number_of_observed_links,
      .groups = "drop"
    )
}

make_checks <- function(dataset, link_table, exact_n_summary, ref_weights){
  main_strata <- exact_n_summary %>%
    filter(eligible_exact_n) %>%
    group_by(dataset, full_n) %>%
    summarise(min_links = min(number_of_links), .groups = "drop")

  data.frame(
    dataset = dataset,
    all_links_K_ge_1 = all(link_table$full_K >= 1),
    all_links_K_le_n = all(link_table$full_K <= link_table$full_n),
    all_main_exact_n_have_min_5_each_group = ifelse(nrow(main_strata) == 0, TRUE,
                                                    all(main_strata$min_links >= min_links_per_group_per_n)),
    n1_excluded_from_main = !any(link_table$full_n == 1 & link_table$eligible_for_exact_n_main_analysis),
    reference_weights_sum_to_1 = ifelse(nrow(ref_weights) == 0, TRUE,
                                        all(abs(ref_weights %>% group_by(dataset) %>% summarise(s = sum(reference_weight), .groups = "drop") %>% pull(s) - 1) < 1e-8)),
    bootstrap_is_consumer_level = TRUE,
    no_site_removal_simulation_model_or_null_model_used = TRUE,
    n_links = nrow(link_table),
    n_eligible_exact_n = n_distinct(link_table$full_n[link_table$eligible_for_exact_n_main_analysis]),
    n_bootstrap_replicates = n_boot
  )
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tibble)
  })

  message("Running script 25 degree-conditional conversion: ", dataset)

  dataset_dir <- file.path(sep_out, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  built <- build_link_table(dataset)
  link_table <- built$link_table
  degree_groups <- built$degree_groups

  exact_n_summary <- summarise_exact_n(link_table)

  eligible_n <- exact_n_summary %>%
    filter(eligible_exact_n) %>%
    distinct(dataset, full_n)

  link_table <- link_table %>%
    left_join(eligible_n %>% mutate(eligible_for_exact_n_main_analysis = TRUE),
              by = c("dataset", "full_n")) %>%
    mutate(
      eligible_for_exact_n_main_analysis = replace_na(eligible_for_exact_n_main_analysis, FALSE),
      eligible_for_exact_n_main_analysis = eligible_for_exact_n_main_analysis & full_n > 1
    )

  exact_n_summary <- exact_n_summary %>%
    mutate(eligible_exact_n = full_n %in% eligible_n$full_n)

  ref_weights <- make_reference_weights(link_table)

  est_link <- standardized_conversion_link_weighted(link_table, ref_weights)
  boot_link <- bootstrap_standardized(link_table, ref_weights, "Link-weighted")
  std_link <- make_standardized_rows(est_link, link_table, ref_weights,
                                     "Link-weighted", boot_link)

  est_cons <- standardized_conversion_consumer_weighted(link_table, ref_weights)
  boot_cons <- bootstrap_standardized(link_table, ref_weights,
                                      "Consumer-weighted sensitivity")
  std_cons <- make_standardized_rows(est_cons, link_table, ref_weights,
                                     "Consumer-weighted sensitivity", boot_cons)

  std_all <- bind_rows(std_link, std_cons)

  n1_table <- make_n1_table(link_table)
  checks <- make_checks(dataset, link_table, exact_n_summary, ref_weights)

  write.csv2(link_table,
             file.path(dataset_dir, paste0(dataset, "_25_link_level_conversion_data.csv")),
             row.names = FALSE)
  write.csv2(exact_n_summary,
             file.path(dataset_dir, paste0(dataset, "_25_conversion_by_exact_n_and_degree.csv")),
             row.names = FALSE)
  write.csv2(std_all,
             file.path(dataset_dir, paste0(dataset, "_25_standardized_conversion_by_degree.csv")),
             row.names = FALSE)
  write.csv2(n1_table,
             file.path(dataset_dir, paste0(dataset, "_25_n1_link_structure_by_degree.csv")),
             row.names = FALSE)
  write.csv2(degree_groups,
             file.path(dataset_dir, paste0(dataset, "_25_initial_consumer_degree_groups.csv")),
             row.names = FALSE)

  list(
    link_table = link_table,
    exact_n_summary = exact_n_summary,
    standardized = std_all,
    n1_table = n1_table,
    degree_groups = degree_groups,
    checks = checks,
    ref_weights = ref_weights
  )
}

## ---------------------------
## Run
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)

names(all_outputs) <- all_dataset_names

link_level_all <- bind_rows(lapply(all_outputs, `[[`, "link_table"))
exact_n_all <- bind_rows(lapply(all_outputs, `[[`, "exact_n_summary"))
standardized_all <- bind_rows(lapply(all_outputs, `[[`, "standardized"))
n1_all <- bind_rows(lapply(all_outputs, `[[`, "n1_table"))
degree_groups_all <- bind_rows(lapply(all_outputs, `[[`, "degree_groups"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))
ref_weights_all <- bind_rows(lapply(all_outputs, `[[`, "ref_weights"))

## ---------------------------
## Save combined tables
## ---------------------------

write.csv2(link_level_all,
           file.path(combined_out, "25_link_level_conversion_data.csv"),
           row.names = FALSE)

write.csv2(exact_n_all,
           file.path(combined_out, "25_conversion_by_exact_n_and_degree.csv"),
           row.names = FALSE)

write.csv2(standardized_all,
           file.path(combined_out, "25_standardized_conversion_by_degree.csv"),
           row.names = FALSE)

write.csv2(n1_all,
           file.path(combined_out, "25_n1_link_structure_by_degree.csv"),
           row.names = FALSE)

write.csv2(degree_groups_all,
           file.path(combined_out, "25_initial_consumer_degree_groups.csv"),
           row.names = FALSE)

write.csv2(checks_all,
           file.path(combined_out, "25_degree_conditional_conversion_checks.csv"),
           row.names = FALSE)

## ---------------------------
## Figures
## ---------------------------

degree_cols <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#7570b3",
  "Higher initial degree" = "#d95f02"
)

support_cols_n <- c(
  "Recorded together in 1 site" = "#66c2a5",
  "Recorded together in 2 sites" = "#fc8d62",
  "Recorded together in 3-4 sites" = "#8da0cb",
  "Recorded together in 5 or more sites" = "#e78ac3"
)

## Figure 1: conversion by exact n and degree group
plot_exact <- exact_n_all %>%
  filter(eligible_exact_n, full_n > 1) %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = names(degree_cols))
  )

p1 <- ggplot(plot_exact,
             aes(x = full_n,
                 y = mean_conversion,
                 colour = `consumer-degree group`,
                 fill = `consumer-degree group`)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  ## These are IQR intervals from the observed link distribution, not fitted curves.
  geom_ribbon(aes(ymin = q25_conversion, ymax = q75_conversion),
              alpha = 0.12, colour = NA, na.rm = TRUE) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  scale_fill_manual(values = degree_cols, drop = FALSE) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom") +
  xlab("Number of sites where pair was recorded together") +
  ylab("Fraction of those sites with observed interaction") +
  ggtitle("Conversion given exact co-occurrence support",
          subtitle = "Comparisons shown only where every degree group has at least 5 links at the same co-occurrence support")

ggsave(file.path(combined_out, "25_conversion_given_exact_cooccurrence_support.png"),
       p1, width = 14, height = 7, dpi = 300)

## Figure 2: standardized conversion and higher-minus-lower contrast
std_plot_top <- standardized_all %>%
  filter(weighting == "Link-weighted",
         result_type == "Degree group standardized conversion") %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = names(degree_cols))
  )

std_plot_bottom <- standardized_all %>%
  filter(weighting == "Link-weighted",
         result_type == "Higher minus lower") %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

p2a <- ggplot(std_plot_top,
              aes(x = dataset,
                  y = standardized_conversion_or_contrast,
                  colour = `consumer-degree group`)) +
  geom_pointrange(aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
                  position = position_dodge(width = 0.55), linewidth = 0.45,
                  na.rm = TRUE) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  ylab("Standardized conversion") +
  ggtitle("Standardized conversion by initial consumer degree")

p2b <- ggplot(std_plot_bottom,
              aes(x = dataset,
                  y = standardized_conversion_or_contrast)) +
  geom_hline(yintercept = 0, colour = "grey65") +
  geom_pointrange(aes(ymin = bootstrap_q025, ymax = bootstrap_q975),
                  colour = "#333333", linewidth = 0.45, na.rm = TRUE) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Higher minus lower") +
  ggtitle("Difference in standardized conversion")

if(requireNamespace("patchwork", quietly = TRUE)){
  p2 <- p2a / p2b
  ggsave(file.path(combined_out, "25_standardized_conversion_by_initial_degree.png"),
         p2, width = 12, height = 8, dpi = 300)
} else {
  ggsave(file.path(combined_out, "25_standardized_conversion_by_initial_degree_top.png"),
         p2a, width = 12, height = 5, dpi = 300)
  ggsave(file.path(combined_out, "25_standardized_conversion_by_initial_degree_bottom.png"),
         p2b, width = 12, height = 5, dpi = 300)
}

## Figure 3: n = 1 structure
n1_plot <- n1_all %>%
  select(dataset, `consumer-degree group`, proportion_with_full_n_1, proportion_with_full_n_ge_2) %>%
  pivot_longer(cols = c(proportion_with_full_n_1, proportion_with_full_n_ge_2),
               names_to = "n_class", values_to = "proportion") %>%
  mutate(
    n_class = recode(n_class,
                     proportion_with_full_n_1 = "n = 1",
                     proportion_with_full_n_ge_2 = "n >= 2"),
    dataset = factor(dataset, levels = all_dataset_names),
    `consumer-degree group` = factor(`consumer-degree group`, levels = names(degree_cols))
  )

p3 <- ggplot(n1_plot,
             aes(x = `consumer-degree group`, y = proportion, fill = n_class)) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_fill_manual(values = c("n = 1" = "#bdbdbd", "n >= 2" = "#636363")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  ylab("Proportion of observed links") +
  labs(fill = "Co-occurrence support") +
  ggtitle("One-site opportunities by initial consumer degree")

ggsave(file.path(combined_out, "25_n1_links_by_initial_degree.png"),
       p3, width = 14, height = 7, dpi = 300)

## ---------------------------
## Console summary and notes
## ---------------------------

message("\nValidation summary:")
print(checks_all %>%
        select(dataset,
               all_links_K_ge_1,
               all_links_K_le_n,
               all_main_exact_n_have_min_5_each_group,
               n1_excluded_from_main,
               reference_weights_sum_to_1,
               bootstrap_is_consumer_level,
               no_site_removal_simulation_model_or_null_model_used,
               n_links,
               n_eligible_exact_n),
      n = Inf)

notes <- c(
  "Script 25: degree-conditional interaction conversion",
  "",
  "This is a purely empirical full-network analysis.",
  "No site-removal simulations, probabilistic models, null models, regressions, GAMs, power laws, or beta-binomial models were fitted or used.",
  "The main comparison excludes n = 1 because every observed link with n = 1 has K/n = 1 by construction.",
  "Exact-n strata are included only where lower-, middle-, and higher-degree consumer groups each have at least five observed links.",
  "Bootstrap intervals are generated by resampling consumers with replacement within degree groups and retaining all sampled links of each consumer.",
  "The bootstrap is used only for uncertainty display and no p-values are reported."
)
writeLines(notes, con = file.path(combined_out, "25_interpretation_notes.txt"))

future::plan(future::sequential)

message("Finished script 25 degree-conditional interaction conversion.")
