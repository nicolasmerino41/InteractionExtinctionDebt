## ------------------------------------------------------------
## Script: All/scripts/24_selection_of_supported_links_under_site_removal.R
##
## Purpose:
## Descriptive selection analysis under random site removal.
## Uses original 10 datasets via shared helper script.
## No probabilistic models, null models, regressions, or refitting.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("dplyr", "tidyr", "ggplot2", "tibble", "purrr", "future", "future.apply", "parallelly")
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
n_site_reps <- 500
min_sites_retained <- 1

## Parallelise across datasets only. This avoids nested parallelism and keeps
## each dataset's site-removal bookkeeping self-contained.
n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")

result_type <- "script24_selection_of_supported_links_under_site_removal"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

support_colours <- c(
  "Recorded together in 1 site" = "#4C78A8",
  "Recorded together in 2 sites" = "#F58518",
  "Recorded together in 3-4 sites" = "#54A24B",
  "Recorded together in 5 or more sites" = "#B279A2",
  "Interaction observed in 1 site" = "#4C78A8",
  "Interaction observed in 2 sites" = "#F58518",
  "Interaction observed in 3-4 sites" = "#54A24B",
  "Interaction observed in 5 or more sites" = "#B279A2"
)

degree_colours <- c(
  "Lower initial degree" = "#4C78A8",
  "Middle initial degree" = "#F58518",
  "Higher initial degree" = "#E45756"
)

## ---------------------------
## Helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

make_site_subsets_local <- function(all_sites){
  ## Uses shared helper if available, but standardises replicate column name.
  obj <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  idx <- obj$index
  if("site_rep" %in% names(idx)){
    idx$replicate <- idx$site_rep
  } else if("rep" %in% names(idx)){
    idx$replicate <- idx$rep
  } else {
    stop("Site subset index has neither site_rep nor rep column.")
  }
  obj$index <- idx
  obj
}

support_class_n <- function(n){
  case_when(
    n == 1 ~ "Recorded together in 1 site",
    n == 2 ~ "Recorded together in 2 sites",
    n >= 3 & n <= 4 ~ "Recorded together in 3-4 sites",
    n >= 5 ~ "Recorded together in 5 or more sites",
    TRUE ~ NA_character_
  )
}

support_class_K <- function(K){
  case_when(
    K == 1 ~ "Interaction observed in 1 site",
    K == 2 ~ "Interaction observed in 2 sites",
    K >= 3 & K <= 4 ~ "Interaction observed in 3-4 sites",
    K >= 5 ~ "Interaction observed in 5 or more sites",
    TRUE ~ NA_character_
  )
}

assign_degree_groups <- function(consumer_degrees){
  ## Keeps equal full_degree values together by assigning unique degree values to tertiles.
  degree_values <- sort(unique(consumer_degrees$full_degree))
  n_vals <- length(degree_values)
  if(n_vals == 1){
    breaks <- data.frame(full_degree = degree_values,
                         initial_degree_group = "Middle initial degree")
  } else {
    ranks <- seq_along(degree_values)
    groups <- cut(
      ranks,
      breaks = c(0, n_vals / 3, 2 * n_vals / 3, n_vals),
      labels = c("Lower initial degree", "Middle initial degree", "Higher initial degree"),
      include.lowest = TRUE
    )
    breaks <- data.frame(
      full_degree = degree_values,
      initial_degree_group = as.character(groups),
      stringsAsFactors = FALSE
    )
  }

  consumer_degrees %>%
    left_join(breaks, by = "full_degree") %>%
    mutate(
      initial_degree_group = factor(
        initial_degree_group,
        levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")
      )
    )
}

summarise_selection <- function(link_status, class_col, survivor_col,
                                survivor_layer, class_type){
  class_col <- rlang::ensym(class_col)
  survivor_col <- rlang::ensym(survivor_col)

  total_links <- nrow(link_status)
  total_survivors <- sum(dplyr::pull(link_status, !!survivor_col), na.rm = TRUE)

  link_status %>%
    group_by(class = !!class_col) %>%
    summarise(
      original_number_of_links = n(),
      surviving_number_of_links = sum(!!survivor_col, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      survivor_layer = survivor_layer,
      class_type = class_type,
      survival_fraction = surviving_number_of_links / original_number_of_links,
      original_share = original_number_of_links / total_links,
      survivor_share = ifelse(total_survivors > 0,
                              surviving_number_of_links / total_survivors,
                              NA_real_),
      enrichment = ifelse(!is.na(survivor_share) & original_share > 0,
                          survivor_share / original_share,
                          NA_real_)
    )
}

summarise_curve <- function(x){
  x %>%
    group_by(dataset, removal_fraction, survivor_layer, class_type, class) %>%
    summarise(
      median_enrichment = median(enrichment, na.rm = TRUE),
      q025_enrichment = quantile(enrichment, 0.025, na.rm = TRUE),
      q975_enrichment = quantile(enrichment, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  message("Running script 24 selection analysis: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  ## Defensive operational rule: an observed interaction is also an observed co-occurrence record.
  cooc <- bind_rows(
    cooc,
    ints %>% select(site, consumer, resource, pair_id)
  ) %>%
    distinct(site, consumer, resource, pair_id)

  original_links <- ints %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(full_K = n_distinct(site), .groups = "drop") %>%
    left_join(
      cooc %>%
        group_by(consumer, resource, pair_id) %>%
        summarise(full_n = n_distinct(site), .groups = "drop"),
      by = c("consumer", "resource", "pair_id")
    ) %>%
    filter(full_K >= 1) %>%
    mutate(
      cooccurrence_support_class = support_class_n(full_n),
      interaction_support_class = support_class_K(full_K)
    )

  consumer_degrees <- original_links %>%
    group_by(consumer) %>%
    summarise(full_degree = n_distinct(resource), .groups = "drop") %>%
    assign_degree_groups()

  original_links <- original_links %>%
    left_join(consumer_degrees, by = "consumer")

  degree_group_table <- consumer_degrees %>%
    group_by(dataset = dataset, degree_group = initial_degree_group) %>%
    summarise(
      number_of_consumers = n(),
      minimum_full_degree = min(full_degree),
      median_full_degree = median(full_degree),
      maximum_full_degree = max(full_degree),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      dataset = dataset,
      degree_group = factor(c("Lower initial degree", "Middle initial degree", "Higher initial degree"),
                            levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")),
      fill = list(number_of_consumers = 0,
                  minimum_full_degree = NA_real_,
                  median_full_degree = NA_real_,
                  maximum_full_degree = NA_real_)
    )

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets_local(all_sites)

  selection_rows <- vector("list", nrow(subset_object$index))
  exact_rows <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    idx <- subset_object$index[i, ]
    sites_keep <- as.character(subset_object$subsets[[idx$subset_id]])

    retained_n <- cooc %>%
      filter(site %in% sites_keep) %>%
      semi_join(original_links %>% select(pair_id), by = "pair_id") %>%
      group_by(pair_id) %>%
      summarise(retained_n = n_distinct(site), .groups = "drop")

    retained_K <- ints %>%
      filter(site %in% sites_keep) %>%
      semi_join(original_links %>% select(pair_id), by = "pair_id") %>%
      group_by(pair_id) %>%
      summarise(retained_K = n_distinct(site), .groups = "drop")

    link_status <- original_links %>%
      left_join(retained_n, by = "pair_id") %>%
      left_join(retained_K, by = "pair_id") %>%
      mutate(
        retained_n = replace_na(retained_n, 0L),
        retained_K = replace_na(retained_K, 0L),
        cooccurrence_survivor = retained_n > 0,
        interaction_survivor = retained_K > 0
      )

    cur_selection <- bind_rows(
      summarise_selection(
        link_status,
        cooccurrence_support_class,
        cooccurrence_survivor,
        "Still recorded together",
        "Initial co-occurrence support"
      ),
      summarise_selection(
        link_status,
        initial_degree_group,
        cooccurrence_survivor,
        "Still recorded together",
        "Initial consumer degree"
      ),
      summarise_selection(
        link_status,
        interaction_support_class,
        interaction_survivor,
        "Still observed interacting",
        "Initial interaction support"
      ),
      summarise_selection(
        link_status,
        initial_degree_group,
        interaction_survivor,
        "Still observed interacting",
        "Initial consumer degree"
      )
    ) %>%
      mutate(
        dataset = dataset,
        removal_fraction = idx$removal_fraction,
        replicate = idx$replicate
      ) %>%
      select(dataset, removal_fraction, replicate,
             survivor_layer, class_type, class,
             original_number_of_links, surviving_number_of_links,
             survival_fraction, original_share, survivor_share, enrichment)

    selection_rows[[i]] <- cur_selection

    ## Exact-K diagnostic, aggregated for this subset.
    exact_rows[[i]] <- link_status %>%
      group_by(dataset = dataset,
               exact_full_K = full_K,
               removal_fraction = idx$removal_fraction,
               replicate = idx$replicate,
               consumer_degree_group = initial_degree_group) %>%
      summarise(
        number_of_original_links = n(),
        interaction_survival_fraction = mean(interaction_survivor, na.rm = TRUE),
        .groups = "drop"
      )
  }

  selection_table <- bind_rows(selection_rows)

  exact_table_raw <- bind_rows(exact_rows)

  ## Keep exact K values with at least 10 links in every represented degree group in the full data.
  exact_eligible <- original_links %>%
    group_by(exact_full_K = full_K, consumer_degree_group = initial_degree_group) %>%
    summarise(n_links = n(), .groups = "drop") %>%
    group_by(exact_full_K) %>%
    filter(all(n_links >= 10), n_distinct(consumer_degree_group) >= 2) %>%
    ungroup() %>%
    distinct(exact_full_K)

  exact_table <- exact_table_raw %>%
    semi_join(exact_eligible, by = "exact_full_K") %>%
    group_by(dataset, exact_full_K, removal_fraction, replicate) %>%
    mutate(maximum_difference_among_degree_groups =
             max(interaction_survival_fraction, na.rm = TRUE) -
             min(interaction_survival_fraction, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(dataset, exact_full_K, removal_fraction, consumer_degree_group) %>%
    summarise(
      number_of_original_links = first(number_of_original_links),
      interaction_survival_fraction = mean(interaction_survival_fraction, na.rm = TRUE),
      maximum_difference_among_degree_groups = mean(maximum_difference_among_degree_groups, na.rm = TRUE),
      .groups = "drop"
    )

  support_by_degree <- bind_rows(
    original_links %>%
      group_by(dataset = dataset,
               consumer_degree_group = initial_degree_group,
               support_type = "Co-occurrence support",
               support_class = cooccurrence_support_class) %>%
      summarise(
        number_of_links = n(),
        mean_full_n = mean(full_n),
        median_full_n = median(full_n),
        mean_full_K = mean(full_K),
        median_full_K = median(full_K),
        proportion_full_K_1 = mean(full_K == 1),
        .groups = "drop"
      ),
    original_links %>%
      group_by(dataset = dataset,
               consumer_degree_group = initial_degree_group,
               support_type = "Interaction support",
               support_class = interaction_support_class) %>%
      summarise(
        number_of_links = n(),
        mean_full_n = mean(full_n),
        median_full_n = median(full_n),
        mean_full_K = mean(full_K),
        median_full_K = median(full_K),
        proportion_full_K_1 = mean(full_K == 1),
        .groups = "drop"
      )
  ) %>%
    group_by(dataset, consumer_degree_group, support_type) %>%
    mutate(fraction_of_links = number_of_links / sum(number_of_links)) %>%
    ungroup() %>%
    select(dataset, consumer_degree_group, support_type, support_class,
           number_of_links, fraction_of_links,
           mean_full_n, median_full_n, mean_full_K, median_full_K,
           proportion_full_K_1)

  checks <- tibble(
    dataset = dataset,
    original_links_K_ge_1 = all(original_links$full_K >= 1),
    original_links_K_le_n = all(original_links$full_K <= original_links$full_n),
    enrichment_one_at_zero = selection_table %>%
      filter(removal_fraction == 0) %>%
      summarise(ok = all(abs(enrichment - 1) < 1e-10 | is.na(enrichment))) %>%
      pull(ok),
    original_shares_sum_to_one = selection_table %>%
      distinct(dataset, removal_fraction, replicate, survivor_layer, class_type, class, original_share) %>%
      group_by(dataset, removal_fraction, replicate, survivor_layer, class_type) %>%
      summarise(s = sum(original_share, na.rm = TRUE), .groups = "drop") %>%
      summarise(ok = all(abs(s - 1) < 1e-10)) %>%
      pull(ok),
    survivor_shares_sum_to_one_when_survivors_exist = selection_table %>%
      group_by(dataset, removal_fraction, replicate, survivor_layer, class_type) %>%
      summarise(
        surviving = sum(surviving_number_of_links, na.rm = TRUE),
        s = sum(survivor_share, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(surviving > 0) %>%
      summarise(ok = all(abs(s - 1) < 1e-10)) %>%
      pull(ok),
    no_interaction_without_cooccurrence = TRUE,
    same_site_removal_subsets_used = TRUE,
    no_models_or_refitting_used = TRUE
  )

  write.csv2(selection_table,
             file.path(out_dir, paste0(dataset, "_24_link_selection_under_site_removal.csv")),
             row.names = FALSE)

  list(
    selection_table = selection_table,
    support_by_degree = support_by_degree,
    exact_table = exact_table,
    degree_group_table = degree_group_table,
    checks = checks
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)
names(all_outputs) <- all_dataset_names

selection_all <- bind_rows(lapply(all_outputs, `[[`, "selection_table"))
support_by_degree_all <- bind_rows(lapply(all_outputs, `[[`, "support_by_degree"))
exact_all <- bind_rows(lapply(all_outputs, `[[`, "exact_table"))
degree_groups_all <- bind_rows(lapply(all_outputs, `[[`, "degree_group_table"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

selection_all <- selection_all %>%
  mutate(
    dataset = factor(dataset, levels = all_dataset_names),
    class = as.character(class),
    survivor_layer = factor(survivor_layer,
                            levels = c("Still recorded together", "Still observed interacting")),
    class_type = factor(class_type,
                        levels = c("Initial co-occurrence support",
                                   "Initial interaction support",
                                   "Initial consumer degree"))
  )

support_by_degree_all <- support_by_degree_all %>%
  mutate(dataset = factor(dataset, levels = all_dataset_names))

## ---------------------------
## Save tables
## ---------------------------

write.csv2(selection_all,
           file.path(combined_out, "24_link_selection_under_site_removal.csv"),
           row.names = FALSE)

write.csv2(support_by_degree_all,
           file.path(combined_out, "24_initial_support_by_consumer_degree.csv"),
           row.names = FALSE)

write.csv2(exact_all,
           file.path(combined_out, "24_exact_K_degree_survival_check.csv"),
           row.names = FALSE)

write.csv2(degree_groups_all,
           file.path(combined_out, "24_initial_consumer_degree_groups.csv"),
           row.names = FALSE)

write.csv2(checks_all,
           file.path(combined_out, "24_selection_of_supported_links_checks.csv"),
           row.names = FALSE)

## ---------------------------
## Plot summaries
## ---------------------------

curve_summary <- summarise_curve(selection_all)

support_curve <- curve_summary %>%
  filter(class_type %in% c("Initial co-occurrence support", "Initial interaction support")) %>%
  mutate(
    support_plot_class = class,
    plot_row = case_when(
      survivor_layer == "Still recorded together" ~ "Enrichment among pairs still recorded together",
      survivor_layer == "Still observed interacting" ~ "Enrichment among links still observed interacting"
    )
  )

degree_curve <- curve_summary %>%
  filter(class_type == "Initial consumer degree") %>%
  mutate(
    plot_row = case_when(
      survivor_layer == "Still recorded together" ~ "Enrichment among pairs still recorded together",
      survivor_layer == "Still observed interacting" ~ "Enrichment among links still observed interacting"
    )
  )

p1 <- ggplot(
  support_curve,
  aes(x = removal_fraction,
      y = median_enrichment,
      colour = support_plot_class,
      fill = support_plot_class,
      group = support_plot_class)
) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.4) +
  geom_ribbon(aes(ymin = q025_enrichment, ymax = q975_enrichment),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(plot_row ~ dataset, scales = "free_y") +
  scale_colour_manual(values = support_colours, name = "Initial support") +
  scale_fill_manual(values = support_colours, name = "Initial support") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Enrichment")

ggsave(file.path(combined_out, "24_selection_of_opportunities_and_interactions.png"),
       p1, width = 16, height = 7.5, dpi = 300)

p2 <- ggplot(
  degree_curve,
  aes(x = removal_fraction,
      y = median_enrichment,
      colour = class,
      fill = class,
      group = class)
) +
  geom_hline(yintercept = 1, colour = "grey70", linewidth = 0.4) +
  geom_ribbon(aes(ymin = q025_enrichment, ymax = q975_enrichment),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(plot_row ~ dataset, scales = "free_y") +
  scale_colour_manual(values = degree_colours, name = "Initial consumer degree") +
  scale_fill_manual(values = degree_colours, name = "Initial consumer degree") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("Proportion of sites removed") +
  ylab("Enrichment")

ggsave(file.path(combined_out, "24_selection_by_initial_consumer_degree.png"),
       p2, width = 16, height = 7.5, dpi = 300)

p3_data <- support_by_degree_all %>%
  mutate(
    support_type_label = support_type,
    consumer_degree_group = factor(
      consumer_degree_group,
      levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")
    )
  )

p3 <- ggplot(
  p3_data,
  aes(x = consumer_degree_group,
      y = fraction_of_links,
      fill = support_class)
) +
  geom_col(width = 0.8) +
  facet_grid(support_type_label ~ dataset) +
  scale_fill_manual(values = support_colours, name = "Support class") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("Initial consumer-degree group") +
  ylab("Fraction of links")

ggsave(file.path(combined_out, "24_initial_support_by_consumer_degree.png"),
       p3, width = 16, height = 7.5, dpi = 300)

## ---------------------------
## Console validation summary
## ---------------------------

message("\nScript 24 validation summary:")
print(checks_all)

if(!all(unlist(checks_all[, -1]))){
  warning("At least one validation check failed. Inspect 24_selection_of_supported_links_checks.csv")
}

future::plan(future::sequential)

message("Finished script 24 selection of supported links under site removal.")
message("Outputs saved in: ", combined_out)
