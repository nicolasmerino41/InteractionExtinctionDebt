## ------------------------------------------------------------
## Script: All/scripts/22_degree_dependent_opportunity_and_interaction_loss.R
##
## Purpose:
## Descriptive site-removal analysis of consumer disappearance,
## retained co-occurrence opportunity with original partners, retained
## interactions with original partners, and interaction loss despite
## retained opportunity, grouped by initial consumer degree.
##
## Uses the 10 original datasets through the shared helper script.
## Does not fit or use any probabilistic/null/beta-binomial model.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

set.seed(123)

## ---------------------------
## Settings
## ---------------------------

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500
min_sites_retained <- 1

result_type <- "script22_degree_dependent_opportunity_and_interaction_loss"
dirs <- make_output_dirs(result_type)
sep_out <- dirs$separated
combined_out <- dirs$combined

## ---------------------------
## Helpers
## ---------------------------

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

make_site_subsets_22 <- function(all_sites){
  ## Uses shared helper when available. Kept as wrapper so this script is explicit.
  make_site_subsets(all_sites, removal_levels, n_site_reps)
}

assign_degree_groups <- function(consumer_degree){
  consumer_degree <- consumer_degree %>% arrange(full_degree, consumer)
  n_cons <- nrow(consumer_degree)

  if(n_cons == 0){
    consumer_degree$initial_degree_group <- character()
    return(consumer_degree)
  }

  ## Keep tied degrees together whenever possible by ranking unique degree values.
  degree_values <- sort(unique(consumer_degree$full_degree))
  n_unique <- length(degree_values)

  if(n_unique == 1){
    lower_vals <- degree_values
    middle_vals <- numeric(0)
    higher_vals <- numeric(0)
  } else {
    degree_rank <- seq_along(degree_values)
    tercile <- ceiling(3 * degree_rank / n_unique)
    lower_vals <- degree_values[tercile == 1]
    middle_vals <- degree_values[tercile == 2]
    higher_vals <- degree_values[tercile == 3]
  }

  consumer_degree %>%
    mutate(
      initial_degree_group = case_when(
        full_degree %in% lower_vals ~ "Lower initial degree",
        full_degree %in% middle_vals ~ "Middle initial degree",
        full_degree %in% higher_vals ~ "Higher initial degree",
        TRUE ~ "Middle initial degree"
      ),
      initial_degree_group = factor(
        initial_degree_group,
        levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")
      )
    )
}

summarise_group_definitions <- function(consumer_degree){
  consumer_degree %>%
    group_by(dataset, initial_degree_group) %>%
    summarise(
      number_of_consumers = n(),
      minimum_full_degree = min(full_degree, na.rm = TRUE),
      maximum_full_degree = max(full_degree, na.rm = TRUE),
      median_full_degree = median(full_degree, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(
      dataset,
      initial_degree_group = factor(
        c("Lower initial degree", "Middle initial degree", "Higher initial degree"),
        levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")
      ),
      fill = list(number_of_consumers = 0)
    )
}

run_subset_for_dataset <- function(dataset, consumer_degree, original_links,
                                   consumer_presence_table,
                                   cooc_triples, empirical_site_interactions,
                                   sites_keep, removal_fraction, site_rep){

  retained_consumers <- consumer_presence_table %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer) %>%
    mutate(consumer_present = 1L)

  retained_n_by_link <- cooc_triples %>%
    filter(site %in% sites_keep) %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource)) %>%
    group_by(pair_id) %>%
    summarise(retained_n = n_distinct(site), .groups = "drop")

  retained_K_by_link <- empirical_site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(site, consumer, resource) %>%
    mutate(pair_id = make_pair_id(consumer, resource)) %>%
    group_by(pair_id) %>%
    summarise(retained_K = n_distinct(site), .groups = "drop")

  link_status <- original_links %>%
    select(dataset, consumer, resource, pair_id, full_degree) %>%
    left_join(retained_n_by_link, by = "pair_id") %>%
    left_join(retained_K_by_link, by = "pair_id") %>%
    mutate(
      retained_n = replace_na(retained_n, 0L),
      retained_K = replace_na(retained_K, 0L),
      opportunity_retained = retained_n > 0,
      interaction_retained = retained_K > 0,
      missing_despite_opportunity = retained_n > 0 & retained_K == 0
    )

  by_consumer <- link_status %>%
    group_by(dataset, consumer) %>%
    summarise(
      opportunity_retention_all_consumers = sum(opportunity_retained) / first(full_degree),
      interaction_retention_all_consumers = sum(interaction_retained) / first(full_degree),
      interaction_missing_despite_opportunity_all_consumers = sum(missing_despite_opportunity) / first(full_degree),
      .groups = "drop"
    ) %>%
    right_join(consumer_degree, by = c("dataset", "consumer")) %>%
    left_join(retained_consumers, by = "consumer") %>%
    mutate(
      consumer_present = replace_na(consumer_present, 0L),
      opportunity_retention_all_consumers = replace_na(opportunity_retention_all_consumers, 0),
      interaction_retention_all_consumers = replace_na(interaction_retention_all_consumers, 0),
      interaction_missing_despite_opportunity_all_consumers = replace_na(interaction_missing_despite_opportunity_all_consumers, 0),
      opportunity_retention_if_present = ifelse(consumer_present == 1L,
                                                opportunity_retention_all_consumers, NA_real_),
      interaction_retention_if_present = ifelse(consumer_present == 1L,
                                                interaction_retention_all_consumers, NA_real_),
      interaction_missing_despite_opportunity_if_present = ifelse(consumer_present == 1L,
                                                                  interaction_missing_despite_opportunity_all_consumers, NA_real_),
      removal_fraction = removal_fraction,
      site_rep = site_rep
    ) %>%
    select(dataset, consumer, initial_degree_group, full_degree,
           removal_fraction, site_rep, consumer_present,
           opportunity_retention_all_consumers,
           interaction_retention_all_consumers,
           interaction_missing_despite_opportunity_all_consumers,
           opportunity_retention_if_present,
           interaction_retention_if_present,
           interaction_missing_despite_opportunity_if_present)

  by_consumer
}

make_group_summary <- function(consumer_data){
  base <- consumer_data %>%
    group_by(dataset, removal_fraction, site_rep, initial_degree_group) %>%
    summarise(
      number_of_original_consumers_in_group = n(),
      number_of_consumers_still_present = sum(consumer_present == 1L, na.rm = TRUE),
      fraction_consumers_still_present = mean(consumer_present == 1L, na.rm = TRUE),
      mean_opportunity_retention_among_present = mean(opportunity_retention_if_present, na.rm = TRUE),
      mean_interaction_retention_among_present = mean(interaction_retention_if_present, na.rm = TRUE),
      mean_interaction_missing_despite_opportunity_among_present = mean(interaction_missing_despite_opportunity_if_present, na.rm = TRUE),
      .groups = "drop"
    )

  base %>%
    pivot_longer(
      cols = c(fraction_consumers_still_present,
               mean_opportunity_retention_among_present,
               mean_interaction_retention_among_present,
               mean_interaction_missing_despite_opportunity_among_present),
      names_to = "metric",
      values_to = "group_mean"
    ) %>%
    mutate(
      metric = recode(
        metric,
        fraction_consumers_still_present = "Fraction of consumers still present",
        mean_opportunity_retention_among_present = "Opportunity retention among present consumers",
        mean_interaction_retention_among_present = "Interaction retention among present consumers",
        mean_interaction_missing_despite_opportunity_among_present = "Interaction missing despite opportunity among present consumers"
      )
    )
}

make_aggregated_group_summary <- function(group_summary){
  group_summary %>%
    group_by(dataset, removal_fraction, initial_degree_group, metric) %>%
    summarise(
      number_of_original_consumers_in_group = first(number_of_original_consumers_in_group),
      median = median(group_mean, na.rm = TRUE),
      q025 = quantile(group_mean, 0.025, na.rm = TRUE),
      q975 = quantile(group_mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_correlations <- function(consumer_data){
  consumer_data %>%
    filter(removal_fraction > 0) %>%
    group_by(dataset, removal_fraction, consumer, full_degree) %>%
    summarise(
      probability_still_present = mean(consumer_present == 1L, na.rm = TRUE),
      mean_opportunity_retention_if_present = mean(opportunity_retention_if_present, na.rm = TRUE),
      mean_interaction_retention_if_present = mean(interaction_retention_if_present, na.rm = TRUE),
      mean_interaction_missing_despite_opportunity_if_present = mean(interaction_missing_despite_opportunity_if_present, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      cor_presence = safe_spearman(full_degree, probability_still_present),
      n_presence = sum(is.finite(full_degree) & is.finite(probability_still_present)),
      cor_opportunity = safe_spearman(full_degree, mean_opportunity_retention_if_present),
      n_opportunity = sum(is.finite(full_degree) & is.finite(mean_opportunity_retention_if_present)),
      cor_interaction = safe_spearman(full_degree, mean_interaction_retention_if_present),
      n_interaction = sum(is.finite(full_degree) & is.finite(mean_interaction_retention_if_present)),
      cor_missing = safe_spearman(full_degree, mean_interaction_missing_despite_opportunity_if_present),
      n_missing = sum(is.finite(full_degree) & is.finite(mean_interaction_missing_despite_opportunity_if_present)),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(cor_presence, cor_opportunity, cor_interaction, cor_missing),
      names_to = "test_code",
      values_to = "spearman_correlation"
    ) %>%
    mutate(
      test_name = recode(
        test_code,
        cor_presence = "Degree vs consumer presence",
        cor_opportunity = "Degree vs opportunity retention",
        cor_interaction = "Degree vs interaction retention",
        cor_missing = "Degree vs interaction missing despite opportunity"
      ),
      number_of_consumers_included = case_when(
        test_code == "cor_presence" ~ n_presence,
        test_code == "cor_opportunity" ~ n_opportunity,
        test_code == "cor_interaction" ~ n_interaction,
        test_code == "cor_missing" ~ n_missing
      )
    ) %>%
    select(dataset, removal_fraction, test_name, spearman_correlation, number_of_consumers_included)
}

make_checks <- function(dataset, consumer_degree, consumer_data, original_links,
                        cooc_triples, empirical_site_interactions){
  zero_data <- consumer_data %>% filter(removal_fraction == 0)

  interaction_cells <- empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(in_interactions = TRUE)

  cooc_cells <- cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(in_cooc = TRUE)

  outside <- interaction_cells %>%
    left_join(cooc_cells, by = c("site", "consumer", "resource")) %>%
    filter(is.na(in_cooc))

  decomposition <- consumer_data %>%
    filter(consumer_present == 1L) %>%
    mutate(
      opportunity_loss = 1 - opportunity_retention_if_present,
      decomposition_sum = interaction_retention_if_present +
        interaction_missing_despite_opportunity_if_present + opportunity_loss,
      decomposition_ok = abs(decomposition_sum - 1) < 1e-9
    )

  data.frame(
    dataset = dataset,
    all_original_consumers_have_full_degree_ge_1 = all(consumer_degree$full_degree >= 1),
    consumer_presence_equals_1_at_zero_removal = all(zero_data$consumer_present == 1L),
    opportunity_retention_equals_1_at_zero_removal = all(abs(zero_data$opportunity_retention_all_consumers - 1) < 1e-9),
    interaction_retention_equals_1_at_zero_removal = all(abs(zero_data$interaction_retention_all_consumers - 1) < 1e-9),
    interaction_retention_never_exceeds_opportunity_retention = all(
      consumer_data$interaction_retention_all_consumers <= consumer_data$opportunity_retention_all_consumers + 1e-9,
      na.rm = TRUE
    ),
    decomposition_sums_to_1_for_present_consumers = all(decomposition$decomposition_ok, na.rm = TRUE),
    no_interaction_outside_cooccurrence_cell = nrow(outside) == 0,
    no_probabilistic_or_null_model_used = TRUE
  )
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  message("Running script 22 degree-dependent opportunity and interaction loss: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)

  cooc_triples <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  empirical_site_interactions <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource),
      pair_id = make_pair_id(consumer, resource)
    )

  ## Consumers are considered present when recorded as consumers in occupancy
  ## if available; otherwise this falls back to co-occurrence records.
  if("occupancy" %in% names(site_tables) && nrow(site_tables$occupancy) > 0){
    consumer_presence_table <- site_tables$occupancy %>%
      filter(trophic_level == "consumer") %>%
      transmute(site = as.character(site), consumer = as.character(species)) %>%
      distinct(site, consumer)
  } else {
    consumer_presence_table <- cooc_triples %>%
      distinct(site, consumer)
  }

  original_links <- empirical_site_interactions %>%
    distinct(dataset = dataset, consumer, resource, pair_id)

  consumer_degree <- original_links %>%
    group_by(dataset, consumer) %>%
    summarise(full_degree = n_distinct(resource), .groups = "drop") %>%
    arrange(consumer) %>%
    assign_degree_groups()

  original_links <- original_links %>%
    left_join(consumer_degree %>% select(dataset, consumer, full_degree, initial_degree_group),
              by = c("dataset", "consumer"))

  degree_group_definitions <- summarise_group_definitions(consumer_degree)

  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets_22(all_sites)

  consumer_rows <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]

    consumer_rows[[i]] <- run_subset_for_dataset(
      dataset = dataset,
      consumer_degree = consumer_degree,
      original_links = original_links,
      consumer_presence_table = consumer_presence_table,
      cooc_triples = cooc_triples,
      empirical_site_interactions = empirical_site_interactions,
      sites_keep = sites_keep,
      removal_fraction = subset_object$index$removal_fraction[i],
      site_rep = subset_object$index$site_rep[i]
    )
  }

  consumer_data <- bind_rows(consumer_rows)
  group_summary <- make_group_summary(consumer_data)
  group_summary_aggregated <- make_aggregated_group_summary(group_summary)
  correlations <- make_correlations(consumer_data)
  checks <- make_checks(dataset, consumer_degree, consumer_data, original_links,
                        cooc_triples, empirical_site_interactions)

  write.csv2(consumer_data,
             file.path(out_dir, paste0(dataset, "_22_degree_dependent_loss_consumer_data.csv")),
             row.names = FALSE)
  write.csv2(group_summary,
             file.path(out_dir, paste0(dataset, "_22_degree_dependent_loss_group_summary.csv")),
             row.names = FALSE)
  write.csv2(group_summary_aggregated,
             file.path(out_dir, paste0(dataset, "_22_degree_dependent_loss_group_summary_aggregated.csv")),
             row.names = FALSE)
  write.csv2(correlations,
             file.path(out_dir, paste0(dataset, "_22_degree_dependent_loss_correlations.csv")),
             row.names = FALSE)
  write.csv2(degree_group_definitions,
             file.path(out_dir, paste0(dataset, "_22_initial_degree_group_definition.csv")),
             row.names = FALSE)
  write.csv2(checks,
             file.path(out_dir, paste0(dataset, "_22_degree_dependent_loss_checks.csv")),
             row.names = FALSE)

  list(
    consumer_data = consumer_data,
    group_summary = group_summary,
    group_summary_aggregated = group_summary_aggregated,
    correlations = correlations,
    degree_group_definitions = degree_group_definitions,
    checks = checks
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- lapply(all_dataset_names, run_one_dataset)
names(all_outputs) <- all_dataset_names

consumer_data_all <- bind_rows(lapply(all_outputs, `[[`, "consumer_data"))
group_summary_all <- bind_rows(lapply(all_outputs, `[[`, "group_summary"))
group_summary_aggregated_all <- bind_rows(lapply(all_outputs, `[[`, "group_summary_aggregated"))
correlations_all <- bind_rows(lapply(all_outputs, `[[`, "correlations"))
degree_group_definitions_all <- bind_rows(lapply(all_outputs, `[[`, "degree_group_definitions"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

consumer_data_all$dataset <- factor(consumer_data_all$dataset, levels = all_dataset_names)
group_summary_all$dataset <- factor(group_summary_all$dataset, levels = all_dataset_names)
group_summary_aggregated_all$dataset <- factor(group_summary_aggregated_all$dataset, levels = all_dataset_names)
correlations_all$dataset <- factor(correlations_all$dataset, levels = all_dataset_names)
degree_group_definitions_all$dataset <- factor(degree_group_definitions_all$dataset, levels = all_dataset_names)
checks_all$dataset <- factor(checks_all$dataset, levels = all_dataset_names)

write.csv2(consumer_data_all,
           file.path(combined_out, "22_degree_dependent_loss_consumer_data.csv"),
           row.names = FALSE)
write.csv2(group_summary_all,
           file.path(combined_out, "22_degree_dependent_loss_group_summary.csv"),
           row.names = FALSE)
write.csv2(group_summary_aggregated_all,
           file.path(combined_out, "22_degree_dependent_loss_group_summary_aggregated.csv"),
           row.names = FALSE)
write.csv2(correlations_all,
           file.path(combined_out, "22_degree_dependent_loss_correlations.csv"),
           row.names = FALSE)
write.csv2(degree_group_definitions_all,
           file.path(combined_out, "22_initial_degree_group_definition.csv"),
           row.names = FALSE)
write.csv2(checks_all,
           file.path(combined_out, "22_degree_dependent_loss_checks.csv"),
           row.names = FALSE)

message("Validation checks:")
print(checks_all)

if(any(!checks_all$all_original_consumers_have_full_degree_ge_1) ||
   any(!checks_all$consumer_presence_equals_1_at_zero_removal) ||
   any(!checks_all$opportunity_retention_equals_1_at_zero_removal) ||
   any(!checks_all$interaction_retention_equals_1_at_zero_removal) ||
   any(!checks_all$interaction_retention_never_exceeds_opportunity_retention) ||
   any(!checks_all$decomposition_sums_to_1_for_present_consumers) ||
   any(!checks_all$no_interaction_outside_cooccurrence_cell) ||
   any(!checks_all$no_probabilistic_or_null_model_used)){
  stop("At least one script 22 validation check failed. Inspect 22_degree_dependent_loss_checks.csv")
}

## ---------------------------
## Figures
## ---------------------------

plot_group_definitions <- degree_group_definitions_all %>%
  mutate(plot_group = number_of_consumers >= 3) %>%
  select(dataset, initial_degree_group, plot_group)

plot_data <- group_summary_aggregated_all %>%
  left_join(plot_group_definitions,
            by = c("dataset", "initial_degree_group")) %>%
  filter(plot_group)

group_colours <- c(
  "Lower initial degree" = "#1b9e77",
  "Middle initial degree" = "#7570b3",
  "Higher initial degree" = "#d95f02"
)

## Figure 1: consumer presence
p1_data <- plot_data %>%
  filter(metric == "Fraction of consumers still present")

p1 <- ggplot(p1_data,
             aes(x = removal_fraction,
                 y = median,
                 ymin = q025,
                 ymax = q975,
                 colour = initial_degree_group,
                 fill = initial_degree_group)) +
  geom_ribbon(alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = group_colours, name = "Initial-degree group") +
  scale_fill_manual(values = group_colours, name = "Initial-degree group") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Fraction of consumers still present")

ggsave(file.path(combined_out, "22_consumer_presence_by_initial_degree.png"),
       p1, width = 14, height = 7, dpi = 300)

## Figure 2: opportunity and interaction retention
p2_data <- plot_data %>%
  filter(metric %in% c("Opportunity retention among present consumers",
                       "Interaction retention among present consumers")) %>%
  mutate(
    metric = recode(
      metric,
      "Opportunity retention among present consumers" =
        "Fraction of original interaction partners still recorded together",
      "Interaction retention among present consumers" =
        "Fraction of original interaction partners still observed interacting"
    ),
    metric = factor(metric, levels = c(
      "Fraction of original interaction partners still recorded together",
      "Fraction of original interaction partners still observed interacting"
    ))
  )

p2 <- ggplot(p2_data,
             aes(x = removal_fraction,
                 y = median,
                 ymin = q025,
                 ymax = q975,
                 colour = initial_degree_group,
                 fill = initial_degree_group)) +
  geom_ribbon(alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(metric ~ dataset) +
  scale_colour_manual(values = group_colours, name = "Initial-degree group") +
  scale_fill_manual(values = group_colours, name = "Initial-degree group") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8)) +
  xlab("Proportion of sites removed") +
  ylab("")

ggsave(file.path(combined_out, "22_opportunity_and_interaction_retention_by_initial_degree.png"),
       p2, width = 18, height = 6.5, dpi = 300)

## Figure 3: missing despite opportunity
p3_data <- plot_data %>%
  filter(metric == "Interaction missing despite opportunity among present consumers")

p3 <- ggplot(p3_data,
             aes(x = removal_fraction,
                 y = median,
                 ymin = q025,
                 ymax = q975,
                 colour = initial_degree_group,
                 fill = initial_degree_group)) +
  geom_ribbon(alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ dataset, ncol = 5) +
  scale_colour_manual(values = group_colours, name = "Initial-degree group") +
  scale_fill_manual(values = group_colours, name = "Initial-degree group") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Fraction of original interaction partners still recorded together but no longer interacting")

ggsave(file.path(combined_out, "22_interaction_missing_despite_opportunity_by_initial_degree.png"),
       p3, width = 14, height = 7, dpi = 300)

## Figure 4: correlation summary
p4_data <- correlations_all %>%
  filter(removal_fraction > 0) %>%
  mutate(
    test_name = factor(test_name, levels = c(
      "Degree vs consumer presence",
      "Degree vs opportunity retention",
      "Degree vs interaction retention",
      "Degree vs interaction missing despite opportunity"
    )),
    removal_fraction_label = paste0("Removed ", removal_fraction)
  )

p4 <- ggplot(p4_data,
             aes(x = dataset,
                 y = spearman_correlation,
                 colour = factor(removal_fraction))) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(size = 1.9, alpha = 0.9,
             position = position_jitter(width = 0.12, height = 0)) +
  facet_wrap(~ test_name, ncol = 1) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Spearman correlation") +
  labs(colour = "Proportion removed")

ggsave(file.path(combined_out, "22_degree_dependent_loss_correlation_summary.png"),
       p4, width = 12, height = 10, dpi = 300)

message("Finished script 22 degree-dependent opportunity and interaction loss.")
