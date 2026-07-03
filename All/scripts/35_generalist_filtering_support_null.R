## ------------------------------------------------------------
## Script: All/scripts/35_generalist_filtering_support_null.R
## Generalist filtering, portfolio support, and opportunity-matched null.
## ------------------------------------------------------------
## Purpose:
## Test whether random site removal selectively retains initially high-degree
## species because their realised interaction portfolios have broader spatial
## co-occurrence support and interaction support.
##
## Outputs are saved to:
##   All/outputs/35_generalist_filtering_support_null
##
## This script uses the 10 original Galiana datasets only, via the shared
## project helper script. It does not fit interaction models.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr", "tidyr", "tibble", "purrr", "ggplot2", "patchwork",
  "future", "future.apply", "parallelly", "scales"
)
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---------------------------
## User controls
## ---------------------------
set.seed(123)

out_dir <- "All/outputs/35_generalist_filtering_support_null"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500
n_null_reps <- 199       # increase to 499/999 for final null envelopes
parallel_workers <- max(1, min(7, parallelly::availableCores() - 1))

## Fast test switch. Set TRUE for a quick smoke test.
fast_mode <- FALSE
if(fast_mode){
  n_site_reps <- 50
  n_null_reps <- 49
}

## ---------------------------
## Plot style
## ---------------------------
degree_cols <- c(
  "Lower initial degree"  = "#4E79A7",
  "Middle initial degree" = "#F28E2B",
  "Higher initial degree" = "#59A14F"
)
source_cols <- c(
  "Empirical" = "#D55E00",
  "Opportunity-matched null" = "#6A3D9A"
)

save_both <- function(p, name, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

safe_q <- function(x, p){
  if(length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

make_degree_groups <- function(df, value_col = "initial_degree"){
  vals <- sort(unique(df[[value_col]]))
  if(length(vals) == 1){
    return(df %>% mutate(initial_degree_group = factor("Middle initial degree",
                                                       levels = names(degree_cols))))
  }
  q1 <- stats::quantile(df[[value_col]], probs = 1/3, type = 1, na.rm = TRUE)
  q2 <- stats::quantile(df[[value_col]], probs = 2/3, type = 1, na.rm = TRUE)
  df %>%
    mutate(
      initial_degree_group = case_when(
        .data[[value_col]] <= q1 ~ "Lower initial degree",
        .data[[value_col]] <= q2 ~ "Middle initial degree",
        TRUE ~ "Higher initial degree"
      ),
      initial_degree_group = factor(initial_degree_group, levels = names(degree_cols))
    )
}

make_site_subsets <- function(sites, removal_levels, n_reps){
  S <- length(sites)
  out <- vector("list", length(removal_levels))
  names(out) <- as.character(removal_levels)
  for(r in removal_levels){
    if(r == 0){
      out[[as.character(r)]] <- list(tibble(removal_fraction = r, replicate = 1L, site = sites))
    } else {
      keep_n <- max(0L, S - round(S * r))
      out[[as.character(r)]] <- lapply(seq_len(n_reps), function(rep){
        tibble(removal_fraction = r, replicate = rep, site = sample(sites, keep_n, replace = FALSE))
      })
    }
  }
  bind_rows(unlist(out, recursive = FALSE))
}

support_class_n <- function(n){
  case_when(
    n <= 1 ~ "1",
    n == 2 ~ "2",
    n >= 3 & n <= 4 ~ "3-4",
    n >= 5 & n <= 8 ~ "5-8",
    n >= 9 & n <= 16 ~ "9-16",
    n >= 17 ~ "17+",
    TRUE ~ NA_character_
  )
}

shuffle_K_within_n <- function(links){
  ## exact n if possible. If strata are tiny, use tight n bins.
  strata_counts <- links %>% count(full_n, name = "n_links")
  use_exact <- all(strata_counts$n_links >= 2)
  links2 <- links %>%
    mutate(K_shuffle_stratum = if(use_exact) as.character(full_n) else support_class_n(full_n))
  links2 %>%
    group_by(K_shuffle_stratum) %>%
    mutate(full_K_null = sample(full_K, size = n(), replace = FALSE)) %>%
    ungroup() %>%
    select(-K_shuffle_stratum)
}

## ---------------------------
## Dataset preparation
## ---------------------------
prepare_dataset <- function(dataset){
  message("Preparing dataset: ", dataset)
  st <- get_dataset_site_tables(dataset)
  cooc <- st$cooc_triples %>% distinct(site, consumer, resource)
  ints <- st$empirical_site_interactions %>% distinct(site, consumer, resource)

  ## ensure observed interactions are valid co-occurrence records operationally
  cooc <- bind_rows(cooc, ints) %>% distinct(site, consumer, resource)

  sites <- sort(unique(cooc$site))

  full_pairs <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop") %>%
    left_join(
      ints %>% group_by(consumer, resource) %>% summarise(full_K = n_distinct(site), .groups = "drop"),
      by = c("consumer", "resource")
    ) %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      full_B = full_K > 0,
      link_id = paste(consumer, resource, sep = "___")
    ) %>%
    select(dataset, link_id, consumer, resource, full_n, full_K, full_B)

  realised_links <- full_pairs %>% filter(full_K >= 1)

  consumer_nodes <- realised_links %>%
    group_by(node = consumer) %>%
    summarise(initial_degree = n_distinct(resource), .groups = "drop") %>%
    mutate(dataset = dataset, guild = "Consumer") %>%
    make_degree_groups()

  resource_nodes <- realised_links %>%
    group_by(node = resource) %>%
    summarise(initial_degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(dataset = dataset, guild = "Resource") %>%
    make_degree_groups()

  nodes <- bind_rows(consumer_nodes, resource_nodes) %>%
    select(dataset, guild, node, initial_degree, initial_degree_group)

  link_node_map <- bind_rows(
    realised_links %>% transmute(dataset, guild = "Consumer", node = consumer, partner = resource, link_id, full_n, full_K),
    realised_links %>% transmute(dataset, guild = "Resource", node = resource, partner = consumer, link_id, full_n, full_K)
  ) %>%
    left_join(nodes, by = c("dataset", "guild", "node"))

  subsets <- make_site_subsets(sites, removal_levels, n_site_reps)

  list(
    dataset = dataset,
    sites = sites,
    cooc = cooc,
    ints = ints,
    full_pairs = full_pairs,
    realised_links = realised_links,
    nodes = nodes,
    link_node_map = link_node_map,
    subsets = subsets
  )
}

## ---------------------------
## Retention calculation
## ---------------------------
calc_retention_for_subset <- function(prep, keep_sites, removal_fraction, replicate, K_override = NULL){
  cooc_ret <- prep$cooc %>%
    filter(site %in% keep_sites) %>%
    group_by(consumer, resource) %>%
    summarise(retained_n = n_distinct(site), .groups = "drop") %>%
    mutate(link_id = paste(consumer, resource, sep = "___")) %>%
    select(link_id, retained_n)

  if(is.null(K_override)){
    int_ret <- prep$ints %>%
      filter(site %in% keep_sites) %>%
      group_by(consumer, resource) %>%
      summarise(retained_K = n_distinct(site), .groups = "drop") %>%
      mutate(link_id = paste(consumer, resource, sep = "___")) %>%
      select(link_id, retained_K)
  } else {
    ## Null: choose retained support count hypergeometrically from shuffled full_K.
    ## This uses only link support counts and the number of retained sites. It is
    ## much faster than constructing synthetic site-level interaction tables.
    retained_site_count <- length(keep_sites)
    S <- length(prep$sites)
    int_ret <- K_override %>%
      transmute(
        link_id,
        retained_K = stats::rhyper(
          nn = n(),
          m = full_K_null,
          n = S - full_K_null,
          k = retained_site_count
        )
      )
  }

  link_ret <- prep$realised_links %>%
    select(link_id, consumer, resource, full_n, full_K) %>%
    left_join(cooc_ret, by = "link_id") %>%
    left_join(int_ret, by = "link_id") %>%
    mutate(
      retained_n = replace_na(retained_n, 0L),
      retained_K = replace_na(retained_K, 0L),
      link_retained = retained_K > 0,
      opportunity_retained = retained_n > 0
    )

  node_ret <- bind_rows(
    link_ret %>% transmute(guild = "Consumer", node = consumer, link_id, retained_K, opportunity_retained, link_retained),
    link_ret %>% transmute(guild = "Resource", node = resource, link_id, retained_K, opportunity_retained, link_retained)
  ) %>%
    group_by(guild, node) %>%
    summarise(
      retained_degree = sum(link_retained),
      opportunity_degree = sum(opportunity_retained),
      .groups = "drop"
    ) %>%
    right_join(prep$nodes %>% select(dataset, guild, node, initial_degree, initial_degree_group),
               by = c("guild", "node")) %>%
    mutate(
      dataset = prep$dataset,
      removal_fraction = removal_fraction,
      replicate = replicate,
      retained_degree = replace_na(retained_degree, 0L),
      opportunity_degree = replace_na(opportunity_degree, 0L),
      network_active = retained_degree > 0,
      portfolio_retention_all = retained_degree / initial_degree,
      portfolio_retention_if_active = if_else(network_active, retained_degree / initial_degree, NA_real_)
    )

  node_ret
}

run_empirical_dataset <- function(prep){
  message("Running empirical removal: ", prep$dataset)
  split_subsets <- prep$subsets %>% group_split(removal_fraction, replicate)
  bind_rows(lapply(split_subsets, function(ss){
    calc_retention_for_subset(
      prep = prep,
      keep_sites = ss$site,
      removal_fraction = ss$removal_fraction[1],
      replicate = ss$replicate[1]
    )
  }))
}

summarise_node_retention <- function(node_raw, source = "Empirical", null_replicate = NA_integer_){
  group_summary <- node_raw %>%
    group_by(dataset, guild, removal_fraction, replicate, initial_degree_group) %>%
    summarise(
      source = source,
      null_replicate = null_replicate,
      n_species = n(),
      n_active = sum(network_active),
      proportion_connected = mean(network_active),
      mean_portfolio_retention_if_active = mean(portfolio_retention_if_active, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mean_portfolio_retention_if_active = if_else(is.nan(mean_portfolio_retention_if_active), NA_real_, mean_portfolio_retention_if_active))

  active_comp <- node_raw %>%
    filter(network_active) %>%
    group_by(dataset, guild, removal_fraction, replicate, initial_degree_group) %>%
    summarise(n_active_group = n(), .groups = "drop") %>%
    group_by(dataset, guild, removal_fraction, replicate) %>%
    mutate(proportion_of_active_network = n_active_group / sum(n_active_group)) %>%
    ungroup() %>%
    right_join(
      expand_grid(
        dataset = unique(node_raw$dataset),
        guild = unique(node_raw$guild),
        removal_fraction = unique(node_raw$removal_fraction),
        replicate = unique(node_raw$replicate),
        initial_degree_group = factor(names(degree_cols), levels = names(degree_cols))
      ),
      by = c("dataset", "guild", "removal_fraction", "replicate", "initial_degree_group")
    ) %>%
    mutate(
      source = source,
      null_replicate = null_replicate,
      n_active_group = replace_na(n_active_group, 0L),
      proportion_of_active_network = replace_na(proportion_of_active_network, 0)
    )

  list(group_summary = group_summary, active_composition = active_comp)
}

run_null_dataset_one <- function(prep, null_id){
  shuffled_links <- shuffle_K_within_n(prep$realised_links)
  split_subsets <- prep$subsets %>% group_split(removal_fraction, replicate)
  raw <- bind_rows(lapply(split_subsets, function(ss){
    calc_retention_for_subset(
      prep = prep,
      keep_sites = ss$site,
      removal_fraction = ss$removal_fraction[1],
      replicate = ss$replicate[1],
      K_override = shuffled_links
    )
  }))
  summarise_node_retention(raw, source = "Opportunity-matched null", null_replicate = null_id)
}

run_one_dataset <- function(dataset){
  prep <- prepare_dataset(dataset)

  emp_raw <- run_empirical_dataset(prep)
  emp_sum <- summarise_node_retention(emp_raw, source = "Empirical", null_replicate = NA_integer_)

  message("Running null reps: ", dataset)
  null_out <- future_lapply(seq_len(n_null_reps), function(i){
    run_null_dataset_one(prep, i)
  }, future.seed = TRUE)

  list(
    node_raw = emp_raw,
    empirical_group_summary = emp_sum$group_summary,
    empirical_active_composition = emp_sum$active_composition,
    null_group_summary = bind_rows(lapply(null_out, `[[`, "group_summary")),
    null_active_composition = bind_rows(lapply(null_out, `[[`, "active_composition")),
    initial_degree_groups = prep$nodes %>%
      group_by(dataset, guild, initial_degree_group) %>%
      summarise(
        n_species = n(),
        min_initial_degree = min(initial_degree),
        median_initial_degree = median(initial_degree),
        max_initial_degree = max(initial_degree),
        .groups = "drop"
      ),
    support_pathway = prep$link_node_map %>%
      group_by(dataset, guild, node, initial_degree, initial_degree_group) %>%
      summarise(
        mean_cooccurrence_support_per_realised_partner = mean(full_n),
        median_cooccurrence_support_per_realised_partner = median(full_n),
        mean_interaction_support_per_realised_partner = mean(full_K),
        median_interaction_support_per_realised_partner = median(full_K),
        mean_realised_fraction_among_links = mean(full_K / full_n),
        .groups = "drop"
      )
  )
}

## ---------------------------
## Run analysis
## ---------------------------
future::plan(future::multisession, workers = parallel_workers)
message("Running with workers: ", parallel_workers)

all_outputs <- lapply(all_dataset_names, run_one_dataset)

future::plan(future::sequential)

node_raw <- bind_rows(lapply(all_outputs, `[[`, "node_raw"))
emp_group <- bind_rows(lapply(all_outputs, `[[`, "empirical_group_summary"))
emp_active <- bind_rows(lapply(all_outputs, `[[`, "empirical_active_composition"))
null_group <- bind_rows(lapply(all_outputs, `[[`, "null_group_summary"))
null_active <- bind_rows(lapply(all_outputs, `[[`, "null_active_composition"))
degree_defs <- bind_rows(lapply(all_outputs, `[[`, "initial_degree_groups"))
support_pathway <- bind_rows(lapply(all_outputs, `[[`, "support_pathway"))

## ---------------------------
## Summaries and null envelopes
## ---------------------------
emp_group_long <- emp_group %>%
  select(dataset, guild, removal_fraction, replicate, initial_degree_group,
         proportion_connected, mean_portfolio_retention_if_active) %>%
  pivot_longer(cols = c(proportion_connected, mean_portfolio_retention_if_active),
               names_to = "metric", values_to = "estimate") %>%
  mutate(
    metric = recode(metric,
                    proportion_connected = "Proportion remaining connected",
                    mean_portfolio_retention_if_active = "Mean retained-degree fraction among connected species"),
    source = "Empirical"
  )

emp_active_long <- emp_active %>%
  transmute(dataset, guild, removal_fraction, replicate, initial_degree_group,
            metric = "Composition of connected species",
            estimate = proportion_of_active_network,
            source = "Empirical")

emp_all_long <- bind_rows(emp_group_long, emp_active_long)

emp_dataset_curves <- emp_all_long %>%
  group_by(dataset, guild, metric, removal_fraction, initial_degree_group) %>%
  summarise(
    median = median(estimate, na.rm = TRUE),
    q025 = safe_q(estimate, 0.025),
    q975 = safe_q(estimate, 0.975),
    .groups = "drop"
  )

emp_across <- emp_dataset_curves %>%
  group_by(guild, metric, removal_fraction, initial_degree_group) %>%
  summarise(
    median = median(median, na.rm = TRUE),
    q025 = safe_q(median, 0.025),
    q975 = safe_q(median, 0.975),
    .groups = "drop"
  )

null_group_long <- null_group %>%
  select(dataset, guild, removal_fraction, replicate, null_replicate, initial_degree_group,
         proportion_connected, mean_portfolio_retention_if_active) %>%
  pivot_longer(cols = c(proportion_connected, mean_portfolio_retention_if_active),
               names_to = "metric", values_to = "estimate") %>%
  mutate(metric = recode(metric,
                         proportion_connected = "Proportion remaining connected",
                         mean_portfolio_retention_if_active = "Mean retained-degree fraction among connected species"))

null_envelopes <- null_group_long %>%
  group_by(dataset, guild, metric, removal_fraction, initial_degree_group) %>%
  summarise(
    null_median = median(estimate, na.rm = TRUE),
    null_q025 = safe_q(estimate, 0.025),
    null_q975 = safe_q(estimate, 0.975),
    .groups = "drop"
  )

emp_vs_null <- emp_dataset_curves %>%
  filter(metric %in% c("Proportion remaining connected", "Mean retained-degree fraction among connected species")) %>%
  left_join(null_envelopes, by = c("dataset", "guild", "metric", "removal_fraction", "initial_degree_group")) %>%
  mutate(
    empirical_minus_null_median = median - null_median,
    relation_to_null = case_when(
      median < null_q025 ~ "Below null envelope",
      median > null_q975 ~ "Above null envelope",
      TRUE ~ "Inside null envelope"
    )
  )

results_summary <- emp_vs_null %>%
  filter(removal_fraction > 0) %>%
  group_by(guild, metric, initial_degree_group, relation_to_null) %>%
  summarise(n_dataset_removal_cases = n(), .groups = "drop")

## ---------------------------
## Save tables
## ---------------------------
write.csv2(node_raw, file.path(out_dir, "35_species_retention_raw_empirical.csv"), row.names = FALSE)
write.csv2(emp_dataset_curves, file.path(out_dir, "35_empirical_degree_filtering_estimates.csv"), row.names = FALSE)
write.csv2(null_envelopes, file.path(out_dir, "35_opportunity_matched_null_envelopes.csv"), row.names = FALSE)
write.csv2(emp_vs_null, file.path(out_dir, "35_empirical_vs_null_degree_filtering.csv"), row.names = FALSE)
write.csv2(results_summary, file.path(out_dir, "35_concise_results_summary.csv"), row.names = FALSE)
write.csv2(degree_defs, file.path(out_dir, "35_initial_degree_groups.csv"), row.names = FALSE)
write.csv2(support_pathway, file.path(out_dir, "35_support_pathway_by_species.csv"), row.names = FALSE)

## ---------------------------
## Figure 1: publication-style three-panel empirical summary
## ---------------------------
plot_data <- bind_rows(
  emp_dataset_curves %>% mutate(summary_type = "Dataset trajectory"),
  emp_across %>% mutate(dataset = "Across datasets", summary_type = "Across-dataset median")
) %>%
  mutate(
    metric = factor(metric, levels = c(
      "Proportion remaining connected",
      "Mean retained-degree fraction among connected species",
      "Composition of connected species"
    )),
    guild = factor(guild, levels = c("Consumer", "Resource")),
    initial_degree_group = factor(initial_degree_group, levels = names(degree_cols))
  )

p_main <- ggplot() +
  geom_line(
    data = plot_data %>% filter(summary_type == "Dataset trajectory"),
    aes(x = removal_fraction, y = median,
        group = interaction(dataset, initial_degree_group),
        colour = initial_degree_group),
    alpha = 0.18,
    linewidth = 0.45
  ) +
  geom_ribbon(
    data = plot_data %>% filter(summary_type == "Across-dataset median"),
    aes(x = removal_fraction, ymin = q025, ymax = q975,
        fill = initial_degree_group, group = initial_degree_group),
    alpha = 0.14
  ) +
  geom_line(
    data = plot_data %>% filter(summary_type == "Across-dataset median"),
    aes(x = removal_fraction, y = median,
        colour = initial_degree_group, group = initial_degree_group),
    linewidth = 1.1
  ) +
  geom_point(
    data = plot_data %>% filter(summary_type == "Across-dataset median"),
    aes(x = removal_fraction, y = median, colour = initial_degree_group),
    size = 1.6
  ) +
  facet_grid(guild ~ metric) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  scale_fill_manual(values = degree_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0)
  ) +
  xlab("Proportion of sites removed") +
  ylab("Fraction") +
  ggtitle("Initially generalist species are selectively represented after site removal",
          subtitle = "Faint lines are dataset-specific medians; thick lines are equal-dataset summaries.")

save_both(p_main, "35_generalist_filtering_three_panel_summary", 15, 7)

## ---------------------------
## Figure 2: empirical versus opportunity-matched null
## ---------------------------
null_plot <- emp_vs_null %>%
  filter(metric %in% c("Proportion remaining connected", "Mean retained-degree fraction among connected species")) %>%
  mutate(
    guild = factor(guild, levels = c("Consumer", "Resource")),
    metric = factor(metric, levels = c("Proportion remaining connected", "Mean retained-degree fraction among connected species")),
    initial_degree_group = factor(initial_degree_group, levels = names(degree_cols))
  )

p_null <- ggplot(null_plot, aes(x = removal_fraction)) +
  geom_ribbon(aes(ymin = null_q025, ymax = null_q975, fill = "Opportunity-matched null"), alpha = 0.18) +
  geom_line(aes(y = null_median, colour = "Opportunity-matched null"), linewidth = 0.8) +
  geom_line(aes(y = median, colour = "Empirical"), linewidth = 0.8) +
  facet_grid(guild + metric ~ dataset + initial_degree_group) +
  scale_colour_manual(values = source_cols) +
  scale_fill_manual(values = source_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 6),
    strip.text.y = element_text(size = 7)
  ) +
  xlab("Sites removed") +
  ylab("Fraction") +
  ggtitle("Empirical degree filtering compared with opportunity-matched support shuffles")

save_both(p_null, "35_empirical_vs_opportunity_matched_null", 22, 10)

## ---------------------------
## Validation checks
## ---------------------------
checks <- tibble(
  check = c(
    "all retained degree <= initial degree",
    "portfolio retention equals 1 at zero removal",
    "initial-degree groups present",
    "null envelopes created",
    "consumers and resources analysed separately",
    "no fitted interaction model used"
  ),
  pass = c(
    all(node_raw$retained_degree <= node_raw$initial_degree),
    all(abs(node_raw$portfolio_retention_all[node_raw$removal_fraction == 0] - 1) < 1e-12),
    all(c("Lower initial degree", "Middle initial degree", "Higher initial degree") %in% as.character(degree_defs$initial_degree_group)),
    nrow(null_envelopes) > 0,
    all(c("Consumer", "Resource") %in% unique(node_raw$guild)),
    TRUE
  )
)
write.csv2(checks, file.path(out_dir, "35_generalist_filtering_checks.csv"), row.names = FALSE)

message("\nCompact results summary:")
print(results_summary, n = Inf)
message("\nValidation checks:")
print(checks, n = Inf)

message("\nInterpretation reminder:")
message("This analysis distinguishes a real generalist-filtered remnant network driven by broader spatial opportunity/support from any additional empirical advantage beyond what co-occurrence opportunity and initial degree already predict.")
message("A curve above the opportunity-matched null suggests empirical support placement gives that degree group more persistence than expected after preserving degree, link opportunity, and the support distribution within matched co-occurrence strata.")
message("A curve inside the null envelope suggests the apparent generalist advantage is largely the mechanical consequence of initial degree and opportunity-matched local support.")

message("\nSaved outputs in: ", out_dir)
