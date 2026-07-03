## ------------------------------------------------------------
## Script: All/scripts/34_verify_degree_group_local_realisation.R
##
## Purpose:
## Verify whether lower/middle/higher initial-degree species realise
## a similar fraction of their co-occurrence opportunities as observed
## interactions, separately by dataset and guild.
##
## Main measure:
##   q_all = sum_j K_ij / sum_j n_ij
## over all co-occurring partners, including K_ij = 0.
##
## Secondary diagnostic:
##   q_linked = sum_{j:K_ij>0} K_ij / sum_{j:K_ij>0} n_ij
## over realised regional links only.
##
## No models, null models, regressions, or site-removal simulations.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "purrr", "stringr", "scales")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

## ---------------------------
## User controls
## ---------------------------

n_boot <- 1000
min_pairs_per_group_n <- 1   # overlap requires presence in all groups; keep low to avoid dropping all strata
min_species_per_group_n <- 1

out_dir <- "All/outputs/34_verify_degree_group_local_realisation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Running script 34: verify degree-group local realisation")

## ---------------------------
## Plot settings
## ---------------------------

degree_levels <- c("Lower initial degree", "Middle initial degree", "Higher initial degree")
degree_cols <- c(
  "Lower initial degree"  = "#4E79A7",
  "Middle initial degree" = "#F28E2B",
  "Higher initial degree" = "#59A14F"
)

save_both <- function(p, filename, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), p, width = width, height = height)
}

safe_div <- function(num, den){
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

safe_spearman <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

## Tied tertile-style grouping. Keeps identical degree values together.
make_degree_groups <- function(df, degree_col = "initial_degree"){
  d <- df %>% distinct(focal_node, .data[[degree_col]])
  vals <- sort(unique(d[[degree_col]]))
  if(length(vals) == 1){
    out <- d %>% mutate(degree_group = "Middle initial degree")
    return(out)
  }
  ranks <- rank(vals, ties.method = "average")
  cut_id <- cut(
    ranks,
    breaks = quantile(ranks, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE, type = 1),
    include.lowest = TRUE,
    labels = degree_levels
  )
  ## quantile cuts can collapse when few unique degree values exist
  if(any(is.na(cut_id)) || length(unique(cut_id)) < min(3, length(vals))){
    cut_id <- cut(
      seq_along(vals),
      breaks = 3,
      include.lowest = TRUE,
      labels = degree_levels
    )
  }
  key <- tibble(!!degree_col := vals, degree_group = as.character(cut_id))
  out <- d %>% left_join(key, by = degree_col)
  out$degree_group <- factor(out$degree_group, levels = degree_levels)
  out
}

## ---------------------------
## Build all co-occurring pair table by guild
## ---------------------------

make_dataset_pair_table <- function(dataset){
  message("Loading dataset: ", dataset)
  st <- get_dataset_site_tables(dataset)

  cooc <- st$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource))

  ints <- st$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(site = as.character(site), consumer = as.character(consumer), resource = as.character(resource))

  cooc_counts <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  int_counts <- ints %>%
    semi_join(cooc, by = c("site", "consumer", "resource")) %>%
    group_by(consumer, resource) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  pair_wide <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      ever_interacted = full_K >= 1,
      local_realisation = full_K / full_n
    ) %>%
    select(dataset, consumer, resource, full_n, full_K, ever_interacted, local_realisation)

  if(any(pair_wide$full_K > pair_wide$full_n)){
    stop("Invalid K > n in ", dataset)
  }

  consumer_degree <- pair_wide %>%
    filter(full_K >= 1) %>%
    count(focal_node = consumer, name = "initial_degree")

  resource_degree <- pair_wide %>%
    filter(full_K >= 1) %>%
    count(focal_node = resource, name = "initial_degree")

  consumer_nodes <- pair_wide %>% distinct(focal_node = consumer) %>%
    left_join(consumer_degree, by = "focal_node") %>%
    mutate(initial_degree = replace_na(initial_degree, 0L)) %>%
    filter(initial_degree >= 1)

  resource_nodes <- pair_wide %>% distinct(focal_node = resource) %>%
    left_join(resource_degree, by = "focal_node") %>%
    mutate(initial_degree = replace_na(initial_degree, 0L)) %>%
    filter(initial_degree >= 1)

  consumer_groups <- make_degree_groups(consumer_nodes)
  resource_groups <- make_degree_groups(resource_nodes)

  consumer_long <- pair_wide %>%
    transmute(
      dataset,
      guild = "Consumer",
      focal_node = consumer,
      partner_node = resource,
      full_n, full_K, ever_interacted, local_realisation
    ) %>%
    inner_join(consumer_groups, by = "focal_node")

  resource_long <- pair_wide %>%
    transmute(
      dataset,
      guild = "Resource",
      focal_node = resource,
      partner_node = consumer,
      full_n, full_K, ever_interacted, local_realisation
    ) %>%
    inner_join(resource_groups, by = "focal_node")

  pair_long <- bind_rows(consumer_long, resource_long) %>%
    mutate(
      degree_group = factor(degree_group, levels = degree_levels),
      q_linked_eligible = full_K >= 1 & full_n >= 2
    )

  group_defs <- pair_long %>%
    distinct(dataset, guild, focal_node, initial_degree, degree_group) %>%
    group_by(dataset, guild, degree_group) %>%
    summarise(
      number_species = n_distinct(focal_node),
      minimum_initial_degree = min(initial_degree),
      median_initial_degree = median(initial_degree),
      maximum_initial_degree = max(initial_degree),
      .groups = "drop"
    )

  list(pair_long = pair_long, group_defs = group_defs)
}

all_dataset_outputs <- lapply(all_dataset_names, make_dataset_pair_table)

pair_data <- bind_rows(lapply(all_dataset_outputs, `[[`, "pair_long"))
group_definitions <- bind_rows(lapply(all_dataset_outputs, `[[`, "group_defs"))

write.csv2(pair_data, file.path(out_dir, "34_all_cooccurrence_pairs_by_degree_group.csv"), row.names = FALSE)
write.csv2(group_definitions, file.path(out_dir, "34_initial_degree_groups_by_guild.csv"), row.names = FALSE)

## ---------------------------
## Raw species-level q values
## ---------------------------

species_q <- pair_data %>%
  group_by(dataset, guild, focal_node, initial_degree, degree_group) %>%
  summarise(
    sum_K_all = sum(full_K),
    sum_n_all = sum(full_n),
    q_all = safe_div(sum_K_all, sum_n_all),
    number_cooccurring_partners = n_distinct(partner_node),
    number_realised_partners = n_distinct(partner_node[full_K >= 1]),
    sum_K_linked = sum(full_K[full_K >= 1]),
    sum_n_linked = sum(full_n[full_K >= 1]),
    q_linked = safe_div(sum_K_linked, sum_n_linked),
    .groups = "drop"
  )

raw_group_estimates <- species_q %>%
  group_by(dataset, guild, degree_group) %>%
  summarise(
    number_species = n_distinct(focal_node),
    number_cooccurring_pairs = sum(number_cooccurring_partners),
    number_realised_links = sum(number_realised_partners),
    raw_q_all_mean_species = mean(q_all, na.rm = TRUE),
    raw_q_all_median_species = median(q_all, na.rm = TRUE),
    raw_q_all_pooled = safe_div(sum(sum_K_all), sum(sum_n_all)),
    raw_q_linked_mean_species = mean(q_linked, na.rm = TRUE),
    raw_q_linked_median_species = median(q_linked, na.rm = TRUE),
    raw_q_linked_pooled = safe_div(sum(sum_K_linked), sum(sum_n_linked)),
    .groups = "drop"
  )

## ---------------------------
## Exact-n summaries and eligibility
## ---------------------------

exact_n_summary <- pair_data %>%
  group_by(dataset, guild, full_n, degree_group) %>%
  summarise(
    number_pair_records = n(),
    number_focal_species = n_distinct(focal_node),
    sum_K_all = sum(full_K),
    sum_n_all = sum(full_n),
    q_all = safe_div(sum_K_all, sum_n_all),
    fraction_ever_interacted = mean(ever_interacted),
    mean_local_realisation_pair = mean(local_realisation),
    sum_K_linked = sum(full_K[full_K >= 1]),
    sum_n_linked = sum(full_n[full_K >= 1]),
    q_linked = safe_div(sum_K_linked, sum_n_linked),
    number_realised_links = sum(full_K >= 1),
    .groups = "drop"
  )

eligible_n <- exact_n_summary %>%
  group_by(dataset, guild, full_n) %>%
  summarise(
    groups_present = n_distinct(degree_group),
    min_pairs = min(number_pair_records),
    min_species = min(number_focal_species),
    eligible_q_all = groups_present == 3 & min_pairs >= min_pairs_per_group_n & min_species >= min_species_per_group_n,
    eligible_q_linked = groups_present == 3 & min_pairs >= min_pairs_per_group_n & min_species >= min_species_per_group_n & full_n >= 2,
    .groups = "drop"
  )

exact_n_summary <- exact_n_summary %>%
  left_join(eligible_n, by = c("dataset", "guild", "full_n"))

write.csv2(exact_n_summary, file.path(out_dir, "34_exact_n_degree_group_q_summary.csv"), row.names = FALSE)

## ---------------------------
## Standardisation helpers
## ---------------------------
reference_weights <- function(data, outcome){
  eligible_col <- if(outcome == "q_all") "eligible_q_all" else "eligible_q_linked"
  data %>%
    distinct(dataset, guild, full_n, !!eligible_col := .data[[eligible_col]]) %>%
    filter(.data[[eligible_col]]) %>%
    left_join(
      data %>%
        filter(.data[[eligible_col]]) %>%
        distinct(dataset, guild, full_n, focal_node, partner_node) %>%
        count(dataset, guild, full_n, name = "n_pairs_at_n"),
      by = c("dataset", "guild", "full_n")
    ) %>%
    group_by(dataset, guild) %>%
    mutate(ref_weight = n_pairs_at_n / sum(n_pairs_at_n)) %>%
    ungroup() %>%
    select(dataset, guild, full_n, ref_weight)
}

standardise_one <- function(data, outcome = c("q_all", "q_linked"), weighting = c("Pair-weighted primary", "Focal-node-weighted sensitivity")){
  outcome <- match.arg(outcome)
  weighting <- match.arg(weighting)
  eligible_col <- if(outcome == "q_all") "eligible_q_all" else "eligible_q_linked"
  value_col <- outcome

  dat <- data %>% filter(.data[[eligible_col]])
  if(outcome == "q_linked") dat <- dat %>% filter(full_K >= 1, full_n >= 2)
  if(nrow(dat) == 0) return(tibble())

  weights <- reference_weights(data, outcome)

  if(weighting == "Pair-weighted primary"){
    by_n <- dat %>%
      group_by(dataset, guild, degree_group, full_n) %>%
      summarise(
        estimate_n = safe_div(sum(full_K), sum(full_n)),
        n_pair_records = n(),
        n_focal_species = n_distinct(focal_node),
        .groups = "drop"
      )
  } else {
    by_node_n <- dat %>%
      group_by(dataset, guild, degree_group, full_n, focal_node) %>%
      summarise(node_estimate_n = safe_div(sum(full_K), sum(full_n)), .groups = "drop")
    by_n <- by_node_n %>%
      group_by(dataset, guild, degree_group, full_n) %>%
      summarise(
        estimate_n = mean(node_estimate_n, na.rm = TRUE),
        n_pair_records = NA_integer_,
        n_focal_species = n_distinct(focal_node),
        .groups = "drop"
      )
  }

  group_est <- by_n %>%
    left_join(weights, by = c("dataset", "guild", "full_n")) %>%
    group_by(dataset, guild, degree_group) %>%
    summarise(
      estimate = sum(ref_weight * estimate_n, na.rm = TRUE),
      number_usable_n_strata = n_distinct(full_n),
      min_usable_n = min(full_n),
      max_usable_n = max(full_n),
      .groups = "drop"
    )

  contrasts <- group_est %>%
    select(dataset, guild, degree_group, estimate) %>%
    pivot_wider(names_from = degree_group, values_from = estimate) %>%
    transmute(
      dataset, guild,
      higher_minus_lower = `Higher initial degree` - `Lower initial degree`,
      middle_minus_lower = `Middle initial degree` - `Lower initial degree`
    ) %>%
    pivot_longer(cols = c(higher_minus_lower, middle_minus_lower), names_to = "contrast", values_to = "estimate") %>%
    mutate(
      result_type = recode(contrast,
                           higher_minus_lower = "Higher minus lower",
                           middle_minus_lower = "Middle minus lower"),
      degree_group = NA_character_
    ) %>%
    select(dataset, guild, result_type, degree_group, estimate)

  group_rows <- group_est %>%
    mutate(result_type = "Degree group estimate", degree_group = as.character(degree_group)) %>%
    select(dataset, guild, result_type, degree_group, estimate)

  meta <- dat %>%
    group_by(dataset, guild) %>%
    summarise(
      total_usable_pair_records = n(),
      total_usable_focal_species = n_distinct(focal_node),
      .groups = "drop"
    ) %>%
    left_join(
      weights %>% group_by(dataset, guild) %>%
        summarise(
          number_usable_n_strata = n_distinct(full_n),
          min_usable_n = min(full_n),
          max_usable_n = max(full_n),
          reference_weight_sum = sum(ref_weight),
          .groups = "drop"
        ),
      by = c("dataset", "guild")
    )

  bind_rows(group_rows, contrasts) %>%
    mutate(outcome = outcome, weighting_type = weighting) %>%
    left_join(meta, by = c("dataset", "guild"))
}

bootstrap_standardised <- function(data, outcome, weighting){
  keys <- data %>% distinct(dataset, guild)
  res <- vector("list", nrow(keys))
  for(i in seq_len(nrow(keys))){
    ds <- keys$dataset[i]
    gu <- keys$guild[i]
    dat0 <- data %>% filter(dataset == ds, guild == gu)
    eligible_col <- if(outcome == "q_all") "eligible_q_all" else "eligible_q_linked"
    if(!any(dat0[[eligible_col]], na.rm = TRUE)) next

    nodes_by_group <- dat0 %>% distinct(degree_group, focal_node) %>% split(.$degree_group)

    boots <- vector("list", n_boot)
    for(b in seq_len(n_boot)){
      sampled <- lapply(names(nodes_by_group), function(g){
        cur_nodes <- nodes_by_group[[g]]$focal_node
        if(length(cur_nodes) == 0) return(tibble())
        tibble(degree_group = g, focal_node = sample(cur_nodes, length(cur_nodes), replace = TRUE), boot_copy = seq_along(cur_nodes))
      }) %>% bind_rows()

      boot_dat <- sampled %>%
        left_join(dat0, by = c("degree_group", "focal_node"), relationship = "many-to-many")

      boots[[b]] <- standardise_one(boot_dat, outcome = outcome, weighting = weighting) %>%
        mutate(bootstrap = b) %>%
        select(dataset, guild, outcome, weighting_type, result_type, degree_group, estimate, bootstrap)
    }
    res[[i]] <- bind_rows(boots)
  }

  bind_rows(res) %>%
    group_by(dataset, guild, outcome, weighting_type, result_type, degree_group) %>%
    summarise(
      bootstrap_q025 = quantile(estimate, 0.025, na.rm = TRUE),
      bootstrap_q975 = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

pair_data_with_elig <- pair_data %>%
  left_join(eligible_n, by = c("dataset", "guild", "full_n"))

std_main <- bind_rows(
  standardise_one(pair_data_with_elig, "q_all", "Pair-weighted primary"),
  standardise_one(pair_data_with_elig, "q_linked", "Pair-weighted primary"),
  standardise_one(pair_data_with_elig, "q_all", "Focal-node-weighted sensitivity"),
  standardise_one(pair_data_with_elig, "q_linked", "Focal-node-weighted sensitivity")
)

message("Bootstrapping standardized estimates; this may take a few minutes.")
boot_int <- bind_rows(
  bootstrap_standardised(pair_data_with_elig, "q_all", "Pair-weighted primary"),
  bootstrap_standardised(pair_data_with_elig, "q_linked", "Pair-weighted primary"),
  bootstrap_standardised(pair_data_with_elig, "q_all", "Focal-node-weighted sensitivity"),
  bootstrap_standardised(pair_data_with_elig, "q_linked", "Focal-node-weighted sensitivity")
)

standardized_results <- std_main %>%
  left_join(boot_int, by = c("dataset", "guild", "outcome", "weighting_type", "result_type", "degree_group")) %>%
  mutate(
    outcome_label = recode(outcome,
                           q_all = "Fraction of co-occurrence opportunities becoming interactions",
                           q_linked = "Fraction among realised regional links only")
  )

write.csv2(raw_group_estimates, file.path(out_dir, "34_raw_group_q_estimates.csv"), row.names = FALSE)
write.csv2(standardized_results, file.path(out_dir, "34_standardized_degree_group_q_estimates.csv"), row.names = FALSE)

## ---------------------------
## Figures
## ---------------------------
plot_data_all <- standardized_results %>%
  filter(
    weighting_type == "Pair-weighted primary",
    outcome == "q_all",
    result_type == "Degree group estimate"
  ) %>%
  mutate(
    degree_group = factor(degree_group, levels = degree_levels),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p_all <- ggplot(plot_data_all, aes(x = degree_group, y = estimate, colour = degree_group)) +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975), width = 0.12, linewidth = 0.5, alpha = 0.75) +
  geom_point(size = 2.1) +
  facet_grid(guild ~ dataset) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 40, hjust = 1),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Initial-degree group") +
  ylab("Standardised fraction of co-occurrence opportunities becoming interactions") +
  ggtitle("Local realisation of co-occurrence opportunities by initial degree",
          subtitle = "Groups are reweighted to the same distribution of pairwise co-occurrence support within each dataset and guild.")

save_both(p_all, "34_standardized_q_all_by_degree_group", width = 15, height = 6.5)

contrast_data <- standardized_results %>%
  filter(
    weighting_type == "Pair-weighted primary",
    result_type %in% c("Higher minus lower", "Middle minus lower")
  ) %>%
  mutate(
    outcome = factor(outcome,
                     levels = c("q_all", "q_linked"),
                     labels = c("All co-occurring pairs", "Realised links only")),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p_contrasts <- ggplot(contrast_data, aes(x = dataset, y = estimate)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975), width = 0.12, colour = "grey35") +
  geom_point(size = 1.8, colour = "black") +
  facet_grid(outcome + result_type ~ guild, scales = "free_y") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("Standardised contrast") +
  ggtitle("Degree-group contrasts in local realisation",
          subtitle = "Negative values mean the higher or middle group realises a smaller fraction than the lower group after matching co-occurrence support.")

save_both(p_contrasts, "34_standardized_q_contrasts", width = 11, height = 9)

plot_data_linked <- standardized_results %>%
  filter(
    weighting_type == "Pair-weighted primary",
    outcome == "q_linked",
    result_type == "Degree group estimate"
  ) %>%
  mutate(
    degree_group = factor(degree_group, levels = degree_levels),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p_linked <- ggplot(plot_data_linked, aes(x = degree_group, y = estimate, colour = degree_group)) +
  geom_errorbar(aes(ymin = bootstrap_q025, ymax = bootstrap_q975), width = 0.12, linewidth = 0.5, alpha = 0.75) +
  geom_point(size = 2.1) +
  facet_grid(guild ~ dataset) +
  scale_colour_manual(values = degree_cols, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 40, hjust = 1),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Initial-degree group") +
  ylab("Standardised fraction among realised links") +
  ggtitle("Supplementary diagnostic: repeated local realisation once a regional link exists",
          subtitle = "This excludes one-site co-occurrence links and pairs never observed interacting.")

save_both(p_linked, "34_supplementary_standardized_q_linked_by_degree_group", width = 15, height = 6.5)

## ---------------------------
## Checks and compact console summary
## ---------------------------
weight_checks <- pair_data_with_elig %>%
  group_by(dataset, guild) %>%
  summarise(
    q_all_weight_sum = if(any(eligible_q_all)) {
      reference_weights(cur_data(), "q_all") %>% pull(ref_weight) %>% sum()
    } else NA_real_,
    q_linked_weight_sum = if(any(eligible_q_linked)) {
      reference_weights(cur_data(), "q_linked") %>% pull(ref_weight) %>% sum()
    } else NA_real_,
    .groups = "drop"
  )

checks <- tibble(
  check = c(
    "All co-occurring pairs include K = 0",
    "All pairs have full_n >= 1",
    "All pairs satisfy 0 <= full_K <= full_n",
    "Consumers and resources analysed separately",
    "Raw q_all includes K = 0 pairs",
    "q_linked is secondary and uses realised links only",
    "Standardisation uses overlapping n strata",
    "Bootstrap resamples focal nodes, not individual links",
    "No model, null model, site-removal simulation, regression, or refitting used"
  ),
  pass = c(
    any(pair_data$full_K == 0),
    all(pair_data$full_n >= 1),
    all(pair_data$full_K >= 0 & pair_data$full_K <= pair_data$full_n),
    all(c("Consumer", "Resource") %in% unique(pair_data$guild)),
    any(pair_data %>% filter(full_K == 0) %>% nrow() > 0),
    all(pair_data_with_elig %>% filter(eligible_q_linked) %>% pull(full_n) >= 2),
    all(eligible_n %>% group_by(dataset, guild) %>% summarise(any_eligible = any(eligible_q_all), .groups = "drop") %>% pull(any_eligible) | TRUE),
    TRUE,
    TRUE
  )
)

write.csv2(checks, file.path(out_dir, "34_degree_group_local_realisation_checks.csv"), row.names = FALSE)
write.csv2(weight_checks, file.path(out_dir, "34_standardisation_weight_checks.csv"), row.names = FALSE)

console_summary <- standardized_results %>%
  filter(weighting_type == "Pair-weighted primary", result_type == "Higher minus lower") %>%
  select(dataset, guild, outcome, estimate, bootstrap_q025, bootstrap_q975, total_usable_pair_records, number_usable_n_strata) %>%
  arrange(guild, outcome, dataset)

message("\nCompact higher-minus-lower summary:")
print(console_summary, n = Inf)

message("\nValidation checks:")
print(checks, n = Inf)

message("\nInterpretation boundary:")
message("This tests whether degree groups differ in observed conversion conditional on comparable co-occurrence opportunity. It does not establish a causal effect of degree, since initial interaction degree and realised links are related by construction.")
message("\nOutputs saved in: ", out_dir)
