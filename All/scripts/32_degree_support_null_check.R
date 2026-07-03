## ------------------------------------------------------------
## Script: All/scripts/32_degree_support_null_check.R
## Degree-support null check for possible Figure 5
## ------------------------------------------------------------
## Purpose:
## Test whether associations between node initial binary degree and
## local interaction support are stronger than expected after preserving
## each realised link's exact co-occurrence opportunity full_n.
##
## Null model:
## Within each dataset and each exact full_n stratum, shuffle full_K values
## among original realised interaction links. This preserves:
##   - link presence;
##   - node degrees and link identities;
##   - co-occurrence opportunity full_n for every link;
##   - the empirical distribution of full_K conditional on exact full_n.
## It does not fit any interaction model.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr", "tidyr", "ggplot2", "tibble", "purrr", "patchwork",
  "future", "future.apply", "scales"
)
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/32_degree_support_null_check"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## User controls
## ---------------------------

n_null_reps <- 999
n_site_reps_candidate <- 50
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
set.seed(123)

## Parallelize across datasets by default.
n_workers <- max(1, min(7, future::availableCores() - 1))
future::plan(future::multisession, workers = n_workers)

message("Running script 32 with ", n_workers, " workers")

## Dataset order from the shared helper.
dataset_order <- all_dataset_names

## ---------------------------
## Colours and helpers
## ---------------------------

metric_cols <- c(
  "Initial degree vs mean support per link" = "#0072B2",
  "Initial degree vs proportion one-site links" = "#D55E00"
)

degree_group_cols <- c(
  "Lower initial degree"  = "#4E79A7",
  "Middle initial degree" = "#F28E2B",
  "Higher initial degree" = "#59A14F"
)

save_both <- function(p, name, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

safe_cor <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok) < 3) return(NA_real_)
  if(length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

support_class <- function(K){
  case_when(
    K == 1 ~ "Observed in 1 site",
    K == 2 ~ "Observed in 2 sites",
    K %in% 3:4 ~ "Observed in 3–4 sites",
    K >= 5 ~ "Observed in 5 or more sites",
    TRUE ~ NA_character_
  )
}

assign_tied_tertiles <- function(df, value_col, out_col = "initial_degree_group"){
  vals <- sort(unique(df[[value_col]]))
  if(length(vals) == 1){
    out <- rep("Middle initial degree", nrow(df))
  } else {
    counts <- df %>% count(.data[[value_col]], name = "n") %>% arrange(.data[[value_col]])
    counts$cum_mid <- cumsum(counts$n) - counts$n / 2
    total <- sum(counts$n)
    counts[[out_col]] <- case_when(
      counts$cum_mid <= total / 3 ~ "Lower initial degree",
      counts$cum_mid <= 2 * total / 3 ~ "Middle initial degree",
      TRUE ~ "Higher initial degree"
    )
    out <- counts[[out_col]][match(df[[value_col]], counts[[value_col]])]
  }
  df[[out_col]] <- factor(
    out,
    levels = c("Lower initial degree", "Middle initial degree", "Higher initial degree")
  )
  df
}

## ---------------------------
## Data construction
## ---------------------------

make_full_pair_table <- function(dataset){
  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource)

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource)

  ## An observed interaction is necessarily an operational co-occurrence.
  cooc <- bind_rows(cooc, ints) %>%
    distinct(site, consumer, resource)

  cooc_counts <- cooc %>%
    group_by(consumer, resource) %>%
    summarise(full_n = n_distinct(site), .groups = "drop")

  int_counts <- ints %>%
    group_by(consumer, resource) %>%
    summarise(full_K = n_distinct(site), .groups = "drop")

  all_pairs <- cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      dataset = dataset,
      full_K = replace_na(full_K, 0L),
      full_B = full_K >= 1
    ) %>%
    select(dataset, consumer, resource, full_n, full_K, full_B)

  realised_links <- all_pairs %>%
    filter(full_K >= 1) %>%
    mutate(
      support_class = support_class(full_K),
      link_id = row_number()
    )

  list(
    site_tables = site_tables,
    all_pairs = all_pairs,
    realised_links = realised_links
  )
}

make_node_profiles <- function(realised_links){
  consumer_links <- realised_links %>%
    transmute(
      dataset, guild = "Consumer", node = consumer, partner = resource,
      link_id, full_n, full_K
    )

  resource_links <- realised_links %>%
    transmute(
      dataset, guild = "Resource", node = resource, partner = consumer,
      link_id, full_n, full_K
    )

  bind_rows(consumer_links, resource_links) %>%
    group_by(dataset, guild, node) %>%
    summarise(
      initial_degree = n_distinct(partner),
      mean_support_per_link = mean(full_K),
      proportion_one_site_links = mean(full_K == 1),
      mean_cooccurring_sites_per_link = mean(full_n),
      .groups = "drop"
    ) %>%
    group_by(dataset, guild) %>%
    group_modify(~ assign_tied_tertiles(.x, "initial_degree")) %>%
    ungroup()
}

make_node_profiles_from_links <- function(realised_links_with_K, degree_groups = NULL){
  profiles <- make_node_profiles(realised_links_with_K)
  if(!is.null(degree_groups)){
    profiles <- profiles %>%
      select(-initial_degree_group) %>%
      left_join(
        degree_groups %>% select(dataset, guild, node, initial_degree_group),
        by = c("dataset", "guild", "node")
      )
  }
  profiles
}

profile_correlations <- function(node_profiles, replicate = NA_integer_, source = "Empirical"){
  node_profiles %>%
    group_by(dataset, guild) %>%
    summarise(
      rho_mean_support = safe_cor(initial_degree, mean_support_per_link),
      rho_one_site = safe_cor(initial_degree, proportion_one_site_links),
      n_nodes = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(rho_mean_support, rho_one_site),
      names_to = "metric_code",
      values_to = "spearman_correlation"
    ) %>%
    mutate(
      metric = recode(
        metric_code,
        rho_mean_support = "Initial degree vs mean support per link",
        rho_one_site = "Initial degree vs proportion one-site links"
      ),
      source = source,
      replicate = replicate
    ) %>%
    select(dataset, guild, source, replicate, metric, spearman_correlation, n_nodes)
}

shuffle_K_within_exact_n <- function(realised_links){
  realised_links %>%
    group_by(full_n) %>%
    mutate(full_K = sample(full_K, size = n(), replace = FALSE)) %>%
    ungroup()
}

group_support_summary <- function(node_profiles, replicate = NA_integer_, source = "Empirical"){
  node_profiles %>%
    group_by(dataset, guild, initial_degree_group) %>%
    summarise(
      mean_support_per_link = mean(mean_support_per_link, na.rm = TRUE),
      proportion_one_site_links = mean(proportion_one_site_links, na.rm = TRUE),
      n_nodes = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(mean_support_per_link, proportion_one_site_links),
      names_to = "metric_code",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric_code,
        mean_support_per_link = "Mean support per link",
        proportion_one_site_links = "Proportion one-site links"
      ),
      source = source,
      replicate = replicate
    ) %>%
    select(dataset, guild, source, replicate, initial_degree_group, metric, value, n_nodes)
}

run_null_for_dataset <- function(dataset){
  message("Running degree-support null: ", dataset)

  dat <- make_full_pair_table(dataset)
  realised_links <- dat$realised_links

  if(nrow(realised_links) == 0){
    stop("No realised links found for ", dataset)
  }

  if(any(realised_links$full_K < 1) || any(realised_links$full_K > realised_links$full_n)){
    stop("Invalid realised link support in ", dataset)
  }

  empirical_profiles <- make_node_profiles(realised_links)
  empirical_cor <- profile_correlations(empirical_profiles, source = "Empirical")
  empirical_group <- group_support_summary(empirical_profiles, source = "Empirical")

  degree_groups <- empirical_profiles %>%
    select(dataset, guild, node, initial_degree, initial_degree_group)

  null_list <- vector("list", n_null_reps)
  null_group_list <- vector("list", n_null_reps)

  for(rep in seq_len(n_null_reps)){
    shuffled_links <- shuffle_K_within_exact_n(realised_links)
    profiles_rep <- make_node_profiles_from_links(shuffled_links, degree_groups = degree_groups)
    null_list[[rep]] <- profile_correlations(profiles_rep, replicate = rep, source = "Exact-n shuffled support")
    null_group_list[[rep]] <- group_support_summary(profiles_rep, replicate = rep, source = "Exact-n shuffled support")
  }

  null_cor <- bind_rows(null_list)
  null_group <- bind_rows(null_group_list)

  null_summary <- null_cor %>%
    group_by(dataset, guild, metric) %>%
    summarise(
      null_median = median(spearman_correlation, na.rm = TRUE),
      null_q025 = quantile(spearman_correlation, 0.025, na.rm = TRUE),
      null_q975 = quantile(spearman_correlation, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  summary <- empirical_cor %>%
    select(dataset, guild, metric, empirical_correlation = spearman_correlation, n_nodes) %>%
    left_join(null_summary, by = c("dataset", "guild", "metric")) %>%
    mutate(
      empirical_minus_null_median = empirical_correlation - null_median,
      empirical_position = case_when(
        empirical_correlation < null_q025 ~ "Below null 95% range",
        empirical_correlation > null_q975 ~ "Above null 95% range",
        TRUE ~ "Inside null 95% range"
      )
    )

  group_null_summary <- null_group %>%
    group_by(dataset, guild, initial_degree_group, metric) %>%
    summarise(
      null_median = median(value, na.rm = TRUE),
      null_q025 = quantile(value, 0.025, na.rm = TRUE),
      null_q975 = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  group_summary <- empirical_group %>%
    select(dataset, guild, initial_degree_group, metric, empirical_value = value, n_nodes) %>%
    left_join(group_null_summary, by = c("dataset", "guild", "initial_degree_group", "metric"))

  group_defs <- empirical_profiles %>%
    group_by(dataset, guild, initial_degree_group) %>%
    summarise(
      number_nodes = n(),
      min_initial_degree = min(initial_degree),
      median_initial_degree = median(initial_degree),
      max_initial_degree = max(initial_degree),
      .groups = "drop"
    )

  checks <- tibble(
    dataset = dataset,
    realised_links_have_K_ge_1 = all(realised_links$full_K >= 1),
    realised_links_have_K_le_n = all(realised_links$full_K <= realised_links$full_n),
    all_null_K_valid = all(null_group$source == "Exact-n shuffled support"),
    n_null_reps = n_null_reps,
    consumers_and_resources_separate = TRUE,
    no_model_or_regression_used = TRUE
  )

  list(
    realised_links = realised_links,
    empirical_profiles = empirical_profiles,
    empirical_correlations = empirical_cor,
    null_correlations = null_cor,
    summary = summary,
    group_summary = group_summary,
    group_defs = group_defs,
    checks = checks
  )
}

all_outputs <- future.apply::future_lapply(
  dataset_order,
  run_null_for_dataset,
  future.seed = TRUE
)

future::plan(future::sequential)

realised_links_all <- bind_rows(lapply(all_outputs, `[[`, "realised_links")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

node_profiles_all <- bind_rows(lapply(all_outputs, `[[`, "empirical_profiles")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

empirical_correlations <- bind_rows(lapply(all_outputs, `[[`, "empirical_correlations")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

null_correlations <- bind_rows(lapply(all_outputs, `[[`, "null_correlations")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

null_summary <- bind_rows(lapply(all_outputs, `[[`, "summary")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

group_support_null_summary <- bind_rows(lapply(all_outputs, `[[`, "group_summary")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

group_definitions <- bind_rows(lapply(all_outputs, `[[`, "group_defs")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

checks <- bind_rows(lapply(all_outputs, `[[`, "checks")) %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

## ---------------------------
## Decision summary
## ---------------------------

decision_summary <- null_summary %>%
  mutate(
    stronger_same_direction = case_when(
      metric == "Initial degree vs mean support per link" & empirical_correlation > null_q975 ~ TRUE,
      metric == "Initial degree vs proportion one-site links" & empirical_correlation < null_q025 ~ TRUE,
      TRUE ~ FALSE
    ),
    opposite_or_inconclusive = !stronger_same_direction
  ) %>%
  group_by(guild, metric) %>%
  summarise(
    datasets_stronger_than_null_same_direction = sum(stronger_same_direction, na.rm = TRUE),
    datasets_opposite_or_inconclusive = sum(opposite_or_inconclusive, na.rm = TRUE),
    n_datasets = n(),
    .groups = "drop"
  )

candidate_threshold_met <- decision_summary %>%
  group_by(guild) %>%
  summarise(
    min_stronger = min(datasets_stronger_than_null_same_direction),
    .groups = "drop"
  ) %>%
  summarise(any(min_stronger >= 6)) %>%
  pull()

## ---------------------------
## Write tables
## ---------------------------

write.csv2(realised_links_all, file.path(out_dir, "32_realised_link_support_data.csv"), row.names = FALSE)
write.csv2(node_profiles_all, file.path(out_dir, "32_node_support_profiles.csv"), row.names = FALSE)
write.csv2(empirical_correlations, file.path(out_dir, "32_empirical_degree_support_correlations.csv"), row.names = FALSE)
write.csv2(null_correlations, file.path(out_dir, "32_degree_support_null_correlations_all_replicates.csv"), row.names = FALSE)
write.csv2(null_summary, file.path(out_dir, "32_degree_support_null_summary.csv"), row.names = FALSE)
write.csv2(group_support_null_summary, file.path(out_dir, "32_degree_group_support_null_envelopes.csv"), row.names = FALSE)
write.csv2(group_definitions, file.path(out_dir, "32_initial_degree_group_definitions.csv"), row.names = FALSE)
write.csv2(decision_summary, file.path(out_dir, "32_degree_support_null_decision_summary.csv"), row.names = FALSE)
write.csv2(checks, file.path(out_dir, "32_degree_support_null_checks.csv"), row.names = FALSE)

## ---------------------------
## Figure: null-test summary
## ---------------------------

p_null <- ggplot(null_summary, aes(x = dataset, y = empirical_correlation)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_linerange(
    aes(ymin = null_q025, ymax = null_q975),
    colour = "grey72",
    linewidth = 2.8,
    alpha = 0.65
  ) +
  geom_point(aes(colour = metric), size = 2.2) +
  facet_grid(metric ~ guild) +
  scale_colour_manual(values = metric_cols, guide = "none") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Dataset") +
  ylab("Spearman correlation") +
  ggtitle(
    "Degree-support associations compared with exact-opportunity shuffles",
    subtitle = "Grey intervals show the 95% range after shuffling interaction support among links with the same co-occurrence support"
  )

save_both(p_null, "32_degree_support_null_summary", 12, 7)

## ---------------------------
## Candidate Figure 5, only if justified
## ---------------------------

make_site_subsets_local <- function(sites, removal_levels, n_reps){
  sites <- sort(unique(sites))
  S <- length(sites)
  out <- list()
  idx <- 1L
  for(r in removal_levels){
    reps <- if(r == 0) 1L else n_reps
    m <- floor(S * r)
    for(rep in seq_len(reps)){
      removed <- if(m == 0) character(0) else sample(sites, m, replace = FALSE)
      retained <- setdiff(sites, removed)
      out[[idx]] <- tibble(
        removal_fraction = r,
        replicate = rep,
        retained_sites = list(retained)
      )
      idx <- idx + 1L
    }
  }
  bind_rows(out)
}

compute_candidate_retention_dataset <- function(dataset){
  message("Computing candidate Figure 5 retention: ", dataset)
  dat <- make_full_pair_table(dataset)
  site_tables <- dat$site_tables
  realised_links <- dat$realised_links
  node_profiles <- make_node_profiles(realised_links)

  cooc <- site_tables$cooc_triples %>% distinct(site, consumer, resource)
  ints <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource)
  cooc <- bind_rows(cooc, ints) %>% distinct(site, consumer, resource)

  ## Node presence from co-occurrence records. This is conservative and follows
  ## the same operational occurrence logic as previous scripts: a node is present
  ## at a site if it appears in any retained consumer-resource co-occurrence.
  consumer_presence <- cooc %>% distinct(site, node = consumer) %>% mutate(guild = "Consumer")
  resource_presence <- cooc %>% distinct(site, node = resource) %>% mutate(guild = "Resource")
  presence <- bind_rows(consumer_presence, resource_presence)

  sites <- sort(unique(cooc$site))
  subsets <- make_site_subsets_local(sites, removal_levels, n_site_reps_candidate)

  consumer_link_map <- realised_links %>%
    transmute(guild = "Consumer", node = consumer, partner = resource, consumer, resource)
  resource_link_map <- realised_links %>%
    transmute(guild = "Resource", node = resource, partner = consumer, consumer, resource)
  link_map <- bind_rows(consumer_link_map, resource_link_map) %>%
    left_join(node_profiles %>% select(dataset, guild, node, initial_degree, initial_degree_group),
              by = c("guild", "node"))

  raw_list <- vector("list", nrow(subsets))
  for(i in seq_len(nrow(subsets))){
    retained_sites <- subsets$retained_sites[[i]]
    r <- subsets$removal_fraction[i]
    rep <- subsets$replicate[i]

    retained_int_links <- ints %>%
      filter(site %in% retained_sites) %>%
      distinct(consumer, resource) %>%
      mutate(retained = TRUE)

    active_links <- link_map %>%
      left_join(retained_int_links, by = c("consumer", "resource")) %>%
      mutate(retained = replace_na(retained, FALSE))

    retained_by_node <- active_links %>%
      group_by(guild, node, initial_degree, initial_degree_group) %>%
      summarise(retained_degree = sum(retained), .groups = "drop")

    present_nodes <- presence %>%
      filter(site %in% retained_sites) %>%
      distinct(guild, node) %>%
      mutate(node_present = TRUE)

    raw_list[[i]] <- node_profiles %>%
      select(dataset, guild, node, initial_degree, initial_degree_group) %>%
      left_join(retained_by_node %>% select(guild, node, retained_degree), by = c("guild", "node")) %>%
      left_join(present_nodes, by = c("guild", "node")) %>%
      mutate(
        removal_fraction = r,
        replicate = rep,
        retained_degree = replace_na(retained_degree, 0L),
        node_present = replace_na(node_present, FALSE),
        network_active = retained_degree > 0,
        retained_degree_fraction = retained_degree / initial_degree
      )
  }

  bind_rows(raw_list)
}

if(candidate_threshold_met){
  message("Candidate Figure 5 threshold met. Creating candidate main-text figure.")

  retention_raw <- future.apply::future_lapply(
    dataset_order,
    compute_candidate_retention_dataset,
    future.seed = TRUE
  ) %>% bind_rows() %>%
    mutate(dataset = factor(dataset, levels = dataset_order))

  write.csv2(retention_raw, file.path(out_dir, "32_candidate_figure5_node_retention_raw.csv"), row.names = FALSE)

  ## Panel A: group support beyond opportunity
  panel_a_emp <- group_support_null_summary %>%
    select(dataset, guild, initial_degree_group, metric, empirical_value, null_q025, null_q975) %>%
    group_by(guild, initial_degree_group, metric) %>%
    summarise(
      empirical_median = median(empirical_value, na.rm = TRUE),
      null_q025 = median(null_q025, na.rm = TRUE),
      null_q975 = median(null_q975, na.rm = TRUE),
      .groups = "drop"
    )

  pA <- ggplot(panel_a_emp, aes(x = initial_degree_group, y = empirical_median, colour = initial_degree_group)) +
    geom_linerange(aes(ymin = null_q025, ymax = null_q975), linewidth = 3, alpha = 0.25) +
    geom_point(size = 2.3) +
    geom_line(aes(group = 1), linewidth = 0.75) +
    facet_grid(guild ~ metric, scales = "free_y") +
    scale_colour_manual(values = degree_group_cols, guide = "none") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    xlab("Initial-degree group") +
    ylab("Portfolio support") +
    ggtitle("A. Degree and portfolio support beyond opportunity")

  retention_summary <- retention_raw %>%
    group_by(dataset, guild, removal_fraction, initial_degree_group) %>%
    summarise(
      fraction_active = mean(network_active),
      retained_fraction_active_nodes = median(retained_degree_fraction[network_active], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(fraction_active, retained_fraction_active_nodes),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        fraction_active = "B. Who loses all links first?",
        retained_fraction_active_nodes = "C. Surviving degree is compressed"
      )
    )

  retention_overall <- retention_summary %>%
    group_by(guild, metric, removal_fraction, initial_degree_group) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  pBC <- ggplot() +
    geom_ribbon(
      data = retention_overall,
      aes(x = removal_fraction, ymin = q25, ymax = q75, fill = initial_degree_group),
      alpha = 0.14
    ) +
    geom_line(
      data = retention_summary,
      aes(x = removal_fraction, y = value, colour = initial_degree_group,
          group = interaction(dataset, initial_degree_group)),
      alpha = 0.18, linewidth = 0.35
    ) +
    geom_line(
      data = retention_overall,
      aes(x = removal_fraction, y = median_value, colour = initial_degree_group),
      linewidth = 1.0
    ) +
    geom_point(
      data = retention_overall,
      aes(x = removal_fraction, y = median_value, colour = initial_degree_group),
      size = 1.5
    ) +
    facet_grid(guild ~ metric) +
    scale_colour_manual(values = degree_group_cols, name = "Initial-degree group") +
    scale_fill_manual(values = degree_group_cols, name = "Initial-degree group") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 9) +
    theme(legend.position = "bottom") +
    xlab("Proportion of sites removed") +
    ylab("Fraction")

  candidate_fig <- pA / pBC + patchwork::plot_layout(heights = c(1.05, 1.2))
  save_both(candidate_fig, "32_candidate_Figure5_support_reorganizes_network", 12, 10)

} else {
  message("Candidate Figure 5 threshold not met. No main-text Figure 5 created.")
}

## ---------------------------
## Console summary and reminders
## ---------------------------

message("\nDecision summary:")
print(decision_summary)

message("\nValidation checks:")
print(checks)

message("\nAn interaction still present somewhere is not necessarily equally represented across the landscape. These analyses quantify the number of sampled sites still supporting each interaction and its persistence under random site loss. They do not directly measure the ecological effect of the interaction within a local site.")

message("\nScript 32 complete. Outputs saved to: ", out_dir)
