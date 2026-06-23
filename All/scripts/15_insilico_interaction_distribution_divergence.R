## ------------------------------------------------------------
## Script: All/scripts/15_insilico_interaction_distribution_divergence.R
##
## Purpose:
## Small in silico diagnostic: isolate whether redistributing the same
## observed local interaction records (K) across the same observed regional
## links, within exact co-occurrence-frequency n groups, changes divergence
## from the homogeneous M0 site-removal expectation.
##
## This script does NOT add a new ecological model.
## It preserves, for each original dataset:
##   - the observed regional-link list;
##   - each link's n = number of co-occurring sites;
##   - the total number of observed regional links;
##   - total sum(K);
##   - within each exact n, the number of links and total K.
##
## Outputs:
##   15_insilico_divergence_curves.png
##   15_insilico_divergence_summary.png
##   15_insilico_checks.csv
##
## Source:
##   All/scripts/00_dataset_loaders_and_helpers_all.R
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages_extra <- c("dplyr", "tidyr", "tibble", "ggplot2", "purrr")
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
n_site_reps <- 100
n_artificial_sims <- 200
min_sites_retained <- 2

## If TRUE, artificial K allocations are randomly permuted among links within
## each exact n group in each simulation. This avoids making arbitrary row order
## decide which pair gets the high K values in the concentrated scenario.
randomise_K_allocation_within_n_each_sim <- TRUE

dirs <- make_output_dirs("script15_insilico_interaction_distribution_divergence")
sep_out <- dirs$separated
combined_out <- dirs$combined

if(length(all_dataset_names) != 10){
  warning("all_dataset_names has length ", length(all_dataset_names),
          ". This script assumes the original 10 Galiana datasets only.")
}

## ---------------------------
## Helpers
## ---------------------------
make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

plink_M0 <- function(n, p){
  1 - (1 - p)^n
}

expected_links_M0 <- function(n_vec, p){
  sum(plink_M0(n_vec, p), na.rm = TRUE)
}

fit_p_M0 <- function(n_vec, observed_links){
  if(length(n_vec) == 0 || is.na(observed_links)) return(NA_real_)
  if(observed_links <= 0) return(0)
  if(observed_links >= length(n_vec)) return(1)
  f <- function(p){ expected_links_M0(n_vec, p) - observed_links }
  tryCatch(uniroot(f, interval = c(0, 1), tol = 1e-10)$root,
           error = function(e) NA_real_)
}

make_dataset_tables <- function(dataset){
  site_tables <- get_dataset_site_tables(dataset)

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      pair_id = make_pair_id(consumer, resource)
    )

  interactions <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(
      site = as.character(site),
      pair_id = make_pair_id(consumer, resource)
    )

  cooc_counts <- cooc %>%
    group_by(consumer, resource, pair_id) %>%
    summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

  int_counts <- interactions %>%
    group_by(pair_id) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop")

  pair_table_all <- cooc_counts %>%
    left_join(int_counts, by = "pair_id") %>%
    mutate(
      dataset = dataset,
      n_interacting_sites = tidyr::replace_na(n_interacting_sites, 0L),
      is_realised_link = n_interacting_sites > 0,
      empirical_repeatability = n_interacting_sites / n_cooccurring_sites
    ) %>%
    select(dataset, consumer, resource, pair_id,
           n_cooccurring_sites, n_interacting_sites,
           is_realised_link, empirical_repeatability)

  observed_links <- pair_table_all %>%
    filter(is_realised_link) %>%
    arrange(n_cooccurring_sites, pair_id)

  cooc_sites_by_observed_link <- cooc %>%
    semi_join(observed_links %>% select(pair_id), by = "pair_id") %>%
    group_by(pair_id) %>%
    summarise(cooc_sites = list(sort(unique(site))), .groups = "drop")

  observed_links <- observed_links %>%
    left_join(cooc_sites_by_observed_link, by = "pair_id")

  list(
    cooc = cooc,
    interactions = interactions,
    pair_table_all = pair_table_all,
    observed_links = observed_links
  )
}

## Scenario A: as even as possible, all K >= 1 and K <= n.
allocate_K_even <- function(n, m, total_K){
  stopifnot(total_K >= m, total_K <= m * n)
  base <- floor(total_K / m)
  remainder <- total_K %% m
  out <- rep(base, m)
  if(remainder > 0){
    out[seq_len(remainder)] <- out[seq_len(remainder)] + 1L
  }
  stopifnot(all(out >= 1), all(out <= n), sum(out) == total_K)
  out
}

## Scenario B: as concentrated as possible, while all links keep K >= 1.
allocate_K_concentrated <- function(n, m, total_K){
  stopifnot(total_K >= m, total_K <= m * n)
  out <- rep(1L, m)
  remaining <- total_K - m
  i <- 1L
  while(remaining > 0){
    add <- min(n - out[i], remaining)
    out[i] <- out[i] + add
    remaining <- remaining - add
    i <- i + 1L
  }
  stopifnot(all(out >= 1), all(out <= n), sum(out) == total_K)
  out
}

make_artificial_link_table <- function(observed_links, scenario){
  allocator <- switch(
    scenario,
    "Scenario_A_even" = allocate_K_even,
    "Scenario_B_concentrated" = allocate_K_concentrated,
    stop("Unknown artificial scenario: ", scenario)
  )

  observed_links %>%
    group_by(n_cooccurring_sites) %>%
    group_modify(function(.x, .y){
      n <- .y$n_cooccurring_sites[[1]]
      m <- nrow(.x)
      total_K <- sum(.x$n_interacting_sites)
      K_new <- allocator(n = n, m = m, total_K = total_K)
      if(randomise_K_allocation_within_n_each_sim){
        K_new <- sample(K_new, length(K_new), replace = FALSE)
      }
      .x$n_interacting_sites_artificial <- as.integer(K_new)
      .x
    }) %>%
    ungroup()
}

check_K_preservation <- function(empirical_links, artificial_links, scenario_name){
  
  total_ok <- sum(empirical_links$n_interacting_sites) ==
    sum(artificial_links$n_interacting_sites_artificial)
  
  by_n <- empirical_links %>%
    group_by(n_cooccurring_sites) %>%
    summarise(K_empirical = sum(n_interacting_sites), .groups = "drop") %>%
    left_join(
      artificial_links %>%
        group_by(n_cooccurring_sites) %>%
        summarise(K_artificial = sum(n_interacting_sites_artificial), .groups = "drop"),
      by = "n_cooccurring_sites"
    ) %>%
    mutate(ok = K_empirical == K_artificial)
  
  if(!total_ok || any(!by_n$ok)){
    print(by_n)
    stop("K preservation failed for ", scenario_name)
  }
  
  invisible(TRUE)
}

assign_interactions_to_allowed_sites <- function(artificial_links){
  bind_rows(lapply(seq_len(nrow(artificial_links)), function(i){
    row <- artificial_links[i, ]
    allowed_sites <- row$cooc_sites[[1]]
    K <- row$n_interacting_sites_artificial[[1]]

    if(K < 1 || K > length(allowed_sites)){
      stop("Invalid K assignment for pair ", row$pair_id[[1]],
           ": K = ", K, ", available co-occurring sites = ", length(allowed_sites))
    }

    data.frame(
      site = sample(allowed_sites, K, replace = FALSE),
      consumer = row$consumer[[1]],
      resource = row$resource[[1]],
      pair_id = row$pair_id[[1]],
      stringsAsFactors = FALSE
    )
  }))
}

make_site_subsets_local <- function(sites){
  sites <- sort(unique(as.character(sites)))
  n_sites <- length(sites)
  subset_list <- list()
  subset_index <- list()
  counter <- 1L

  for(removal in removal_levels){
    n_keep <- max(min_sites_retained, round(n_sites * (1 - removal)))
    n_keep <- min(n_keep, n_sites)
    reps_here <- ifelse(removal == 0, 1, n_site_reps)

    for(r in seq_len(reps_here)){
      keep <- if(removal == 0) sites else sample(sites, n_keep, replace = FALSE)
      subset_list[[counter]] <- sort(keep)
      subset_index[[counter]] <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        actual_removal_fraction = 1 - length(keep) / n_sites,
        site_rep = r,
        n_sites_kept = length(keep)
      )
      counter <- counter + 1L
    }
  }

  list(subsets = subset_list, index = bind_rows(subset_index))
}

retained_links_given_subsets <- function(interaction_table, subset_object,
                                         expected_M0_by_subset,
                                         observed_full_links,
                                         dataset, scenario, sim_id){
  bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    sites_keep <- subset_object$subsets[[i]]
    idx <- subset_object$index[i, ]

    retained <- interaction_table %>%
      filter(site %in% sites_keep) %>%
      distinct(pair_id) %>%
      nrow()

    expected_M0 <- expected_M0_by_subset$expected_M0_retained_links[i]

    data.frame(
      dataset = dataset,
      scenario = scenario,
      sim_id = sim_id,
      subset_id = idx$subset_id,
      removal_fraction = idx$removal_fraction,
      actual_removal_fraction = idx$actual_removal_fraction,
      site_rep = idx$site_rep,
      n_sites_kept = idx$n_sites_kept,
      retained_regional_links = retained,
      expected_M0_retained_links = expected_M0,
      D_excess_retained_links = (retained - expected_M0) / observed_full_links
    )
  }))
}
check_K_preservation <- function(empirical_links, artificial_links, scenario_name){

  total_ok <- sum(empirical_links$n_interacting_sites) ==
    sum(artificial_links$n_interacting_sites_artificial)

  by_n <- empirical_links %>%
    group_by(n_cooccurring_sites) %>%
    summarise(K_empirical = sum(n_interacting_sites), .groups = "drop") %>%
    left_join(
      artificial_links %>%
        group_by(n_cooccurring_sites) %>%
        summarise(K_artificial = sum(n_interacting_sites_artificial), .groups = "drop"),
      by = "n_cooccurring_sites"
    ) %>%
    mutate(ok = K_empirical == K_artificial)

  if(!total_ok || any(!by_n$ok)){
    print(by_n)
    stop("K preservation failed for ", scenario_name)
  }

  invisible(TRUE)
}
make_check_row <- function(dataset, scenario, link_table_for_check,
                           empirical_group_totals, observed_full_links,
                           total_K_empirical, p_M0){
  group_check <- link_table_for_check %>%
    group_by(n_cooccurring_sites) %>%
    summarise(total_K_scenario = sum(K_check), .groups = "drop") %>%
    left_join(empirical_group_totals, by = "n_cooccurring_sites") %>%
    mutate(group_K_identical = total_K_scenario == total_K_empirical_n)

  data.frame(
    dataset = dataset,
    scenario = scenario,
    n_regional_links = nrow(link_table_for_check),
    total_sum_K = sum(link_table_for_check$K_check),
    same_number_of_links_as_empirical = nrow(link_table_for_check) == observed_full_links,
    same_total_K_as_empirical = sum(link_table_for_check$K_check) == total_K_empirical,
    total_K_identical_within_every_exact_n_group = all(group_check$group_K_identical),
    calibrated_M0_p = p_M0,
    proportion_links_with_K_eq_1 = mean(link_table_for_check$K_check == 1),
    mean_K_among_observed_links = mean(link_table_for_check$K_check),
    variance_K_among_observed_links = ifelse(nrow(link_table_for_check) > 1,
                                             var(link_table_for_check$K_check),
                                             NA_real_),
    n_exact_n_groups = nrow(group_check)
  )
}

manual_inspection <- function(dataset, observed_links, art_A_links, art_B_links){
  empirical_group_totals <- observed_links %>%
    group_by(n_cooccurring_sites) %>%
    summarise(total_K_empirical_n = sum(n_interacting_sites), .groups = "drop")

  check_group <- function(x){
    x %>%
      group_by(n_cooccurring_sites) %>%
      summarise(total_K = sum(n_interacting_sites_artificial), .groups = "drop") %>%
      left_join(empirical_group_totals, by = "n_cooccurring_sites") %>%
      summarise(ok = all(total_K == total_K_empirical_n)) %>%
      pull(ok)
  }

  cat("\nManual inspection dataset:", dataset, "\n")
  cat("  Empirical links:", nrow(observed_links),
      " total K:", sum(observed_links$n_interacting_sites), "\n")
  cat("  Scenario A links:", nrow(art_A_links),
      " total K:", sum(art_A_links$n_interacting_sites_artificial),
      " prop K=1:", round(mean(art_A_links$n_interacting_sites_artificial == 1), 4), "\n")
  cat("  Scenario B links:", nrow(art_B_links),
      " total K:", sum(art_B_links$n_interacting_sites_artificial),
      " prop K=1:", round(mean(art_B_links$n_interacting_sites_artificial == 1), 4), "\n")
  cat("  A has fewer one-site links than B:",
      mean(art_A_links$n_interacting_sites_artificial == 1) <
        mean(art_B_links$n_interacting_sites_artificial == 1), "\n")
  cat("  A preserves total K within exact n groups:", check_group(art_A_links), "\n")
  cat("  B preserves total K within exact n groups:", check_group(art_B_links), "\n")
  cat("  Site reassignment constraint: each artificial interaction is sampled without replacement from that pair's own co-occurrence-site list.\n")
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  message("Running script 15 in silico diagnostic: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tabs <- make_dataset_tables(dataset)
  cooc <- tabs$cooc
  interactions <- tabs$interactions
  pair_table_all <- tabs$pair_table_all
  observed_links <- tabs$observed_links

  observed_full_links <- nrow(observed_links)
  total_K_empirical <- sum(observed_links$n_interacting_sites)

  ## M0 is calibrated from the original full dataset using all co-occurring pairs,
  ## exactly as in previous scripts; it is not refitted after site removal.
  p_M0 <- fit_p_M0(
    n_vec = pair_table_all$n_cooccurring_sites,
    observed_links = observed_full_links
  )

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets_local(all_sites)

  ## Precompute M0 expectation for each site subset from all co-occurring pairs,
  ## including pairs that were not observed as regional links.
  expected_M0_by_subset <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){
    sites_keep <- subset_object$subsets[[i]]
    cooc_retained_counts <- cooc %>%
      filter(site %in% sites_keep) %>%
      group_by(pair_id) %>%
      summarise(n_retained = n_distinct(site), .groups = "drop")

    data.frame(
      subset_id = subset_object$index$subset_id[i],
      expected_M0_retained_links = expected_links_M0(cooc_retained_counts$n_retained, p_M0)
    )
  }))

  ## Empirical curve.
  empirical_rows <- retained_links_given_subsets(
    interaction_table = interactions,
    subset_object = subset_object,
    expected_M0_by_subset = expected_M0_by_subset,
    observed_full_links = observed_full_links,
    dataset = dataset,
    scenario = "Empirical",
    sim_id = 0L
  )

  empirical_group_totals <- observed_links %>%
    group_by(n_cooccurring_sites) %>%
    summarise(total_K_empirical_n = sum(n_interacting_sites), .groups = "drop")

  ## Checks for empirical.
  checks <- list()
  checks[["Empirical"]] <- make_check_row(
    dataset = dataset,
    scenario = "Empirical",
    link_table_for_check = observed_links %>% mutate(K_check = n_interacting_sites),
    empirical_group_totals = empirical_group_totals,
    observed_full_links = observed_full_links,
    total_K_empirical = total_K_empirical,
    p_M0 = p_M0
  )

  artificial_rows <- list()
  artificial_check_collect <- list()

  for(scenario in c("Scenario_A_even", "Scenario_B_concentrated")){
    scenario_rows <- vector("list", n_artificial_sims)
    scenario_check_rows <- vector("list", n_artificial_sims)

    for(sim in seq_len(n_artificial_sims)){
      artificial_links <- make_artificial_link_table(observed_links, scenario = scenario)
      
      check_K_preservation(
        empirical_links = observed_links,
        artificial_links = artificial_links,
        scenario_name = paste(dataset, scenario, "sim", sim)
      )
      
      artificial_interactions <- assign_interactions_to_allowed_sites(artificial_links)

      ## Hard safety check: artificial interaction sites must be allowed co-occurrence sites.
      assignment_ok <- artificial_interactions %>%
        left_join(observed_links %>% select(pair_id, cooc_sites), by = "pair_id") %>%
        rowwise() %>%
        mutate(site_allowed = site %in% cooc_sites) %>%
        ungroup() %>%
        summarise(ok = all(site_allowed)) %>%
        pull(ok)

      if(!isTRUE(assignment_ok)){
        stop("Artificial interaction assigned outside pair-specific co-occurrence sites in ",
             dataset, " / ", scenario, " / sim ", sim)
      }

      scenario_rows[[sim]] <- retained_links_given_subsets(
        interaction_table = artificial_interactions,
        subset_object = subset_object,
        expected_M0_by_subset = expected_M0_by_subset,
        observed_full_links = observed_full_links,
        dataset = dataset,
        scenario = scenario,
        sim_id = sim
      )

      scenario_check_rows[[sim]] <- artificial_links %>%
        transmute(dataset = dataset,
                  scenario = scenario,
                  sim_id = sim,
                  pair_id = pair_id,
                  n_cooccurring_sites = n_cooccurring_sites,
                  K_check = n_interacting_sites_artificial)
    }

    artificial_rows[[scenario]] <- bind_rows(scenario_rows)
    artificial_check_collect[[scenario]] <- bind_rows(scenario_check_rows)

    ## Check rows are summarized over all simulations; because K allocation can be
    ## permuted among links within n, variance may vary slightly across simulations.
    scenario_check_summary <- artificial_check_collect[[scenario]] %>%
      group_by(sim_id) %>%
      group_modify(function(.x, .y){
        make_check_row(
          dataset = dataset,
          scenario = scenario,
          link_table_for_check = .x,
          empirical_group_totals = empirical_group_totals,
          observed_full_links = observed_full_links,
          total_K_empirical = total_K_empirical,
          p_M0 = p_M0
        )
      }) %>%
      ungroup() %>%
      summarise(
        dataset = dataset,
        scenario = scenario,
        n_regional_links = unique(n_regional_links)[1],
        total_sum_K = unique(total_sum_K)[1],
        same_number_of_links_as_empirical = all(same_number_of_links_as_empirical),
        same_total_K_as_empirical = all(same_total_K_as_empirical),
        total_K_identical_within_every_exact_n_group = all(total_K_identical_within_every_exact_n_group),
        calibrated_M0_p = unique(calibrated_M0_p)[1],
        proportion_links_with_K_eq_1 = median(proportion_links_with_K_eq_1),
        mean_K_among_observed_links = median(mean_K_among_observed_links),
        variance_K_among_observed_links = median(variance_K_among_observed_links),
        n_exact_n_groups = unique(n_exact_n_groups)[1],
        .groups = "drop"
      )

    checks[[scenario]] <- scenario_check_summary
  }

  ## Manual inspection once, using a fresh representative allocation for A and B.
  if(dataset == all_dataset_names[[1]]){
    art_A_links <- make_artificial_link_table(observed_links, "Scenario_A_even")
    art_B_links <- make_artificial_link_table(observed_links, "Scenario_B_concentrated")
    manual_inspection(dataset, observed_links, art_A_links, art_B_links)
  }

  all_rows <- bind_rows(empirical_rows, bind_rows(artificial_rows))
  check_rows <- bind_rows(checks)

  write.csv2(all_rows,
             file.path(out_dir, paste0(dataset, "_15_insilico_divergence_by_subset.csv")),
             row.names = FALSE)
  write.csv2(check_rows,
             file.path(out_dir, paste0(dataset, "_15_insilico_checks.csv")),
             row.names = FALSE)

  list(rows = all_rows, checks = check_rows)
}

## ---------------------------
## Run
## ---------------------------

all_outputs <- lapply(all_dataset_names, run_one_dataset)
names(all_outputs) <- all_dataset_names

rows_all <- bind_rows(lapply(all_outputs, `[[`, "rows"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

write.csv2(
  rows_all,
  file.path(combined_out, "15_insilico_divergence_by_subset.csv"),
  row.names = FALSE
)

write.csv2(
  checks_all,
  file.path(combined_out, "15_insilico_checks.csv"),
  row.names = FALSE
)

cat("\nCompact check table:\n")
print(
  checks_all %>%
    select(dataset, scenario, n_regional_links, total_sum_K,
           total_K_identical_within_every_exact_n_group,
           calibrated_M0_p, proportion_links_with_K_eq_1,
           mean_K_among_observed_links, variance_K_among_observed_links),
  n = Inf
)

## ---------------------------
## Plot summaries
## ---------------------------

scenario_labels <- c(
  "Empirical" = "empirical data",
  "Scenario_A_even" = "Scenario A, interactions spread evenly",
  "Scenario_B_concentrated" = "Scenario B, interactions concentrated"
)

scenario_colours <- c(
  "empirical data" = "black",
  "Scenario A, interactions spread evenly" = "#1f78b4",
  "Scenario B, interactions concentrated" = "#e31a1c"
)

## Empirical: summarize across site-removal reps.
empirical_curve <- rows_all %>%
  filter(scenario == "Empirical") %>%
  group_by(dataset, scenario, actual_removal_fraction) %>%
  summarise(
    D_median = median(D_excess_retained_links, na.rm = TRUE),
    D_low = NA_real_,
    D_high = NA_real_,
    .groups = "drop"
  )

## Artificial: first average across site-removal reps within each simulation and
## removal level, then take median and central 95% interval across simulations.
artificial_curve <- rows_all %>%
  filter(scenario != "Empirical") %>%
  group_by(dataset, scenario, sim_id, actual_removal_fraction) %>%
  summarise(D_sim = median(D_excess_retained_links, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, scenario, actual_removal_fraction) %>%
  summarise(
    D_median = median(D_sim, na.rm = TRUE),
    D_low = quantile(D_sim, 0.025, na.rm = TRUE),
    D_high = quantile(D_sim, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

curve_plot_data <- bind_rows(empirical_curve, artificial_curve) %>%
  mutate(
    scenario_label = scenario_labels[scenario],
    scenario_label = factor(scenario_label, levels = scenario_labels),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p1 <- ggplot(curve_plot_data,
             aes(x = actual_removal_fraction,
                 y = D_median,
                 colour = scenario_label,
                 fill = scenario_label)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_ribbon(
    data = curve_plot_data %>% filter(scenario != "Empirical"),
    aes(ymin = D_low, ymax = D_high),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = scenario_colours, name = NULL) +
  scale_fill_manual(values = scenario_colours, name = NULL) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Excess retained links relative to homogeneous model, D(r)") +
  ggtitle("In silico redistribution of observed interaction records",
          subtitle = "Artificial scenarios preserve link list, n values, total K, and total K within exact n groups")

ggsave(
  file.path(combined_out, "15_insilico_divergence_curves.png"),
  p1,
  width = 14,
  height = 7,
  dpi = 300
)

## Summary statistic: average D(r) over all nonzero removal levels.
empirical_summary <- rows_all %>%
  filter(scenario == "Empirical", actual_removal_fraction > 0) %>%
  group_by(dataset, scenario, actual_removal_fraction) %>%
  summarise(D_removal = median(D_excess_retained_links, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, scenario) %>%
  summarise(mean_D_nonzero_removal = mean(D_removal, na.rm = TRUE), .groups = "drop")

artificial_summary <- rows_all %>%
  filter(scenario != "Empirical", actual_removal_fraction > 0) %>%
  group_by(dataset, scenario, sim_id, actual_removal_fraction) %>%
  summarise(D_removal = median(D_excess_retained_links, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, scenario, sim_id) %>%
  summarise(mean_D_sim = mean(D_removal, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, scenario) %>%
  summarise(mean_D_nonzero_removal = median(mean_D_sim, na.rm = TRUE), .groups = "drop")

summary_plot_data <- bind_rows(empirical_summary, artificial_summary) %>%
  mutate(
    scenario_label = scenario_labels[scenario],
    scenario_label = factor(scenario_label, levels = scenario_labels),
    dataset = factor(dataset, levels = all_dataset_names)
  )

p2 <- ggplot(summary_plot_data,
             aes(x = dataset,
                 y = mean_D_nonzero_removal,
                 colour = scenario_label)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_point(size = 2.5, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = scenario_colours, name = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Average D(r) across nonzero removal levels") +
  ggtitle("Average excess retained links under interaction-record redistribution")

ggsave(
  file.path(combined_out, "15_insilico_divergence_summary.png"),
  p2,
  width = 12,
  height = 5.5,
  dpi = 300
)

message("Finished script 15 in silico interaction-distribution diagnostic.")
message("Outputs written to: ", combined_out)
