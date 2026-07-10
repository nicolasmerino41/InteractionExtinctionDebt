## ------------------------------------------------------------
## Script:
## All/scripts/36b_Figure4_generalist_specialist_support_tests.R
##
## Purpose:
## Test whether Specialists and Generalists differ in:
##
##   K  = mean number of sites supporting each realised interaction
##   q  = mean K_ij / n_ij among realised interaction partners
##   CO = mean co-occurrence support n_ij / total dataset sites
##        among realised interaction partners
##
## Primary inferential unit:
##   Dataset.
##
## For every dataset and guild, calculate:
##   contrast = Generalists - Specialists
##
## Then test whether the equal-dataset mean contrast differs from zero
## using a two-sided sign-flip permutation test.
##
## Consumers and resources are tested separately.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "purrr",
  "scales"
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

n_boot <- 10000
n_permutations_large_n <- 99999

out_dir <- "All/outputs/36b_generalist_specialist_support_tests"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

degree_levels <- c(
  "Specialists",
  "Generalists"
)

degree_cols <- c(
  "Specialists" = "#4E79A7",
  "Generalists" = "#59A14F"
)

outcome_labels <- c(
  K = "K: local interaction support",
  q = "q: relative local interaction realisation",
  CO = "CO: co-occurrence support"
)

theme_pub <- function(base_size = 10){
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

save_png <- function(p, filename, width = 11, height = 8){
  ggplot2::ggsave(
    filename = file.path(out_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

replace_na_int <- function(x, value = 0L){
  x[is.na(x)] <- value
  x
}

## ---------------------------
## Degree grouping
## ---------------------------

## This matches the two-group logic used in Script 36:
## unique degree values are ranked within each dataset and guild.
## The lower half of unique degree values is called Specialists,
## and the upper half is called Generalists.
##
## Identical degree values always remain in the same group.

make_degree_groups <- function(df, value_col = "initial_degree"){
  
  vals <- sort(unique(df[[value_col]]))
  
  ## A contrast cannot be estimated if all nodes have the same degree.
  if(length(vals) < 2){
    return(
      df %>%
        dplyr::mutate(
          initial_degree_group = factor(
            NA_character_,
            levels = degree_levels
          )
        )
    )
  }
  
  degree_key <- tibble::tibble(
    !!value_col := vals,
    degree_rank_fraction =
      rank(vals, ties.method = "average") / length(vals)
  ) %>%
    dplyr::mutate(
      initial_degree_group = dplyr::case_when(
        degree_rank_fraction <= 0.5 ~ "Specialists",
        TRUE ~ "Generalists"
      )
    ) %>%
    dplyr::select(
      -degree_rank_fraction
    )
  
  df %>%
    dplyr::left_join(
      degree_key,
      by = value_col
    ) %>%
    dplyr::mutate(
      initial_degree_group = factor(
        initial_degree_group,
        levels = degree_levels
      )
    )
}

## ---------------------------
## Build species-level outcomes
## ---------------------------

make_guild_pair_table <- function(full_links, guild_name){
  
  if(guild_name == "Consumer"){
    
    full_links %>%
      dplyr::transmute(
        guild = "Consumer",
        focal_node = consumer,
        partner_node = resource,
        full_n,
        full_K
      )
    
  } else {
    
    full_links %>%
      dplyr::transmute(
        guild = "Resource",
        focal_node = resource,
        partner_node = consumer,
        full_n,
        full_K
      )
  }
}

run_one_dataset <- function(dataset){
  
  message("Generalist-specialist tests: ", dataset)
  
  st <- get_dataset_site_tables(dataset)
  
  cooc <- st$cooc_triples %>%
    dplyr::distinct(
      site,
      consumer,
      resource
    ) %>%
    dplyr::mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )
  
  ints <- st$empirical_site_interactions %>%
    dplyr::distinct(
      site,
      consumer,
      resource
    ) %>%
    dplyr::mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )
  
  ## An observed interaction implies that both species were together.
  cooc <- dplyr::bind_rows(
    cooc,
    ints
  ) %>%
    dplyr::distinct(
      site,
      consumer,
      resource
    )
  
  number_sites <- dplyr::n_distinct(cooc$site)
  
  pair_table <- cooc %>%
    dplyr::count(
      consumer,
      resource,
      name = "full_n"
    ) %>%
    dplyr::left_join(
      ints %>%
        dplyr::count(
          consumer,
          resource,
          name = "full_K"
        ),
      by = c("consumer", "resource")
    ) %>%
    dplyr::mutate(
      full_K = replace_na_int(full_K, 0L)
    )
  
  if(any(pair_table$full_n < 1)){
    stop("Invalid full_n < 1 in ", dataset)
  }
  
  if(any(pair_table$full_K > pair_table$full_n)){
    stop("Invalid full_K > full_n in ", dataset)
  }
  
  ## All three requested outcomes are conditional on the regional
  ## interaction existing at least once.
  full_links <- pair_table %>%
    dplyr::filter(full_K > 0)
  
  if(nrow(full_links) == 0){
    warning("No realised regional links in ", dataset)
    return(tibble::tibble())
  }
  
  guild_pairs <- dplyr::bind_rows(
    make_guild_pair_table(full_links, "Consumer"),
    make_guild_pair_table(full_links, "Resource")
  )
  
  species_values <- guild_pairs %>%
    dplyr::group_by(
      guild,
      focal_node
    ) %>%
    dplyr::summarise(
      initial_degree = dplyr::n_distinct(partner_node),
      
      ## K_i = arithmetic mean of K_ij over realised partners.
      K = mean(
        full_K,
        na.rm = TRUE
      ),
      
      ## q_i = arithmetic mean of K_ij / n_ij over realised partners.
      ## This is conditional local realisation, not an all-opportunity
      ## conversion rate.
      q = mean(
        full_K / full_n,
        na.rm = TRUE
      ),
      
      ## CO_i = arithmetic mean co-occurrence support per realised
      ## partner, expressed as a fraction of all sites in the dataset.
      CO = mean(
        full_n / number_sites,
        na.rm = TRUE
      ),
      
      number_realised_partners = dplyr::n_distinct(partner_node),
      
      .groups = "drop"
    ) %>%
    dplyr::group_by(guild) %>%
    dplyr::group_modify(
      ~ make_degree_groups(.x, "initial_degree")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dataset = dataset,
      number_dataset_sites = number_sites
    ) %>%
    dplyr::select(
      dataset,
      guild,
      focal_node,
      initial_degree,
      initial_degree_group,
      number_dataset_sites,
      number_realised_partners,
      K,
      q,
      CO
    )
  
  species_values
}

## ---------------------------
## Run datasets
## ---------------------------

species_values <- dplyr::bind_rows(
  lapply(
    all_dataset_names,
    run_one_dataset
  )
)

species_values <- species_values %>%
  dplyr::mutate(
    initial_degree_group = factor(
      initial_degree_group,
      levels = degree_levels
    )
  )

## ---------------------------
## Dataset-level group means
## ---------------------------

make_long_outcomes <- function(dat){
  
  dplyr::bind_rows(
    dat %>%
      dplyr::transmute(
        dataset,
        guild,
        focal_node,
        initial_degree,
        initial_degree_group,
        outcome = "K",
        value = K
      ),
    
    dat %>%
      dplyr::transmute(
        dataset,
        guild,
        focal_node,
        initial_degree,
        initial_degree_group,
        outcome = "q",
        value = q
      ),
    
    dat %>%
      dplyr::transmute(
        dataset,
        guild,
        focal_node,
        initial_degree,
        initial_degree_group,
        outcome = "CO",
        value = CO
      )
  ) %>%
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = c("K", "q", "CO")
      )
    )
}

outcome_long <- make_long_outcomes(species_values)

dataset_group_means <- outcome_long %>%
  dplyr::filter(
    !is.na(initial_degree_group),
    is.finite(value)
  ) %>%
  dplyr::group_by(
    dataset,
    guild,
    outcome,
    initial_degree_group
  ) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    median_value = median(value, na.rm = TRUE),
    number_species = dplyr::n_distinct(focal_node),
    .groups = "drop"
  )

## Build one row per dataset × guild × outcome.
specialist_means <- dataset_group_means %>%
  dplyr::filter(initial_degree_group == "Specialists") %>%
  dplyr::transmute(
    dataset,
    guild,
    outcome,
    specialist_mean = mean_value,
    n_specialists = number_species
  )

generalist_means <- dataset_group_means %>%
  dplyr::filter(initial_degree_group == "Generalists") %>%
  dplyr::transmute(
    dataset,
    guild,
    outcome,
    generalist_mean = mean_value,
    n_generalists = number_species
  )

dataset_contrasts <- specialist_means %>%
  dplyr::inner_join(
    generalist_means,
    by = c(
      "dataset",
      "guild",
      "outcome"
    )
  ) %>%
  dplyr::mutate(
    contrast = generalist_mean - specialist_mean
  )

## ---------------------------
## Sign-flip test
## ---------------------------

sign_flip_test <- function(x, n_large = 99999){
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  if(n < 2){
    return(tibble::tibble(
      n_datasets = n,
      mean_contrast = ifelse(n == 1, mean(x), NA_real_),
      permutation_p = NA_real_,
      permutation_type = NA_character_
    ))
  }
  
  observed <- mean(x)
  
  if(n <= 20){
    
    ## Exact test over all 2^n possible sign assignments.
    sign_grid <- expand.grid(
      rep(
        list(c(-1, 1)),
        n
      )
    )
    
    sign_matrix <- as.matrix(sign_grid)
    
    permuted_means <- rowMeans(
      sweep(
        sign_matrix,
        MARGIN = 2,
        STATS = x,
        FUN = "*"
      )
    )
    
    p_value <- mean(
      abs(permuted_means) >= abs(observed) - 1e-12
    )
    
    test_type <- paste0(
      "Exact sign-flip, ",
      2^n,
      " permutations"
    )
    
  } else {
    
    ## Monte Carlo approximation if there are too many datasets
    ## for exact enumeration.
    permuted_means <- replicate(
      n_large,
      mean(
        x * sample(
          c(-1, 1),
          size = n,
          replace = TRUE
        )
      )
    )
    
    p_value <- (
      1 +
        sum(
          abs(permuted_means) >= abs(observed)
        )
    ) / (
      n_large + 1
    )
    
    test_type <- paste0(
      "Monte Carlo sign-flip, ",
      n_large,
      " permutations"
    )
  }
  
  tibble::tibble(
    n_datasets = n,
    mean_contrast = observed,
    permutation_p = p_value,
    permutation_type = test_type
  )
}

## ---------------------------
## Dataset-bootstrap interval
## ---------------------------

bootstrap_dataset_contrast <- function(x, n_boot = 10000){
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  if(n < 2){
    return(tibble::tibble(
      bootstrap_q025 = NA_real_,
      bootstrap_q975 = NA_real_
    ))
  }
  
  boot_means <- replicate(
    n_boot,
    mean(
      sample(
        x,
        size = n,
        replace = TRUE
      )
    )
  )
  
  tibble::tibble(
    bootstrap_q025 = stats::quantile(
      boot_means,
      0.025,
      na.rm = TRUE
    ),
    bootstrap_q975 = stats::quantile(
      boot_means,
      0.975,
      na.rm = TRUE
    )
  )
}

## ---------------------------
## Sensitivity tests
## ---------------------------

sensitivity_tests <- function(x){
  
  x <- x[is.finite(x)]
  
  if(length(x) < 2){
    return(tibble::tibble(
      paired_dataset_t_p = NA_real_,
      wilcoxon_signed_rank_p = NA_real_
    ))
  }
  
  t_result <- tryCatch(
    stats::t.test(
      x,
      mu = 0,
      alternative = "two.sided"
    ),
    error = function(e) NULL
  )
  
  wilcox_result <- tryCatch(
    stats::wilcox.test(
      x,
      mu = 0,
      alternative = "two.sided",
      exact = FALSE,
      conf.int = FALSE
    ),
    error = function(e) NULL
  )
  
  tibble::tibble(
    paired_dataset_t_p = if(is.null(t_result)){
      NA_real_
    } else {
      t_result$p.value
    },
    
    wilcoxon_signed_rank_p = if(is.null(wilcox_result)){
      NA_real_
    } else {
      wilcox_result$p.value
    }
  )
}

## ---------------------------
## Run inference
## ---------------------------

test_results <- dataset_contrasts %>%
  dplyr::group_by(
    guild,
    outcome
  ) %>%
  dplyr::group_modify(function(.x, .y){
    
    x <- .x$contrast
    
    dplyr::bind_cols(
      sign_flip_test(
        x,
        n_large = n_permutations_large_n
      ),
      bootstrap_dataset_contrast(
        x,
        n_boot = n_boot
      ),
      sensitivity_tests(x)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    outcome_label = unname(
      outcome_labels[as.character(outcome)]
    ),
    
    direction = dplyr::case_when(
      mean_contrast > 0 ~ "Generalists higher",
      mean_contrast < 0 ~ "Specialists higher",
      TRUE ~ "No mean difference"
    ),
    
    significant_0.05 = !is.na(permutation_p) &
      permutation_p < 0.05,
    
    permutation_p_BH = stats::p.adjust(
      permutation_p,
      method = "BH"
    ),
    
    significant_BH_0.05 = !is.na(permutation_p_BH) &
      permutation_p_BH < 0.05
  )

## ---------------------------
## Plot data
## ---------------------------

dataset_contrasts_plot <- dataset_contrasts %>%
  dplyr::mutate(
    outcome_label = factor(
      unname(outcome_labels[as.character(outcome)]),
      levels = unname(outcome_labels)
    )
  )

test_results_plot <- test_results %>%
  dplyr::mutate(
    outcome_label = factor(
      outcome_label,
      levels = unname(outcome_labels)
    )
  )

contrast_plot <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = dataset_contrasts_plot,
    ggplot2::aes(
      x = dataset,
      y = contrast,
      group = dataset
    ),
    colour = "grey45",
    alpha = 0.75,
    size = 1.7
  ) +
  ggplot2::geom_errorbar(
    data = test_results_plot,
    ggplot2::aes(
      x = "Equal-dataset mean",
      ymin = bootstrap_q025,
      ymax = bootstrap_q975
    ),
    width = 0.16,
    linewidth = 0.8,
    colour = "#222222"
  ) +
  ggplot2::geom_point(
    data = test_results_plot,
    ggplot2::aes(
      x = "Equal-dataset mean",
      y = mean_contrast
    ),
    size = 3,
    colour = "#222222"
  ) +
  ggplot2::facet_grid(
    outcome_label ~ guild,
    scales = "free_y"
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Generalists minus Specialists",
    title = "Dataset-level contrasts between Generalists and Specialists",
    subtitle = paste0(
      "Large point and interval: equal-dataset mean contrast and ",
      "dataset-bootstrap 95% interval"
    )
  )

## ---------------------------
## Group-value plot
## ---------------------------

group_overall <- dataset_group_means %>%
  dplyr::group_by(
    guild,
    outcome,
    initial_degree_group
  ) %>%
  dplyr::summarise(
    mean_value = mean(mean_value, na.rm = TRUE),
    q25 = stats::quantile(mean_value, 0.25, na.rm = TRUE),
    q75 = stats::quantile(mean_value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    outcome_label = factor(
      unname(outcome_labels[as.character(outcome)]),
      levels = unname(outcome_labels)
    )
  )

group_plot_data <- dataset_group_means %>%
  dplyr::mutate(
    outcome_label = factor(
      unname(outcome_labels[as.character(outcome)]),
      levels = unname(outcome_labels)
    )
  )

group_plot <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = group_plot_data,
    ggplot2::aes(
      x = initial_degree_group,
      y = mean_value,
      group = dataset,
      colour = initial_degree_group
    ),
    alpha = 0.25,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = group_plot_data,
    ggplot2::aes(
      x = initial_degree_group,
      y = mean_value,
      colour = initial_degree_group
    ),
    alpha = 0.25,
    size = 1.3
  ) +
  ggplot2::geom_errorbar(
    data = group_overall,
    ggplot2::aes(
      x = initial_degree_group,
      ymin = q25,
      ymax = q75,
      colour = initial_degree_group
    ),
    width = 0.14,
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = group_overall,
    ggplot2::aes(
      x = initial_degree_group,
      y = mean_value,
      colour = initial_degree_group
    ),
    size = 3
  ) +
  ggplot2::facet_grid(
    outcome_label ~ guild,
    scales = "free_y"
  ) +
  ggplot2::scale_colour_manual(
    values = degree_cols,
    drop = FALSE
  ) +
  theme_pub(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1
    )
  ) +
  ggplot2::labs(
    x = "Initial-degree group",
    y = "Dataset-level group mean",
    colour = NULL,
    title = "Generalist and Specialist values across datasets"
  )

## ---------------------------
## Validation checks
## ---------------------------

checks <- tibble::tibble(
  check = c(
    "All analysed species have initial degree greater than zero",
    "All K values are positive",
    "All q values are between zero and one",
    "All CO values are between zero and one",
    "Every tested dataset has both degree groups",
    "Primary test uses one contrast per dataset",
    "Consumers and resources are tested separately"
  ),
  pass = c(
    all(species_values$initial_degree > 0, na.rm = TRUE),
    all(species_values$K > 0, na.rm = TRUE),
    all(
      species_values$q >= 0 &
        species_values$q <= 1,
      na.rm = TRUE
    ),
    all(
      species_values$CO >= 0 &
        species_values$CO <= 1,
      na.rm = TRUE
    ),
    all(
      dataset_contrasts$n_specialists > 0 &
        dataset_contrasts$n_generalists > 0
    ),
    TRUE,
    all(
      c("Consumer", "Resource") %in%
        unique(dataset_contrasts$guild)
    )
  )
)

## ---------------------------
## Save outputs
## ---------------------------

write.csv2(
  species_values,
  file.path(
    out_dir,
    "36b_species_level_K_q_CO_by_degree_group.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_group_means,
  file.path(
    out_dir,
    "36b_dataset_degree_group_means_K_q_CO.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_contrasts,
  file.path(
    out_dir,
    "36b_dataset_generalist_minus_specialist_contrasts.csv"
  ),
  row.names = FALSE
)

write.csv2(
  test_results,
  file.path(
    out_dir,
    "36b_generalist_specialist_significance_tests.csv"
  ),
  row.names = FALSE
)

write.csv2(
  checks,
  file.path(
    out_dir,
    "36b_generalist_specialist_test_checks.csv"
  ),
  row.names = FALSE
)

save_png(
  contrast_plot,
  "36b_dataset_contrasts_K_q_CO.png",
  width = 12,
  height = 9
)

save_png(
  group_plot,
  "36b_group_values_K_q_CO.png",
  width = 10,
  height = 9
)

## ---------------------------
## Console summary
## ---------------------------

message("\nPrimary equal-dataset significance tests:")
print(
  test_results %>%
    dplyr::select(
      guild,
      outcome,
      n_datasets,
      mean_contrast,
      bootstrap_q025,
      bootstrap_q975,
      permutation_p,
      permutation_p_BH,
      direction,
      significant_0.05,
      significant_BH_0.05
    ),
  n = Inf
)

message("\nValidation checks:")
print(checks, n = Inf)

message("\nInterpretation:")
message(
  "Positive contrasts mean that Generalists have higher values than Specialists. ",
  "The sign-flip permutation test evaluates whether dataset-specific contrasts ",
  "are systematically displaced from zero while giving every dataset equal weight."
)

message("\nSaved outputs in: ", out_dir)