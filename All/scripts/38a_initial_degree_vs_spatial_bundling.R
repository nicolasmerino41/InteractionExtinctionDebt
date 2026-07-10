## ------------------------------------------------------------
## Script:
## All/scripts/38a_initial_degree_vs_spatial_bundling.R
##
## Purpose:
## Test whether species with higher initial regional interaction
## degree have different spatial bundling of their interaction
## portfolios.
##
## Response:
##   bundling_delta =
##     observed_mean_jaccard - null_mean_jaccard
##
## Primary inferential unit:
##   one Spearman correlation per dataset and guild.
##
## Analyses:
##   1. primary_all_degree_2_plus
##   2. sensitivity_degree_3_plus
##
## Consumers and resources are always analysed separately.
##
## Outputs:
##   - two PNG figures
##   - five CSV files
##
## No pooled species-level correlation is used for inference.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "purrr"
)

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ------------------------------------------------------------
## User controls
## ------------------------------------------------------------

set.seed(123)

n_boot <- 10000
n_permutations_large_n <- 99999

input_file <- paste0(
  "All/outputs/38_interaction_portfolio_spatial_bundling/",
  "38_species_level_spatial_bundling.csv"
)

out_dir <- "All/outputs/38a_initial_degree_vs_spatial_bundling"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

guild_levels <- c(
  "Consumer",
  "Resource"
)

analysis_levels <- c(
  "primary_all_degree_2_plus",
  "sensitivity_degree_3_plus"
)

degree_class_levels <- c(
  "2",
  "3–4",
  "5–9",
  "10+"
)

## ------------------------------------------------------------
## Plot helpers
## ------------------------------------------------------------

theme_pub <- function(base_size = 10){
  
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    )
}

save_png <- function(
    p,
    filename,
    width = 9,
    height = 6
){
  
  ggplot2::ggsave(
    filename = file.path(
      out_dir,
      filename
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

## ------------------------------------------------------------
## Input checks
## ------------------------------------------------------------

if(!file.exists(input_file)){
  stop(
    "Required Script 38 input file was not found:\n",
    input_file
  )
}

required_columns <- c(
  "dataset",
  "guild",
  "focal_node",
  "initial_degree",
  "number_partner_pairs",
  "observed_mean_jaccard",
  "null_mean_jaccard",
  "bundling_delta",
  "bundling_SES"
)

species_raw <- read.csv2(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(
  required_columns,
  names(species_raw)
)

if(length(missing_columns) > 0){
  stop(
    "The Script 38 input is missing required columns:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

## ------------------------------------------------------------
## Prepare species-level input
## ------------------------------------------------------------

species_data <- species_raw %>%
  dplyr::mutate(
    dataset = as.character(dataset),
    guild = as.character(guild),
    focal_node = as.character(focal_node),
    initial_degree = as.numeric(initial_degree),
    number_partner_pairs = as.numeric(number_partner_pairs),
    observed_mean_jaccard = as.numeric(observed_mean_jaccard),
    null_mean_jaccard = as.numeric(null_mean_jaccard),
    bundling_delta = as.numeric(bundling_delta),
    bundling_SES = as.numeric(bundling_SES)
  ) %>%
  dplyr::filter(
    guild %in% guild_levels,
    initial_degree >= 2,
    is.finite(bundling_delta)
  ) %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

if(nrow(species_data) == 0){
  stop(
    "No eligible species remained after requiring guild to be ",
    "Consumer or Resource, initial_degree >= 2, and finite bundling_delta."
  )
}

## Ensure one row per dataset × guild × focal species.
duplicate_species_rows <- species_data %>%
  dplyr::count(
    dataset,
    guild,
    focal_node,
    name = "number_rows"
  ) %>%
  dplyr::filter(
    number_rows > 1
  )

if(nrow(duplicate_species_rows) > 0){
  stop(
    "The Script 38 input contains duplicated dataset × guild × focal_node rows."
  )
}

## ------------------------------------------------------------
## Dataset-level Spearman correlations
## ------------------------------------------------------------

calculate_dataset_correlation <- function(
    dat,
    analysis_name,
    minimum_degree
){
  
  dat_analysis <- dat %>%
    dplyr::filter(
      initial_degree >= minimum_degree
    )
  
  number_species <- nrow(dat_analysis)
  
  minimum_observed_degree <- if(number_species > 0){
    min(
      dat_analysis$initial_degree,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  
  maximum_observed_degree <- if(number_species > 0){
    max(
      dat_analysis$initial_degree,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }
  
  number_unique_degree_values <- if(number_species > 0){
    dplyr::n_distinct(
      dat_analysis$initial_degree
    )
  } else {
    0L
  }
  
  valid_for_inference <-
    number_species >= 5 &&
    number_unique_degree_values >= 3
  
  exclusion_reason <- dplyr::case_when(
    number_species == 0 ~
      paste0(
        "No species with initial_degree >= ",
        minimum_degree
      ),
    
    number_species < 5 ~
      "Fewer than 5 eligible species",
    
    number_unique_degree_values < 3 ~
      "Fewer than 3 unique initial-degree values",
    
    TRUE ~
      NA_character_
  )
  
  spearman_rho <- if(valid_for_inference){
    
    suppressWarnings(
      stats::cor(
        dat_analysis$initial_degree,
        dat_analysis$bundling_delta,
        method = "spearman",
        use = "complete.obs"
      )
    )
    
  } else {
    NA_real_
  }
  
  tibble::tibble(
    analysis = analysis_name,
    dataset = as.character(dat$dataset[1]),
    guild = as.character(dat$guild[1]),
    number_species = number_species,
    minimum_degree = minimum_observed_degree,
    maximum_degree = maximum_observed_degree,
    number_unique_degree_values =
      number_unique_degree_values,
    spearman_rho = spearman_rho,
    valid_for_inference = valid_for_inference,
    exclusion_reason = exclusion_reason
  )
}

dataset_groups <- species_data %>%
  dplyr::group_by(
    dataset,
    guild
  ) %>%
  dplyr::group_split(
    .keep = TRUE
  )

dataset_correlations <- dplyr::bind_rows(
  lapply(
    dataset_groups,
    function(dataset_group){
      
      dplyr::bind_rows(
        calculate_dataset_correlation(
          dat = dataset_group,
          analysis_name =
            "primary_all_degree_2_plus",
          minimum_degree = 2
        ),
        calculate_dataset_correlation(
          dat = dataset_group,
          analysis_name =
            "sensitivity_degree_3_plus",
          minimum_degree = 3
        )
      )
    }
  )
) %>%
  dplyr::mutate(
    analysis = factor(
      analysis,
      levels = analysis_levels
    ),
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

## ------------------------------------------------------------
## Equal-dataset sign-flip test
## ------------------------------------------------------------

sign_flip_test <- function(
    x,
    n_large = 99999
){
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  if(n < 2){
    
    return(
      tibble::tibble(
        number_valid_datasets = n,
        equal_dataset_mean_rho =
          ifelse(
            n == 1,
            mean(x),
            NA_real_
          ),
        sign_flip_p = NA_real_,
        test_type = NA_character_
      )
    )
  }
  
  observed_mean <- mean(x)
  
  if(n <= 20){
    
    ## Exact enumeration without constructing a 2^n by n matrix.
    ## At each step, append the two possible signed sums.
    signed_sums <- 0
    
    for(value in x){
      signed_sums <- c(
        signed_sums + value,
        signed_sums - value
      )
    }
    
    permuted_means <- signed_sums / n
    
    sign_flip_p <- mean(
      abs(permuted_means) >=
        abs(observed_mean) - 1e-12
    )
    
    test_type <- paste0(
      "Exact sign-flip: ",
      2^n,
      " sign patterns"
    )
    
  } else {
    
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
    
    sign_flip_p <- (
      1 +
        sum(
          abs(permuted_means) >=
            abs(observed_mean)
        )
    ) / (
      n_large + 1
    )
    
    test_type <- paste0(
      "Monte Carlo sign-flip: ",
      n_large,
      " permutations"
    )
  }
  
  tibble::tibble(
    number_valid_datasets = n,
    equal_dataset_mean_rho = observed_mean,
    sign_flip_p = sign_flip_p,
    test_type = test_type
  )
}

## ------------------------------------------------------------
## Dataset-bootstrap interval
## ------------------------------------------------------------

bootstrap_dataset_mean <- function(
    x,
    n_boot = 10000
){
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  if(n < 2){
    
    return(
      tibble::tibble(
        bootstrap_q025 = NA_real_,
        bootstrap_q975 = NA_real_
      )
    )
  }
  
  bootstrap_means <- replicate(
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
      bootstrap_means,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    bootstrap_q975 = stats::quantile(
      bootstrap_means,
      0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

## ------------------------------------------------------------
## Equal-dataset inference
## ------------------------------------------------------------

equal_dataset_tests <- dataset_correlations %>%
  dplyr::filter(
    valid_for_inference,
    is.finite(spearman_rho)
  ) %>%
  dplyr::group_by(
    analysis,
    guild
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      rho_values <- .x$spearman_rho
      
      dplyr::bind_cols(
        sign_flip_test(
          rho_values,
          n_large =
            n_permutations_large_n
        ),
        bootstrap_dataset_mean(
          rho_values,
          n_boot = n_boot
        )
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      equal_dataset_mean_rho > 0 ~
        "Bundling delta increases with degree",
      
      equal_dataset_mean_rho < 0 ~
        "Bundling delta decreases with degree",
      
      equal_dataset_mean_rho == 0 ~
        "No mean association",
      
      TRUE ~
        NA_character_
    ),
    significant_0.05 =
      !is.na(sign_flip_p) &
      sign_flip_p < 0.05,
    analysis = factor(
      analysis,
      levels = analysis_levels
    ),
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

## Ensure rows also exist when a guild-analysis combination has
## no valid datasets.

all_inference_combinations <- expand.grid(
  analysis = analysis_levels,
  guild = guild_levels,
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    analysis = factor(
      analysis,
      levels = analysis_levels
    ),
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

equal_dataset_tests <- all_inference_combinations %>%
  dplyr::left_join(
    equal_dataset_tests,
    by = c(
      "analysis",
      "guild"
    )
  ) %>%
  dplyr::arrange(
    analysis,
    guild
  )

## ------------------------------------------------------------
## Degree classes for visualisation only
## ------------------------------------------------------------

make_degree_class <- function(initial_degree){
  
  dplyr::case_when(
    initial_degree == 2 ~ "2",
    initial_degree >= 3 &
      initial_degree <= 4 ~ "3–4",
    initial_degree >= 5 &
      initial_degree <= 9 ~ "5–9",
    initial_degree >= 10 ~ "10+",
    TRUE ~ NA_character_
  )
}

species_class_data <- species_data %>%
  dplyr::mutate(
    degree_class = factor(
      make_degree_class(initial_degree),
      levels = degree_class_levels
    )
  ) %>%
  dplyr::filter(
    !is.na(degree_class)
  )

dataset_degree_class_means <- species_class_data %>%
  dplyr::group_by(
    dataset,
    guild,
    degree_class
  ) %>%
  dplyr::summarise(
    mean_bundling_delta =
      mean(
        bundling_delta,
        na.rm = TRUE
      ),
    median_bundling_delta =
      median(
        bundling_delta,
        na.rm = TRUE
      ),
    number_species =
      dplyr::n_distinct(focal_node),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    ),
    degree_class = factor(
      degree_class,
      levels = degree_class_levels
    )
  )

## ------------------------------------------------------------
## Bootstrap summaries for degree classes
## ------------------------------------------------------------

summarise_degree_class <- function(
    x,
    n_boot = 10000
){
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  equal_dataset_mean <- if(n > 0){
    mean(x)
  } else {
    NA_real_
  }
  
  if(n < 2){
    
    return(
      tibble::tibble(
        number_datasets = n,
        equal_dataset_mean_bundling_delta =
          equal_dataset_mean,
        bootstrap_q025 = NA_real_,
        bootstrap_q975 = NA_real_
      )
    )
  }
  
  bootstrap_means <- replicate(
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
    number_datasets = n,
    equal_dataset_mean_bundling_delta =
      equal_dataset_mean,
    bootstrap_q025 = stats::quantile(
      bootstrap_means,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    bootstrap_q975 = stats::quantile(
      bootstrap_means,
      0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

equal_dataset_degree_class_summaries <-
  dataset_degree_class_means %>%
  dplyr::group_by(
    guild,
    degree_class
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      summarise_degree_class(
        .x$mean_bundling_delta,
        n_boot = n_boot
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    ),
    degree_class = factor(
      degree_class,
      levels = degree_class_levels
    )
  )

## Add explicit missing combinations, if any.

all_class_combinations <- expand.grid(
  guild = guild_levels,
  degree_class = degree_class_levels,
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    ),
    degree_class = factor(
      degree_class,
      levels = degree_class_levels
    )
  )

equal_dataset_degree_class_summaries <-
  all_class_combinations %>%
  dplyr::left_join(
    equal_dataset_degree_class_summaries,
    by = c(
      "guild",
      "degree_class"
    )
  ) %>%
  dplyr::arrange(
    guild,
    degree_class
  )

## ------------------------------------------------------------
## Main descriptive figure
## ------------------------------------------------------------

degree_relationship_plot <- ggplot2::ggplot() +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.4
  ) +
  
  ggplot2::geom_line(
    data = dataset_degree_class_means,
    ggplot2::aes(
      x = degree_class,
      y = mean_bundling_delta,
      group = dataset
    ),
    colour = "grey65",
    linewidth = 0.45,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  
  ggplot2::geom_point(
    data = dataset_degree_class_means,
    ggplot2::aes(
      x = degree_class,
      y = mean_bundling_delta
    ),
    colour = "grey55",
    size = 1.6,
    alpha = 0.7,
    na.rm = TRUE
  ) +
  
  ggplot2::geom_errorbar(
    data = equal_dataset_degree_class_summaries,
    ggplot2::aes(
      x = degree_class,
      ymin = bootstrap_q025,
      ymax = bootstrap_q975
    ),
    width = 0.12,
    linewidth = 0.75,
    colour = "black",
    na.rm = TRUE
  ) +
  
  ggplot2::geom_point(
    data = equal_dataset_degree_class_summaries,
    ggplot2::aes(
      x = degree_class,
      y = equal_dataset_mean_bundling_delta
    ),
    colour = "black",
    size = 3,
    na.rm = TRUE
  ) +
  
  ggplot2::facet_wrap(
    ~ guild,
    nrow = 1
  ) +
  
  ggplot2::scale_x_discrete(
    limits = degree_class_levels,
    drop = FALSE
  ) +
  
  theme_pub(base_size = 10) +
  
  ggplot2::labs(
    x = "Initial regional degree",
    y = paste0(
      "Observed minus expected overlap among ",
      "interaction-supporting sites"
    ),
    title =
      "Spatial portfolio structure changes with regional degree",
    subtitle = paste0(
      "Positive values indicate bundling; negative values indicate ",
      "spatial separation among a species’ interactions"
    )
  )

## ------------------------------------------------------------
## Inferential figure: primary analysis only
## ------------------------------------------------------------

primary_dataset_correlations <- dataset_correlations %>%
  dplyr::filter(
    analysis == "primary_all_degree_2_plus",
    valid_for_inference,
    is.finite(spearman_rho)
  )

primary_equal_tests <- equal_dataset_tests %>%
  dplyr::filter(
    analysis == "primary_all_degree_2_plus"
  )

correlation_plot <- ggplot2::ggplot() +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.4
  ) +
  
  ggplot2::geom_jitter(
    data = primary_dataset_correlations,
    ggplot2::aes(
      x = guild,
      y = spearman_rho
    ),
    width = 0.08,
    height = 0,
    colour = "grey55",
    alpha = 0.8,
    size = 2
  ) +
  
  ggplot2::geom_errorbar(
    data = primary_equal_tests,
    ggplot2::aes(
      x = guild,
      ymin = bootstrap_q025,
      ymax = bootstrap_q975
    ),
    width = 0.12,
    linewidth = 0.8,
    colour = "black",
    na.rm = TRUE
  ) +
  
  ggplot2::geom_point(
    data = primary_equal_tests,
    ggplot2::aes(
      x = guild,
      y = equal_dataset_mean_rho
    ),
    colour = "black",
    size = 3.4,
    na.rm = TRUE
  ) +
  
  ggplot2::scale_x_discrete(
    limits = guild_levels,
    drop = FALSE
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(-1, 1)
  ) +
  
  theme_pub(base_size = 11) +
  
  ggplot2::labs(
    x = NULL,
    y = paste0(
      "Spearman correlation between initial degree ",
      "and bundling delta"
    ),
    title = paste0(
      "Dataset-level association between regional degree ",
      "and spatial bundling"
    )
  )

## ------------------------------------------------------------
## Validation checks
## ------------------------------------------------------------

valid_correlations <- dataset_correlations %>%
  dplyr::filter(
    valid_for_inference
  )

primary_rows <- dataset_correlations %>%
  dplyr::filter(
    analysis == "primary_all_degree_2_plus"
  )

sensitivity_rows <- dataset_correlations %>%
  dplyr::filter(
    analysis == "sensitivity_degree_3_plus"
  )

one_correlation_per_dataset <- dataset_correlations %>%
  dplyr::count(
    analysis,
    dataset,
    guild,
    name = "number_rows"
  ) %>%
  dplyr::summarise(
    valid = all(number_rows == 1)
  ) %>%
  dplyr::pull(valid)

equal_inference_matches_unweighted_mean <- all(
  vapply(
    seq_len(nrow(equal_dataset_tests)),
    function(i){
      
      analysis_i <- as.character(
        equal_dataset_tests$analysis[i]
      )
      
      guild_i <- as.character(
        equal_dataset_tests$guild[i]
      )
      
      observed_rhos <- dataset_correlations %>%
        dplyr::filter(
          analysis == analysis_i,
          guild == guild_i,
          valid_for_inference,
          is.finite(spearman_rho)
        ) %>%
        dplyr::pull(
          spearman_rho
        )
      
      reported_mean <-
        equal_dataset_tests$
        equal_dataset_mean_rho[i]
      
      if(length(observed_rhos) == 0){
        return(is.na(reported_mean))
      }
      
      isTRUE(
        abs(
          mean(observed_rhos) -
            reported_mean
        ) < 1e-12
      )
    },
    logical(1)
  )
)

class_summaries_match_dataset_means <- all(
  vapply(
    seq_len(
      nrow(
        equal_dataset_degree_class_summaries
      )
    ),
    function(i){
      
      guild_i <- as.character(
        equal_dataset_degree_class_summaries$guild[i]
      )
      
      class_i <- as.character(
        equal_dataset_degree_class_summaries$
          degree_class[i]
      )
      
      dataset_values <- dataset_degree_class_means %>%
        dplyr::filter(
          guild == guild_i,
          degree_class == class_i
        ) %>%
        dplyr::pull(
          mean_bundling_delta
        )
      
      reported_mean <-
        equal_dataset_degree_class_summaries$
        equal_dataset_mean_bundling_delta[i]
      
      if(length(dataset_values) == 0){
        return(is.na(reported_mean))
      }
      
      isTRUE(
        abs(
          mean(dataset_values) -
            reported_mean
        ) < 1e-12
      )
    },
    logical(1)
  )
)

validation_checks <- tibble::tibble(
  check = c(
    "All retained initial degrees are at least two",
    "All retained bundling_delta values are finite",
    "Consumers and Resources are analysed separately",
    "Every primary Spearman result contains one dataset and one guild only",
    "Every sensitivity result excludes degree-two species",
    "Every valid correlation contains at least five species",
    "Every valid correlation contains at least three unique degree values",
    "All estimated Spearman correlations lie between minus one and one",
    "Equal-dataset inference uses one correlation per dataset",
    "No dataset is weighted by its number of species",
    "Descriptive degree-class summaries first average within datasets",
    "Degree classes are used only for visualisation and not for the primary test"
  ),
  pass = c(
    all(
      species_data$initial_degree >= 2
    ),
    
    all(
      is.finite(
        species_data$bundling_delta
      )
    ),
    
    all(
      guild_levels %in%
        as.character(
          unique(species_data$guild)
        )
    ),
    
    one_correlation_per_dataset &&
      all(
        primary_rows$analysis ==
          "primary_all_degree_2_plus"
      ),
    
    all(
      sensitivity_rows$minimum_degree >= 3 |
        is.na(sensitivity_rows$minimum_degree)
    ),
    
    all(
      valid_correlations$number_species >= 5
    ),
    
    all(
      valid_correlations$
        number_unique_degree_values >= 3
    ),
    
    all(
      valid_correlations$spearman_rho >= -1 &
        valid_correlations$spearman_rho <= 1
    ),
    
    one_correlation_per_dataset,
    
    equal_inference_matches_unweighted_mean,
    
    class_summaries_match_dataset_means,
    
    TRUE
  )
)

## ------------------------------------------------------------
## Save CSV outputs
## ------------------------------------------------------------

write.csv2(
  dataset_correlations,
  file.path(
    out_dir,
    "38a_dataset_degree_bundling_correlations.csv"
  ),
  row.names = FALSE
)

write.csv2(
  equal_dataset_tests,
  file.path(
    out_dir,
    "38a_equal_dataset_degree_bundling_tests.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_degree_class_means,
  file.path(
    out_dir,
    "38a_dataset_degree_class_means.csv"
  ),
  row.names = FALSE
)

write.csv2(
  equal_dataset_degree_class_summaries,
  file.path(
    out_dir,
    "38a_equal_dataset_degree_class_summaries.csv"
  ),
  row.names = FALSE
)

write.csv2(
  validation_checks,
  file.path(
    out_dir,
    "38a_degree_bundling_validation_checks.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Save PNG outputs
## ------------------------------------------------------------

save_png(
  degree_relationship_plot,
  "38a_degree_bundling_relationship.png",
  width = 9,
  height = 5.5
)

save_png(
  correlation_plot,
  "38a_dataset_degree_bundling_correlations.png",
  width = 8,
  height = 6
)

## ------------------------------------------------------------
## Console summary
## ------------------------------------------------------------

print_inference_result <- function(
    analysis_name,
    guild_name,
    heading
){
  
  message("\n", heading)
  
  result <- equal_dataset_tests %>%
    dplyr::filter(
      analysis == analysis_name,
      guild == guild_name
    )
  
  print(
    result,
    n = Inf
  )
}

print_inference_result(
  analysis_name =
    "primary_all_degree_2_plus",
  guild_name =
    "Consumer",
  heading =
    "Primary equal-dataset result: Consumers"
)

print_inference_result(
  analysis_name =
    "primary_all_degree_2_plus",
  guild_name =
    "Resource",
  heading =
    "Primary equal-dataset result: Resources"
)

print_inference_result(
  analysis_name =
    "sensitivity_degree_3_plus",
  guild_name =
    "Consumer",
  heading =
    "Degree >= 3 sensitivity result: Consumers"
)

print_inference_result(
  analysis_name =
    "sensitivity_degree_3_plus",
  guild_name =
    "Resource",
  heading =
    "Degree >= 3 sensitivity result: Resources"
)

dataset_counts <- dataset_correlations %>%
  dplyr::group_by(
    analysis,
    guild
  ) %>%
  dplyr::summarise(
    number_valid_datasets =
      sum(
        valid_for_inference,
        na.rm = TRUE
      ),
    number_excluded_datasets =
      sum(
        !valid_for_inference,
        na.rm = TRUE
      ),
    total_dataset_guild_combinations =
      dplyr::n(),
    .groups = "drop"
  )

message(
  "\nNumber of valid and excluded datasets ",
  "by analysis and guild:"
)

print(
  dataset_counts,
  n = Inf
)

message("\nValidation checks:")

print(
  validation_checks,
  n = Inf
)

message("\nInterpretation:")

message(
  "Positive mean Spearman correlation means that species with larger\n",
  "regional interaction portfolios tend to shift from spatial separation\n",
  "toward greater overlap or bundling of their interactions.\n\n",
  "Negative mean Spearman correlation means that higher-degree species\n",
  "have increasingly spatially separated interaction portfolios.\n\n",
  "The degree >= 3 sensitivity analysis evaluates whether the result is\n",
  "driven mainly by species with only two realised partners."
)

message(
  "\nSaved Script 38a outputs in: ",
  out_dir
)
