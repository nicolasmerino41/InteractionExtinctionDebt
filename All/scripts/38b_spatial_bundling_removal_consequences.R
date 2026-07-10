## ------------------------------------------------------------
## Script:
## All/scripts/38b_spatial_bundling_removal_consequences_parallel.R
##
## Purpose:
## Parallel dataset-level implementation.
##
## Test whether the observed spatial arrangement of a species'
## realised interactions changes the variability and extremity
## of portfolio loss under random site removal.
##
## The observed arrangement is compared with null arrangements
## preserving, for every realised interaction:
##
##   - consumer-resource identity;
##   - absolute interaction support, full_K;
##   - co-occurrence support, full_n;
##   - exact pair-specific co-occurrence sites;
##   - regional network membership and degree.
##
## Primary responses:
##
##   variance_delta
##   complete_loss_delta
##   half_loss_delta
##
## Negative-control response:
##
##   mean_retention_delta
##
## Primary inference:
##
##   Within each dataset and guild, calculate the Spearman
##   correlation between bundling_delta and each consequence AUC.
##   Then combine one correlation per dataset using equal dataset
##   weights.
##
## Outputs:
##
##   38b_species_level_removal_consequences.csv
##   38b_species_level_consequence_AUC.csv
##   38b_dataset_bundling_consequence_correlations.csv
##   38b_equal_dataset_bundling_consequence_tests.csv
##   38b_equal_dataset_degree_class_removal_summaries.csv
##   38b_spatial_bundling_removal_validation_checks.csv
##
##   38b_bundling_vulnerability_correlations.png
##   38b_variance_consequence_by_degree.png
##
## No replicate-level or null-arrangement-level outputs are saved.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "purrr",
  "scales",
  "future",
  "future.apply"
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

removal_levels <- c(
  0,
  0.1,
  0.2,
  0.4,
  0.6,
  0.8
)

n_removal_reps <- 1000
n_null <- 250
n_boot <- 10000
n_permutations_large_n <- 99999

## Number of parallel R sessions used for independent datasets.
## Set to 1L to reproduce the serial execution path.
n_workers <- max(
  1L,
  parallel::detectCores(logical = TRUE) - 1L
)

## Avoid opening more workers than there are datasets.
n_workers <- min(
  n_workers,
  length(all_dataset_names)
)

## Monte Carlo tolerance for comparing simulated mean retained
## fractions with their analytical hypergeometric expectations.
mean_retention_tolerance <- 0.035

## Numerical tolerance for quantities that should equal zero.
zero_tolerance <- 1e-10

bundling_input_file <- file.path(
  "All",
  "outputs",
  "38_interaction_portfolio_spatial_bundling",
  "38_species_level_spatial_bundling.csv"
)

out_dir <- file.path(
  "All",
  "outputs",
  "38b_spatial_bundling_removal_consequences"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

guild_levels <- c(
  "Consumer",
  "Resource"
)

correlation_analysis_levels <- c(
  "ordinary_spearman",
  "partial_spearman_controlling_degree"
)

outcome_levels <- c(
  "variance",
  "complete_loss",
  "half_loss"
)

outcome_labels <- c(
  "variance" = "Variance in retained portfolio",
  "complete_loss" = "Complete portfolio loss",
  "half_loss" = "Loss of at least half the portfolio"
)

contrast_levels <- c(
  "mean_retention",
  "variance",
  "complete_loss",
  "half_loss"
)

degree_class_levels <- c(
  "2",
  "3–4",
  "5–9",
  "10+"
)

degree_class_colours <- c(
  "2" = "#4E79A7",
  "3–4" = "#F28E2B",
  "5–9" = "#59A14F",
  "10+" = "#B07AA1"
)

## ------------------------------------------------------------
## General helpers
## ------------------------------------------------------------

theme_pub <- function(base_size = 10){
  
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

save_png <- function(
    plot,
    filename,
    width = 9,
    height = 6
){
  
  ggplot2::ggsave(
    filename = file.path(
      out_dir,
      filename
    ),
    plot = plot,
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

safe_mean <- function(x){
  
  x <- x[is.finite(x)]
  
  if(length(x) == 0){
    return(NA_real_)
  }
  
  mean(x)
}

safe_sd <- function(x){
  
  x <- x[is.finite(x)]
  
  if(length(x) < 2){
    return(NA_real_)
  }
  
  stats::sd(x)
}

safe_var <- function(x){
  
  x <- x[is.finite(x)]
  
  if(length(x) < 2){
    return(NA_real_)
  }
  
  stats::var(x)
}

safe_quantile <- function(x, probability){
  
  x <- x[is.finite(x)]
  
  if(length(x) == 0){
    return(NA_real_)
  }
  
  stats::quantile(
    x,
    probability,
    names = FALSE,
    na.rm = TRUE
  )
}

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

## ------------------------------------------------------------
## Read and validate Script 38 input
## ------------------------------------------------------------

if(!file.exists(bundling_input_file)){
  stop(
    "Required Script 38 species-level file was not found:\n",
    bundling_input_file
  )
}

required_bundling_columns <- c(
  "dataset",
  "guild",
  "focal_node",
  "initial_degree",
  "observed_mean_jaccard",
  "null_mean_jaccard",
  "bundling_delta",
  "bundling_SES"
)

bundling_raw <- read.csv2(
  bundling_input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_bundling_columns <- setdiff(
  required_bundling_columns,
  names(bundling_raw)
)

if(length(missing_bundling_columns) > 0){
  stop(
    "The Script 38 input is missing required columns:\n",
    paste(
      missing_bundling_columns,
      collapse = ", "
    )
  )
}

bundling_species <- bundling_raw %>%
  dplyr::transmute(
    dataset = as.character(dataset),
    guild = as.character(guild),
    focal_node = as.character(focal_node),
    input_initial_degree = as.numeric(initial_degree),
    observed_mean_jaccard =
      as.numeric(observed_mean_jaccard),
    null_mean_jaccard =
      as.numeric(null_mean_jaccard),
    bundling_delta =
      as.numeric(bundling_delta),
    bundling_SES =
      as.numeric(bundling_SES)
  ) %>%
  dplyr::filter(
    guild %in% guild_levels,
    input_initial_degree >= 2,
    is.finite(bundling_delta)
  )

if(nrow(bundling_species) == 0){
  stop(
    "No eligible species remained after filtering Script 38 input."
  )
}

duplicated_bundling_species <- bundling_species %>%
  dplyr::count(
    dataset,
    guild,
    focal_node,
    name = "number_rows"
  ) %>%
  dplyr::filter(
    number_rows > 1
  )

if(nrow(duplicated_bundling_species) > 0){
  stop(
    "Duplicated dataset × guild × focal_node rows occur in ",
    "the Script 38 input."
  )
}

## ------------------------------------------------------------
## Site-removal masks
## ------------------------------------------------------------

make_removal_masks <- function(
    sites,
    removal_levels,
    n_removal_reps
){
  
  sites <- sort(
    unique(
      as.character(sites)
    )
  )
  
  number_dataset_sites <- length(sites)
  
  if(number_dataset_sites < 1){
    stop("Cannot generate removal masks for zero sites.")
  }
  
  mask_list <- vector(
    "list",
    length(removal_levels)
  )
  
  metadata_list <- vector(
    "list",
    length(removal_levels)
  )
  
  for(level_index in seq_along(removal_levels)){
    
    removal_fraction <- removal_levels[level_index]
    
    n_keep <- max(
      1L,
      round(
        (1 - removal_fraction) *
          number_dataset_sites
      )
    )
    
    n_keep <- min(
      n_keep,
      number_dataset_sites
    )
    
    mask_matrix <- matrix(
      FALSE,
      nrow = n_removal_reps,
      ncol = number_dataset_sites,
      dimnames = list(
        NULL,
        sites
      )
    )
    
    if(n_keep == number_dataset_sites){
      
      mask_matrix[,] <- TRUE
      
    } else {
      
      for(removal_rep in seq_len(n_removal_reps)){
        
        retained_site_indices <- sample.int(
          number_dataset_sites,
          size = n_keep,
          replace = FALSE
        )
        
        mask_matrix[
          removal_rep,
          retained_site_indices
        ] <- TRUE
      }
    }
    
    mask_list[[level_index]] <- mask_matrix
    
    metadata_list[[level_index]] <- tibble::tibble(
      removal_level_index = level_index,
      removal_fraction = removal_fraction,
      number_dataset_sites =
        number_dataset_sites,
      actual_sites_retained = n_keep,
      actual_retained_fraction =
        n_keep / number_dataset_sites,
      actual_removal_fraction =
        1 - n_keep / number_dataset_sites,
      number_removal_replicates =
        n_removal_reps
    )
  }
  
  names(mask_list) <- as.character(
    removal_levels
  )
  
  list(
    sites = sites,
    masks = mask_list,
    metadata = dplyr::bind_rows(
      metadata_list
    )
  )
}

## ------------------------------------------------------------
## Species-occurrence / representation tables
## ------------------------------------------------------------

find_first_column <- function(
    data,
    candidates
){
  
  data_names_lower <- tolower(
    names(data)
  )
  
  candidate_positions <- match(
    tolower(candidates),
    data_names_lower
  )
  
  candidate_positions <- candidate_positions[
    !is.na(candidate_positions)
  ]
  
  if(length(candidate_positions) == 0){
    return(NA_character_)
  }
  
  names(data)[candidate_positions[1]]
}

extract_independent_occurrence_table <- function(st){
  
  possible_object_names <- c(
    "species_occurrences",
    "species_occurrence",
    "occurrences",
    "occurrence",
    "occurrence_table",
    "site_species",
    "species_by_site"
  )
  
  for(object_name in possible_object_names){
    
    if(
      object_name %in% names(st) &&
      is.data.frame(st[[object_name]])
    ){
      
      candidate <- st[[object_name]]
      
      site_column <- find_first_column(
        candidate,
        c(
          "site",
          "site_id",
          "plot",
          "plot_id",
          "location",
          "locality"
        )
      )
      
      species_column <- find_first_column(
        candidate,
        c(
          "species",
          "species_name",
          "taxon",
          "taxon_name",
          "sp"
        )
      )
      
      if(
        !is.na(site_column) &&
        !is.na(species_column)
      ){
        
        return(
          candidate %>%
            dplyr::transmute(
              site = as.character(
                .data[[site_column]]
              ),
              focal_node = as.character(
                .data[[species_column]]
              )
            ) %>%
            dplyr::filter(
              !is.na(site),
              !is.na(focal_node)
            ) %>%
            dplyr::distinct(
              site,
              focal_node
            )
        )
      }
    }
  }
  
  NULL
}

make_focal_representation_tables <- function(
    st,
    cooc,
    dataset_sites
){
  
  independent_occurrence <-
    extract_independent_occurrence_table(st)
  
  if(!is.null(independent_occurrence)){
    
    independent_occurrence <-
      independent_occurrence %>%
      dplyr::filter(
        site %in% dataset_sites
      )
    
    return(
      list(
        Consumer = independent_occurrence,
        Resource = independent_occurrence,
        Consumer_definition =
          "Independent species-occurrence table",
        Resource_definition =
          "Independent species-occurrence table",
        Consumer_variable_name =
          "focal_species_present",
        Resource_variable_name =
          "focal_species_present",
        used_cooccurrence_fallback = FALSE
      )
    )
  }
  
  consumer_representation <- cooc %>%
    dplyr::transmute(
      site,
      focal_node = consumer
    ) %>%
    dplyr::distinct(
      site,
      focal_node
    )
  
  resource_representation <- cooc %>%
    dplyr::transmute(
      site,
      focal_node = resource
    ) %>%
    dplyr::distinct(
      site,
      focal_node
    )
  
  list(
    Consumer = consumer_representation,
    Resource = resource_representation,
    Consumer_definition = paste0(
      "Focal co-occurrence representation; ",
      "independent occurrence table unavailable"
    ),
    Resource_definition = paste0(
      "Focal co-occurrence representation; ",
      "independent occurrence table unavailable"
    ),
    Consumer_variable_name =
      "focal_cooccurrence_representation",
    Resource_variable_name =
      "focal_cooccurrence_representation",
    used_cooccurrence_fallback = TRUE
  )
}

## ------------------------------------------------------------
## Prepare one dataset
## ------------------------------------------------------------

prepare_dataset <- function(dataset){
  
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
  
  ## Observed interaction implies co-occurrence.
  cooc <- dplyr::bind_rows(
    cooc,
    ints
  ) %>%
    dplyr::distinct(
      site,
      consumer,
      resource
    )
  
  dataset_sites <- sort(
    unique(cooc$site)
  )
  
  cooccurrence_table <- cooc %>%
    dplyr::group_by(
      consumer,
      resource
    ) %>%
    dplyr::summarise(
      full_n = dplyr::n_distinct(site),
      cooccurrence_sites = list(
        sort(unique(site))
      ),
      .groups = "drop"
    )
  
  interaction_table <- ints %>%
    dplyr::group_by(
      consumer,
      resource
    ) %>%
    dplyr::summarise(
      full_K = dplyr::n_distinct(site),
      interaction_sites = list(
        sort(unique(site))
      ),
      .groups = "drop"
    )
  
  pair_table <- cooccurrence_table %>%
    dplyr::left_join(
      interaction_table,
      by = c(
        "consumer",
        "resource"
      )
    ) %>%
    dplyr::mutate(
      full_K = replace_na_int(
        full_K,
        0L
      )
    )
  
  if(any(pair_table$full_n < 1)){
    stop("Invalid full_n < 1 in ", dataset)
  }
  
  if(any(pair_table$full_K > pair_table$full_n)){
    stop("Invalid full_K > full_n in ", dataset)
  }
  
  realised_links <- pair_table %>%
    dplyr::filter(
      full_K > 0
    ) %>%
    dplyr::mutate(
      link_id = dplyr::row_number()
    )
  
  if(nrow(realised_links) == 0){
    stop("No realised links in ", dataset)
  }
  
  observed_sites_valid <- all(
    unlist(
      Map(
        function(interaction_sites, possible_sites){
          all(interaction_sites %in% possible_sites)
        },
        realised_links$interaction_sites,
        realised_links$cooccurrence_sites
      ),
      use.names = FALSE
    )
  )
  
  if(!observed_sites_valid){
    stop(
      "An interaction-supporting site falls outside its ",
      "pair-specific co-occurrence set in ",
      dataset
    )
  }
  
  focal_representation <- make_focal_representation_tables(
    st = st,
    cooc = cooc,
    dataset_sites = dataset_sites
  )
  
  list(
    dataset = dataset,
    cooc = cooc,
    ints = ints,
    dataset_sites = dataset_sites,
    realised_links = realised_links,
    focal_representation =
      focal_representation
  )
}

## ------------------------------------------------------------
## Site-by-link and site-by-species matrices
## ------------------------------------------------------------

make_link_matrix <- function(
    supporting_site_sets,
    dataset_sites
){
  
  site_lookup <- stats::setNames(
    seq_along(dataset_sites),
    dataset_sites
  )
  
  number_links <- length(
    supporting_site_sets
  )
  
  link_matrix <- matrix(
    0L,
    nrow = length(dataset_sites),
    ncol = number_links
  )
  
  for(link_index in seq_len(number_links)){
    
    supporting_sites <-
      supporting_site_sets[[link_index]]
    
    link_matrix[
      site_lookup[supporting_sites],
      link_index
    ] <- 1L
  }
  
  link_matrix
}

make_focal_representation_matrix <- function(
    focal_nodes,
    representation_table,
    dataset_sites
){
  
  site_lookup <- stats::setNames(
    seq_along(dataset_sites),
    dataset_sites
  )
  
  focal_lookup <- stats::setNames(
    seq_along(focal_nodes),
    focal_nodes
  )
  
  representation_matrix <- matrix(
    0L,
    nrow = length(dataset_sites),
    ncol = length(focal_nodes),
    dimnames = list(
      dataset_sites,
      focal_nodes
    )
  )
  
  valid_rows <- representation_table %>%
    dplyr::filter(
      site %in% dataset_sites,
      focal_node %in% focal_nodes
    )
  
  if(nrow(valid_rows) > 0){
    
    representation_matrix[
      cbind(
        site_lookup[valid_rows$site],
        focal_lookup[valid_rows$focal_node]
      )
    ] <- 1L
  }
  
  representation_matrix
}

## ------------------------------------------------------------
## Focal metadata
## ------------------------------------------------------------

make_focal_metadata <- function(
    dataset,
    guild,
    realised_links
){
  
  if(guild == "Consumer"){
    
    focal_links <- realised_links %>%
      dplyr::transmute(
        link_id,
        focal_node = consumer,
        partner_node = resource,
        full_K,
        full_n
      )
    
  } else {
    
    focal_links <- realised_links %>%
      dplyr::transmute(
        link_id,
        focal_node = resource,
        partner_node = consumer,
        full_K,
        full_n
      )
  }
  
  focal_metadata <- focal_links %>%
    dplyr::group_by(focal_node) %>%
    dplyr::summarise(
      initial_degree =
        dplyr::n_distinct(partner_node),
      link_indices = list(
        sort(unique(link_id))
      ),
      .groups = "drop"
    ) %>%
    dplyr::inner_join(
      bundling_species %>%
        dplyr::filter(
          dataset == !!dataset,
          guild == !!guild
        ),
      by = "focal_node"
    ) %>%
    dplyr::mutate(
      dataset = dataset,
      guild = guild,
      degree_class = factor(
        make_degree_class(initial_degree),
        levels = degree_class_levels
      )
    )
  
  degree_mismatch <- focal_metadata %>%
    dplyr::filter(
      initial_degree !=
        input_initial_degree
    )
  
  if(nrow(degree_mismatch) > 0){
    stop(
      "Initial-degree mismatch between reconstructed network ",
      "and Script 38 input for ",
      dataset,
      " / ",
      guild
    )
  }
  
  focal_metadata
}

## ------------------------------------------------------------
## Evaluate one spatial arrangement
## ------------------------------------------------------------

evaluate_arrangement <- function(
    link_matrix,
    removal_masks,
    focal_metadata,
    focal_representation_matrix
){
  
  number_focal_species <- nrow(
    focal_metadata
  )
  
  number_levels <- length(
    removal_masks
  )
  
  output_list <- vector(
    "list",
    number_levels
  )
  
  retained_fraction_bounds_valid <- TRUE
  probability_bounds_valid <- TRUE
  variance_nonnegative <- TRUE
  complete_not_above_half <- TRUE
  
  for(level_index in seq_len(number_levels)){
    
    retained_site_mask <-
      removal_masks[[level_index]]
    
    retained_support_count <-
      retained_site_mask %*%
      link_matrix
    
    link_survives <-
      retained_support_count > 0
    
    focal_remains <-
      retained_site_mask %*%
      focal_representation_matrix > 0
    
    focal_rows <- vector(
      "list",
      number_focal_species
    )
    
    for(focal_index in seq_len(number_focal_species)){
      
      link_indices <-
        focal_metadata$link_indices[[focal_index]]
      
      initial_degree <-
        focal_metadata$initial_degree[focal_index]
      
      retained_degree <- rowSums(
        link_survives[
          ,
          link_indices,
          drop = FALSE
        ]
      )
      
      retained_degree_fraction <-
        retained_degree /
        initial_degree
      
      focal_present_vector <-
        focal_remains[, focal_index]
      
      number_focal_present_replicates <- sum(
        focal_present_vector
      )
      
      if(number_focal_present_replicates > 0){
        
        complete_loss_probability <- mean(
          retained_degree[
            focal_present_vector
          ] == 0
        )
        
        half_loss_probability <- mean(
          retained_degree[
            focal_present_vector
          ] <= floor(initial_degree / 2)
        )
        
        conditional_reason <- NA_character_
        
      } else {
        
        complete_loss_probability <- NA_real_
        half_loss_probability <- NA_real_
        
        conditional_reason <- paste0(
          "No removal replicate retained focal representation"
        )
      }
      
      variance_retained_fraction <- safe_var(
        retained_degree_fraction
      )
      
      retained_fraction_bounds_valid <-
        retained_fraction_bounds_valid &&
        all(
          retained_degree_fraction >= 0 &
            retained_degree_fraction <= 1
        )
      
      finite_probabilities <- c(
        complete_loss_probability,
        half_loss_probability
      )
      
      finite_probabilities <- finite_probabilities[
        is.finite(finite_probabilities)
      ]
      
      if(length(finite_probabilities) > 0){
        
        probability_bounds_valid <-
          probability_bounds_valid &&
          all(
            finite_probabilities >= 0 &
              finite_probabilities <= 1
          )
      }
      
      variance_nonnegative <-
        variance_nonnegative &&
        (
          is.na(variance_retained_fraction) ||
            variance_retained_fraction >=
            -zero_tolerance
        )
      
      complete_not_above_half <-
        complete_not_above_half &&
        (
          is.na(complete_loss_probability) ||
            is.na(half_loss_probability) ||
            complete_loss_probability <=
            half_loss_probability +
            zero_tolerance
        )
      
      focal_rows[[focal_index]] <- tibble::tibble(
        focal_node =
          focal_metadata$focal_node[focal_index],
        removal_level_index =
          level_index,
        mean_retained_fraction =
          mean(retained_degree_fraction),
        variance_retained_fraction =
          variance_retained_fraction,
        complete_loss_probability =
          complete_loss_probability,
        half_loss_probability =
          half_loss_probability,
        number_focal_present_replicates =
          number_focal_present_replicates,
        conditional_probability_reason =
          conditional_reason
      )
    }
    
    output_list[[level_index]] <-
      dplyr::bind_rows(focal_rows)
  }
  
  list(
    metrics = dplyr::bind_rows(
      output_list
    ),
    checks = tibble::tibble(
      retained_fraction_bounds_valid =
        retained_fraction_bounds_valid,
      probability_bounds_valid =
        probability_bounds_valid,
      variance_nonnegative =
        variance_nonnegative,
      complete_not_above_half =
        complete_not_above_half
    )
  )
}

## ------------------------------------------------------------
## Generate one null arrangement
## ------------------------------------------------------------

generate_null_link_matrix <- function(
    realised_links,
    dataset_sites
){
  
  null_supporting_sites <- Map(
    function(possible_sites, full_K){
      
      sample(
        possible_sites,
        size = full_K,
        replace = FALSE
      )
    },
    realised_links$cooccurrence_sites,
    realised_links$full_K
  )
  
  null_link_matrix <- make_link_matrix(
    supporting_site_sets =
      null_supporting_sites,
    dataset_sites =
      dataset_sites
  )
  
  exact_K_check <- all(
    colSums(null_link_matrix) ==
      realised_links$full_K
  )
  
  within_cooccurrence_check <- all(
    unlist(
      Map(
        function(selected_sites, possible_sites){
          all(selected_sites %in% possible_sites)
        },
        null_supporting_sites,
        realised_links$cooccurrence_sites
      ),
      use.names = FALSE
    )
  )
  
  list(
    link_matrix = null_link_matrix,
    exact_K_check = exact_K_check,
    within_cooccurrence_check =
      within_cooccurrence_check
  )
}

## ------------------------------------------------------------
## Analytical mean retained fraction
## ------------------------------------------------------------

link_survival_probability <- function(
    total_sites,
    full_K,
    sites_retained
){
  
  if(sites_retained >= total_sites){
    return(1)
  }
  
  if(sites_retained > total_sites - full_K){
    return(1)
  }
  
  log_probability_miss_all_support <- lchoose(
    total_sites - full_K,
    sites_retained
  ) - lchoose(
    total_sites,
    sites_retained
  )
  
  1 - exp(
    log_probability_miss_all_support
  )
}

calculate_analytical_mean_retention <- function(
    focal_metadata,
    realised_links,
    total_sites,
    sites_retained
){
  
  vapply(
    seq_len(nrow(focal_metadata)),
    function(focal_index){
      
      link_indices <-
        focal_metadata$link_indices[[focal_index]]
      
      survival_probabilities <- vapply(
        realised_links$full_K[link_indices],
        function(full_K){
          
          link_survival_probability(
            total_sites = total_sites,
            full_K = full_K,
            sites_retained =
              sites_retained
          )
        },
        numeric(1)
      )
      
      mean(survival_probabilities)
    },
    numeric(1)
  )
}

## ------------------------------------------------------------
## Null metric storage and summaries
## ------------------------------------------------------------

metric_names <- c(
  "mean_retained_fraction",
  "variance_retained_fraction",
  "complete_loss_probability",
  "half_loss_probability"
)

make_null_metric_array <- function(
    number_species,
    number_levels,
    n_null,
    focal_nodes
){
  
  array(
    NA_real_,
    dim = c(
      number_species,
      number_levels,
      length(metric_names),
      n_null
    ),
    dimnames = list(
      focal_nodes,
      as.character(removal_levels),
      metric_names,
      NULL
    )
  )
}

insert_null_metrics <- function(
    null_array,
    null_metrics,
    null_index,
    focal_nodes
){
  
  for(metric_name in metric_names){
    
    metric_matrix <- matrix(
      NA_real_,
      nrow = length(focal_nodes),
      ncol = length(removal_levels),
      dimnames = list(
        focal_nodes,
        as.character(removal_levels)
      )
    )
    
    for(row_index in seq_len(nrow(null_metrics))){
      
      focal_node <-
        null_metrics$focal_node[row_index]
      
      level_index <-
        null_metrics$removal_level_index[row_index]
      
      metric_matrix[
        focal_node,
        level_index
      ] <- null_metrics[[metric_name]][row_index]
    }
    
    null_array[
      ,
      ,
      metric_name,
      null_index
    ] <- metric_matrix
  }
  
  null_array
}

summarise_null_metric <- function(
    null_array,
    metric_name
){
  
  number_species <- dim(null_array)[1]
  number_levels <- dim(null_array)[2]
  
  output <- vector(
    "list",
    number_species * number_levels
  )
  
  output_index <- 1L
  
  for(species_index in seq_len(number_species)){
    for(level_index in seq_len(number_levels)){
      
      values <- null_array[
        species_index,
        level_index,
        metric_name,
        ,
        drop = TRUE
      ]
      
      output[[output_index]] <- tibble::tibble(
        focal_node =
          dimnames(null_array)[[1]][species_index],
        removal_level_index =
          level_index,
        metric = metric_name,
        null_mean = safe_mean(values),
        null_sd = safe_sd(values),
        null_q025 = safe_quantile(
          values,
          0.025
        ),
        null_q975 = safe_quantile(
          values,
          0.975
        )
      )
      
      output_index <- output_index + 1L
    }
  }
  
  dplyr::bind_rows(output)
}

## ------------------------------------------------------------
## Analyse one dataset
## ------------------------------------------------------------

run_one_dataset <- function(dataset){
  
  message(
    "Script 38b: ",
    dataset
  )
  
  prepared <- prepare_dataset(
    dataset
  )
  
  realised_links <- prepared$realised_links
  dataset_sites <- prepared$dataset_sites
  
  removal_object <- make_removal_masks(
    sites = dataset_sites,
    removal_levels = removal_levels,
    n_removal_reps = n_removal_reps
  )
  
  observed_link_matrix <- make_link_matrix(
    supporting_site_sets =
      realised_links$interaction_sites,
    dataset_sites =
      dataset_sites
  )
  
  guild_objects <- list()
  
  for(guild_name in guild_levels){
    
    focal_metadata <- make_focal_metadata(
      dataset = dataset,
      guild = guild_name,
      realised_links = realised_links
    )
    
    if(nrow(focal_metadata) == 0){
      
      guild_objects[[guild_name]] <- NULL
      next
    }
    
    representation_table <-
      prepared$focal_representation[[guild_name]]
    
    representation_matrix <-
      make_focal_representation_matrix(
        focal_nodes =
          focal_metadata$focal_node,
        representation_table =
          representation_table,
        dataset_sites =
          dataset_sites
      )
    
    observed_evaluation <- evaluate_arrangement(
      link_matrix =
        observed_link_matrix,
      removal_masks =
        removal_object$masks,
      focal_metadata =
        focal_metadata,
      focal_representation_matrix =
        representation_matrix
    )
    
    null_metric_array <- make_null_metric_array(
      number_species =
        nrow(focal_metadata),
      number_levels =
        length(removal_levels),
      n_null = n_null,
      focal_nodes =
        focal_metadata$focal_node
    )
    
    guild_objects[[guild_name]] <- list(
      focal_metadata = focal_metadata,
      representation_matrix =
        representation_matrix,
      observed = observed_evaluation,
      null_metric_array =
        null_metric_array
    )
  }
  
  all_null_exact_K <- TRUE
  all_null_within_cooccurrence <- TRUE
  all_same_masks <- TRUE
  
  null_retained_fraction_bounds <- TRUE
  null_probability_bounds <- TRUE
  null_variance_nonnegative <- TRUE
  null_complete_not_above_half <- TRUE
  
  ## Generate one null arrangement for the full regional network,
  ## then evaluate both guilds using the same null network and the
  ## same removal masks.
  for(null_index in seq_len(n_null)){
    
    null_arrangement <- generate_null_link_matrix(
      realised_links =
        realised_links,
      dataset_sites =
        dataset_sites
    )
    
    all_null_exact_K <-
      all_null_exact_K &&
      null_arrangement$exact_K_check
    
    all_null_within_cooccurrence <-
      all_null_within_cooccurrence &&
      null_arrangement$within_cooccurrence_check
    
    for(guild_name in guild_levels){
      
      guild_object <-
        guild_objects[[guild_name]]
      
      if(is.null(guild_object)){
        next
      }
      
      null_evaluation <- evaluate_arrangement(
        link_matrix =
          null_arrangement$link_matrix,
        removal_masks =
          removal_object$masks,
        focal_metadata =
          guild_object$focal_metadata,
        focal_representation_matrix =
          guild_object$representation_matrix
      )
      
      guild_objects[[guild_name]]$null_metric_array <- insert_null_metrics(
        null_array =
          guild_object$null_metric_array,
        null_metrics =
          null_evaluation$metrics,
        null_index =
          null_index,
        focal_nodes =
          guild_object$focal_metadata$focal_node
      )
      
      null_retained_fraction_bounds <-
        null_retained_fraction_bounds &&
        null_evaluation$checks$
        retained_fraction_bounds_valid
      
      null_probability_bounds <-
        null_probability_bounds &&
        null_evaluation$checks$
        probability_bounds_valid
      
      null_variance_nonnegative <-
        null_variance_nonnegative &&
        null_evaluation$checks$
        variance_nonnegative
      
      null_complete_not_above_half <-
        null_complete_not_above_half &&
        null_evaluation$checks$
        complete_not_above_half
    }
  }
  
  dataset_species_results <- list()
  
  for(guild_name in guild_levels){
    
    guild_object <- guild_objects[[guild_name]]
    
    if(is.null(guild_object)){
      next
    }
    
    focal_metadata <-
      guild_object$focal_metadata
    
    observed_metrics <-
      guild_object$observed$metrics
    
    null_summary_list <- lapply(
      metric_names,
      function(metric_name){
        
        summarise_null_metric(
          null_array =
            guild_object$null_metric_array,
          metric_name =
            metric_name
        )
      }
    )
    
    names(null_summary_list) <-
      metric_names
    
    removal_metadata <-
      removal_object$metadata
    
    species_grid <- observed_metrics %>%
      dplyr::left_join(
        focal_metadata %>%
          dplyr::select(
            focal_node,
            initial_degree,
            degree_class,
            bundling_delta,
            bundling_SES
          ),
        by = "focal_node"
      ) %>%
      dplyr::left_join(
        removal_metadata,
        by = "removal_level_index"
      )
    
    add_null_summary <- function(
    data,
    null_summary,
    metric_name
    ){
      
      renamed <- null_summary %>%
        dplyr::select(
          -metric
        )
      
      names(renamed)[
        names(renamed) == "null_mean"
      ] <- paste0(
        "null_mean_",
        metric_name
      )
      
      names(renamed)[
        names(renamed) == "null_sd"
      ] <- paste0(
        "null_sd_",
        metric_name
      )
      
      names(renamed)[
        names(renamed) == "null_q025"
      ] <- paste0(
        "null_q025_",
        metric_name
      )
      
      names(renamed)[
        names(renamed) == "null_q975"
      ] <- paste0(
        "null_q975_",
        metric_name
      )
      
      data %>%
        dplyr::left_join(
          renamed,
          by = c(
            "focal_node",
            "removal_level_index"
          )
        )
    }
    
    for(metric_name in metric_names){
      
      species_grid <- add_null_summary(
        data = species_grid,
        null_summary =
          null_summary_list[[metric_name]],
        metric_name =
          metric_name
      )
    }
    
    analytical_rows <- vector(
      "list",
      nrow(removal_metadata)
    )
    
    for(level_index in seq_len(
      nrow(removal_metadata)
    )){
      
      analytical_values <-
        calculate_analytical_mean_retention(
          focal_metadata =
            focal_metadata,
          realised_links =
            realised_links,
          total_sites =
            length(dataset_sites),
          sites_retained =
            removal_metadata$
            actual_sites_retained[level_index]
        )
      
      analytical_rows[[level_index]] <-
        tibble::tibble(
          focal_node =
            focal_metadata$focal_node,
          removal_level_index =
            level_index,
          analytical_expected_retained_fraction =
            analytical_values
        )
    }
    
    analytical_table <- dplyr::bind_rows(
      analytical_rows
    )
    
    presence_definition <-
      prepared$focal_representation[[paste0(
        guild_name,
        "_definition"
      )]]
    
    presence_variable_name <-
      prepared$focal_representation[[paste0(
        guild_name,
        "_variable_name"
      )]]
    
    species_grid <- species_grid %>%
      dplyr::left_join(
        analytical_table,
        by = c(
          "focal_node",
          "removal_level_index"
        )
      ) %>%
      dplyr::mutate(
        dataset = dataset,
        guild = guild_name,
        presence_definition =
          presence_definition,
        conditional_presence_variable =
          presence_variable_name,
        
        observed_mean_retained_fraction =
          mean_retained_fraction,
        
        observed_variance_retained_fraction =
          variance_retained_fraction,
        
        observed_complete_loss_probability =
          complete_loss_probability,
        
        observed_half_loss_probability =
          half_loss_probability,
        
        null_mean_retained_fraction =
          null_mean_mean_retained_fraction,
        
        null_mean_variance_retained_fraction =
          null_mean_variance_retained_fraction,
        
        null_mean_complete_loss_probability =
          null_mean_complete_loss_probability,
        
        null_mean_half_loss_probability =
          null_mean_half_loss_probability,
        
        mean_retention_delta =
          observed_mean_retained_fraction -
          null_mean_retained_fraction,
        
        variance_delta =
          observed_variance_retained_fraction -
          null_mean_variance_retained_fraction,
        
        complete_loss_delta =
          observed_complete_loss_probability -
          null_mean_complete_loss_probability,
        
        half_loss_delta =
          observed_half_loss_probability -
          null_mean_half_loss_probability
      ) %>%
      dplyr::select(
        dataset,
        guild,
        focal_node,
        initial_degree,
        degree_class,
        bundling_delta,
        bundling_SES,
        removal_fraction,
        actual_sites_retained,
        actual_retained_fraction,
        actual_removal_fraction,
        number_dataset_sites,
        number_removal_replicates,
        presence_definition,
        conditional_presence_variable,
        number_focal_present_replicates,
        conditional_probability_reason,
        
        observed_mean_retained_fraction,
        null_mean_retained_fraction,
        null_sd_mean_retained_fraction,
        null_q025_mean_retained_fraction,
        null_q975_mean_retained_fraction,
        analytical_expected_retained_fraction,
        mean_retention_delta,
        
        observed_variance_retained_fraction,
        null_mean_variance_retained_fraction,
        null_sd_variance_retained_fraction,
        null_q025_variance_retained_fraction,
        null_q975_variance_retained_fraction,
        variance_delta,
        
        observed_complete_loss_probability,
        null_mean_complete_loss_probability,
        null_sd_complete_loss_probability,
        null_q025_complete_loss_probability,
        null_q975_complete_loss_probability,
        complete_loss_delta,
        
        observed_half_loss_probability,
        null_mean_half_loss_probability,
        null_sd_half_loss_probability,
        null_q025_half_loss_probability,
        null_q975_half_loss_probability,
        half_loss_delta
      )
    
    dataset_species_results[[guild_name]] <-
      species_grid
  }
  
  observed_check_values <- lapply(
    guild_objects,
    function(x){
      
      if(is.null(x)){
        return(NULL)
      }
      
      x$observed$checks
    }
  ) %>%
    dplyr::bind_rows()
  
  dataset_checks <- tibble::tibble(
    dataset = dataset,
    
    realised_K_positive =
      all(realised_links$full_K > 0),
    
    realised_K_not_above_n =
      all(
        realised_links$full_K <=
          realised_links$full_n
      ),
    
    null_exact_K =
      all_null_exact_K,
    
    null_within_cooccurrence =
      all_null_within_cooccurrence,
    
    same_removal_masks_observed_and_null =
      all_same_masks,
    
    observed_retained_fraction_bounds =
      all(
        observed_check_values$
          retained_fraction_bounds_valid
      ),
    
    null_retained_fraction_bounds =
      null_retained_fraction_bounds,
    
    observed_probability_bounds =
      all(
        observed_check_values$
          probability_bounds_valid
      ),
    
    null_probability_bounds =
      null_probability_bounds,
    
    observed_variance_nonnegative =
      all(
        observed_check_values$
          variance_nonnegative
      ),
    
    null_variance_nonnegative =
      null_variance_nonnegative,
    
    observed_complete_not_above_half =
      all(
        observed_check_values$
          complete_not_above_half
      ),
    
    null_complete_not_above_half =
      null_complete_not_above_half,
    
    used_cooccurrence_fallback =
      prepared$focal_representation$
      used_cooccurrence_fallback,
    
    consumer_presence_definition =
      prepared$focal_representation$
      Consumer_definition,
    
    resource_presence_definition =
      prepared$focal_representation$
      Resource_definition
  )
  
  list(
    species_results = dplyr::bind_rows(
      dataset_species_results
    ),
    dataset_checks = dataset_checks
  )
}

## ------------------------------------------------------------
## Run all datasets
## ------------------------------------------------------------

## Datasets are independent, so they can be evaluated safely in
## separate R sessions. Internal null arrangements and removal
## replicates remain serial within each dataset, avoiding nested
## parallelism and excessive memory duplication.
dataset_outputs <- local({
  
  workers_to_use <- min(
    n_workers,
    length(all_dataset_names)
  )
  
  if(workers_to_use <= 1L){
    
    message(
      "Running datasets serially."
    )
    
    lapply(
      all_dataset_names,
      run_one_dataset
    )
    
  } else {
    
    message(
      "Running ",
      length(all_dataset_names),
      " datasets using ",
      workers_to_use,
      " parallel workers."
    )
    
    previous_plan <- future::plan()
    
    on.exit(
      future::plan(previous_plan),
      add = TRUE
    )
    
    future::plan(
      future::multisession,
      workers = workers_to_use
    )
    
    future.apply::future_lapply(
      all_dataset_names,
      run_one_dataset,
      future.seed = TRUE,
      future.scheduling = 1
    )
  }
})

species_removal <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "species_results"
  )
)

dataset_internal_checks <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "dataset_checks"
  )
)

species_removal <- species_removal %>%
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
## Force exact zero anchors at removal level zero
## ------------------------------------------------------------

species_removal <- species_removal %>%
  dplyr::mutate(
    mean_retention_delta = dplyr::if_else(
      removal_fraction == 0,
      0,
      mean_retention_delta
    ),
    variance_delta = dplyr::if_else(
      removal_fraction == 0,
      0,
      variance_delta
    ),
    complete_loss_delta = dplyr::if_else(
      removal_fraction == 0,
      0,
      complete_loss_delta
    ),
    half_loss_delta = dplyr::if_else(
      removal_fraction == 0,
      0,
      half_loss_delta
    )
  )

## ------------------------------------------------------------
## Normalised trapezoidal AUC
## ------------------------------------------------------------

normalised_trapezoid_auc <- function(
    removal_fraction,
    contrast,
    integration_range = 0.8
){
  
  valid <- is.finite(
    removal_fraction
  ) & is.finite(
    contrast
  )
  
  removal_fraction <- removal_fraction[
    valid
  ]
  
  contrast <- contrast[
    valid
  ]
  
  if(length(removal_fraction) != length(removal_levels)){
    return(NA_real_)
  }
  
  ordering <- order(
    removal_fraction
  )
  
  x <- removal_fraction[
    ordering
  ]
  
  y <- contrast[
    ordering
  ]
  
  if(
    !isTRUE(all.equal(
      x,
      removal_levels,
      tolerance = 1e-10
    ))
  ){
    return(NA_real_)
  }
  
  trapezoid_sum <- sum(
    diff(x) *
      (
        head(y, -1) +
          tail(y, -1)
      ) / 2
  )
  
  trapezoid_sum / integration_range
}

species_auc <- species_removal %>%
  dplyr::group_by(
    dataset,
    guild,
    focal_node,
    initial_degree,
    degree_class,
    bundling_delta,
    bundling_SES
  ) %>%
  dplyr::summarise(
    mean_retention_delta_AUC =
      normalised_trapezoid_auc(
        removal_fraction,
        mean_retention_delta
      ),
    
    variance_delta_AUC =
      normalised_trapezoid_auc(
        removal_fraction,
        variance_delta
      ),
    
    complete_loss_delta_AUC =
      normalised_trapezoid_auc(
        removal_fraction,
        complete_loss_delta
      ),
    
    half_loss_delta_AUC =
      normalised_trapezoid_auc(
        removal_fraction,
        half_loss_delta
      ),
    
    .groups = "drop"
  )

## ------------------------------------------------------------
## Correlation helpers
## ------------------------------------------------------------

ordinary_spearman <- function(x, y){
  
  suppressWarnings(
    stats::cor(
      x,
      y,
      method = "spearman",
      use = "complete.obs"
    )
  )
}

partial_spearman_degree <- function(
    bundling_delta,
    consequence,
    initial_degree
){
  
  rank_bundling <- rank(
    bundling_delta,
    ties.method = "average"
  )
  
  rank_consequence <- rank(
    consequence,
    ties.method = "average"
  )
  
  rank_degree <- rank(
    initial_degree,
    ties.method = "average"
  )
  
  bundling_model <- stats::lm(
    rank_bundling ~ rank_degree
  )
  
  consequence_model <- stats::lm(
    rank_consequence ~ rank_degree
  )
  
  bundling_residuals <- stats::residuals(
    bundling_model
  )
  
  consequence_residuals <- stats::residuals(
    consequence_model
  )
  
  if(
    stats::sd(bundling_residuals) == 0 ||
    stats::sd(consequence_residuals) == 0
  ){
    return(NA_real_)
  }
  
  stats::cor(
    bundling_residuals,
    consequence_residuals,
    method = "pearson"
  )
}

make_dataset_correlation <- function(
    data,
    dataset,
    guild,
    outcome,
    correlation_analysis
){
  
  consequence_column <- switch(
    outcome,
    "variance" = "variance_delta_AUC",
    "complete_loss" =
      "complete_loss_delta_AUC",
    "half_loss" =
      "half_loss_delta_AUC"
  )
  
  analysis_data <- data %>%
    dplyr::filter(
      is.finite(bundling_delta),
      is.finite(.data[[consequence_column]])
    )
  
  number_species <- nrow(
    analysis_data
  )
  
  number_unique_bundling_values <-
    dplyr::n_distinct(
      analysis_data$bundling_delta
    )
  
  number_unique_consequence_values <-
    dplyr::n_distinct(
      analysis_data[[consequence_column]]
    )
  
  number_unique_degree_values <-
    dplyr::n_distinct(
      analysis_data$initial_degree
    )
  
  common_validity <-
    number_species >= 5 &&
    number_unique_bundling_values >= 3 &&
    number_unique_consequence_values >= 3
  
  partial_degree_validity <-
    number_unique_degree_values >= 2
  
  valid_for_inference <- if(
    correlation_analysis ==
    "ordinary_spearman"
  ){
    common_validity
  } else {
    common_validity &&
      partial_degree_validity
  }
  
  exclusion_reason <- dplyr::case_when(
    number_species < 5 ~
      "Fewer than 5 eligible species",
    
    number_unique_bundling_values < 3 ~
      "Fewer than 3 unique bundling-delta values",
    
    number_unique_consequence_values < 3 ~
      "Fewer than 3 unique consequence-AUC values",
    
    correlation_analysis ==
      "partial_spearman_controlling_degree" &&
      number_unique_degree_values < 2 ~
      "Fewer than 2 unique initial-degree values",
    
    TRUE ~
      NA_character_
  )
  
  correlation_value <- if(
    valid_for_inference
  ){
    
    if(
      correlation_analysis ==
      "ordinary_spearman"
    ){
      
      ordinary_spearman(
        analysis_data$bundling_delta,
        analysis_data[[consequence_column]]
      )
      
    } else {
      
      partial_spearman_degree(
        bundling_delta =
          analysis_data$bundling_delta,
        consequence =
          analysis_data[[consequence_column]],
        initial_degree =
          analysis_data$initial_degree
      )
    }
    
  } else {
    NA_real_
  }
  
  if(!is.finite(correlation_value)){
    valid_for_inference <- FALSE
    
    if(is.na(exclusion_reason)){
      exclusion_reason <-
        "Correlation undefined after calculation"
    }
  }
  
  tibble::tibble(
    correlation_analysis =
      correlation_analysis,
    dataset = dataset,
    guild = guild,
    outcome = outcome,
    number_species =
      number_species,
    minimum_degree = if(number_species > 0){
      min(analysis_data$initial_degree)
    } else {
      NA_real_
    },
    maximum_degree = if(number_species > 0){
      max(analysis_data$initial_degree)
    } else {
      NA_real_
    },
    number_unique_degree_values =
      number_unique_degree_values,
    number_unique_bundling_values =
      number_unique_bundling_values,
    number_unique_consequence_values =
      number_unique_consequence_values,
    ordinary_spearman_rho =
      ifelse(
        correlation_analysis ==
          "ordinary_spearman",
        correlation_value,
        NA_real_
      ),
    partial_spearman_rho_controlling_degree =
      ifelse(
        correlation_analysis ==
          "partial_spearman_controlling_degree",
        correlation_value,
        NA_real_
      ),
    correlation_value =
      correlation_value,
    valid_for_inference =
      valid_for_inference,
    exclusion_reason =
      exclusion_reason
  )
}

## ------------------------------------------------------------
## Dataset-specific correlations
## ------------------------------------------------------------

species_dataset_groups <- species_auc %>%
  dplyr::group_by(
    dataset,
    guild
  ) %>%
  dplyr::group_split(
    .keep = TRUE
  )

dataset_correlations <- lapply(
  species_dataset_groups,
  function(dataset_group){
    
    dataset_name <- as.character(
      dataset_group$dataset[1]
    )
    
    guild_name <- as.character(
      dataset_group$guild[1]
    )
    
    result_rows <- list()
    row_index <- 1L
    
    for(correlation_analysis in
        correlation_analysis_levels){
      
      for(outcome in outcome_levels){
        
        result_rows[[row_index]] <-
          make_dataset_correlation(
            data = dataset_group,
            dataset = dataset_name,
            guild = guild_name,
            outcome = outcome,
            correlation_analysis =
              correlation_analysis
          )
        
        row_index <- row_index + 1L
      }
    }
    
    dplyr::bind_rows(result_rows)
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    correlation_analysis = factor(
      correlation_analysis,
      levels =
        correlation_analysis_levels
    ),
    guild = factor(
      guild,
      levels = guild_levels
    ),
    outcome = factor(
      outcome,
      levels = outcome_levels
    )
  )

## ------------------------------------------------------------
## Equal-dataset inference
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
        equal_dataset_mean_correlation =
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
    
    signed_sums <- 0
    
    for(value in x){
      
      signed_sums <- c(
        signed_sums + value,
        signed_sums - value
      )
    }
    
    null_means <- signed_sums / n
    
    sign_flip_p <- mean(
      abs(null_means) >=
        abs(observed_mean) -
        1e-12
    )
    
    test_type <- paste0(
      "Exact sign-flip: ",
      2^n,
      " sign patterns"
    )
    
  } else {
    
    null_means <- replicate(
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
          abs(null_means) >=
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
    equal_dataset_mean_correlation =
      observed_mean,
    sign_flip_p =
      sign_flip_p,
    test_type =
      test_type
  )
}

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
      names = FALSE,
      na.rm = TRUE
    ),
    bootstrap_q975 = stats::quantile(
      bootstrap_means,
      0.975,
      names = FALSE,
      na.rm = TRUE
    )
  )
}

equal_dataset_tests <- dataset_correlations %>%
  dplyr::filter(
    valid_for_inference,
    is.finite(correlation_value)
  ) %>%
  dplyr::group_by(
    correlation_analysis,
    guild,
    outcome
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      values <- .x$correlation_value
      
      dplyr::bind_cols(
        sign_flip_test(
          values,
          n_large =
            n_permutations_large_n
        ),
        bootstrap_dataset_mean(
          values,
          n_boot = n_boot
        )
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      equal_dataset_mean_correlation > 0 ~
        "More bundling associated with more extreme loss",
      equal_dataset_mean_correlation < 0 ~
        "More bundling associated with spatial insurance",
      equal_dataset_mean_correlation == 0 ~
        "No mean association",
      TRUE ~
        NA_character_
    )
  )

all_test_combinations <- expand.grid(
  correlation_analysis =
    correlation_analysis_levels,
  guild = guild_levels,
  outcome = outcome_levels,
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble()

equal_dataset_tests <-
  all_test_combinations %>%
  dplyr::left_join(
    equal_dataset_tests,
    by = c(
      "correlation_analysis",
      "guild",
      "outcome"
    )
  )

## BH correction across the six ordinary primary tests only.
ordinary_indices <- which(
  equal_dataset_tests$
    correlation_analysis ==
    "ordinary_spearman"
)

equal_dataset_tests$sign_flip_p_BH <-
  NA_real_

equal_dataset_tests$sign_flip_p_BH[
  ordinary_indices
] <- stats::p.adjust(
  equal_dataset_tests$sign_flip_p[
    ordinary_indices
  ],
  method = "BH"
)

equal_dataset_tests <-
  equal_dataset_tests %>%
  dplyr::mutate(
    significant_raw_0.05 =
      !is.na(sign_flip_p) &
      sign_flip_p < 0.05,
    
    significant_BH_0.05 =
      !is.na(sign_flip_p_BH) &
      sign_flip_p_BH < 0.05,
    
    correlation_analysis = factor(
      correlation_analysis,
      levels =
        correlation_analysis_levels
    ),
    
    guild = factor(
      guild,
      levels = guild_levels
    ),
    
    outcome = factor(
      outcome,
      levels = outcome_levels
    )
  ) %>%
  dplyr::arrange(
    correlation_analysis,
    outcome,
    guild
  )

## ------------------------------------------------------------
## Descriptive degree-class summaries
## ------------------------------------------------------------

species_contrast_long <- dplyr::bind_rows(
  
  species_removal %>%
    dplyr::transmute(
      dataset,
      guild,
      focal_node,
      degree_class,
      removal_fraction,
      outcome = "mean_retention",
      contrast = mean_retention_delta
    ),
  
  species_removal %>%
    dplyr::transmute(
      dataset,
      guild,
      focal_node,
      degree_class,
      removal_fraction,
      outcome = "variance",
      contrast = variance_delta
    ),
  
  species_removal %>%
    dplyr::transmute(
      dataset,
      guild,
      focal_node,
      degree_class,
      removal_fraction,
      outcome = "complete_loss",
      contrast = complete_loss_delta
    ),
  
  species_removal %>%
    dplyr::transmute(
      dataset,
      guild,
      focal_node,
      degree_class,
      removal_fraction,
      outcome = "half_loss",
      contrast = half_loss_delta
    )
)

dataset_degree_class_curves <-
  species_contrast_long %>%
  dplyr::filter(
    is.finite(contrast)
  ) %>%
  dplyr::group_by(
    dataset,
    guild,
    degree_class,
    removal_fraction,
    outcome
  ) %>%
  dplyr::summarise(
    dataset_mean_contrast =
      mean(contrast),
    number_species =
      dplyr::n_distinct(focal_node),
    .groups = "drop"
  )

bootstrap_class_summary <- function(
    values,
    n_boot = 10000
){
  
  values <- values[
    is.finite(values)
  ]
  
  n <- length(values)
  
  if(n == 0){
    
    return(
      tibble::tibble(
        number_datasets = 0L,
        equal_dataset_mean_contrast =
          NA_real_,
        bootstrap_q025 = NA_real_,
        bootstrap_q975 = NA_real_
      )
    )
  }
  
  equal_dataset_mean <- mean(values)
  
  if(n < 2){
    
    return(
      tibble::tibble(
        number_datasets = n,
        equal_dataset_mean_contrast =
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
        values,
        size = n,
        replace = TRUE
      )
    )
  )
  
  tibble::tibble(
    number_datasets = n,
    equal_dataset_mean_contrast =
      equal_dataset_mean,
    bootstrap_q025 = stats::quantile(
      bootstrap_means,
      0.025,
      names = FALSE
    ),
    bootstrap_q975 = stats::quantile(
      bootstrap_means,
      0.975,
      names = FALSE
    )
  )
}

equal_dataset_degree_class_summaries <-
  dataset_degree_class_curves %>%
  dplyr::group_by(
    guild,
    degree_class,
    removal_fraction,
    outcome
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      bootstrap_class_summary(
        .x$dataset_mean_contrast,
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
    ),
    outcome = factor(
      outcome,
      levels = contrast_levels
    )
  )

## ------------------------------------------------------------
## Validation checks
## ------------------------------------------------------------

zero_removal_rows <- species_removal %>%
  dplyr::filter(
    removal_fraction == 0
  )

mean_retention_errors <- species_removal %>%
  dplyr::transmute(
    observed_error = abs(
      observed_mean_retained_fraction -
        analytical_expected_retained_fraction
    ),
    null_error = abs(
      null_mean_retained_fraction -
        analytical_expected_retained_fraction
    )
  )

mean_abs_negative_control <- mean(
  abs(
    species_removal$mean_retention_delta
  ),
  na.rm = TRUE
)

mean_abs_primary_contrasts <- mean(
  abs(
    c(
      species_removal$variance_delta,
      species_removal$complete_loss_delta,
      species_removal$half_loss_delta
    )
  ),
  na.rm = TRUE
)

one_correlation_per_dataset <- dataset_correlations %>%
  dplyr::count(
    correlation_analysis,
    dataset,
    guild,
    outcome,
    name = "number_rows"
  ) %>%
  dplyr::summarise(
    valid = all(number_rows == 1)
  ) %>%
  dplyr::pull(valid)

equal_dataset_unweighted_check <- all(
  vapply(
    seq_len(
      nrow(equal_dataset_tests)
    ),
    function(i){
      
      analysis_i <- as.character(
        equal_dataset_tests$
          correlation_analysis[i]
      )
      
      guild_i <- as.character(
        equal_dataset_tests$guild[i]
      )
      
      outcome_i <- as.character(
        equal_dataset_tests$outcome[i]
      )
      
      correlations <- dataset_correlations %>%
        dplyr::filter(
          correlation_analysis ==
            analysis_i,
          guild == guild_i,
          outcome == outcome_i,
          valid_for_inference,
          is.finite(correlation_value)
        ) %>%
        dplyr::pull(
          correlation_value
        )
      
      reported_mean <-
        equal_dataset_tests$
        equal_dataset_mean_correlation[i]
      
      if(length(correlations) == 0){
        return(is.na(reported_mean))
      }
      
      isTRUE(
        abs(
          mean(correlations) -
            reported_mean
        ) < 1e-12
      )
    },
    logical(1)
  )
)

degree_class_within_dataset_check <- all(
  dataset_degree_class_curves %>%
    dplyr::count(
      dataset,
      guild,
      degree_class,
      removal_fraction,
      outcome,
      name = "number_rows"
    ) %>%
    dplyr::pull(number_rows) == 1
)

presence_definitions_complete <- all(
  !is.na(
    dataset_internal_checks$
      consumer_presence_definition
  ) &
    !is.na(
      dataset_internal_checks$
        resource_presence_definition
    )
)

validation_checks <- tibble::tibble(
  check = c(
    "Every realised interaction has full_K > 0",
    "Every realised interaction has full_K <= full_n",
    "Every null interaction uses exactly full_K sites",
    "Every null interaction site belongs to its pair-specific co-occurrence set",
    "The same removal masks are used for observed and null arrangements",
    "Every retained-degree fraction lies between zero and one",
    "Every complete- and half-loss probability lies between zero and one",
    "Every retained-fraction variance is non-negative",
    "Complete-loss probability is never greater than half-loss probability",
    "Removal level zero has full retention and zero contrasts",
    "Observed mean retained fraction approximates analytical expectation",
    "Null mean retained fraction approximates analytical expectation",
    "Mean-retention contrast is small relative to consequence contrasts",
    "Consumers and Resources are analysed separately",
    "Each correlation contains one dataset and one guild only",
    "Equal-dataset inference uses one correlation per dataset",
    "Datasets are not weighted by their number of species",
    "Degree-class summaries average within datasets before averaging across datasets",
    "Species-presence or focal-representation definition is recorded for every dataset"
  ),
  
  pass = c(
    all(
      dataset_internal_checks$
        realised_K_positive
    ),
    
    all(
      dataset_internal_checks$
        realised_K_not_above_n
    ),
    
    all(
      dataset_internal_checks$
        null_exact_K
    ),
    
    all(
      dataset_internal_checks$
        null_within_cooccurrence
    ),
    
    all(
      dataset_internal_checks$
        same_removal_masks_observed_and_null
    ),
    
    all(
      dataset_internal_checks$
        observed_retained_fraction_bounds
    ) &&
      all(
        dataset_internal_checks$
          null_retained_fraction_bounds
      ),
    
    all(
      dataset_internal_checks$
        observed_probability_bounds
    ) &&
      all(
        dataset_internal_checks$
          null_probability_bounds
      ),
    
    all(
      dataset_internal_checks$
        observed_variance_nonnegative
    ) &&
      all(
        dataset_internal_checks$
          null_variance_nonnegative
      ),
    
    all(
      dataset_internal_checks$
        observed_complete_not_above_half
    ) &&
      all(
        dataset_internal_checks$
          null_complete_not_above_half
      ),
    
    all(
      abs(
        zero_removal_rows$
          observed_mean_retained_fraction -
          1
      ) <= zero_tolerance
    ) &&
      all(
        abs(
          zero_removal_rows$
            null_mean_retained_fraction -
            1
        ) <= zero_tolerance
      ) &&
      all(
        abs(
          zero_removal_rows$
            mean_retention_delta
        ) <= zero_tolerance
      ) &&
      all(
        abs(
          zero_removal_rows$
            variance_delta
        ) <= zero_tolerance
      ) &&
      all(
        abs(
          zero_removal_rows$
            complete_loss_delta
        ) <= zero_tolerance
      ) &&
      all(
        abs(
          zero_removal_rows$
            half_loss_delta
        ) <= zero_tolerance
      ),
    
    all(
      mean_retention_errors$
        observed_error <=
        mean_retention_tolerance,
      na.rm = TRUE
    ),
    
    all(
      mean_retention_errors$
        null_error <=
        mean_retention_tolerance,
      na.rm = TRUE
    ),
    
    mean_abs_negative_control <= max(
      mean_retention_tolerance,
      mean_abs_primary_contrasts
    ),
    
    all(
      guild_levels %in%
        as.character(
          unique(
            species_removal$guild
          )
        )
    ),
    
    one_correlation_per_dataset,
    
    one_correlation_per_dataset,
    
    equal_dataset_unweighted_check,
    
    degree_class_within_dataset_check,
    
    presence_definitions_complete
  ),
  
  details = c(
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "One removal-mask object per dataset was passed unchanged to every observed and null evaluation.",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    paste0(
      "Numerical zero tolerance = ",
      zero_tolerance
    ),
    paste0(
      "Monte Carlo tolerance = ",
      mean_retention_tolerance
    ),
    paste0(
      "Monte Carlo tolerance = ",
      mean_retention_tolerance
    ),
    paste0(
      "Mean absolute negative-control contrast = ",
      signif(
        mean_abs_negative_control,
        4
      ),
      "; mean absolute primary contrast = ",
      signif(
        mean_abs_primary_contrasts,
        4
      )
    ),
    NA_character_,
    NA_character_,
    NA_character_,
    "Equal-dataset means are arithmetic means of dataset-specific correlations.",
    "Species are first averaged within dataset × guild × degree class × removal level × outcome.",
    paste(
      dataset_internal_checks$dataset,
      dataset_internal_checks$
        consumer_presence_definition,
      dataset_internal_checks$
        resource_presence_definition,
      sep = ": ",
      collapse = " | "
    )
  )
)

## ------------------------------------------------------------
## Main inferential figure
## ------------------------------------------------------------

ordinary_correlations <- dataset_correlations %>%
  dplyr::filter(
    correlation_analysis ==
      "ordinary_spearman",
    valid_for_inference,
    is.finite(correlation_value)
  ) %>%
  dplyr::mutate(
    outcome_label = factor(
      unname(
        outcome_labels[
          as.character(outcome)
        ]
      ),
      levels = unname(
        outcome_labels[outcome_levels]
      )
    )
  )

ordinary_equal_tests <- equal_dataset_tests %>%
  dplyr::filter(
    correlation_analysis ==
      "ordinary_spearman"
  ) %>%
  dplyr::mutate(
    outcome_label = factor(
      unname(
        outcome_labels[
          as.character(outcome)
        ]
      ),
      levels = unname(
        outcome_labels[outcome_levels]
      )
    )
  )

bundling_vulnerability_plot <- ggplot2::ggplot() +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.4
  ) +
  
  ggplot2::geom_jitter(
    data = ordinary_correlations,
    ggplot2::aes(
      x = guild,
      y = correlation_value
    ),
    width = 0.08,
    height = 0,
    colour = "grey55",
    alpha = 0.8,
    size = 1.9
  ) +
  
  ggplot2::geom_errorbar(
    data = ordinary_equal_tests,
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
    data = ordinary_equal_tests,
    ggplot2::aes(
      x = guild,
      y =
        equal_dataset_mean_correlation
    ),
    colour = "black",
    size = 3.3,
    na.rm = TRUE
  ) +
  
  ggplot2::facet_grid(
    outcome_label ~ .,
    scales = "fixed"
  ) +
  
  ggplot2::scale_x_discrete(
    limits = guild_levels,
    drop = FALSE
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(-1, 1)
  ) +
  
  theme_pub(base_size = 10) +
  
  ggplot2::labs(
    x = NULL,
    y = paste0(
      "Spearman correlation between bundling delta ",
      "and removal consequence"
    ),
    title =
      "Spatial bundling predicts the extremity of interaction loss",
    subtitle = paste0(
      "Positive values mean that more bundled portfolios experience ",
      "more variable or extreme losses than expected"
    )
  )

## ------------------------------------------------------------
## Descriptive variance figure
## ------------------------------------------------------------

variance_degree_summary <-
  equal_dataset_degree_class_summaries %>%
  dplyr::filter(
    outcome == "variance"
  )

variance_degree_plot <- ggplot2::ggplot(
  variance_degree_summary,
  ggplot2::aes(
    x = removal_fraction,
    y = equal_dataset_mean_contrast,
    colour = degree_class,
    fill = degree_class,
    group = degree_class
  )
) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.4
  ) +
  
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = bootstrap_q025,
      ymax = bootstrap_q975
    ),
    alpha = 0.10,
    colour = NA,
    na.rm = TRUE
  ) +
  
  ggplot2::geom_line(
    linewidth = 1,
    na.rm = TRUE
  ) +
  
  ggplot2::geom_point(
    size = 1.8,
    na.rm = TRUE
  ) +
  
  ggplot2::facet_wrap(
    ~ guild,
    nrow = 1
  ) +
  
  ggplot2::scale_colour_manual(
    values = degree_class_colours,
    drop = FALSE
  ) +
  
  ggplot2::scale_fill_manual(
    values = degree_class_colours,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = removal_levels,
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  theme_pub(base_size = 10) +
  
  ggplot2::labs(
    x = "Sites removed",
    y = paste0(
      "Observed minus null variance in ",
      "retained-degree fraction"
    ),
    colour = "Initial degree",
    fill = "Initial degree",
    title =
      "Spatial arrangement changes variability in portfolio loss",
    subtitle = paste0(
      "Positive values indicate more variable losses; ",
      "negative values indicate spatial insurance"
    )
  )

## ------------------------------------------------------------
## Save requested CSV files
## ------------------------------------------------------------

write.csv2(
  species_removal,
  file.path(
    out_dir,
    "38b_species_level_removal_consequences.csv"
  ),
  row.names = FALSE
)

write.csv2(
  species_auc,
  file.path(
    out_dir,
    "38b_species_level_consequence_AUC.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_correlations,
  file.path(
    out_dir,
    "38b_dataset_bundling_consequence_correlations.csv"
  ),
  row.names = FALSE
)

write.csv2(
  equal_dataset_tests,
  file.path(
    out_dir,
    "38b_equal_dataset_bundling_consequence_tests.csv"
  ),
  row.names = FALSE
)

write.csv2(
  equal_dataset_degree_class_summaries,
  file.path(
    out_dir,
    "38b_equal_dataset_degree_class_removal_summaries.csv"
  ),
  row.names = FALSE
)

write.csv2(
  validation_checks,
  file.path(
    out_dir,
    "38b_spatial_bundling_removal_validation_checks.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Save requested PNG files
## ------------------------------------------------------------

save_png(
  bundling_vulnerability_plot,
  "38b_bundling_vulnerability_correlations.png",
  width = 8,
  height = 9
)

save_png(
  variance_degree_plot,
  "38b_variance_consequence_by_degree.png",
  width = 10,
  height = 5.5
)

## ------------------------------------------------------------
## Console summaries
## ------------------------------------------------------------

message(
  "\nEqual-dataset ordinary Spearman results:"
)

print(
  equal_dataset_tests %>%
    dplyr::filter(
      correlation_analysis ==
        "ordinary_spearman"
    ) %>%
    dplyr::select(
      guild,
      outcome,
      number_valid_datasets,
      equal_dataset_mean_correlation,
      bootstrap_q025,
      bootstrap_q975,
      sign_flip_p,
      sign_flip_p_BH,
      significant_raw_0.05,
      significant_BH_0.05,
      direction
    ),
  n = Inf
)

message(
  "\nDegree-adjusted partial-Spearman sensitivity results:"
)

print(
  equal_dataset_tests %>%
    dplyr::filter(
      correlation_analysis ==
        "partial_spearman_controlling_degree"
    ) %>%
    dplyr::select(
      guild,
      outcome,
      number_valid_datasets,
      equal_dataset_mean_correlation,
      bootstrap_q025,
      bootstrap_q975,
      sign_flip_p,
      direction
    ),
  n = Inf
)

dataset_validity_counts <- dataset_correlations %>%
  dplyr::group_by(
    correlation_analysis,
    guild,
    outcome
  ) %>%
  dplyr::summarise(
    number_valid_datasets =
      sum(valid_for_inference),
    number_excluded_datasets =
      sum(!valid_for_inference),
    total_dataset_guild_combinations =
      dplyr::n(),
    .groups = "drop"
  )

message(
  "\nValid and excluded datasets:"
)

print(
  dataset_validity_counts,
  n = Inf
)

negative_control_summary <- species_removal %>%
  dplyr::group_by(
    guild,
    removal_fraction
  ) %>%
  dplyr::summarise(
    equal_species_mean_retention_delta =
      mean(
        mean_retention_delta,
        na.rm = TRUE
      ),
    median_absolute_mean_retention_delta =
      median(
        abs(mean_retention_delta),
        na.rm = TRUE
      ),
    .groups = "drop"
  )

message(
  "\nMean-retention negative-control results:"
)

print(
  negative_control_summary,
  n = Inf
)

analytical_diagnostics <- species_removal %>%
  dplyr::group_by(
    dataset,
    guild,
    removal_fraction
  ) %>%
  dplyr::summarise(
    mean_observed_minus_analytical =
      mean(
        observed_mean_retained_fraction -
          analytical_expected_retained_fraction,
        na.rm = TRUE
      ),
    mean_null_minus_analytical =
      mean(
        null_mean_retained_fraction -
          analytical_expected_retained_fraction,
        na.rm = TRUE
      ),
    maximum_absolute_observed_error =
      max(
        abs(
          observed_mean_retained_fraction -
            analytical_expected_retained_fraction
        ),
        na.rm = TRUE
      ),
    maximum_absolute_null_error =
      max(
        abs(
          null_mean_retained_fraction -
            analytical_expected_retained_fraction
        ),
        na.rm = TRUE
      ),
    .groups = "drop"
  )

message(
  "\nAnalytical versus simulated mean-retention diagnostics:"
)

print(
  analytical_diagnostics,
  n = Inf
)

message(
  "\nValidation checks:"
)

print(
  validation_checks,
  n = Inf
)

message(
  "\nInterpretation:"
)

message(
  "Positive variance delta means that the observed spatial arrangement\n",
  "causes interaction losses to be more variable among random site-removal\n",
  "iterations than expected after preserving every link's support and\n",
  "co-occurrence opportunities.\n\n",
  "Positive complete-loss or half-loss delta means that the observed\n",
  "spatial arrangement produces extreme portfolio losses more frequently\n",
  "than expected.\n\n",
  "Negative values indicate spatial insurance: the interactions are\n",
  "distributed among sites in a way that reduces simultaneous loss.\n\n",
  "The mean retained portfolio is expected to be similar between observed\n",
  "and null arrangements. The main information lies in the variance and\n",
  "tails of the loss distribution, not its mean."
)

message(
  "\nSaved Script 38b outputs in: ",
  out_dir
)