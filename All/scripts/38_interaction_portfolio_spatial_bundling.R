## ------------------------------------------------------------
## Script:
## All/scripts/38_interaction_portfolio_spatial_bundling_parallel.R
##
## Purpose:
## Test whether the realised interactions belonging to the same
## focal species are spatially bundled into the same local sites,
## beyond expectations based on:
##
##   1. each interaction's observed local support, full_K; and
##   2. the exact pair-specific sites where the interaction could occur.
##
## Consumers and resources are analysed separately.
##
## Primary species-level response:
##
##   bundling_delta =
##     observed_mean_jaccard - null_mean_jaccard
##
## Primary inferential unit:
##   one dataset-level mean bundling delta per dataset and guild.
##
## Parallelisation:
##   datasets are processed in parallel using future.apply.
##   Null replicates within focal species remain sequential.
##
## Outputs:
##   - two PNG figures
##   - four CSV files
##
## No site removal, degree-group comparison, degree-retention analysis,
## portfolio-collapse analysis, or alternative bundling metrics.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "purrr",
  "future",
  "future.apply",
  "parallelly"
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

n_null <- 1000
n_boot <- 10000
n_permutations_large_n <- 99999

## Parallelise across datasets.
use_parallel <- TRUE

n_workers <- max(
  1,
  min(
    length(all_dataset_names),
    parallelly::availableCores() - 1
  )
)

out_dir <- "All/outputs/38_interaction_portfolio_spatial_bundling"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

guild_levels <- c(
  "Consumer",
  "Resource"
)

## ------------------------------------------------------------
## Parallel plan
## ------------------------------------------------------------

if(use_parallel && n_workers > 1){
  
  message(
    "Using parallel processing across datasets with ",
    n_workers,
    " workers."
  )
  
  future::plan(
    future::multisession,
    workers = n_workers
  )
  
} else {
  
  message("Using sequential processing.")
  
  future::plan(
    future::sequential
  )
}

## ------------------------------------------------------------
## Plotting helpers
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

replace_na_int <- function(
    x,
    value = 0L
){
  
  x[is.na(x)] <- value
  x
}

## ------------------------------------------------------------
## Jaccard helpers
## ------------------------------------------------------------

jaccard_sets <- function(
    set_a,
    set_b
){
  
  set_a <- unique(
    as.character(set_a)
  )
  
  set_b <- unique(
    as.character(set_b)
  )
  
  union_size <- length(
    union(
      set_a,
      set_b
    )
  )
  
  if(union_size == 0){
    return(NA_real_)
  }
  
  intersection_size <- length(
    intersect(
      set_a,
      set_b
    )
  )
  
  intersection_size / union_size
}

mean_pairwise_jaccard <- function(site_sets){
  
  number_sets <- length(site_sets)
  
  if(number_sets < 2){
    return(NA_real_)
  }
  
  partner_pairs <- utils::combn(
    seq_len(number_sets),
    2
  )
  
  pairwise_values <- apply(
    partner_pairs,
    2,
    function(pair_index){
      
      jaccard_sets(
        site_sets[[pair_index[1]]],
        site_sets[[pair_index[2]]]
      )
    }
  )
  
  mean(
    pairwise_values,
    na.rm = TRUE
  )
}

## ------------------------------------------------------------
## Null model for one focal species
## ------------------------------------------------------------

analyse_focal_species <- function(
    focal_data,
    dataset,
    guild,
    focal_node,
    focal_species_cooccurrence_occupancy,
    n_null
){
  
  focal_data <- focal_data %>%
    dplyr::arrange(partner_node)
  
  initial_degree <- dplyr::n_distinct(
    focal_data$partner_node
  )
  
  if(initial_degree < 2){
    
    stop(
      "analyse_focal_species() received a focal species with degree < 2: ",
      dataset,
      " / ",
      guild,
      " / ",
      focal_node
    )
  }
  
  cooccurrence_sets <- focal_data$cooccurrence_sites
  interaction_sets <- focal_data$interaction_sites
  interaction_support <- focal_data$full_K
  
  number_partner_pairs <- choose(
    initial_degree,
    2
  )
  
  observed_mean_jaccard <- mean_pairwise_jaccard(
    interaction_sets
  )
  
  interaction_footprint_sites <- length(
    unique(
      unlist(
        interaction_sets,
        use.names = FALSE
      )
    )
  )
  
  total_local_interaction_records <- sum(
    interaction_support,
    na.rm = TRUE
  )
  
  null_values <- numeric(n_null)
  
  all_null_allocations_preserve_K <- TRUE
  all_null_sites_inside_pair_cooccurrence <- TRUE
  
  for(null_rep in seq_len(n_null)){
    
    null_interaction_sets <- Map(
      function(possible_sites, observed_K){
        
        selected_sites <- sample(
          possible_sites,
          size = observed_K,
          replace = FALSE
        )
        
        if(length(selected_sites) != observed_K){
          all_null_allocations_preserve_K <<- FALSE
        }
        
        if(!all(selected_sites %in% possible_sites)){
          all_null_sites_inside_pair_cooccurrence <<- FALSE
        }
        
        selected_sites
      },
      cooccurrence_sets,
      interaction_support
    )
    
    null_values[null_rep] <- mean_pairwise_jaccard(
      null_interaction_sets
    )
  }
  
  null_mean_jaccard <- mean(
    null_values,
    na.rm = TRUE
  )
  
  null_sd_jaccard <- stats::sd(
    null_values,
    na.rm = TRUE
  )
  
  null_q025_jaccard <- stats::quantile(
    null_values,
    0.025,
    na.rm = TRUE,
    names = FALSE
  )
  
  null_q975_jaccard <- stats::quantile(
    null_values,
    0.975,
    na.rm = TRUE,
    names = FALSE
  )
  
  bundling_delta <-
    observed_mean_jaccard -
    null_mean_jaccard
  
  bundling_SES <- if(
    is.na(null_sd_jaccard) ||
    null_sd_jaccard == 0
  ){
    NA_real_
  } else {
    bundling_delta / null_sd_jaccard
  }
  
  observed_jaccard_valid <-
    is.finite(observed_mean_jaccard) &&
    observed_mean_jaccard >= 0 &&
    observed_mean_jaccard <= 1
  
  null_jaccards_valid <- all(
    is.finite(null_values) &
      null_values >= 0 &
      null_values <= 1
  )
  
  tibble::tibble(
    dataset = dataset,
    guild = guild,
    focal_node = focal_node,
    initial_degree = initial_degree,
    focal_species_cooccurrence_occupancy =
      focal_species_cooccurrence_occupancy,
    interaction_footprint_sites =
      interaction_footprint_sites,
    total_local_interaction_records =
      total_local_interaction_records,
    number_partner_pairs =
      number_partner_pairs,
    observed_mean_jaccard =
      observed_mean_jaccard,
    null_mean_jaccard =
      null_mean_jaccard,
    null_sd_jaccard =
      null_sd_jaccard,
    null_q025_jaccard =
      null_q025_jaccard,
    null_q975_jaccard =
      null_q975_jaccard,
    bundling_delta =
      bundling_delta,
    bundling_SES =
      bundling_SES,
    
    ## Internal validation fields.
    observed_jaccard_valid =
      observed_jaccard_valid,
    null_jaccards_valid =
      null_jaccards_valid,
    null_allocations_preserve_K =
      all_null_allocations_preserve_K,
    null_sites_inside_pair_cooccurrence =
      all_null_sites_inside_pair_cooccurrence
  )
}

## ------------------------------------------------------------
## Guild-specific table construction
## ------------------------------------------------------------

make_guild_link_table <- function(
    pair_table,
    guild_name
){
  
  if(guild_name == "Consumer"){
    
    pair_table %>%
      dplyr::transmute(
        guild = "Consumer",
        focal_node = consumer,
        partner_node = resource,
        full_n,
        full_K,
        cooccurrence_sites,
        interaction_sites
      )
    
  } else {
    
    pair_table %>%
      dplyr::transmute(
        guild = "Resource",
        focal_node = resource,
        partner_node = consumer,
        full_n,
        full_K,
        cooccurrence_sites,
        interaction_sites
      )
  }
}

make_focal_occupancy_table <- function(
    cooc,
    guild_name
){
  
  if(guild_name == "Consumer"){
    
    cooc %>%
      dplyr::group_by(
        focal_node = consumer
      ) %>%
      dplyr::summarise(
        focal_species_cooccurrence_occupancy =
          dplyr::n_distinct(site),
        .groups = "drop"
      )
    
  } else {
    
    cooc %>%
      dplyr::group_by(
        focal_node = resource
      ) %>%
      dplyr::summarise(
        focal_species_cooccurrence_occupancy =
          dplyr::n_distinct(site),
        .groups = "drop"
      )
  }
}

## ------------------------------------------------------------
## Analyse one dataset
## ------------------------------------------------------------

run_one_dataset <- function(dataset){
  
  message("Spatial bundling analysis: ", dataset)
  
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
  
  ## An observed interaction necessarily implies co-occurrence.
  cooc <- dplyr::bind_rows(
    cooc,
    ints
  ) %>%
    dplyr::distinct(
      site,
      consumer,
      resource
    )
  
  cooccurrence_counts <- cooc %>%
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
  
  interaction_counts <- ints %>%
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
  
  pair_table_all <- cooccurrence_counts %>%
    dplyr::left_join(
      interaction_counts,
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
  
  if(any(pair_table_all$full_n < 1)){
    stop("Invalid full_n < 1 in ", dataset)
  }
  
  if(any(pair_table_all$full_K > pair_table_all$full_n)){
    stop("Invalid full_K > full_n in ", dataset)
  }
  
  ## Retain only realised regional interactions.
  pair_table <- pair_table_all %>%
    dplyr::filter(full_K > 0)
  
  if(nrow(pair_table) == 0){
    
    warning(
      "No realised regional interactions in ",
      dataset
    )
    
    return(
      list(
        species = tibble::tibble(),
        exclusions = tibble::tibble(),
        pair_checks = tibble::tibble()
      )
    )
  }
  
  ## All realised interactions must have an interaction-site list.
  if(any(
    vapply(
      pair_table$interaction_sites,
      is.null,
      logical(1)
    )
  )){
    
    stop(
      "Missing interaction-site list for a realised link in ",
      dataset
    )
  }
  
  pair_lengths_valid <- all(
    vapply(
      pair_table$cooccurrence_sites,
      length,
      integer(1)
    ) == pair_table$full_n
  ) &&
    all(
      vapply(
        pair_table$interaction_sites,
        length,
        integer(1)
      ) == pair_table$full_K
    )
  
  observed_interactions_inside_cooccurrence <- all(
    unlist(
      Map(
        function(
    interaction_sites,
    cooccurrence_sites
        ){
          
          all(
            interaction_sites %in%
              cooccurrence_sites
          )
        },
    pair_table$interaction_sites,
    pair_table$cooccurrence_sites
      ),
    use.names = FALSE
    )
  )
  
  guild_outputs <- lapply(
    guild_levels,
    function(guild_name){
      
      guild_links <- make_guild_link_table(
        pair_table,
        guild_name
      )
      
      focal_occupancy <- make_focal_occupancy_table(
        cooc,
        guild_name
      )
      
      focal_degree <- guild_links %>%
        dplyr::group_by(focal_node) %>%
        dplyr::summarise(
          initial_degree =
            dplyr::n_distinct(partner_node),
          .groups = "drop"
        ) %>%
        dplyr::left_join(
          focal_occupancy,
          by = "focal_node"
        )
      
      exclusions <- tibble::tibble(
        dataset = dataset,
        guild = guild_name,
        number_focal_species_with_realised_links =
          nrow(focal_degree),
        number_eligible_species =
          sum(focal_degree$initial_degree >= 2),
        number_excluded_degree_below_two =
          sum(focal_degree$initial_degree < 2)
      )
      
      eligible_nodes <- focal_degree %>%
        dplyr::filter(
          initial_degree >= 2
        )
      
      if(nrow(eligible_nodes) == 0){
        
        warning(
          "No eligible focal species with degree >= 2 in ",
          dataset,
          " / ",
          guild_name
        )
        
        return(
          list(
            species = tibble::tibble(),
            exclusions = exclusions
          )
        )
      }
      
      species_results <- lapply(
        seq_len(nrow(eligible_nodes)),
        function(i){
          
          focal_node_i <-
            eligible_nodes$focal_node[i]
          
          occupancy_i <-
            eligible_nodes$
            focal_species_cooccurrence_occupancy[i]
          
          focal_data_i <- guild_links %>%
            dplyr::filter(
              focal_node == focal_node_i
            )
          
          analyse_focal_species(
            focal_data = focal_data_i,
            dataset = dataset,
            guild = guild_name,
            focal_node = focal_node_i,
            focal_species_cooccurrence_occupancy =
              occupancy_i,
            n_null = n_null
          )
        }
      ) %>%
        dplyr::bind_rows()
      
      list(
        species = species_results,
        exclusions = exclusions
      )
    }
  )
  
  species_results <- dplyr::bind_rows(
    lapply(
      guild_outputs,
      `[[`,
      "species"
    )
  )
  
  exclusions <- dplyr::bind_rows(
    lapply(
      guild_outputs,
      `[[`,
      "exclusions"
    )
  )
  
  pair_checks <- tibble::tibble(
    dataset = dataset,
    all_realised_links_have_positive_K =
      all(pair_table$full_K > 0),
    all_realised_links_have_K_le_n =
      all(
        pair_table$full_K <=
          pair_table$full_n
      ),
    pair_site_list_lengths_match_counts =
      pair_lengths_valid,
    observed_interaction_sites_inside_cooccurrence =
      observed_interactions_inside_cooccurrence
  )
  
  list(
    species = species_results,
    exclusions = exclusions,
    pair_checks = pair_checks
  )
}

## ------------------------------------------------------------
## Run all datasets in parallel
## ------------------------------------------------------------

dataset_outputs <- if(use_parallel && n_workers > 1){
  
  future.apply::future_lapply(
    all_dataset_names,
    run_one_dataset,
    future.seed = TRUE,
    future.scheduling = 1
  )
  
} else {
  
  lapply(
    all_dataset_names,
    run_one_dataset
  )
}

## Return to sequential execution after dataset analysis.
future::plan(
  future::sequential
)

## ------------------------------------------------------------
## Combine dataset outputs
## ------------------------------------------------------------

species_internal <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "species"
  )
)

exclusions <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "exclusions"
  )
)

pair_checks <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "pair_checks"
  )
)

species_internal <- species_internal %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

## ------------------------------------------------------------
## Dataset-level summaries
## ------------------------------------------------------------

dataset_summary <- species_internal %>%
  dplyr::group_by(
    dataset,
    guild
  ) %>%
  dplyr::summarise(
    dataset_mean_bundling_delta =
      mean(
        bundling_delta,
        na.rm = TRUE
      ),
    dataset_median_bundling_delta =
      median(
        bundling_delta,
        na.rm = TRUE
      ),
    dataset_mean_bundling_SES =
      mean(
        bundling_SES,
        na.rm = TRUE
      ),
    proportion_species_positive_delta =
      mean(
        bundling_delta > 0,
        na.rm = TRUE
      ),
    number_eligible_species =
      dplyr::n_distinct(focal_node),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    exclusions %>%
      dplyr::select(
        dataset,
        guild,
        number_excluded_degree_below_two
      ),
    by = c(
      "dataset",
      "guild"
    )
  ) %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
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
        n_datasets = n,
        equal_dataset_mean_bundling_delta =
          ifelse(
            n == 1,
            mean(x),
            NA_real_
          ),
        permutation_p = NA_real_,
        permutation_type = NA_character_
      )
    )
  }
  
  observed_mean <- mean(x)
  
  if(n <= 20){
    
    sign_grid <- expand.grid(
      rep(
        list(c(-1, 1)),
        n
      )
    )
    
    sign_matrix <- as.matrix(
      sign_grid
    )
    
    permuted_means <- rowMeans(
      sweep(
        sign_matrix,
        MARGIN = 2,
        STATS = x,
        FUN = "*"
      )
    )
    
    permutation_p <- mean(
      abs(permuted_means) >=
        abs(observed_mean) - 1e-12
    )
    
    permutation_type <- paste0(
      "Exact sign-flip, ",
      2^n,
      " permutations"
    )
    
  } else {
    
    permuted_means <- replicate(
      n_large,
      mean(
        x *
          sample(
            c(-1, 1),
            size = n,
            replace = TRUE
          )
      )
    )
    
    permutation_p <- (
      1 +
        sum(
          abs(permuted_means) >=
            abs(observed_mean)
        )
    ) / (
      n_large + 1
    )
    
    permutation_type <- paste0(
      "Monte Carlo sign-flip, ",
      n_large,
      " permutations"
    )
  }
  
  tibble::tibble(
    n_datasets = n,
    equal_dataset_mean_bundling_delta =
      observed_mean,
    permutation_p =
      permutation_p,
    permutation_type =
      permutation_type
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
    bootstrap_q025 =
      stats::quantile(
        boot_means,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),
    bootstrap_q975 =
      stats::quantile(
        boot_means,
        0.975,
        na.rm = TRUE,
        names = FALSE
      )
  )
}

equal_dataset_tests <- dataset_summary %>%
  dplyr::group_by(guild) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      x <- .x$dataset_mean_bundling_delta
      
      dplyr::bind_cols(
        sign_flip_test(
          x,
          n_large =
            n_permutations_large_n
        ),
        bootstrap_dataset_mean(
          x,
          n_boot = n_boot
        )
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      equal_dataset_mean_bundling_delta > 0 ~
        "Spatial bundling",
      equal_dataset_mean_bundling_delta < 0 ~
        "Spatial separation",
      TRUE ~
        "No mean departure"
    ),
    significant_0.05 =
      !is.na(permutation_p) &
      permutation_p < 0.05,
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

## ------------------------------------------------------------
## Validation checks
## ------------------------------------------------------------

species_checks <- tibble::tibble(
  check = c(
    "Every analysed interaction has full_K > 0",
    "Every analysed interaction has full_K <= full_n",
    "Every eligible focal species has at least two realised partners",
    "Every observed and null Jaccard value lies between zero and one",
    "Every null allocation uses exactly full_K sites per interaction",
    "Every null interaction site belongs to its pair-specific co-occurrence set",
    "Consumers and resources are analysed separately",
    "Primary inference uses one dataset mean per dataset and guild",
    "Equal-dataset means give every dataset equal weight"
  ),
  pass = c(
    all(
      pair_checks$
        all_realised_links_have_positive_K
    ),
    all(
      pair_checks$
        all_realised_links_have_K_le_n
    ),
    all(
      species_internal$initial_degree >= 2
    ),
    all(
      species_internal$observed_jaccard_valid &
        species_internal$null_jaccards_valid
    ),
    all(
      species_internal$
        null_allocations_preserve_K
    ),
    all(
      species_internal$
        null_sites_inside_pair_cooccurrence
    ),
    all(
      guild_levels %in%
        as.character(
          unique(species_internal$guild)
        )
    ),
    dataset_summary %>%
      dplyr::count(
        dataset,
        guild
      ) %>%
      dplyr::summarise(
        valid = all(n == 1)
      ) %>%
      dplyr::pull(valid),
    all(
      vapply(
        guild_levels,
        function(guild_i){
          
          dataset_values <- dataset_summary %>%
            dplyr::filter(
              guild == guild_i
            ) %>%
            dplyr::pull(
              dataset_mean_bundling_delta
            )
          
          reported_mean <- equal_dataset_tests %>%
            dplyr::filter(
              guild == guild_i
            ) %>%
            dplyr::pull(
              equal_dataset_mean_bundling_delta
            )
          
          if(
            length(dataset_values) == 0 ||
            length(reported_mean) != 1
          ){
            return(FALSE)
          }
          
          abs(
            mean(
              dataset_values,
              na.rm = TRUE
            ) -
              reported_mean
          ) < 1e-12
        },
        logical(1)
      )
    )
  )
)

additional_checks <- tibble::tibble(
  check = c(
    "Pair site-list lengths match full_n and full_K",
    "Observed interaction sites are contained within pair-specific co-occurrence sites"
  ),
  pass = c(
    all(
      pair_checks$
        pair_site_list_lengths_match_counts
    ),
    all(
      pair_checks$
        observed_interaction_sites_inside_cooccurrence
    )
  )
)

validation_checks <- dplyr::bind_rows(
  species_checks,
  additional_checks
)

## ------------------------------------------------------------
## Remove internal validation columns from species output
## ------------------------------------------------------------

species_output <- species_internal %>%
  dplyr::select(
    dataset,
    guild,
    focal_node,
    initial_degree,
    focal_species_cooccurrence_occupancy,
    interaction_footprint_sites,
    total_local_interaction_records,
    number_partner_pairs,
    observed_mean_jaccard,
    null_mean_jaccard,
    null_sd_jaccard,
    null_q025_jaccard,
    null_q975_jaccard,
    bundling_delta,
    bundling_SES
  )

## ------------------------------------------------------------
## Main figure
## ------------------------------------------------------------

main_plot <- ggplot2::ggplot() +
  
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.45
  ) +
  
  ggplot2::geom_jitter(
    data = dataset_summary,
    ggplot2::aes(
      x = guild,
      y = dataset_mean_bundling_delta
    ),
    width = 0.08,
    height = 0,
    colour = "grey55",
    alpha = 0.8,
    size = 2
  ) +
  
  ggplot2::geom_errorbar(
    data = equal_dataset_tests,
    ggplot2::aes(
      x = guild,
      ymin = bootstrap_q025,
      ymax = bootstrap_q975
    ),
    width = 0.12,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  ggplot2::geom_point(
    data = equal_dataset_tests,
    ggplot2::aes(
      x = guild,
      y = equal_dataset_mean_bundling_delta
    ),
    colour = "black",
    size = 3.5
  ) +
  
  ggplot2::scale_x_discrete(
    limits = guild_levels,
    drop = FALSE
  ) +
  
  theme_pub(base_size = 11) +
  
  ggplot2::labs(
    x = NULL,
    y = paste0(
      "Observed minus expected mean overlap among ",
      "interaction-supporting sites"
    ),
    title =
      "Spatial bundling of species interaction portfolios",
    subtitle = paste0(
      "Positive values indicate that different interactions of the same ",
      "species occur together in the same sites more often than expected"
    )
  )

## ------------------------------------------------------------
## Observed-versus-null diagnostic figure
## ------------------------------------------------------------

diagnostic_plot <- ggplot2::ggplot(
  species_output,
  ggplot2::aes(
    x = null_mean_jaccard,
    y = observed_mean_jaccard
  )
) +
  
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.55
  ) +
  
  ggplot2::geom_point(
    colour = "grey45",
    alpha = 0.35,
    size = 1.5
  ) +
  
  ggplot2::facet_wrap(
    ~ guild,
    nrow = 1
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  
  ggplot2::coord_fixed(
    ratio = 1
  ) +
  
  theme_pub(base_size = 10) +
  
  ggplot2::labs(
    x = "Null mean Jaccard overlap",
    y = "Observed mean Jaccard overlap",
    title = "Observed versus expected spatial overlap among interactions"
  )

## ------------------------------------------------------------
## Save CSV outputs
## ------------------------------------------------------------

write.csv2(
  species_output,
  file.path(
    out_dir,
    "38_species_level_spatial_bundling.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_summary,
  file.path(
    out_dir,
    "38_dataset_level_spatial_bundling.csv"
  ),
  row.names = FALSE
)

write.csv2(
  equal_dataset_tests,
  file.path(
    out_dir,
    "38_equal_dataset_spatial_bundling_tests.csv"
  ),
  row.names = FALSE
)

write.csv2(
  validation_checks,
  file.path(
    out_dir,
    "38_spatial_bundling_validation_checks.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Save PNG outputs
## ------------------------------------------------------------

save_png(
  main_plot,
  "38_spatial_bundling_main.png",
  width = 8,
  height = 6
)

save_png(
  diagnostic_plot,
  "38_spatial_bundling_observed_vs_null.png",
  width = 9,
  height = 5
)

## ------------------------------------------------------------
## Console summary
## ------------------------------------------------------------

message("\nEqual-dataset result: Consumers")

print(
  equal_dataset_tests %>%
    dplyr::filter(
      guild == "Consumer"
    ),
  n = Inf
)

message("\nEqual-dataset result: Resources")

print(
  equal_dataset_tests %>%
    dplyr::filter(
      guild == "Resource"
    ),
  n = Inf
)

analysis_counts <- species_output %>%
  dplyr::group_by(guild) %>%
  dplyr::summarise(
    number_datasets =
      dplyr::n_distinct(dataset),
    number_eligible_species =
      dplyr::n_distinct(
        paste(
          dataset,
          focal_node,
          sep = "___"
        )
      ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    exclusions %>%
      dplyr::group_by(guild) %>%
      dplyr::summarise(
        number_excluded_degree_below_two =
          sum(
            number_excluded_degree_below_two,
            na.rm = TRUE
          ),
        .groups = "drop"
      ),
    by = "guild"
  )

message("\nDatasets and focal species analysed:")

print(
  analysis_counts,
  n = Inf
)

message("\nValidation checks:")

print(
  validation_checks,
  n = Inf
)

message("\nInterpretation:")

message(
  "Positive bundling delta means that the different interactions belonging\n",
  "to the same species are supported by the same local sites more often\n",
  "than expected after preserving the observed support and pair-specific\n",
  "co-occurrence opportunities of every interaction.\n\n",
  "Negative bundling delta means that a species' interactions are more\n",
  "spatially complementary or separated among sites than expected."
)

message(
  "\nSaved Script 38 outputs in: ",
  out_dir
)