## ------------------------------------------------------------
## Script:
## All/scripts/36FINAL_Figure2.R
##
## Purpose:
## Create the final combined Figure 2:
##
##   A. Species, regional links, and local links retained under
##      random site removal.
##
##   B. Mean local interaction support among species and datasets.
##
##   C. State of original regional interactions after site removal.
##
## Layout:
##
##       ┌──────────────────────────────┐
##       │              A               │
##       ├──────────────┬───────────────┤
##       │      B       │       C       │
##       └──────────────┴───────────────┘
##
## Panel A is larger and spans the full top row.
##
## No panel titles or subtitles are included. Only panel letters,
## axes, and legends are shown.
##
## Output:
## All/outputs/36FINAL_Figure2/
## 36FINAL_Figure2.png
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble",
  "purrr",
  "patchwork",
  "scales"
)

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ------------------------------------------------------------
## Controls
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

n_site_reps <- 500L

out_dir <- "All/outputs/36FINAL_Figure2"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  out_dir,
  "36FINAL_Figure2.png"
)

## ------------------------------------------------------------
## Colours and factor ordering
## ------------------------------------------------------------

retention_colours <- c(
  "Species" = "#6A3D9A",
  "Regional links" = "#D55E00",
  "Local links" = "#0072B2"
)

state_levels <- c(
  "Species no longer recorded together",
  "Species co-occurring, not interacting",
  "Interaction still observed"
)

state_colours <- c(
  "Species no longer recorded together" = "#999999",
  "Species co-occurring, not interacting" = "#CC79A7",
  "Interaction still observed" = "#009E73"
)

state_legend_labels <- c(
  "Species no longer recorded together",
  "Species co-occurring, not interacting",
  "Interaction still observed"
)

guild_levels <- c(
  "Consumer",
  "Resource"
)

guild_colours <- c(
  "Consumer" = "#4E79A7",
  "Resource" = "#59A14F"
)

replace_na_int <- function(
    x,
    value = 0L
){
  
  x[is.na(x)] <- value
  x
}

safe_quantile <- function(
    x,
    probability
){
  
  x <- x[
    is.finite(x)
  ]
  
  if(length(x) == 0){
    return(NA_real_)
  }
  
  stats::quantile(
    x,
    probability,
    na.rm = TRUE,
    names = FALSE
  )
}

compact_theme <- function(base_size = 10){
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      
      legend.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      
      legend.box.margin = ggplot2::margin(
        t = -3,
        r = 0,
        b = 0,
        l = 0
      ),
      
      legend.key.height = grid::unit(
        0.35,
        "cm"
      ),
      
      legend.key.width = grid::unit(
        0.65,
        "cm"
      ),
      
      legend.spacing.x = grid::unit(
        0.10,
        "cm"
      ),
      
      strip.background =
        ggplot2::element_blank(),
      
      strip.text =
        ggplot2::element_text(
          face = "bold"
        ),
      
      plot.margin = ggplot2::margin(
        t = 2,
        r = 3,
        b = 2,
        l = 3
      )
    )
}

## ------------------------------------------------------------
## Random retained-site sets
## ------------------------------------------------------------

make_retained_site_sets <- function(
    sites,
    removal_levels,
    n_site_reps
){
  
  sites <- sort(
    unique(
      as.character(sites)
    )
  )
  
  number_sites <- length(sites)
  
  if(number_sites < 1){
    stop("Cannot generate site-removal sets for zero sites.")
  }
  
  retained_rows <- vector(
    "list",
    length(removal_levels)
  )
  
  for(level_index in seq_along(
    removal_levels
  )){
    
    removal_fraction <-
      removal_levels[level_index]
    
    number_removed <- round(
      removal_fraction *
        number_sites
    )
    
    number_removed <- min(
      number_removed,
      number_sites
    )
    
    number_replicates <- if(
      removal_fraction == 0
    ){
      1L
    } else {
      n_site_reps
    }
    
    level_rows <- vector(
      "list",
      number_replicates
    )
    
    for(rep_index in seq_len(
      number_replicates
    )){
      
      removed_sites <- if(
        number_removed == 0
      ){
        character(0)
      } else {
        sample(
          sites,
          size = number_removed,
          replace = FALSE
        )
      }
      
      level_rows[[rep_index]] <-
        tibble::tibble(
          removal_fraction =
            removal_fraction,
          replicate =
            rep_index,
          site = setdiff(
            sites,
            removed_sites
          )
        )
    }
    
    retained_rows[[level_index]] <-
      dplyr::bind_rows(
        level_rows
      )
  }
  
  dplyr::bind_rows(
    retained_rows
  )
}

## ------------------------------------------------------------
## Process one dataset
## ------------------------------------------------------------

run_one_dataset <- function(dataset){
  
  message(
    "Final Figure 2: ",
    dataset
  )
  
  st <- get_dataset_site_tables(
    dataset
  )
  
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
  
  sites <- sort(
    unique(cooc$site)
  )
  
  if(length(sites) == 0){
    stop(
      "No sites found in ",
      dataset
    )
  }
  
  retained_sites <- make_retained_site_sets(
    sites = sites,
    removal_levels = removal_levels,
    n_site_reps = n_site_reps
  )
  
  subset_index <- retained_sites %>%
    dplyr::distinct(
      removal_fraction,
      replicate
    )
  
  ## ----------------------------------------------------------
  ## Full species pool
  ## ----------------------------------------------------------
  
  full_species <- union(
    unique(cooc$consumer),
    unique(cooc$resource)
  )
  
  full_species_richness <-
    length(full_species)
  
  if(full_species_richness == 0){
    stop(
      "No species found in ",
      dataset
    )
  }
  
  ## ----------------------------------------------------------
  ## Full regional interaction links
  ## ----------------------------------------------------------
  
  full_pairs <- cooc %>%
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
  
  if(any(
    full_pairs$full_K >
    full_pairs$full_n
  )){
    stop(
      "Invalid full_K > full_n in ",
      dataset
    )
  }
  
  full_links <- full_pairs %>%
    dplyr::filter(
      full_K > 0
    )
  
  full_regional_link_richness <-
    nrow(full_links)
  
  full_local_link_support <- sum(
    full_links$full_K
  )
  
  if(
    full_regional_link_richness == 0 ||
    full_local_link_support == 0
  ){
    
    warning(
      "No realised regional links in ",
      dataset
    )
    
    return(
      list(
        retention = tibble::tibble(),
        states = tibble::tibble(),
        species_support = tibble::tibble()
      )
    )
  }
  
  ## ----------------------------------------------------------
  ## Retained co-occurrence and interaction support
  ## ----------------------------------------------------------
  
  full_link_grid <- tidyr::crossing(
    subset_index,
    full_links %>%
      dplyr::select(
        consumer,
        resource,
        full_n,
        full_K
      )
  )
  
  retained_cooccurrence <- retained_sites %>%
    dplyr::left_join(
      cooc %>%
        dplyr::semi_join(
          full_links,
          by = c(
            "consumer",
            "resource"
          )
        ) %>%
        dplyr::mutate(
          cooccurrence_here = 1L
        ),
      by = "site",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      !is.na(consumer),
      !is.na(resource)
    ) %>%
    dplyr::group_by(
      removal_fraction,
      replicate,
      consumer,
      resource
    ) %>%
    dplyr::summarise(
      retained_n = sum(
        cooccurrence_here,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  retained_interactions <- retained_sites %>%
    dplyr::left_join(
      ints %>%
        dplyr::semi_join(
          full_links,
          by = c(
            "consumer",
            "resource"
          )
        ) %>%
        dplyr::mutate(
          interaction_here = 1L
        ),
      by = "site",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      !is.na(consumer),
      !is.na(resource)
    ) %>%
    dplyr::group_by(
      removal_fraction,
      replicate,
      consumer,
      resource
    ) %>%
    dplyr::summarise(
      retained_K = sum(
        interaction_here,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  link_replicates <- full_link_grid %>%
    dplyr::left_join(
      retained_cooccurrence,
      by = c(
        "removal_fraction",
        "replicate",
        "consumer",
        "resource"
      )
    ) %>%
    dplyr::left_join(
      retained_interactions,
      by = c(
        "removal_fraction",
        "replicate",
        "consumer",
        "resource"
      )
    ) %>%
    dplyr::mutate(
      retained_n = replace_na_int(
        retained_n,
        0L
      ),
      retained_K = replace_na_int(
        retained_K,
        0L
      )
    )
  
  ## ----------------------------------------------------------
  ## Species still represented in retained sites
  ## ----------------------------------------------------------
  
  retained_species <- retained_sites %>%
    dplyr::left_join(
      cooc,
      by = "site",
      relationship = "many-to-many"
    ) %>%
    dplyr::select(
      removal_fraction,
      replicate,
      consumer,
      resource
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        consumer,
        resource
      ),
      names_to = "guild",
      values_to = "species"
    ) %>%
    dplyr::filter(
      !is.na(species)
    ) %>%
    dplyr::distinct(
      removal_fraction,
      replicate,
      species
    ) %>%
    dplyr::count(
      removal_fraction,
      replicate,
      name =
        "retained_species_richness"
    )
  
  ## ----------------------------------------------------------
  ## Panel A data
  ## ----------------------------------------------------------
  
  retention <- link_replicates %>%
    dplyr::group_by(
      removal_fraction,
      replicate
    ) %>%
    dplyr::summarise(
      dataset = dataset,
      
      local_links_retained =
        sum(
          retained_K,
          na.rm = TRUE
        ) /
        full_local_link_support,
      
      regional_links_retained =
        sum(
          retained_K > 0,
          na.rm = TRUE
        ) /
        full_regional_link_richness,
      
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      retained_species,
      by = c(
        "removal_fraction",
        "replicate"
      )
    ) %>%
    dplyr::mutate(
      retained_species_richness =
        replace_na_int(
          retained_species_richness,
          0L
        ),
      
      species_retained =
        retained_species_richness /
        full_species_richness
    )
  
  ## ----------------------------------------------------------
  ## Panel C data
  ## ----------------------------------------------------------
  
  states <- link_replicates %>%
    dplyr::mutate(
      state = dplyr::case_when(
        retained_K > 0 ~
          "Interaction still observed",
        
        retained_n > 0 ~
          "Species co-occurring, not interacting",
        
        TRUE ~
          "Species no longer recorded together"
      ),
      state = factor(
        state,
        levels = state_levels
      )
    ) %>%
    dplyr::count(
      removal_fraction,
      replicate,
      state,
      name = "number_links"
    ) %>%
    tidyr::complete(
      removal_fraction,
      replicate,
      state = factor(
        state_levels,
        levels = state_levels
      ),
      fill = list(
        number_links = 0L
      )
    ) %>%
    dplyr::group_by(
      removal_fraction,
      replicate
    ) %>%
    dplyr::mutate(
      dataset = dataset,
      fraction =
        number_links /
        sum(number_links)
    ) %>%
    dplyr::ungroup()
  
  ## ----------------------------------------------------------
  ## Panel B data
  ## Mean K per realised partner for each focal species
  ## ----------------------------------------------------------
  
  focal_links <- dplyr::bind_rows(
    
    full_links %>%
      dplyr::transmute(
        dataset = dataset,
        guild = "Consumer",
        focal_node = consumer,
        partner_node = resource,
        full_K
      ),
    
    full_links %>%
      dplyr::transmute(
        dataset = dataset,
        guild = "Resource",
        focal_node = resource,
        partner_node = consumer,
        full_K
      )
  )
  
  species_support <- focal_links %>%
    dplyr::group_by(
      dataset,
      guild,
      focal_node
    ) %>%
    dplyr::summarise(
      initial_degree =
        dplyr::n_distinct(
          partner_node
        ),
      
      mean_K_per_link =
        mean(
          full_K,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  list(
    retention = retention,
    states = states,
    species_support = species_support
  )
}

## ------------------------------------------------------------
## Run all datasets
## ------------------------------------------------------------

dataset_outputs <- lapply(
  all_dataset_names,
  run_one_dataset
)

retention <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "retention"
  )
)

states <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "states"
  )
)

species_support <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    `[[`,
    "species_support"
  )
)

if(nrow(retention) == 0){
  stop(
    "No retention data were produced."
  )
}

if(nrow(states) == 0){
  stop(
    "No interaction-state data were produced."
  )
}

if(nrow(species_support) == 0){
  stop(
    "No species-level local-support data were produced."
  )
}

## ------------------------------------------------------------
## Panel A summaries
## Equal dataset weight
## ------------------------------------------------------------

retention_long <- retention %>%
  dplyr::select(
    dataset,
    removal_fraction,
    replicate,
    species_retained,
    regional_links_retained,
    local_links_retained
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      species_retained,
      regional_links_retained,
      local_links_retained
    ),
    names_to = "layer",
    values_to = "fraction"
  ) %>%
  dplyr::mutate(
    layer = dplyr::recode(
      layer,
      species_retained =
        "Species",
      regional_links_retained =
        "Regional links",
      local_links_retained =
        "Local links"
    ),
    
    layer = factor(
      layer,
      levels = c(
        "Species",
        "Regional links",
        "Local links"
      )
    )
  )

dataset_retention_trajectories <-
  retention_long %>%
  dplyr::group_by(
    dataset,
    layer,
    removal_fraction
  ) %>%
  dplyr::summarise(
    fraction = stats::median(
      fraction,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

equal_dataset_retention <-
  dataset_retention_trajectories %>%
  dplyr::group_by(
    layer,
    removal_fraction
  ) %>%
  dplyr::summarise(
    mean_fraction =
      mean(
        fraction,
        na.rm = TRUE
      ),
    
    q25 =
      safe_quantile(
        fraction,
        0.25
      ),
    
    q75 =
      safe_quantile(
        fraction,
        0.75
      ),
    
    .groups = "drop"
  )

## ------------------------------------------------------------
## Panel C summaries
## Equal dataset weight
## ------------------------------------------------------------

state_mean <- states %>%
  dplyr::mutate(
    state = factor(
      as.character(state),
      levels = state_levels
    )
  ) %>%
  dplyr::group_by(
    dataset,
    removal_fraction,
    state
  ) %>%
  dplyr::summarise(
    fraction =
      mean(
        fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    removal_fraction,
    state
  ) %>%
  dplyr::summarise(
    fraction =
      mean(
        fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state = factor(
      as.character(state),
      levels = state_levels
    )
  )

## ------------------------------------------------------------
## Panel B preparation
## ------------------------------------------------------------
## Desired order written approximately.
desired_dataset_order <- c(
  "Garraf PP ",
  "Garraf PP2 ",
  "Montseny ",
  "Gottin PP ",
  "Nahuel ",
  "Garraf HP ",
  "Quercus ",
  "Olot ",
  "Gottin HP ",
  "Galpar "
)

## Standardise spelling, spaces, punctuation, and capitalisation.
normalise_dataset_name <- function(x){
  
  x <- tolower(x)
  x <- iconv(
    x,
    from = "",
    to = "ASCII//TRANSLIT"
  )
  x <- gsub(
    "[^a-z0-9]",
    "",
    x
  )
  
  x
}

actual_dataset_names <- unique(
  as.character(
    species_support$dataset
  )
)

normalised_actual <- normalise_dataset_name(
  actual_dataset_names
)

match_dataset_name <- function(requested_name){
  
  requested_normalised <- normalise_dataset_name(
    requested_name
  )
  
  ## First try an exact match after normalisation.
  exact_match <- which(
    normalised_actual ==
      requested_normalised
  )
  
  if(length(exact_match) == 1){
    return(
      actual_dataset_names[
        exact_match
      ]
    )
  }
  
  ## Otherwise use approximate spelling matching.
  approximate_match <- agrep(
    pattern = requested_normalised,
    x = normalised_actual,
    max.distance = 0.25,
    value = FALSE
  )
  
  if(length(approximate_match) == 1){
    return(
      actual_dataset_names[
        approximate_match
      ]
    )
  }
  
  if(length(approximate_match) > 1){
    
    distances <- utils::adist(
      requested_normalised,
      normalised_actual[
        approximate_match
      ]
    )
    
    best_match <- approximate_match[
      which.min(distances)
    ]
    
    return(
      actual_dataset_names[
        best_match
      ]
    )
  }
  
  warning(
    "No dataset matched the requested name: ",
    requested_name
  )
  
  NA_character_
}

dataset_order <- vapply(
  desired_dataset_order,
  match_dataset_name,
  character(1)
)

dataset_order <- unique(
  dataset_order[
    !is.na(dataset_order)
  ]
)

## Add any datasets not included in the requested list at the end.
unmatched_datasets <- setdiff(
  actual_dataset_names,
  dataset_order
)

dataset_order <- c(
  dataset_order,
  unmatched_datasets
)

message(
  "\nDataset order used in panel B:"
)

print(dataset_order)

species_support <- species_support %>%
  dplyr::mutate(
    dataset = factor(
      as.character(dataset),
      levels = dataset_order
    ),
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

maximum_species_support <- max(
  species_support$mean_K_per_link,
  na.rm = TRUE
)

upper_support_break <- max(
  5,
  ceiling(
    maximum_species_support / 5
  ) * 5
)

support_axis_breaks <- unique(
  c(
    1:5,
    seq(
      from = 10,
      to = upper_support_break,
      by = 5
    )
  )
)

support_axis_breaks <- support_axis_breaks[
  support_axis_breaks <=
    upper_support_break
]

## ------------------------------------------------------------
## Panel A
## ------------------------------------------------------------
panel_A <- ggplot2::ggplot() +
  
  # ggplot2::geom_ribbon(
  #   data = equal_dataset_retention,
  #   ggplot2::aes(
  #     x = removal_fraction,
  #     ymin = q25,
  #     ymax = q75,
  #   ),
  #   alpha = 0.15,
  #   colour = NA
  # ) +
  
  ggplot2::geom_line(
    data =
      dataset_retention_trajectories,
    ggplot2::aes(
      x = removal_fraction,
      y = fraction,
      colour = layer,
      group = interaction(
        dataset,
        layer
      )
    ),
    alpha = 0.25,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_line(
    data = equal_dataset_retention,
    ggplot2::aes(
      x = removal_fraction,
      y = mean_fraction,
      colour = layer
    ),
    linewidth = 1.35
  ) +
  
  ggplot2::geom_point(
    data = equal_dataset_retention,
    ggplot2::aes(
      x = removal_fraction,
      y = mean_fraction,
      colour = layer
    ),
    size = 2
  ) +
  
  ggplot2::scale_colour_manual(
    values = retention_colours,
    breaks = names(
      retention_colours
    ),
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = removal_levels,
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  ggplot2::scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  
  compact_theme(
    base_size = 11
  ) +
  
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5
      ),
    legend.text = ggplot2::element_text(
      size = 10
    ),
  ) +
  
  ggplot2::labs(
    x = "Sites removed",
    y = "Fraction remaining",
    colour = NULL,
    fill = NULL
  )

## ------------------------------------------------------------
## Panel B
## ------------------------------------------------------------
panel_B <- ggplot2::ggplot(
  species_support,
  ggplot2::aes(
    x = dataset,
    y = mean_K_per_link,
    colour = guild
  )
) +
  
  ggplot2::geom_boxplot(
    ggplot2::aes(
      group = interaction(
        dataset,
        guild
      )
    ),
    position =
      ggplot2::position_dodge(
        width = 0.7
      ),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.12,
    linewidth = 0.45
  ) +
  
  ggplot2::geom_point(
    position =
      ggplot2::position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.7
      ),
    alpha = 0.32,
    size = 1
  ) +
  
  ggplot2::scale_colour_manual(
    values = guild_colours,
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    trans = "log1p",
    breaks = support_axis_breaks,
    labels = support_axis_breaks,
    expand = ggplot2::expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +
  
  compact_theme(
    base_size = 9
  ) +
  
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      size = 8.5,
      margin = ggplot2::margin(b = 0)
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 9
    ),
    
    # Keep the vertical axis title close to the tick labels/axis.
    axis.title.y = ggplot2::element_text(
      vjust = -2,
      size = 10
    ),
    
    axis.title.x = ggplot2::element_text(
      size = 10
    ),
    
    legend.text = ggplot2::element_text(
      size = 9
    ),
    
    legend.key.height = grid::unit(
      0.45,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.80,
      "cm"
    ),
    
    legend.box.margin = ggplot2::margin(
      t = -10,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.margin = ggplot2::margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Mean sites supporting each realised partner",
    colour = NULL
  )

## ------------------------------------------------------------
## Panel C
## ------------------------------------------------------------
panel_C <- ggplot2::ggplot(
  state_mean,
  ggplot2::aes(
    x = removal_fraction,
    y = fraction,
    fill = state
  )
) +
  
  ggplot2::geom_area(
    alpha = 0.95,
    colour = "white",
    linewidth = 0.2
  ) +
  
  ggplot2::scale_fill_manual(
    values = state_colours,
    breaks = state_levels,
    labels = state_legend_labels,
    drop = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = removal_levels,
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  ggplot2::scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  ggplot2::coord_cartesian(
    ylim = c(
      0,
      1
    )
  ) +
  
  compact_theme(
    base_size = 9
  ) +
  
  ggplot2::theme(
    axis.text = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 9,
      margin = ggplot2::margin(t = 1)
    ),
    
    axis.title.x = ggplot2::element_text(
      size = 10,
      margin = ggplot2::margin(
        t = 6,
        b = 3
      )
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 10
    ),
    
    legend.box.margin = ggplot2::margin(
      t = -5,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.margin = ggplot2::margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    
    plot.margin = ggplot2::margin(
      t = 2,
      r = 3,
      b = 2,
      l = 3
    )
  ) +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      ncol = 3,
      byrow = TRUE
    )
  ) +
  
  ggplot2::labs(
    x = "Sites removed",
    y = "Share of original regional interactions",
    fill = NULL
  )

## ------------------------------------------------------------
## Combined figure
## ------------------------------------------------------------
# Release panel C from panel B's unusually deep bottom axis area.
# Panel B's rotated dataset names therefore no longer push panel C's
# x-axis title and legend downward.
bottom_row <- panel_B +
  patchwork::wrap_elements(
    full = panel_C
  ) +
  patchwork::plot_layout(
    widths = c(0.45, 0.55)
  )

final_figure <- panel_A /
  patchwork::plot_spacer() /
  bottom_row +
  
  patchwork::plot_layout(
    heights = c(
      1.42,
      0.07,
      1.15
    )
  ) +
  
  patchwork::plot_annotation(
    # tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      
      plot.tag.position = c(
        0,
        1
      ),
      
      plot.margin = ggplot2::margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1
      )
    )
  ) &
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 2,
      r = 2,
      b = 2,
      l = 2
    )
  )

## ------------------------------------------------------------
## Save
## ------------------------------------------------------------
ggplot2::ggsave(
  filename = output_file,
  plot = final_figure,
  width = 12,
  height = 9.5,
  dpi = 320,
  bg = "white",
  limitsize = FALSE
)

## Optional underlying tables for reproducibility.
write.csv2(
  retention,
  file.path(
    out_dir,
    "36FINAL_retention_data.csv"
  ),
  row.names = FALSE
)

write.csv2(
  states,
  file.path(
    out_dir,
    "36FINAL_interaction_state_data.csv"
  ),
  row.names = FALSE
)

write.csv2(
  species_support,
  file.path(
    out_dir,
    "36FINAL_species_local_support_data.csv"
  ),
  row.names = FALSE
)

message(
  "\nSaved final Figure 2:\n",
  output_file
)