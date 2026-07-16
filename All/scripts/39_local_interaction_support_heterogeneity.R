## ------------------------------------------------------------
## Script:
## All/scripts/39_local_interaction_support_heterogeneity.R
##
## Purpose:
## Describe how heterogeneous local interaction support is:
##
##   1. among realised regional links;
##   2. among species;
##   3. among datasets.
##
## Definitions:
##
##   full_K =
##     number of sites supporting a realised regional interaction
##
##   K_standardised =
##     full_K / number of sites in the dataset
##
## Species-level quantities are calculated separately for Consumers
## and Resources and only across initially realised partners.
##
## This is a descriptive script. It does not perform null modelling,
## site removal, or significance testing.
##
## Outputs:
##   - five separate PNG figures;
##   - one combined 2 × 2 PNG;
##   - three descriptive CSV files.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "scales",
  "gridExtra",
  "grid"
)

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ------------------------------------------------------------
## Output settings
## ------------------------------------------------------------

out_dir <- "All/outputs/39_local_interaction_support_heterogeneity"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_png <- function(
    plot,
    filename,
    width = 10,
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

save_png_grid <- function(
    grob,
    filename,
    width = 13,
    height = 10
){
  
  grDevices::png(
    filename = file.path(
      out_dir,
      filename
    ),
    width = width,
    height = height,
    units = "in",
    res = 320,
    bg = "white"
  )
  
  grid::grid.newpage()
  grid::grid.draw(grob)
  
  grDevices::dev.off()
}

theme_pub <- function(base_size = 10){
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
    )
}

guild_levels <- c(
  "Consumer",
  "Resource"
)

guild_cols <- c(
  "Consumer" = "#4E79A7",
  "Resource" = "#59A14F"
)

## ------------------------------------------------------------
## Small helpers
## ------------------------------------------------------------

replace_na_int <- function(
    x,
    value = 0L
){
  
  x[is.na(x)] <- value
  x
}

safe_cv <- function(x){
  
  x <- x[
    is.finite(x)
  ]
  
  if(
    length(x) < 2 ||
    mean(x) == 0
  ){
    return(NA_real_)
  }
  
  stats::sd(x) / mean(x)
}

safe_gini <- function(x){
  
  x <- x[
    is.finite(x)
  ]
  
  x <- x[
    x >= 0
  ]
  
  n <- length(x)
  
  if(n == 0){
    return(NA_real_)
  }
  
  if(all(x == 0)){
    return(0)
  }
  
  x <- sort(x)
  
  (
    2 * sum(
      seq_len(n) * x
    ) /
      (
        n * sum(x)
      )
  ) -
    (
      n + 1
    ) / n
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

## ------------------------------------------------------------
## Build link-level data for one dataset
## ------------------------------------------------------------

run_one_dataset <- function(dataset){
  
  message(
    "Local-support heterogeneity: ",
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
  
  number_dataset_sites <- dplyr::n_distinct(
    cooc$site
  )
  
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
  
  if(any(pair_table$full_K > pair_table$full_n)){
    stop(
      "Invalid full_K > full_n in ",
      dataset
    )
  }
  
  realised_links <- pair_table %>%
    dplyr::filter(
      full_K > 0
    ) %>%
    dplyr::mutate(
      dataset = dataset,
      number_dataset_sites =
        number_dataset_sites,
      K_standardised =
        full_K /
        number_dataset_sites,
      q =
        full_K /
        full_n,
      link_id = paste(
        consumer,
        resource,
        sep = "___"
      )
    )
  
  if(nrow(realised_links) == 0){
    
    warning(
      "No realised regional links in ",
      dataset
    )
    
    return(
      list(
        links = tibble::tibble(),
        species = tibble::tibble(),
        dataset = tibble::tibble()
      )
    )
  }
  
  ## Convert every regional link into Consumer- and Resource-centred
  ## representations for species-level summaries.
  focal_links <- dplyr::bind_rows(
    
    realised_links %>%
      dplyr::transmute(
        dataset,
        guild = "Consumer",
        focal_node = consumer,
        partner_node = resource,
        link_id,
        number_dataset_sites,
        full_n,
        full_K,
        K_standardised,
        q
      ),
    
    realised_links %>%
      dplyr::transmute(
        dataset,
        guild = "Resource",
        focal_node = resource,
        partner_node = consumer,
        link_id,
        number_dataset_sites,
        full_n,
        full_K,
        K_standardised,
        q
      )
  )
  
  species_summary <- focal_links %>%
    dplyr::group_by(
      dataset,
      guild,
      focal_node
    ) %>%
    dplyr::summarise(
      initial_degree =
        dplyr::n_distinct(partner_node),
      
      mean_K_per_link =
        mean(
          full_K,
          na.rm = TRUE
        ),
      
      median_K_per_link =
        stats::median(
          full_K,
          na.rm = TRUE
        ),
      
      minimum_K =
        min(
          full_K,
          na.rm = TRUE
        ),
      
      maximum_K =
        max(
          full_K,
          na.rm = TRUE
        ),
      
      range_K =
        maximum_K -
        minimum_K,
      
      sd_K =
        if(dplyr::n() >= 2){
          stats::sd(full_K)
        } else {
          NA_real_
        },
      
      cv_K =
        safe_cv(full_K),
      
      gini_K =
        safe_gini(full_K),
      
      mean_K_standardised =
        mean(
          K_standardised,
          na.rm = TRUE
        ),
      
      mean_q =
        mean(
          q,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      guild = factor(
        guild,
        levels = guild_levels
      )
    )
  
  dataset_summary <- realised_links %>%
    dplyr::summarise(
      dataset =
        dplyr::first(dataset),
      
      number_dataset_sites =
        dplyr::first(number_dataset_sites),
      
      number_realised_links =
        dplyr::n(),
      
      mean_K =
        mean(
          full_K,
          na.rm = TRUE
        ),
      
      median_K =
        stats::median(
          full_K,
          na.rm = TRUE
        ),
      
      q25_K =
        safe_quantile(
          full_K,
          0.25
        ),
      
      q75_K =
        safe_quantile(
          full_K,
          0.75
        ),
      
      minimum_K =
        min(
          full_K,
          na.rm = TRUE
        ),
      
      maximum_K =
        max(
          full_K,
          na.rm = TRUE
        ),
      
      cv_K =
        safe_cv(full_K),
      
      gini_K =
        safe_gini(full_K),
      
      mean_K_standardised =
        mean(
          K_standardised,
          na.rm = TRUE
        ),
      
      median_K_standardised =
        stats::median(
          K_standardised,
          na.rm = TRUE
        ),
      
      proportion_single_site_links =
        mean(
          full_K == 1
        ),
      
      .groups = "drop"
    )
  
  return(
    list(
      links = realised_links,
      species = species_summary,
      dataset = dataset_summary
    )
  )
}

## ------------------------------------------------------------
## Run all datasets
## ------------------------------------------------------------
dataset_outputs <- lapply(
  all_dataset_names,
  run_one_dataset
)

names(dataset_outputs) <- all_dataset_names

## Check that every dataset returned the expected object
invalid_outputs <- names(dataset_outputs)[
  !vapply(
    dataset_outputs,
    function(x){
      is.list(x) &&
        all(
          c(
            "links",
            "species",
            "dataset"
          ) %in% names(x)
        ) &&
        is.data.frame(x$links) &&
        is.data.frame(x$species) &&
        is.data.frame(x$dataset)
    },
    logical(1)
  )
]

if(length(invalid_outputs) > 0){
  stop(
    paste0(
      "run_one_dataset() returned an invalid object for: ",
      paste(
        invalid_outputs,
        collapse = ", "
      )
    )
  )
}

link_data <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    function(x){
      x$links
    }
  )
)

species_data <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    function(x){
      x$species
    }
  )
)

dataset_data <- dplyr::bind_rows(
  lapply(
    dataset_outputs,
    function(x){
      x$dataset
    }
  )
)

dataset_order <- dataset_data %>%
  dplyr::arrange(
    mean_K
  ) %>%
  dplyr::distinct(
    dataset
  ) %>%
  dplyr::pull(
    dataset
  )

link_data <- link_data %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = dataset_order
    )
  )

species_data <- species_data %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = dataset_order
    ),
    guild = factor(
      guild,
      levels = guild_levels
    )
  )

dataset_data <- dataset_data %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = dataset_order
    )
  )

## ------------------------------------------------------------
## Figure 1: heterogeneity among links
## ------------------------------------------------------------

plot_links <- ggplot2::ggplot(
  link_data,
  ggplot2::aes(
    x = dataset,
    y = full_K
  )
) +
  ggplot2::geom_violin(
    fill = "grey88",
    colour = "grey55",
    linewidth = 0.35,
    scale = "width",
    trim = TRUE
  ) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.45
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.15,
    size = 0.7,
    colour = "grey30"
  ) +
  ggplot2::scale_y_continuous(
    trans = "log1p",
    breaks = c(
      1,
      2,
      5,
      10,
      20,
      50,
      100,
      200
    )
  ) +
  theme_pub() +
  ggplot2::labs(
    x = NULL,
    y = "Sites supporting each regional interaction, K",
    title = "Local interaction support varies strongly among regional links",
    subtitle = "Each point is one realised consumer–resource interaction"
  )

## ------------------------------------------------------------
## Figure 2: mean support among species
## ------------------------------------------------------------

plot_species_mean <- ggplot2::ggplot(
  species_data,
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
    position = ggplot2::position_dodge(
      width = 0.7
    ),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.12,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.7
    ),
    alpha = 0.32,
    size = 1
  ) +
  ggplot2::scale_colour_manual(
    values = guild_cols,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    trans = "log1p"
  ) +
  theme_pub() +
  ggplot2::labs(
    x = NULL,
    y = "Mean sites supporting each realised partner",
    colour = NULL,
    title = "Species differ in the typical support of their interactions",
    subtitle = "Each point is one focal species, calculated across its realised regional partners"
  )

## ------------------------------------------------------------
## Figure 3: within-species heterogeneity among links
## ------------------------------------------------------------

species_heterogeneity <- species_data %>%
  dplyr::filter(
    initial_degree >= 2,
    is.finite(cv_K)
  )

plot_species_cv <- ggplot2::ggplot(
  species_heterogeneity,
  ggplot2::aes(
    x = dataset,
    y = cv_K,
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
    position = ggplot2::position_dodge(
      width = 0.7
    ),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.12,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.7
    ),
    alpha = 0.32,
    size = 1
  ) +
  ggplot2::scale_colour_manual(
    values = guild_cols,
    drop = FALSE
  ) +
  theme_pub() +
  ggplot2::labs(
    x = NULL,
    y = "Coefficient of variation in K among partners",
    colour = NULL,
    title = "Interactions of the same species differ in local support",
    subtitle = "Only species with at least two realised regional partners are included"
  )

## ------------------------------------------------------------
## Figure 4: heterogeneity among datasets
## ------------------------------------------------------------

plot_datasets <- ggplot2::ggplot(
  dataset_data,
  ggplot2::aes(
    x = mean_K,
    y = gini_K,
    size = number_realised_links
  )
) +
  ggplot2::geom_point(
    alpha = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = dataset
    ),
    nudge_y = 0.015,
    check_overlap = TRUE,
    size = 3
  ) +
  ggplot2::scale_x_continuous(
    trans = "log1p"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.5,
      7
    )
  ) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0
    )
  ) +
  ggplot2::labs(
    x = "Mean local interaction support, K",
    y = "Gini inequality in K among links",
    size = "Regional links",
    title = "Datasets differ in both typical support and support inequality",
    subtitle = "Higher Gini values indicate that local support is concentrated in fewer regional links"
  )

## ------------------------------------------------------------
## Figure 5: standardised support across datasets
## ------------------------------------------------------------

plot_standardised <- ggplot2::ggplot(
  link_data,
  ggplot2::aes(
    x = dataset,
    y = K_standardised
  )
) +
  ggplot2::geom_violin(
    fill = "grey88",
    colour = "grey55",
    linewidth = 0.35,
    scale = "width",
    trim = TRUE
  ) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.45
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.15,
    size = 0.7,
    colour = "grey30"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  theme_pub() +
  ggplot2::labs(
    x = NULL,
    y = "Dataset sites supporting each regional interaction",
    title = "Support heterogeneity remains after accounting for dataset size",
    subtitle = "Local support K divided by the number of sampled sites in each dataset"
  )

## ------------------------------------------------------------
## Combined 2 × 2 diagnostic figure
## ------------------------------------------------------------

combined_title <- grid::textGrob(
  "Heterogeneity in local support across links, species, and datasets",
  gp = grid::gpar(
    fontface = "bold",
    fontsize = 14
  )
)

combined_body <- gridExtra::arrangeGrob(
  plot_links,
  plot_species_mean,
  plot_species_cv,
  plot_datasets,
  ncol = 2
)

combined_figure <- gridExtra::arrangeGrob(
  combined_title,
  combined_body,
  ncol = 1,
  heights = c(
    0.06,
    1
  )
)

## ------------------------------------------------------------
## Validation checks
## ------------------------------------------------------------
validation_checks <- tibble::tibble(
  check = c(
    "Every retained link has full_K greater than zero",
    "Every retained link has full_K no greater than full_n",
    "Every standardised K value is between zero and one",
    "Every species has initial degree greater than zero",
    "Within-species CV is calculated only for species with at least two partners",
    "Consumers and Resources are summarised separately",
    "Every dataset has one dataset-level summary"
  ),
  pass = c(
    all(
      link_data$full_K > 0
    ),
    all(
      link_data$full_K <=
        link_data$full_n
    ),
    all(
      link_data$K_standardised > 0 &
        link_data$K_standardised <= 1
    ),
    all(
      species_data$initial_degree > 0
    ),
    all(
      species_heterogeneity$
        initial_degree >= 2
    ),
    all(
      guild_levels %in%
        as.character(
          unique(species_data$guild)
        )
    ),
    dataset_data %>%
      dplyr::count(
        dataset,
        name = "number_rows"
      ) %>%
      dplyr::summarise(
        valid = all(
          number_rows == 1
        )
      ) %>%
      dplyr::pull(valid)
  )
)

## ------------------------------------------------------------
## Save figures
## ------------------------------------------------------------

save_png(
  plot_links,
  "39A_local_support_among_links_by_dataset.png",
  width = 11,
  height = 6
)

save_png(
  plot_species_mean,
  "39B_mean_local_support_among_species.png",
  width = 11,
  height = 6
)

save_png(
  plot_species_cv,
  "39C_within_species_local_support_heterogeneity.png",
  width = 11,
  height = 6
)

save_png(
  plot_datasets,
  "39D_dataset_support_mean_and_inequality.png",
  width = 8,
  height = 6
)

save_png(
  plot_standardised,
  "39E_standardised_local_support_among_links.png",
  width = 11,
  height = 6
)

save_png_grid(
  combined_figure,
  "39_local_support_heterogeneity_combined.png",
  width = 14,
  height = 10
)

## ------------------------------------------------------------
## Save descriptive tables
## ------------------------------------------------------------

write.csv2(
  link_data,
  file.path(
    out_dir,
    "39_link_level_local_support.csv"
  ),
  row.names = FALSE
)

write.csv2(
  species_data,
  file.path(
    out_dir,
    "39_species_level_local_support_summary.csv"
  ),
  row.names = FALSE
)

write.csv2(
  dataset_data,
  file.path(
    out_dir,
    "39_dataset_level_local_support_summary.csv"
  ),
  row.names = FALSE
)

write.csv2(
  validation_checks,
  file.path(
    out_dir,
    "39_local_support_validation_checks.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Console summary
## ------------------------------------------------------------

message(
  "\nDataset-level support summaries:"
)

print(
  dataset_data %>%
    dplyr::arrange(
      dplyr::desc(gini_K)
    ),
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
  "Figure 39A shows variation in absolute local support among links.\n",
  "Figure 39B shows variation in typical link support among species.\n",
  "Figure 39C shows how unevenly support is distributed among the links\n",
  "of the same species. Figure 39D compares typical support and inequality\n",
  "among datasets. Figure 39E repeats the link-level comparison after\n",
  "standardising K by the number of sampled sites in each dataset."
)

message(
  "\nSaved outputs in: ",
  out_dir
)