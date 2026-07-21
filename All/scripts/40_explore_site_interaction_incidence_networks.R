## ------------------------------------------------------------
## Script:
## All/scripts/40_explore_site_interaction_incidence_networks.R
##
## Goal:
## Explore the hidden site × interaction incidence structure behind
## regional binary interaction networks.
##
## Conceptual framing:
##
## Let M be a site × regional-interaction incidence matrix:
##
##   M[s, l] = 1 when regional interaction link l is observed at site s.
##
## Deleting a site removes one row of M.
##
## Local interaction support declines whenever entries are removed.
## A regional interaction disappears only when its entire column becomes
## empty. For a focal species, overlap among its interaction columns
## determines whether partner losses are spatially separated or bundled.
##
## This script is exploratory. It does not identify "important sites"
## because sites have no independent metadata.
##
## Outputs:
##   - dataset incidence heatmaps;
##   - focal-species portfolio matrices;
##   - incidence summaries;
##   - link- and site-level distribution plots;
##   - degree–bundling summaries;
##   - one real-data separated-versus-bundled schematic;
##   - CSV summary tables.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c(
  "dplyr",
  "ggplot2",
  "tibble",
  "purrr",
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
## User controls
## ------------------------------------------------------------

set.seed(123)

## Use 999 null replicates when the number of eligible focal species
## in a dataset is manageable. Use 250 for larger datasets.
n_null_portfolio_default <- 999L
n_null_portfolio_large <- 250L
maximum_eligible_species_for_999 <- 100L

n_boot <- 5000L

low_support_thresholds <- c(
  "K ≤ 2" = 2L,
  "K = 1" = 1L
)

out_dir <- "All/outputs/40_incidence_networks"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

guild_levels <- c(
  "Consumer",
  "Resource"
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
## Plotting helpers
## ------------------------------------------------------------

theme_pub <- function(base_size = 9){
  
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

theme_matrix <- function(base_size = 8){
  
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1,
        hjust = 0
      ),
      plot.margin = grid::unit(
        c(3, 3, 3, 3),
        "pt"
      )
    )
}

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

save_grob_png <- function(
    grob,
    filename,
    width = 10,
    height = 7
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

sanitize_filename <- function(x){
  
  x <- as.character(x)
  
  x <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- gsub(
    "^_|_$",
    "",
    x
  )
  
  ifelse(
    nchar(x) == 0,
    "unnamed",
    x
  )
}

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
## Matrix ordering
## ------------------------------------------------------------

safe_cluster_order <- function(
    matrix_input,
    margin = c("rows", "columns"),
    fallback_score = NULL
){
  
  margin <- match.arg(margin)
  
  number_objects <- if(
    margin == "rows"
  ){
    nrow(matrix_input)
  } else {
    ncol(matrix_input)
  }
  
  if(number_objects <= 1){
    return(seq_len(number_objects))
  }
  
  clustering_matrix <- if(
    margin == "rows"
  ){
    matrix_input
  } else {
    t(matrix_input)
  }
  
  clustered_order <- tryCatch(
    {
      distance_object <- stats::dist(
        clustering_matrix,
        method = "euclidean"
      )
      
      if(any(!is.finite(distance_object))){
        stop("Non-finite clustering distance.")
      }
      
      stats::hclust(
        distance_object,
        method = "average"
      )$order
    },
    error = function(e){
      NULL
    }
  )
  
  if(!is.null(clustered_order)){
    return(clustered_order)
  }
  
  if(is.null(fallback_score)){
    
    fallback_score <- if(
      margin == "rows"
    ){
      rowSums(matrix_input)
    } else {
      colSums(matrix_input)
    }
  }
  
  order(
    fallback_score,
    decreasing = TRUE,
    na.last = TRUE
  )
}

matrix_to_long <- function(matrix_input){
  
  grid_data <- expand.grid(
    row_index = seq_len(
      nrow(matrix_input)
    ),
    column_index = seq_len(
      ncol(matrix_input)
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  grid_data$value <- as.vector(
    matrix_input
  )
  
  tibble::as_tibble(
    grid_data
  )
}

## ------------------------------------------------------------
## Jaccard and portfolio-null helpers
## ------------------------------------------------------------

mean_pairwise_jaccard_matrix <- function(
    incidence_matrix
){
  
  if(ncol(incidence_matrix) < 2){
    return(NA_real_)
  }
  
  overlap_matrix <- crossprod(
    incidence_matrix
  )
  
  support <- diag(
    overlap_matrix
  )
  
  union_matrix <- outer(
    support,
    support,
    "+"
  ) - overlap_matrix
  
  upper_indices <- upper.tri(
    overlap_matrix
  )
  
  denominator <- union_matrix[
    upper_indices
  ]
  
  numerator <- overlap_matrix[
    upper_indices
  ]
  
  valid <- denominator > 0
  
  if(!any(valid)){
    return(NA_real_)
  }
  
  mean(
    numerator[valid] /
      denominator[valid]
  )
}

make_support_matrix <- function(
    site_sets,
    all_sites
){
  
  site_lookup <- stats::setNames(
    seq_along(all_sites),
    all_sites
  )
  
  out <- matrix(
    0L,
    nrow = length(all_sites),
    ncol = length(site_sets)
  )
  
  for(column_index in seq_along(site_sets)){
    
    selected_sites <- site_sets[[column_index]
    ]
    
    out[
      site_lookup[selected_sites],
      column_index
    ] <- 1L
  }
  
  out
}

simulate_portfolio_null <- function(
    cooccurrence_site_sets,
    support_sizes,
    all_sites,
    n_null
){
  
  null_values <- numeric(
    n_null
  )
  
  exact_support_check <- TRUE
  cooccurrence_check <- TRUE
  
  for(null_index in seq_len(n_null)){
    
    null_sets <- Map(
      function(possible_sites, support_size){
        
        selected <- sample(
          possible_sites,
          size = support_size,
          replace = FALSE
        )
        
        exact_support_check <<-
          exact_support_check &&
          length(selected) == support_size
        
        cooccurrence_check <<-
          cooccurrence_check &&
          all(selected %in% possible_sites)
        
        selected
      },
      cooccurrence_site_sets,
      support_sizes
    )
    
    null_matrix <- make_support_matrix(
      site_sets = null_sets,
      all_sites = all_sites
    )
    
    null_values[null_index] <-
      mean_pairwise_jaccard_matrix(
        null_matrix
      )
  }
  
  list(
    null_values = null_values,
    exact_support_check =
      exact_support_check,
    cooccurrence_check =
      cooccurrence_check
  )
}

## ------------------------------------------------------------
## Data preparation for one dataset
## ------------------------------------------------------------

prepare_one_dataset <- function(dataset){
  
  message(
    "Preparing incidence structure: ",
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
  
  all_sites <- sort(
    unique(cooc$site)
  )
  
  pair_cooccurrence <- cooc %>%
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
  
  pair_interactions <- ints %>%
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
  
  pair_table <- pair_cooccurrence %>%
    dplyr::left_join(
      pair_interactions,
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
      link_id = paste(
        consumer,
        resource,
        sep = "___"
      ),
      link_index = dplyr::row_number()
    )
  
  if(nrow(realised_links) == 0){
    stop(
      "No realised regional interactions in ",
      dataset
    )
  }
  
  interactions_inside_cooccurrence <- all(
    unlist(
      Map(
        function(interaction_sites, possible_sites){
          all(
            interaction_sites %in%
              possible_sites
          )
        },
        realised_links$interaction_sites,
        realised_links$cooccurrence_sites
      ),
      use.names = FALSE
    )
  )
  
  if(!interactions_inside_cooccurrence){
    stop(
      "Observed interaction outside pair-specific co-occurrence set in ",
      dataset
    )
  }
  
  incidence_matrix <- make_support_matrix(
    site_sets =
      realised_links$interaction_sites,
    all_sites =
      all_sites
  )
  
  rownames(incidence_matrix) <-
    all_sites
  
  colnames(incidence_matrix) <-
    realised_links$link_id
  
  site_load <- rowSums(
    incidence_matrix
  )
  
  link_support <- colSums(
    incidence_matrix
  )
  
  if(any(link_support != realised_links$full_K)){
    stop(
      "Incidence-matrix column sums do not equal full_K in ",
      dataset
    )
  }
  
  list(
    dataset = dataset,
    cooc = cooc,
    ints = ints,
    all_sites = all_sites,
    realised_links = realised_links,
    incidence_matrix = incidence_matrix,
    site_load = site_load,
    link_support = link_support
  )
}

## ------------------------------------------------------------
## Dataset incidence heatmap
## ------------------------------------------------------------

make_incidence_heatmap <- function(
    prepared
){
  
  incidence_matrix <-
    prepared$incidence_matrix
  
  row_order <- safe_cluster_order(
    matrix_input = incidence_matrix,
    margin = "rows",
    fallback_score =
      prepared$site_load
  )
  
  column_order <- safe_cluster_order(
    matrix_input = incidence_matrix,
    margin = "columns",
    fallback_score =
      prepared$link_support
  )
  
  ordered_matrix <- incidence_matrix[
    row_order,
    column_order,
    drop = FALSE
  ]
  
  long_matrix <- matrix_to_long(
    ordered_matrix
  )
  
  main_plot <- ggplot2::ggplot(
    long_matrix,
    ggplot2::aes(
      x = column_index,
      y = row_index,
      fill = factor(value)
    )
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(
      values = c(
        "0" = "#F5F5F5",
        "1" = "#303030"
      ),
      guide = "none"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_cartesian(
      expand = FALSE
    ) +
    theme_matrix() +
    ggplot2::labs(
      title = paste0(
        prepared$dataset,
        ": ",
        nrow(incidence_matrix),
        " sites, ",
        ncol(incidence_matrix),
        " regional links, ",
        sum(incidence_matrix),
        " local interaction records"
      )
    )
  
  top_data <- tibble::tibble(
    column_index = seq_len(
      ncol(ordered_matrix)
    ),
    K = colSums(
      ordered_matrix
    )
  )
  
  top_plot <- ggplot2::ggplot(
    top_data,
    ggplot2::aes(
      x = column_index,
      y = K
    )
  ) +
    ggplot2::geom_col(
      fill = "#777777",
      width = 1
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic(
      base_size = 7
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = grid::unit(
        c(1, 2, 1, 2),
        "pt"
      )
    ) +
    ggplot2::labs(
      y = "K"
    )
  
  side_data <- tibble::tibble(
    row_index = seq_len(
      nrow(ordered_matrix)
    ),
    site_load = rowSums(
      ordered_matrix
    )
  )
  
  side_plot <- ggplot2::ggplot(
    side_data,
    ggplot2::aes(
      x = site_load,
      y = row_index
    )
  ) +
    ggplot2::geom_col(
      fill = "#777777",
      width = 1
    ) +
    ggplot2::scale_y_reverse(
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic(
      base_size = 7
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = grid::unit(
        c(2, 1, 2, 1),
        "pt"
      )
    ) +
    ggplot2::labs(
      x = "Site load"
    )
  
  blank_plot <- grid::nullGrob()
  
  gridExtra::arrangeGrob(
    top_plot,
    blank_plot,
    main_plot,
    side_plot,
    ncol = 2,
    widths = c(
      1,
      0.22
    ),
    heights = c(
      0.22,
      1
    )
  )
}

## ------------------------------------------------------------
## Focal-species portfolio preparation
## ------------------------------------------------------------

make_focal_link_table <- function(
    realised_links
){
  
  dplyr::bind_rows(
    
    realised_links %>%
      dplyr::transmute(
        dataset,
        guild = "Consumer",
        species = consumer,
        partner = resource,
        link_index,
        full_n,
        full_K,
        cooccurrence_sites,
        interaction_sites
      ),
    
    realised_links %>%
      dplyr::transmute(
        dataset,
        guild = "Resource",
        species = resource,
        partner = consumer,
        link_index,
        full_n,
        full_K,
        cooccurrence_sites,
        interaction_sites
      )
  ) %>%
    dplyr::mutate(
      guild = factor(
        guild,
        levels = guild_levels
      )
    )
}

analyse_portfolios <- function(
    prepared
){
  
  focal_links <- make_focal_link_table(
    prepared$realised_links
  )
  
  degree_table <- focal_links %>%
    dplyr::group_by(
      dataset,
      guild,
      species
    ) %>%
    dplyr::summarise(
      initial_degree =
        dplyr::n_distinct(partner),
      .groups = "drop"
    )
  
  eligible <- degree_table %>%
    dplyr::filter(
      initial_degree >= 2
    )
  
  number_eligible <- nrow(
    eligible
  )
  
  n_null_used <- if(
    number_eligible <=
    maximum_eligible_species_for_999
  ){
    n_null_portfolio_default
  } else {
    n_null_portfolio_large
  }
  
  message(
    "  Portfolio null replicates for ",
    prepared$dataset,
    ": ",
    n_null_used,
    " (",
    number_eligible,
    " eligible focal species)"
  )
  
  if(number_eligible == 0){
    
    return(
      list(
        summary = tibble::tibble(),
        checks = tibble::tibble(
          dataset = prepared$dataset,
          null_exact_support = TRUE,
          null_within_cooccurrence = TRUE
        )
      )
    )
  }
  
  summary_rows <- vector(
    "list",
    number_eligible
  )
  
  all_exact_support <- TRUE
  all_within_cooccurrence <- TRUE
  
  for(species_index in seq_len(
    number_eligible
  )){
    
    eligible_row <- eligible[
      species_index,
    ]
    
    focal_data <- focal_links %>%
      dplyr::filter(
        guild ==
          eligible_row$guild,
        species ==
          eligible_row$species
      ) %>%
      dplyr::arrange(
        partner
      )
    
    observed_matrix <- make_support_matrix(
      site_sets =
        focal_data$interaction_sites,
      all_sites =
        prepared$all_sites
    )
    
    observed_overlap <-
      mean_pairwise_jaccard_matrix(
        observed_matrix
      )
    
    null_result <- simulate_portfolio_null(
      cooccurrence_site_sets =
        focal_data$cooccurrence_sites,
      support_sizes =
        focal_data$full_K,
      all_sites =
        prepared$all_sites,
      n_null =
        n_null_used
    )
    
    all_exact_support <-
      all_exact_support &&
      null_result$exact_support_check
    
    all_within_cooccurrence <-
      all_within_cooccurrence &&
      null_result$cooccurrence_check
    
    summary_rows[[species_index]] <-
      tibble::tibble(
        dataset = prepared$dataset,
        guild = as.character(
          eligible_row$guild
        ),
        species = as.character(
          eligible_row$species
        ),
        initial_degree =
          eligible_row$initial_degree,
        observed_overlap =
          observed_overlap,
        null_mean_overlap =
          mean(
            null_result$null_values,
            na.rm = TRUE
          ),
        null_q025 =
          safe_quantile(
            null_result$null_values,
            0.025
          ),
        null_q975 =
          safe_quantile(
            null_result$null_values,
            0.975
          ),
        bundling_delta =
          observed_overlap -
          mean(
            null_result$null_values,
            na.rm = TRUE
          ),
        null_replicates =
          n_null_used
      )
  }
  
  list(
    summary = dplyr::bind_rows(
      summary_rows
    ),
    checks = tibble::tibble(
      dataset = prepared$dataset,
      null_exact_support =
        all_exact_support,
      null_within_cooccurrence =
        all_within_cooccurrence
    )
  )
}

## ------------------------------------------------------------
## Representative portfolio selection
## ------------------------------------------------------------

select_representative_portfolios <- function(
    bundling_summary
){
  
  groups <- bundling_summary %>%
    dplyr::group_by(
      dataset,
      guild
    ) %>%
    dplyr::group_split(
      .keep = TRUE
    )
  
  selected <- lapply(
    groups,
    function(group_data){
      
      preferred <- group_data %>%
        dplyr::filter(
          initial_degree >= 3
        )
      
      if(nrow(preferred) == 0){
        
        preferred <- group_data %>%
          dplyr::filter(
            initial_degree >= 2
          )
      }
      
      if(nrow(preferred) == 0){
        return(tibble::tibble())
      }
      
      candidate_rows <- dplyr::bind_rows(
        
        preferred %>%
          dplyr::slice_min(
            bundling_delta,
            n = 1,
            with_ties = FALSE
          ) %>%
          dplyr::mutate(
            portfolio_type =
              "spatially_separated_portfolio"
          ),
        
        preferred %>%
          dplyr::slice_min(
            abs(bundling_delta),
            n = 1,
            with_ties = FALSE
          ) %>%
          dplyr::mutate(
            portfolio_type =
              "null_like_portfolio"
          ),
        
        preferred %>%
          dplyr::slice_max(
            bundling_delta,
            n = 1,
            with_ties = FALSE
          ) %>%
          dplyr::mutate(
            portfolio_type =
              "spatially_bundled_portfolio"
          )
      )
      
      ## A species may be simultaneously the most extreme and closest
      ## to zero in very small groups. Retain it once.
      candidate_rows %>%
        dplyr::distinct(
          dataset,
          guild,
          species,
          .keep_all = TRUE
        )
    }
  )
  
  dplyr::bind_rows(
    selected
  )
}

## ------------------------------------------------------------
## Portfolio matrix plot
## ------------------------------------------------------------

make_portfolio_matrix_plot <- function(
    prepared,
    selection_row
){
  
  focal_links <- make_focal_link_table(
    prepared$realised_links
  ) %>%
    dplyr::filter(
      guild ==
        selection_row$guild,
      species ==
        selection_row$species
    ) %>%
    dplyr::arrange(
      partner
    )
  
  portfolio_matrix <- make_support_matrix(
    site_sets =
      focal_links$interaction_sites,
    all_sites =
      prepared$all_sites
  )
  
  rownames(portfolio_matrix) <-
    prepared$all_sites
  
  colnames(portfolio_matrix) <-
    focal_links$partner
  
  row_order <- safe_cluster_order(
    portfolio_matrix,
    margin = "rows",
    fallback_score =
      rowSums(portfolio_matrix)
  )
  
  column_order <- safe_cluster_order(
    portfolio_matrix,
    margin = "columns",
    fallback_score =
      colSums(portfolio_matrix)
  )
  
  ordered_matrix <- portfolio_matrix[
    row_order,
    column_order,
    drop = FALSE
  ]
  
  plot_data <- matrix_to_long(
    ordered_matrix
  )
  
  show_partner_labels <-
    ncol(ordered_matrix) <= 12
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = column_index,
      y = row_index,
      fill = factor(value)
    )
  ) +
    ggplot2::geom_tile(
      colour = "#E5E5E5",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "0" = "#F5F5F5",
        "1" = "#303030"
      ),
      guide = "none"
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_continuous(
      breaks = if(show_partner_labels){
        seq_len(
          ncol(ordered_matrix)
        )
      } else {
        NULL
      },
      labels = if(show_partner_labels){
        paste0(
          "P",
          seq_len(
            ncol(ordered_matrix)
          )
        )
      } else {
        NULL
      }
    ) +
    ggplot2::coord_fixed(
      ratio = 1
    ) +
    theme_matrix() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 6
      )
    ) +
    ggplot2::labs(
      title = gsub(
        "_",
        " ",
        selection_row$portfolio_type
      ),
      subtitle = paste0(
        selection_row$dataset,
        " | ",
        selection_row$guild,
        " | degree ",
        selection_row$initial_degree,
        "\nObserved overlap = ",
        signif(
          selection_row$observed_overlap,
          3
        ),
        "; null mean = ",
        signif(
          selection_row$null_mean_overlap,
          3
        ),
        "; delta = ",
        signif(
          selection_row$bundling_delta,
          3
        )
      )
    )
}

## ------------------------------------------------------------
## Run data preparation
## ------------------------------------------------------------

prepared_datasets <- lapply(
  all_dataset_names,
  prepare_one_dataset
)

names(prepared_datasets) <-
  all_dataset_names

## ------------------------------------------------------------
## Part 1: dataset incidence heatmaps
## ------------------------------------------------------------

incidence_heatmap_grobs <- list()

for(dataset in all_dataset_names){
  
  prepared <- prepared_datasets[[dataset]
  ]
  
  heatmap_grob <- make_incidence_heatmap(
    prepared
  )
  
  incidence_heatmap_grobs[[dataset]
  ] <- heatmap_grob
  
  save_grob_png(
    grob = heatmap_grob,
    filename = paste0(
      "40_incidence_heatmap_",
      sanitize_filename(dataset),
      ".png"
    ),
    width = 10,
    height = 7
  )
}

if(length(incidence_heatmap_grobs) > 0){
  
  heatmap_columns <- min(
    3,
    length(incidence_heatmap_grobs)
  )
  
  heatmap_rows <- ceiling(
    length(incidence_heatmap_grobs) /
      heatmap_columns
  )
  
  incidence_overview <- do.call(
    gridExtra::arrangeGrob,
    c(
      incidence_heatmap_grobs,
      list(
        ncol = heatmap_columns
      )
    )
  )
  
  save_grob_png(
    grob = incidence_overview,
    filename =
      "40_incidence_heatmaps_all_datasets.png",
    width = 5.2 * heatmap_columns,
    height = 4.2 * heatmap_rows
  )
}

## ------------------------------------------------------------
## Part 2: portfolio-level bundling
## ------------------------------------------------------------

portfolio_outputs <- lapply(
  prepared_datasets,
  analyse_portfolios
)

portfolio_bundling_summary <- dplyr::bind_rows(
  lapply(
    portfolio_outputs,
    `[[`,
    "summary"
  )
) %>%
  dplyr::mutate(
    guild = factor(
      guild,
      levels = guild_levels
    ),
    degree_class = factor(
      make_degree_class(initial_degree),
      levels = degree_class_levels
    )
  )

portfolio_null_checks <- dplyr::bind_rows(
  lapply(
    portfolio_outputs,
    `[[`,
    "checks"
  )
)

write.csv2(
  portfolio_bundling_summary,
  file.path(
    out_dir,
    "40_species_portfolio_bundling_summary.csv"
  ),
  row.names = FALSE
)

representatives <- select_representative_portfolios(
  portfolio_bundling_summary
)

write.csv2(
  representatives,
  file.path(
    out_dir,
    "40_representative_portfolios.csv"
  ),
  row.names = FALSE
)

portfolio_plot_grobs <- list()
portfolio_plots_by_dataset <- list()

if(nrow(representatives) > 0){
  
  for(row_index in seq_len(
    nrow(representatives)
  )){
    
    selection_row <- representatives[
      row_index,
    ]
    
    prepared <- prepared_datasets[[
      as.character(
        selection_row$dataset
      )
    ]
    ]
    
    portfolio_plot <- make_portfolio_matrix_plot(
      prepared = prepared,
      selection_row = selection_row
    )
    
    plot_key <- paste(
      selection_row$dataset,
      selection_row$guild,
      selection_row$species,
      selection_row$portfolio_type,
      sep = "___"
    )
    
    portfolio_plot_grobs[[plot_key]
    ] <- ggplot2::ggplotGrob(
      portfolio_plot
    )
    
    dataset_name <- as.character(
      selection_row$dataset
    )
    
    if(is.null(
      portfolio_plots_by_dataset[[dataset_name]
      ]
    )){
      portfolio_plots_by_dataset[[dataset_name]
      ] <- list()
    }
    
    portfolio_plots_by_dataset[[dataset_name]
    ][[plot_key]
    ] <- ggplot2::ggplotGrob(
      portfolio_plot
    )
    
    save_png(
      plot = portfolio_plot,
      filename = paste0(
        "40_portfolio_matrix_",
        sanitize_filename(
          selection_row$dataset
        ),
        "_",
        sanitize_filename(
          selection_row$guild
        ),
        "_",
        sanitize_filename(
          selection_row$species
        ),
        "_",
        sanitize_filename(
          selection_row$portfolio_type
        ),
        ".png"
      ),
      width = 5.5,
      height = 5
    )
  }
  
  for(dataset_name in names(
    portfolio_plots_by_dataset
  )){
    
    dataset_grobs <-
      portfolio_plots_by_dataset[[dataset_name]
      ]
    
    dataset_panel <- do.call(
      gridExtra::arrangeGrob,
      c(
        dataset_grobs,
        list(
          ncol = min(
            3,
            length(dataset_grobs)
          )
        )
      )
    )
    
    save_grob_png(
      grob = dataset_panel,
      filename = paste0(
        "40_portfolio_examples_",
        sanitize_filename(
          dataset_name
        ),
        ".png"
      ),
      width = 5.2 * min(
        3,
        length(dataset_grobs)
      ),
      height = 5 * ceiling(
        length(dataset_grobs) / 3
      )
    )
  }
  
  contact_columns <- min(
    4,
    length(portfolio_plot_grobs)
  )
  
  contact_rows <- ceiling(
    length(portfolio_plot_grobs) /
      contact_columns
  )
  
  contact_sheet <- do.call(
    gridExtra::arrangeGrob,
    c(
      portfolio_plot_grobs,
      list(
        ncol = contact_columns
      )
    )
  )
  
  save_grob_png(
    grob = contact_sheet,
    filename =
      "40_portfolio_examples_all_datasets.png",
    width = 4.6 * contact_columns,
    height = 4.4 * contact_rows
  )
}

## ------------------------------------------------------------
## Part 3: incidence-structure summaries
## ------------------------------------------------------------

link_level_data <- dplyr::bind_rows(
  lapply(
    prepared_datasets,
    function(prepared){
      
      prepared$realised_links %>%
        dplyr::transmute(
          dataset =
            prepared$dataset,
          link_id,
          consumer,
          resource,
          full_K,
          full_n
        )
    }
  )
)

site_level_data <- dplyr::bind_rows(
  lapply(
    prepared_datasets,
    function(prepared){
      
      incidence <- prepared$incidence_matrix
      
      tibble::tibble(
        dataset = prepared$dataset,
        site = rownames(incidence),
        site_interaction_load =
          rowSums(incidence),
        number_low_support_K_le_2 =
          rowSums(
            incidence[
              ,
              prepared$link_support <= 2,
              drop = FALSE
            ]
          ),
        number_single_site_links =
          rowSums(
            incidence[
              ,
              prepared$link_support == 1,
              drop = FALSE
            ]
          )
      ) %>%
        dplyr::mutate(
          share_low_support_K_le_2 =
            dplyr::if_else(
              site_interaction_load > 0,
              number_low_support_K_le_2 /
                site_interaction_load,
              NA_real_
            ),
          share_single_site_links =
            dplyr::if_else(
              site_interaction_load > 0,
              number_single_site_links /
                site_interaction_load,
              NA_real_
            )
        )
    }
  )
)

dataset_incidence_summary <- dplyr::bind_rows(
  lapply(
    prepared_datasets,
    function(prepared){
      
      tibble::tibble(
        dataset = prepared$dataset,
        number_sites =
          nrow(
            prepared$incidence_matrix
          ),
        number_regional_interaction_links =
          ncol(
            prepared$incidence_matrix
          ),
        total_local_interaction_records =
          sum(
            prepared$incidence_matrix
          ),
        mean_K =
          mean(
            prepared$link_support
          ),
        median_K =
          stats::median(
            prepared$link_support
          ),
        proportion_links_K_equal_1 =
          mean(
            prepared$link_support == 1
          ),
        mean_site_interaction_load =
          mean(
            prepared$site_load
          ),
        median_site_interaction_load =
          stats::median(
            prepared$site_load
          ),
        coefficient_variation_K =
          safe_cv(
            prepared$link_support
          ),
        coefficient_variation_site_load =
          safe_cv(
            prepared$site_load
          )
      )
    }
  )
)

write.csv2(
  dataset_incidence_summary,
  file.path(
    out_dir,
    "40_dataset_incidence_summary.csv"
  ),
  row.names = FALSE
)

write.csv2(
  site_level_data,
  file.path(
    out_dir,
    "40_site_incidence_summary.csv"
  ),
  row.names = FALSE
)

dataset_order <- dataset_incidence_summary %>%
  dplyr::arrange(
    mean_K
  ) %>%
  dplyr::pull(
    dataset
  )

link_level_data <- link_level_data %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = dataset_order
    )
  )

site_level_data <- site_level_data %>%
  dplyr::mutate(
    dataset = factor(
      dataset,
      levels = dataset_order
    )
  )

link_support_plot <- ggplot2::ggplot(
  link_level_data,
  ggplot2::aes(
    x = dataset,
    y = full_K
  )
) +
  ggplot2::geom_violin(
    fill = "#D9D9D9",
    colour = "#666666",
    linewidth = 0.35,
    scale = "width",
    trim = TRUE
  ) +
  ggplot2::geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.4
  ) +
  ggplot2::scale_y_continuous(
    trans = "log1p"
  ) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of sites supporting each interaction, K",
    title =
      "Heterogeneity in site–interaction incidence structure",
    subtitle =
      "Distribution of regional-link support across datasets"
  )

site_load_plot <- ggplot2::ggplot(
  site_level_data,
  ggplot2::aes(
    x = dataset,
    y = site_interaction_load
  )
) +
  ggplot2::geom_violin(
    fill = "#D9D9D9",
    colour = "#666666",
    linewidth = 0.35,
    scale = "width",
    trim = TRUE
  ) +
  ggplot2::geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.4
  ) +
  ggplot2::scale_y_continuous(
    trans = "log1p"
  ) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of interactions observed at each site",
    title =
      "Heterogeneity in site–interaction incidence structure",
    subtitle =
      "Distribution of local interaction load across datasets"
  )

site_threshold_long <- dplyr::bind_rows(
  
  site_level_data %>%
    dplyr::transmute(
      dataset,
      site,
      site_interaction_load,
      threshold = "K ≤ 2",
      low_support_count =
        number_low_support_K_le_2,
      low_support_share =
        share_low_support_K_le_2
    ),
  
  site_level_data %>%
    dplyr::transmute(
      dataset,
      site,
      site_interaction_load,
      threshold = "K = 1",
      low_support_count =
        number_single_site_links,
      low_support_share =
        share_single_site_links
    )
)

site_count_relationship_plot <- ggplot2::ggplot(
  site_threshold_long,
  ggplot2::aes(
    x = site_interaction_load,
    y = low_support_count
  )
) +
  ggplot2::geom_point(
    colour = "#666666",
    alpha = 0.35,
    size = 1.2
  ) +
  ggplot2::geom_smooth(
    method = "loess",
    se = FALSE,
    colour = "black",
    linewidth = 0.8
  ) +
  ggplot2::facet_grid(
    threshold ~ dataset,
    scales = "free"
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0
    )
  ) +
  ggplot2::labs(
    x = "Site interaction load",
    y = "Low-support interactions hosted at site",
    title =
      "Heterogeneity in site–interaction incidence structure",
    subtitle = paste0(
      "Relationship between total local interaction load and ",
      "the number of low-support regional links"
    )
  )

site_share_relationship_plot <- ggplot2::ggplot(
  site_threshold_long,
  ggplot2::aes(
    x = site_interaction_load,
    y = low_support_share
  )
) +
  ggplot2::geom_point(
    colour = "#666666",
    alpha = 0.35,
    size = 1.2
  ) +
  ggplot2::geom_smooth(
    method = "loess",
    se = FALSE,
    colour = "black",
    linewidth = 0.8
  ) +
  ggplot2::facet_grid(
    threshold ~ dataset,
    scales = "free_x"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(
      0,
      1
    )
  ) +
  theme_pub(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0
    )
  ) +
  ggplot2::labs(
    x = "Site interaction load",
    y = "Share of site's interactions with low regional support",
    title =
      "Heterogeneity in site–interaction incidence structure",
    subtitle = paste0(
      "Sensitivity shown for links supported in at most two sites ",
      "and links supported at exactly one site"
    )
  )

save_png(
  link_support_plot,
  "40_link_support_distribution.png",
  width = 11,
  height = 6
)

save_png(
  site_load_plot,
  "40_site_load_distribution.png",
  width = 11,
  height = 6
)

save_png(
  site_count_relationship_plot,
  "40_site_load_vs_low_support_links.png",
  width = 14,
  height = 7
)

save_png(
  site_share_relationship_plot,
  "40_site_load_vs_share_low_support_links.png",
  width = 14,
  height = 7
)

## ------------------------------------------------------------
## Part 4: degree and bundling overview
## ------------------------------------------------------------

dataset_degree_means <- portfolio_bundling_summary %>%
  dplyr::group_by(
    dataset,
    guild,
    initial_degree
  ) %>%
  dplyr::summarise(
    mean_bundling_delta =
      mean(
        bundling_delta,
        na.rm = TRUE
      ),
    number_species =
      dplyr::n(),
    .groups = "drop"
  )

equal_dataset_degree_means <- dataset_degree_means %>%
  dplyr::group_by(
    guild,
    initial_degree
  ) %>%
  dplyr::summarise(
    equal_dataset_mean_delta =
      mean(
        mean_bundling_delta,
        na.rm = TRUE
      ),
    number_datasets =
      dplyr::n_distinct(dataset),
    .groups = "drop"
  )

bundling_degree_plot <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data =
      portfolio_bundling_summary,
    ggplot2::aes(
      x = initial_degree,
      y = bundling_delta
    ),
    colour = "grey60",
    alpha = 0.18,
    size = 1
  ) +
  ggplot2::geom_line(
    data =
      dataset_degree_means,
    ggplot2::aes(
      x = initial_degree,
      y = mean_bundling_delta,
      group = dataset
    ),
    colour = "grey65",
    alpha = 0.35,
    linewidth = 0.45
  ) +
  ggplot2::geom_line(
    data =
      equal_dataset_degree_means,
    ggplot2::aes(
      x = initial_degree,
      y = equal_dataset_mean_delta
    ),
    colour = "black",
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data =
      equal_dataset_degree_means,
    ggplot2::aes(
      x = initial_degree,
      y = equal_dataset_mean_delta
    ),
    colour = "black",
    size = 2
  ) +
  ggplot2::facet_wrap(
    ~ guild,
    nrow = 1,
    scales = "free_x"
  ) +
  ggplot2::scale_x_continuous(
    trans = "log1p"
  ) +
  theme_pub() +
  ggplot2::labs(
    x = "Initial regional degree",
    y = "Observed minus null mean Jaccard overlap",
    title =
      "Spatial portfolio structure across regional degree",
    subtitle = paste0(
      "Grey points are focal species; thin lines are dataset means; ",
      "black lines are equal-dataset means"
    )
  )

dataset_class_means <- portfolio_bundling_summary %>%
  dplyr::filter(
    !is.na(degree_class)
  ) %>%
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
    mean_observed_overlap =
      mean(
        observed_overlap,
        na.rm = TRUE
      ),
    mean_null_overlap =
      mean(
        null_mean_overlap,
        na.rm = TRUE
      ),
    number_species =
      dplyr::n(),
    .groups = "drop"
  )

bootstrap_equal_dataset <- function(
    values,
    n_boot = 5000
){
  
  values <- values[
    is.finite(values)
  ]
  
  n <- length(values)
  
  if(n == 0){
    
    return(
      tibble::tibble(
        mean_value = NA_real_,
        q025 = NA_real_,
        q975 = NA_real_,
        number_datasets = 0L
      )
    )
  }
  
  if(n == 1){
    
    return(
      tibble::tibble(
        mean_value = mean(values),
        q025 = NA_real_,
        q975 = NA_real_,
        number_datasets = 1L
      )
    )
  }
  
  boot_means <- replicate(
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
    mean_value = mean(values),
    q025 = safe_quantile(
      boot_means,
      0.025
    ),
    q975 = safe_quantile(
      boot_means,
      0.975
    ),
    number_datasets = n
  )
}

equal_class_delta <- dataset_class_means %>%
  dplyr::group_by(
    guild,
    degree_class
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      bootstrap_equal_dataset(
        .x$mean_bundling_delta,
        n_boot = n_boot
      )
    }
  ) %>%
  dplyr::ungroup()

bundling_class_plot <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    data = dataset_class_means,
    ggplot2::aes(
      x = degree_class,
      y = mean_bundling_delta,
      group = dataset
    ),
    colour = "grey65",
    alpha = 0.45,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = dataset_class_means,
    ggplot2::aes(
      x = degree_class,
      y = mean_bundling_delta
    ),
    colour = "grey55",
    alpha = 0.7,
    size = 1.5
  ) +
  ggplot2::geom_errorbar(
    data = equal_class_delta,
    ggplot2::aes(
      x = degree_class,
      ymin = q025,
      ymax = q975
    ),
    width = 0.12,
    colour = "black",
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = equal_class_delta,
    ggplot2::aes(
      x = degree_class,
      y = mean_value
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
  theme_pub() +
  ggplot2::labs(
    x = "Initial regional degree class",
    y = "Observed minus null mean Jaccard overlap",
    title =
      "Spatial bundling by regional degree class",
    subtitle = paste0(
      "Black points and intervals are equal-dataset means ",
      "and dataset-bootstrap 95% intervals"
    )
  )

dataset_overlap_long <- dplyr::bind_rows(
  
  dataset_class_means %>%
    dplyr::transmute(
      dataset,
      guild,
      degree_class,
      source = "Observed overlap",
      value =
        mean_observed_overlap
    ),
  
  dataset_class_means %>%
    dplyr::transmute(
      dataset,
      guild,
      degree_class,
      source = "Null mean overlap",
      value =
        mean_null_overlap
    )
)

equal_overlap_class <- dataset_overlap_long %>%
  dplyr::group_by(
    guild,
    degree_class,
    source
  ) %>%
  dplyr::group_modify(
    function(.x, .y){
      
      bootstrap_equal_dataset(
        .x$value,
        n_boot = n_boot
      )
    }
  ) %>%
  dplyr::ungroup()

overlap_class_plot <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = dataset_overlap_long,
    ggplot2::aes(
      x = degree_class,
      y = value,
      group = interaction(
        dataset,
        source
      ),
      linetype = source
    ),
    colour = "grey65",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbar(
    data = equal_overlap_class,
    ggplot2::aes(
      x = degree_class,
      ymin = q025,
      ymax = q975,
      linetype = source
    ),
    position = ggplot2::position_dodge(
      width = 0.22
    ),
    width = 0.08,
    colour = "black",
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = equal_overlap_class,
    ggplot2::aes(
      x = degree_class,
      y = mean_value,
      shape = source
    ),
    position = ggplot2::position_dodge(
      width = 0.22
    ),
    colour = "black",
    size = 2.7,
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
  theme_pub() +
  ggplot2::labs(
    x = "Initial regional degree class",
    y = "Mean pairwise Jaccard overlap",
    linetype = NULL,
    shape = NULL,
    title =
      "Observed and null portfolio overlap by regional degree",
    subtitle = paste0(
      "This separates changes in observed overlap from changes ",
      "in the support-preserving null expectation"
    )
  )

save_png(
  bundling_degree_plot,
  "40_bundling_delta_vs_degree.png",
  width = 10,
  height = 5.5
)

save_png(
  bundling_class_plot,
  "40_bundling_delta_by_degree_class.png",
  width = 9,
  height = 5.5
)

save_png(
  overlap_class_plot,
  "40_observed_vs_null_overlap_by_degree_class.png",
  width = 10,
  height = 5.5
)

## ------------------------------------------------------------
## Part 5: real-data conceptual mini-schematic
## ------------------------------------------------------------

clean_candidates <- portfolio_bundling_summary %>%
  dplyr::filter(
    initial_degree >= 3,
    initial_degree <= 8,
    is.finite(bundling_delta)
  )

if(nrow(clean_candidates) == 0){
  
  clean_candidates <- portfolio_bundling_summary %>%
    dplyr::filter(
      initial_degree >= 2,
      initial_degree <= 10,
      is.finite(bundling_delta)
    )
}

if(nrow(clean_candidates) >= 2){
  
  separated_example <- clean_candidates %>%
    dplyr::slice_min(
      bundling_delta,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(
      portfolio_type =
        "spatially_separated_portfolio"
    )
  
  bundled_example <- clean_candidates %>%
    dplyr::filter(
      !(
        dataset ==
          separated_example$dataset &
          guild ==
          separated_example$guild &
          species ==
          separated_example$species
      )
    ) %>%
    dplyr::slice_max(
      bundling_delta,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(
      portfolio_type =
        "spatially_bundled_portfolio"
    )
  
  if(nrow(bundled_example) == 1){
    
    separated_plot <- make_portfolio_matrix_plot(
      prepared =
        prepared_datasets[[
          as.character(
            separated_example$dataset
          )
        ]
        ],
      selection_row =
        separated_example
    ) +
      ggplot2::labs(
        title =
          "Separated interaction portfolio"
      )
    
    bundled_plot <- make_portfolio_matrix_plot(
      prepared =
        prepared_datasets[[
          as.character(
            bundled_example$dataset
          )
        ]
        ],
      selection_row =
        bundled_example
    ) +
      ggplot2::labs(
        title =
          "Bundled interaction portfolio"
      )
    
    schematic_grob <- gridExtra::arrangeGrob(
      separated_plot,
      bundled_plot,
      ncol = 2,
      top = grid::textGrob(
        paste0(
          "Real examples of separated and bundled ",
          "site × partner incidence structure"
        ),
        gp = grid::gpar(
          fontface = "bold",
          fontsize = 13
        )
      )
    )
    
    save_grob_png(
      grob = schematic_grob,
      filename =
        "40_real_example_separated_vs_bundled_portfolios.png",
      width = 11,
      height = 6
    )
  }
}

## ------------------------------------------------------------
## Validation checks
## ------------------------------------------------------------

validation_checks <- tibble::tibble(
  check = c(
    "Every retained regional link has full_K greater than zero",
    "Every retained regional link satisfies full_K no greater than full_n",
    "Every incidence-matrix column sum equals full_K",
    "Observed interaction sites lie inside pair-specific co-occurrence sets",
    "Every analysed focal species has at least two realised partners",
    "All observed and null mean Jaccard values lie between zero and one",
    "Every portfolio null preserves exact interaction support K",
    "Every portfolio null uses only pair-specific co-occurrence sites",
    "Consumers and Resources are analysed separately",
    "Low-support sensitivity includes both K less than or equal to two and K equal to one"
  ),
  pass = c(
    all(
      link_level_data$full_K > 0
    ),
    
    all(
      link_level_data$full_K <=
        link_level_data$full_n
    ),
    
    all(
      vapply(
        prepared_datasets,
        function(prepared){
          all(
            colSums(
              prepared$incidence_matrix
            ) ==
              prepared$realised_links$full_K
          )
        },
        logical(1)
      )
    ),
    
    all(
      vapply(
        prepared_datasets,
        function(prepared){
          
          all(
            unlist(
              Map(
                function(
    interaction_sites,
    possible_sites
                ){
                  all(
                    interaction_sites %in%
                      possible_sites
                  )
                },
    prepared$realised_links$
      interaction_sites,
    prepared$realised_links$
      cooccurrence_sites
              ),
    use.names = FALSE
            )
          )
        },
    logical(1)
      )
    ),
    
    all(
      portfolio_bundling_summary$
        initial_degree >= 2
    ),
    
    all(
      portfolio_bundling_summary$
        observed_overlap >= 0 &
        portfolio_bundling_summary$
        observed_overlap <= 1 &
        portfolio_bundling_summary$
        null_mean_overlap >= 0 &
        portfolio_bundling_summary$
        null_mean_overlap <= 1
    ),
    
    all(
      portfolio_null_checks$
        null_exact_support
    ),
    
    all(
      portfolio_null_checks$
        null_within_cooccurrence
    ),
    
    all(
      guild_levels %in%
        as.character(
          unique(
            portfolio_bundling_summary$
              guild
          )
        )
    ),
    
    all(
      c(
        "K ≤ 2",
        "K = 1"
      ) %in%
        unique(
          site_threshold_long$
            threshold
        )
    )
  )
)

write.csv2(
  validation_checks,
  file.path(
    out_dir,
    "40_incidence_network_validation_checks.csv"
  ),
  row.names = FALSE
)

## ------------------------------------------------------------
## Console summary
## ------------------------------------------------------------

saved_files <- sort(
  list.files(
    out_dir,
    full.names = FALSE
  )
)

message(
  "\nSaved Script 40 outputs in:\n",
  out_dir
)

message(
  "\nFiles created:"
)

for(filename in saved_files){
  message("  - ", filename)
}

message(
  "\nConceptual interpretation:"
)

message(
  "The site–interaction incidence matrix makes explicit why local support\n",
  "and regional binary links respond differently to site removal. Removing\n",
  "a site deletes one matrix row and immediately reduces the support of all\n",
  "interaction columns represented in that row. A regional interaction is\n",
  "lost only when its entire column becomes empty.\n\n",
  "Within a species portfolio, separated interaction columns distribute\n",
  "partners among different sites, whereas bundled columns place several\n",
  "interactions in the same sites. This column-overlap structure determines\n",
  "whether partner losses are spatially spread out or correlated."
)

message(
  "\nValidation checks:"
)

print(
  validation_checks,
  n = Inf
)