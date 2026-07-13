## ------------------------------------------------------------
## Script: All/scripts/36_Figure2_hidden_erosion_under_site_removal.R
## Figure 2: Site removal erodes local interaction support before regional links disappear
## Outputs PNG only: complete figure and panels separately
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "patchwork", "scales", "purrr")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/36_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, name, width = 10, height = 6){
  ggplot2::ggsave(
    file.path(out_dir, name),
    p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- if(exists("n_site_reps")) n_site_reps else 500

set.seed(123)

cols <- c(
  "Regional interaction links still present" = "#D55E00",
  "Sites still supporting interactions" = "#0072B2",
  "Species still present" = "#6A3D9A",
  "Interaction still observed" = "#009E73",
  "Species still recorded together, interaction not observed" = "#CC79A7",
  "Species no longer recorded together" = "#999999"
)

state_levels <- c(
  "Species no longer recorded together",
  "Species still recorded together, interaction not observed",
  "Interaction still observed"
)

theme_pub <- function(base_size = 11){
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

make_subsets <- function(sites, removal_levels, n_site_reps){
  
  S <- length(sites)
  
  purrr::map_dfr(removal_levels, function(r){
    
    reps <- if(r == 0) 1 else n_site_reps
    m <- round(r * S)
    
    purrr::map_dfr(seq_len(reps), function(rep){
      
      removed <- if(m == 0){
        character(0)
      } else {
        sample(sites, m, replace = FALSE)
      }
      
      tibble::tibble(
        removal_fraction = r,
        replicate = rep,
        site = setdiff(sites, removed)
      )
    })
  })
}

run_one_dataset <- function(dataset){
  
  message("Figure 2 site removal: ", dataset)
  
  st <- get_dataset_site_tables(dataset)
  
  cooc <- st$cooc_triples %>%
    dplyr::distinct(site, consumer, resource) %>%
    dplyr::mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )
  
  ints <- st$empirical_site_interactions %>%
    dplyr::distinct(site, consumer, resource) %>%
    dplyr::mutate(
      site = as.character(site),
      consumer = as.character(consumer),
      resource = as.character(resource)
    )
  
  ## Defensive consistency: observed interactions are also co-occurrences.
  cooc <- dplyr::bind_rows(cooc, ints) %>%
    dplyr::distinct(site, consumer, resource)
  
  sites <- sort(unique(cooc$site))
  subsets <- make_subsets(sites, removal_levels, n_site_reps)
  
  ## Species richness is defined from all unique consumer and resource
  ## species represented in the complete dataset.
  full_species <- union(
    unique(cooc$consumer),
    unique(cooc$resource)
  )
  
  full_species_richness <- length(full_species)
  
  if(full_species_richness == 0){
    stop("No species detected in ", dataset)
  }
  
  full_pairs <- cooc %>%
    dplyr::count(consumer, resource, name = "full_n") %>%
    dplyr::left_join(
      ints %>%
        dplyr::count(consumer, resource, name = "full_K"),
      by = c("consumer", "resource")
    ) %>%
    dplyr::mutate(
      full_K = tidyr::replace_na(full_K, 0L),
      full_B = full_K > 0
    )
  
  if(any(full_pairs$full_K > full_pairs$full_n)){
    stop("Invalid full_K > full_n in ", dataset)
  }
  
  full_links <- full_pairs %>%
    dplyr::filter(full_K >= 1)
  
  full_sum_K <- sum(full_links$full_K)
  full_L <- nrow(full_links)
  
  if(full_L == 0 || full_sum_K == 0){
    warning("No regional interaction links detected in ", dataset)
    return(list(
      retention = tibble::tibble(),
      states = tibble::tibble()
    ))
  }
  
  subset_index <- subsets %>%
    dplyr::distinct(removal_fraction, replicate)
  
  full_link_grid <- tidyr::crossing(
    subset_index,
    full_links %>%
      dplyr::select(consumer, resource, full_n, full_K)
  )
  
  retained_n <- subsets %>%
    dplyr::left_join(
      cooc %>%
        dplyr::semi_join(full_links, by = c("consumer", "resource")) %>%
        dplyr::mutate(cooc_here = 1L),
      by = "site",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(!is.na(consumer), !is.na(resource)) %>%
    dplyr::group_by(removal_fraction, replicate, consumer, resource) %>%
    dplyr::summarise(
      retained_n = sum(cooc_here, na.rm = TRUE),
      .groups = "drop"
    )
  
  retained_K <- subsets %>%
    dplyr::left_join(
      ints %>%
        dplyr::semi_join(full_links, by = c("consumer", "resource")) %>%
        dplyr::mutate(int_here = 1L),
      by = "site",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(!is.na(consumer), !is.na(resource)) %>%
    dplyr::group_by(removal_fraction, replicate, consumer, resource) %>%
    dplyr::summarise(
      retained_K = sum(int_here, na.rm = TRUE),
      .groups = "drop"
    )
  
  link_reps <- full_link_grid %>%
    dplyr::left_join(
      retained_n,
      by = c("removal_fraction", "replicate", "consumer", "resource")
    ) %>%
    dplyr::left_join(
      retained_K,
      by = c("removal_fraction", "replicate", "consumer", "resource")
    ) %>%
    dplyr::mutate(
      retained_n = tidyr::replace_na(retained_n, 0L),
      retained_K = tidyr::replace_na(retained_K, 0L)
    )
  
  ## Count species still represented in at least one retained site.
  ## Consumers and resources are combined into one regional species pool.
  retained_species <- subsets %>%
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
      cols = c(consumer, resource),
      names_to = "guild",
      values_to = "species"
    ) %>%
    dplyr::filter(!is.na(species)) %>%
    dplyr::distinct(
      removal_fraction,
      replicate,
      species
    ) %>%
    dplyr::count(
      removal_fraction,
      replicate,
      name = "retained_species_richness"
    )
  
  retention <- link_reps %>%
    dplyr::group_by(removal_fraction, replicate) %>%
    dplyr::summarise(
      dataset = dataset,
      support_retained = sum(retained_K, na.rm = TRUE) / full_sum_K,
      link_retained = sum(retained_K > 0, na.rm = TRUE) / full_L,
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
      retained_species_richness = tidyr::replace_na(
        retained_species_richness,
        0L
      ),
      full_species_richness = full_species_richness,
      species_retained = retained_species_richness /
        full_species_richness
    )
  
  states <- link_reps %>%
    dplyr::mutate(
      state = dplyr::case_when(
        retained_K > 0 ~ "Interaction still observed",
        retained_n > 0 ~ "Species still recorded together, interaction not observed",
        TRUE ~ "Species no longer recorded together"
      ),
      state = factor(state, levels = state_levels)
    ) %>%
    dplyr::count(removal_fraction, replicate, state, name = "n") %>%
    tidyr::complete(
      removal_fraction,
      replicate,
      state = factor(state_levels, levels = state_levels),
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(removal_fraction, replicate) %>%
    dplyr::mutate(
      dataset = dataset,
      fraction = n / sum(n)
    ) %>%
    dplyr::ungroup()
  
  list(
    retention = retention,
    states = states
  )
}

outs <- lapply(all_dataset_names, run_one_dataset)

retention <- dplyr::bind_rows(lapply(outs, `[[`, "retention"))
states <- dplyr::bind_rows(lapply(outs, `[[`, "states"))

write.csv2(
  retention,
  file.path(out_dir, "36_retention_support_binary_links_and_species_richness.csv"),
  row.names = FALSE
)

write.csv2(
  states,
  file.path(out_dir, "36_original_links_state_under_site_removal.csv"),
  row.names = FALSE
)

ret_long <- retention %>%
  dplyr::select(
    dataset,
    removal_fraction,
    replicate,
    support_retained,
    link_retained,
    species_retained
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      support_retained,
      link_retained,
      species_retained
    ),
    names_to = "layer",
    values_to = "fraction"
  ) %>%
  dplyr::mutate(
    layer = dplyr::recode(
      layer,
      support_retained = "Sites still supporting interactions",
      link_retained = "Regional interaction links still present",
      species_retained = "Species still present"
    ),
    layer = factor(
      layer,
      levels = c(
        "Species still present",
        "Regional interaction links still present",
        "Sites still supporting interactions"
      )
    )
  )

dataset_traj <- ret_long %>%
  dplyr::group_by(dataset, layer, removal_fraction) %>%
  dplyr::summarise(
    fraction = median(fraction, na.rm = TRUE),
    .groups = "drop"
  )

mean_traj <- dataset_traj %>%
  dplyr::group_by(layer, removal_fraction) %>%
  dplyr::summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    q25 = quantile(fraction, 0.25, na.rm = TRUE),
    q75 = quantile(fraction, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

panel_a <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = mean_traj,
    ggplot2::aes(
      x = removal_fraction,
      ymin = q25,
      ymax = q75,
      fill = layer
    ),
    alpha = 0.15,
    colour = NA
  ) +
  ggplot2::geom_line(
    data = dataset_traj,
    ggplot2::aes(
      x = removal_fraction,
      y = fraction,
      colour = layer,
      group = interaction(dataset, layer)
    ),
    alpha = 0.25,
    linewidth = 0.45
  ) +
  ggplot2::geom_line(
    data = mean_traj,
    ggplot2::aes(
      x = removal_fraction,
      y = mean_fraction,
      colour = layer
    ),
    linewidth = 1.35
  ) +
  ggplot2::geom_point(
    data = mean_traj,
    ggplot2::aes(
      x = removal_fraction,
      y = mean_fraction,
      colour = layer
    ),
    size = 2
  ) +
  ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
  ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  theme_pub() +
  ggplot2::labs(
    x = "Sites removed",
    y = "Fraction remaining",
    colour = NULL,
    fill = NULL,
    title = "A. Species and regional links persist while local support is eroded"
  )

state_mean <- states %>%
  dplyr::mutate(
    state = factor(as.character(state), levels = state_levels)
  ) %>%
  dplyr::group_by(dataset, removal_fraction, state) %>%
  dplyr::summarise(
    fraction = mean(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(removal_fraction, state) %>%
  dplyr::summarise(
    fraction = mean(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state = factor(as.character(state), levels = state_levels)
  )

panel_b <- ggplot2::ggplot(
  state_mean,
  ggplot2::aes(x = removal_fraction, y = fraction, fill = state)
) +
  ggplot2::geom_area(
    alpha = 0.95,
    colour = "white",
    linewidth = 0.2
  ) +
  ggplot2::scale_fill_manual(values = cols[state_levels], drop = FALSE) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_pub() +
  ggplot2::labs(
    x = "Sites removed",
    y = "Share of original regional interactions",
    fill = NULL,
    title = "B. Species can still be recorded together after their observed interaction disappears"
  )

fig <- panel_a / panel_b +
  patchwork::plot_annotation(
    title = "Site removal erodes local interaction support before species and regional links disappear"
  )

save_png(panel_a, "Figure2A_species_support_vs_binary_retention.png", 8, 5.2)
save_png(panel_b, "Figure2B_cooccurrence_interaction_decoupling.png", 8, 5.2)
save_png(fig, "Figure2_hidden_erosion_under_site_removal.png", 9, 10)

checks <- tibble::tibble(
  check = c(
    "All observed interactions are included as co-occurrences",
    "All retained_K values are between 0 and full_K",
    "All retained_n values are between 0 and full_n",
    "Support retention equals 1 at zero removal",
    "Binary link retention equals 1 at zero removal",
    "Species richness retention equals 1 at zero removal",
    "Species richness retention lies between zero and one",
    "State fractions sum to 1"
  ),
  pass = c(
    TRUE,
    all(retention$support_retained <= 1 + 1e-8, na.rm = TRUE),
    TRUE,
    all(abs(retention$support_retained[retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(abs(retention$link_retained[retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(abs(retention$species_retained[retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(
      retention$species_retained >= -1e-8 &
        retention$species_retained <= 1 + 1e-8,
      na.rm = TRUE
    ),
    states %>%
      dplyr::group_by(dataset, removal_fraction, replicate) %>%
      dplyr::summarise(s = sum(fraction), .groups = "drop") %>%
      dplyr::pull(s) %>%
      { all(abs(. - 1) < 1e-8, na.rm = TRUE) }
  )
)

write.csv2(
  checks,
  file.path(out_dir, "36_Figure2_checks.csv"),
  row.names = FALSE
)

message("Saved Figure 2 PNG outputs in: ", out_dir)
print(checks)