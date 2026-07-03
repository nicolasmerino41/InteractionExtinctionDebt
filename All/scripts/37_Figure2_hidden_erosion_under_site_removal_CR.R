## ------------------------------------------------------------
## Script: All/scripts/37_Figure2_hidden_erosion_under_site_removal_CR.R
## Figure 2 CR: Site removal erodes local interaction support before regional links disappear
## Consumers and resources analysed separately as focal guilds
## Outputs PNG only: panels and figures separately for consumers and resources
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "ggplot2", "tibble", "scales", "purrr", "gridExtra", "grid")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/37_Figure2_hidden_erosion_under_site_removal_CR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, name, width = 10, height = 6){
  ggplot2::ggsave(
    filename = file.path(out_dir, name),
    plot = p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

save_png_grid <- function(g, name, width = 9, height = 10){
  grDevices::png(
    filename = file.path(out_dir, name),
    width = width,
    height = height,
    units = "in",
    res = 320,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(g)
  grDevices::dev.off()
}

replace_na_int <- function(x, value = 0L){
  x[is.na(x)] <- value
  x
}

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500

set.seed(123)

cols <- c(
  "Regional interaction links still present" = "#D55E00",
  "Sites still supporting interactions" = "#0072B2",
  "Interaction still observed" = "#009E73",
  "Species still recorded together, interaction not observed" = "#CC79A7",
  "Species no longer recorded together" = "#999999"
)

state_levels <- c(
  "Species no longer recorded together",
  "Species still recorded together, interaction not observed",
  "Interaction still observed"
)

guild_levels <- c("Consumer", "Resource")

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

make_full_grid <- function(subset_index, full_links){
  
  subset_index$.tmp_key <- 1L
  full_links$.tmp_key <- 1L
  
  out <- dplyr::left_join(
    subset_index,
    full_links,
    by = ".tmp_key"
  ) %>%
    dplyr::select(-.tmp_key)
  
  subset_index$.tmp_key <- NULL
  full_links$.tmp_key <- NULL
  
  out
}

make_state_grid <- function(subset_index, focal_nodes){
  
  state_tbl <- tibble::tibble(
    state = factor(state_levels, levels = state_levels),
    .tmp_key = 1L
  )
  
  subset_index$.tmp_key <- 1L
  focal_nodes$.tmp_key <- 1L
  
  base <- dplyr::left_join(
    subset_index,
    focal_nodes,
    by = ".tmp_key"
  )
  
  out <- dplyr::left_join(
    base,
    state_tbl,
    by = ".tmp_key"
  ) %>%
    dplyr::select(-.tmp_key)
  
  out
}

make_guild_links <- function(link_reps){
  
  dplyr::bind_rows(
    link_reps %>%
      dplyr::transmute(
        guild = "Consumer",
        focal_node = consumer,
        partner_node = resource,
        removal_fraction,
        replicate,
        full_n,
        full_K,
        retained_n,
        retained_K
      ),
    link_reps %>%
      dplyr::transmute(
        guild = "Resource",
        focal_node = resource,
        partner_node = consumer,
        removal_fraction,
        replicate,
        full_n,
        full_K,
        retained_n,
        retained_K
      )
  ) %>%
    dplyr::mutate(
      guild = factor(guild, levels = guild_levels)
    )
}

run_one_dataset <- function(dataset){
  
  message("Figure 2 CR site removal: ", dataset)
  
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
  
  full_pairs <- cooc %>%
    dplyr::count(consumer, resource, name = "full_n") %>%
    dplyr::left_join(
      ints %>%
        dplyr::count(consumer, resource, name = "full_K"),
      by = c("consumer", "resource")
    ) %>%
    dplyr::mutate(
      full_K = replace_na_int(full_K, 0L),
      full_B = full_K > 0
    )
  
  if(any(full_pairs$full_K > full_pairs$full_n)){
    stop("Invalid full_K > full_n in ", dataset)
  }
  
  full_links <- full_pairs %>%
    dplyr::filter(full_K >= 1) %>%
    dplyr::select(consumer, resource, full_n, full_K)
  
  if(nrow(full_links) == 0 || sum(full_links$full_K) == 0){
    warning("No regional interaction links detected in ", dataset)
    return(list(
      retention = tibble::tibble(),
      states = tibble::tibble(),
      focal_totals = tibble::tibble()
    ))
  }
  
  focal_totals <- dplyr::bind_rows(
    full_links %>%
      dplyr::group_by(focal_node = consumer) %>%
      dplyr::summarise(
        guild = "Consumer",
        full_support = sum(full_K, na.rm = TRUE),
        initial_links = dplyr::n_distinct(resource),
        .groups = "drop"
      ),
    full_links %>%
      dplyr::group_by(focal_node = resource) %>%
      dplyr::summarise(
        guild = "Resource",
        full_support = sum(full_K, na.rm = TRUE),
        initial_links = dplyr::n_distinct(consumer),
        .groups = "drop"
      )
  ) %>%
    dplyr::mutate(
      dataset = dataset,
      guild = factor(guild, levels = guild_levels)
    ) %>%
    dplyr::select(dataset, guild, focal_node, full_support, initial_links)
  
  subset_index <- subsets %>%
    dplyr::distinct(removal_fraction, replicate)
  
  full_link_grid <- make_full_grid(
    subset_index,
    full_links
  )
  
  retained_n <- subsets %>%
    dplyr::left_join(
      cooc %>%
        dplyr::semi_join(full_links, by = c("consumer", "resource")) %>%
        dplyr::mutate(cooc_here = 1L),
      by = "site"
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
      by = "site"
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
      retained_n = replace_na_int(retained_n, 0L),
      retained_K = replace_na_int(retained_K, 0L)
    )
  
  guild_link_reps <- make_guild_links(link_reps)
  
  focal_rep <- guild_link_reps %>%
    dplyr::group_by(guild, focal_node, removal_fraction, replicate) %>%
    dplyr::summarise(
      retained_support = sum(retained_K, na.rm = TRUE),
      retained_links = sum(retained_K > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      focal_totals,
      by = c("guild", "focal_node")
    ) %>%
    dplyr::mutate(
      dataset = dataset,
      support_retained = retained_support / full_support,
      link_retained = retained_links / initial_links
    )
  
  retention <- focal_rep %>%
    dplyr::group_by(dataset, guild, removal_fraction, replicate) %>%
    dplyr::summarise(
      support_retained = mean(support_retained, na.rm = TRUE),
      link_retained = mean(link_retained, na.rm = TRUE),
      number_focal_nodes = dplyr::n_distinct(focal_node),
      .groups = "drop"
    )
  
  state_counts <- guild_link_reps %>%
    dplyr::mutate(
      state = dplyr::case_when(
        retained_K > 0 ~ "Interaction still observed",
        retained_n > 0 ~ "Species still recorded together, interaction not observed",
        TRUE ~ "Species no longer recorded together"
      ),
      state = factor(state, levels = state_levels)
    ) %>%
    dplyr::count(guild, focal_node, removal_fraction, replicate, state, name = "n")
  
  focal_grid <- focal_totals %>%
    dplyr::select(guild, focal_node, initial_links)
  
  states <- make_state_grid(subset_index, focal_grid) %>%
    dplyr::left_join(
      state_counts,
      by = c("guild", "focal_node", "removal_fraction", "replicate", "state")
    ) %>%
    dplyr::mutate(
      n = replace_na_int(n, 0L),
      dataset = dataset,
      fraction_node = n / initial_links
    ) %>%
    dplyr::group_by(dataset, guild, removal_fraction, replicate, state) %>%
    dplyr::summarise(
      fraction = mean(fraction_node, na.rm = TRUE),
      number_focal_nodes = dplyr::n_distinct(focal_node),
      .groups = "drop"
    )
  
  list(
    retention = retention,
    states = states,
    focal_totals = focal_totals
  )
}

outs <- lapply(all_dataset_names, run_one_dataset)

retention <- dplyr::bind_rows(lapply(outs, `[[`, "retention"))
states <- dplyr::bind_rows(lapply(outs, `[[`, "states"))
focal_totals <- dplyr::bind_rows(lapply(outs, `[[`, "focal_totals"))

retention <- retention %>%
  dplyr::mutate(guild = factor(guild, levels = guild_levels))

states <- states %>%
  dplyr::mutate(
    guild = factor(guild, levels = guild_levels),
    state = factor(as.character(state), levels = state_levels)
  )

write.csv2(
  retention,
  file.path(out_dir, "37_CR_retention_support_and_binary_links.csv"),
  row.names = FALSE
)

write.csv2(
  states,
  file.path(out_dir, "37_CR_original_links_state_under_site_removal.csv"),
  row.names = FALSE
)

write.csv2(
  focal_totals,
  file.path(out_dir, "37_CR_focal_node_initial_link_support_totals.csv"),
  row.names = FALSE
)

make_ret_long <- function(dat){
  
  dplyr::bind_rows(
    dat %>%
      dplyr::transmute(
        dataset,
        guild,
        removal_fraction,
        replicate,
        layer = "Sites still supporting interactions",
        fraction = support_retained
      ),
    dat %>%
      dplyr::transmute(
        dataset,
        guild,
        removal_fraction,
        replicate,
        layer = "Regional interaction links still present",
        fraction = link_retained
      )
  ) %>%
    dplyr::mutate(
      layer = factor(
        layer,
        levels = c(
          "Regional interaction links still present",
          "Sites still supporting interactions"
        )
      )
    )
}

plot_one_guild <- function(guild_name){
  
  guild_label <- ifelse(guild_name == "Consumer", "Consumers", "Resources")
  guild_file <- ifelse(guild_name == "Consumer", "consumers", "resources")
  
  ret_guild <- retention %>%
    dplyr::filter(guild == guild_name)
  
  states_guild <- states %>%
    dplyr::filter(guild == guild_name)
  
  ret_long <- make_ret_long(ret_guild)
  
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
      q25 = stats::quantile(fraction, 0.25, na.rm = TRUE),
      q75 = stats::quantile(fraction, 0.75, na.rm = TRUE),
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
      y = "Mean fraction remaining per focal species",
      colour = NULL,
      fill = NULL,
      title = paste0("A. ", guild_label, ": regional links remain while local support is eroded")
    )
  
  state_mean <- states_guild %>%
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
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_pub() +
    ggplot2::labs(
      x = "Sites removed",
      y = "Mean share of original focal links",
      fill = NULL,
      title = paste0("B. ", guild_label, ": co-occurrence can persist after observed interaction disappears")
    )
  
  title_grob <- grid::textGrob(
    paste0(
      guild_label,
      ": site removal erodes local interaction support before regional links disappear"
    ),
    gp = grid::gpar(fontface = "bold", fontsize = 13)
  )
  
  fig_grid <- gridExtra::arrangeGrob(
    title_grob,
    panel_a,
    panel_b,
    ncol = 1,
    heights = c(0.08, 1, 1)
  )
  
  save_png(
    panel_a,
    paste0("37_Figure2A_", guild_file, "_support_vs_binary_retention.png"),
    8,
    5.2
  )
  
  save_png(
    panel_b,
    paste0("37_Figure2B_", guild_file, "_cooccurrence_interaction_decoupling.png"),
    8,
    5.2
  )
  
  save_png_grid(
    fig_grid,
    paste0("37_Figure2_", guild_file, "_hidden_erosion_under_site_removal.png"),
    9,
    10
  )
  
  invisible(list(panel_a = panel_a, panel_b = panel_b, fig = fig_grid))
}

plot_one_guild("Consumer")
plot_one_guild("Resource")

checks <- tibble::tibble(
  check = c(
    "Consumer support retention equals 1 at zero removal",
    "Resource support retention equals 1 at zero removal",
    "Consumer binary link retention equals 1 at zero removal",
    "Resource binary link retention equals 1 at zero removal",
    "State fractions sum to 1 within dataset-guild-replicate-removal"
  ),
  pass = c(
    all(abs(retention$support_retained[retention$guild == "Consumer" & retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(abs(retention$support_retained[retention$guild == "Resource" & retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(abs(retention$link_retained[retention$guild == "Consumer" & retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    all(abs(retention$link_retained[retention$guild == "Resource" & retention$removal_fraction == 0] - 1) < 1e-8, na.rm = TRUE),
    states %>%
      dplyr::group_by(dataset, guild, removal_fraction, replicate) %>%
      dplyr::summarise(s = sum(fraction), .groups = "drop") %>%
      dplyr::pull(s) %>%
      { all(abs(. - 1) < 1e-8, na.rm = TRUE) }
  )
)

write.csv2(
  checks,
  file.path(out_dir, "37_Figure2_CR_checks.csv"),
  row.names = FALSE
)

message("Saved Figure 2 CR PNG outputs in: ", out_dir)
print(checks)