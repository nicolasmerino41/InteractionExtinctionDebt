## ------------------------------------------------------------
## Script: All/scripts/36_Figure4_generalist_filtered_portfolio_compressed.R
## Figure 4: Site loss filters communities toward initially generalist species while compressing their portfolios
## Outputs PNG only: complete figure and panels separately
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "ggplot2", "tibble", "scales", "purrr", "gridExtra", "grid")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/36_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, name, width = 10, height = 7){
  ggplot2::ggsave(
    filename = file.path(out_dir, name),
    plot = p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

save_png_grid <- function(g, name, width = 10, height = 14){
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

degree_cols <- c(
  "Lower initial degree" = "#4E79A7",
  "Middle initial degree" = "#F28E2B",
  "Higher initial degree" = "#59A14F"
)

degree_levels <- names(degree_cols)

theme_pub <- function(base_size = 10){
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

make_degree_groups <- function(df, value_col = "initial_degree"){
  
  vals <- sort(unique(df[[value_col]]))
  
  if(length(vals) == 1){
    return(
      df %>%
        dplyr::mutate(
          initial_degree_group = factor("Middle initial degree", levels = degree_levels)
        )
    )
  }
  
  ranks <- tibble::tibble(
    !!value_col := vals,
    r = rank(vals, ties.method = "average") / length(vals)
  ) %>%
    dplyr::mutate(
      initial_degree_group = dplyr::case_when(
        r <= 1/3 ~ "Lower initial degree",
        r <= 2/3 ~ "Middle initial degree",
        TRUE ~ "Higher initial degree"
      )
    ) %>%
    dplyr::select(-r)
  
  df %>%
    dplyr::left_join(ranks, by = value_col) %>%
    dplyr::mutate(
      initial_degree_group = factor(initial_degree_group, levels = degree_levels)
    )
}

make_subsets <- function(sites){
  
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

make_nodes <- function(links, guild_name){
  
  if(guild_name == "Consumer"){
    links %>%
      dplyr::transmute(
        guild = "Consumer",
        node = consumer,
        partner = resource,
        full_n,
        full_K
      )
  } else {
    links %>%
      dplyr::transmute(
        guild = "Resource",
        node = resource,
        partner = consumer,
        full_n,
        full_K
      )
  }
}

make_full_grid <- function(subset_index, full_links){
  
  subset_index$.tmp_key <- 1L
  full_links$.tmp_key <- 1L
  
  dplyr::left_join(
    subset_index,
    full_links,
    by = ".tmp_key"
  ) %>%
    dplyr::select(-.tmp_key)
}

run_one_dataset <- function(dataset){
  
  message("Figure 4: ", dataset)
  
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
  
  cooc <- dplyr::bind_rows(cooc, ints) %>%
    dplyr::distinct(site, consumer, resource)
  
  sites <- sort(unique(cooc$site))
  S <- length(sites)
  
  full_links <- cooc %>%
    dplyr::count(consumer, resource, name = "full_n") %>%
    dplyr::left_join(
      ints %>%
        dplyr::count(consumer, resource, name = "full_K"),
      by = c("consumer", "resource")
    ) %>%
    dplyr::mutate(
      full_K = replace_na_int(full_K, 0L)
    ) %>%
    dplyr::filter(full_K >= 1)
  
  if(nrow(full_links) == 0){
    warning("No regional links in ", dataset)
    return(list(
      panel_a = tibble::tibble(),
      disappearance = tibble::tibble(),
      compression = tibble::tibble()
    ))
  }
  
  nodes <- dplyr::bind_rows(
    make_nodes(full_links, "Consumer"),
    make_nodes(full_links, "Resource")
  ) %>%
    dplyr::group_by(guild, node) %>%
    dplyr::summarise(
      initial_degree = dplyr::n_distinct(partner),
      mean_n_std = mean(full_n / S, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(guild) %>%
    dplyr::group_modify(~ make_degree_groups(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = dataset)
  
  panel_a <- nodes %>%
    dplyr::group_by(dataset, guild, initial_degree_group) %>%
    dplyr::summarise(
      value = mean(mean_n_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  subsets <- make_subsets(sites)
  
  subset_index <- subsets %>%
    dplyr::distinct(removal_fraction, replicate)
  
  full_link_grid <- make_full_grid(
    subset_index,
    full_links %>%
      dplyr::select(consumer, resource, full_K)
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
  
  rep_links <- full_link_grid %>%
    dplyr::left_join(
      retained_K,
      by = c("removal_fraction", "replicate", "consumer", "resource")
    ) %>%
    dplyr::mutate(
      retained_K = replace_na_int(retained_K, 0L)
    )
  
  node_rep <- dplyr::bind_rows(
    rep_links %>%
      dplyr::transmute(
        guild = "Consumer",
        node = consumer,
        partner = resource,
        removal_fraction,
        replicate,
        retained = retained_K > 0
      ),
    rep_links %>%
      dplyr::transmute(
        guild = "Resource",
        node = resource,
        partner = consumer,
        removal_fraction,
        replicate,
        retained = retained_K > 0
      )
  ) %>%
    dplyr::group_by(guild, node, removal_fraction, replicate) %>%
    dplyr::summarise(
      retained_degree = sum(retained, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      nodes %>%
        dplyr::select(dataset, guild, node, initial_degree, initial_degree_group),
      by = c("guild", "node")
    ) %>%
    dplyr::mutate(
      network_active = retained_degree > 0,
      retained_fraction = retained_degree / initial_degree
    )
  
  disappearance <- node_rep %>%
    dplyr::group_by(dataset, guild, initial_degree_group, removal_fraction, replicate) %>%
    dplyr::summarise(
      value = mean(network_active, na.rm = TRUE),
      .groups = "drop"
    )
  
  compression <- node_rep %>%
    dplyr::filter(network_active) %>%
    dplyr::group_by(dataset, guild, initial_degree_group, removal_fraction, replicate) %>%
    dplyr::summarise(
      value = mean(retained_fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    panel_a = panel_a,
    disappearance = disappearance,
    compression = compression
  )
}

outs <- lapply(all_dataset_names, run_one_dataset)

panel_a_data <- dplyr::bind_rows(lapply(outs, `[[`, "panel_a"))
disappearance <- dplyr::bind_rows(lapply(outs, `[[`, "disappearance"))
compression <- dplyr::bind_rows(lapply(outs, `[[`, "compression"))

panel_a_overall <- panel_a_data %>%
  dplyr::group_by(guild, initial_degree_group) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

panel_a_plot <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = panel_a_data,
    ggplot2::aes(
      x = initial_degree_group,
      y = value,
      group = dataset,
      colour = initial_degree_group
    ),
    alpha = 0.23,
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    data = panel_a_data,
    ggplot2::aes(
      x = initial_degree_group,
      y = value,
      colour = initial_degree_group
    ),
    alpha = 0.23,
    size = 1.2
  ) +
  ggplot2::geom_point(
    data = panel_a_overall,
    ggplot2::aes(
      x = initial_degree_group,
      y = mean_value,
      colour = initial_degree_group
    ),
    size = 3
  ) +
  ggplot2::geom_errorbar(
    data = panel_a_overall,
    ggplot2::aes(
      x = initial_degree_group,
      ymin = q25,
      ymax = q75,
      colour = initial_degree_group
    ),
    width = 0.16,
    linewidth = 0.8
  ) +
  ggplot2::facet_wrap(~ guild, nrow = 1) +
  ggplot2::scale_colour_manual(values = degree_cols, drop = FALSE) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  ) +
  ggplot2::labs(
    x = "Initial-degree group",
    y = "Mean share of sites where each realised pair is recorded together",
    colour = NULL,
    title = "A. Initially generalist species have broader spatial opportunity per realised partner"
  )

summarise_curve <- function(dat){
  
  dataset_traj <- dat %>%
    dplyr::group_by(dataset, guild, initial_degree_group, removal_fraction) %>%
    dplyr::summarise(
      value = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  overall <- dataset_traj %>%
    dplyr::group_by(guild, initial_degree_group, removal_fraction) %>%
    dplyr::summarise(
      mean_value = mean(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(dataset = dataset_traj, overall = overall)
}

dis_s <- summarise_curve(disappearance)
comp_s <- summarise_curve(compression)

make_curve_plot <- function(ss, title, ylab){
  
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ss$overall,
      ggplot2::aes(
        x = removal_fraction,
        ymin = q25,
        ymax = q75,
        fill = initial_degree_group
      ),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = ss$dataset,
      ggplot2::aes(
        x = removal_fraction,
        y = value,
        group = interaction(dataset, initial_degree_group),
        colour = initial_degree_group
      ),
      alpha = 0.22,
      linewidth = 0.45
    ) +
    ggplot2::geom_line(
      data = ss$overall,
      ggplot2::aes(
        x = removal_fraction,
        y = mean_value,
        colour = initial_degree_group
      ),
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = ss$overall,
      ggplot2::aes(
        x = removal_fraction,
        y = mean_value,
        colour = initial_degree_group
      ),
      size = 1.7
    ) +
    ggplot2::facet_wrap(~ guild, nrow = 1) +
    ggplot2::scale_colour_manual(values = degree_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = degree_cols, drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_pub() +
    ggplot2::labs(
      x = "Sites removed",
      y = ylab,
      colour = NULL,
      fill = NULL,
      title = title
    )
}

panel_b_plot <- make_curve_plot(
  dis_s,
  "B. Selective disappearance",
  "Proportion of species remaining connected"
)

panel_c_plot <- make_curve_plot(
  comp_s,
  "C. Portfolio compression among connected species",
  "Retained partners / initial partners"
)

title_grob <- grid::textGrob(
  "Site loss filters communities toward initially generalist species while compressing their portfolios",
  gp = grid::gpar(fontface = "bold", fontsize = 13)
)

fig_grid <- gridExtra::arrangeGrob(
  title_grob,
  panel_a_plot,
  panel_b_plot,
  panel_c_plot,
  ncol = 1,
  heights = c(0.08, 1, 1, 1)
)

save_png(panel_a_plot, "Figure4A_broader_cooccurrence_opportunity.png", 9, 5)
save_png(panel_b_plot, "Figure4B_selective_disappearance.png", 9, 5)
save_png(panel_c_plot, "Figure4C_portfolio_compression.png", 9, 5)
save_png_grid(fig_grid, "Figure4_generalist_filtered_portfolio_compressed.png", 10, 14)

message("Saved Figure 4 PNG outputs in: ", out_dir)