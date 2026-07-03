## ------------------------------------------------------------
## Script: All/scripts/36_Figure3_repeated_support_beyond_homogeneous_opportunity.R
## Figure 3: Repeated local realisation explains persistence beyond homogeneous co-occurrence opportunity
## Baseline: homogeneous UNCONDITIONED per-site interaction probability
## Outputs PNG only: complete figure and panels separately
##
## Panel A update:
##   - x-axis uses discrete co-occurrence-support bins, not numeric midpoints.
##   - empirical values are orange dots.
##   - homogeneous-model values are purple dots with 95% simulation intervals.
##   - no connecting lines.
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

save_png <- function(p, name, width = 12, height = 8){
  ggplot2::ggsave(
    filename = file.path(out_dir, name),
    plot = p,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

save_png_grid <- function(g, name, width = 13, height = 13.5){
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

replace_na_num <- function(x, value = 0){
  x[is.na(x)] <- value
  x
}

theme_pub <- function(base_size = 9){
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 500
n_model_sims_panel_a <- 500
set.seed(123)

cols <- c(
  "Empirical data" = "#D55E00",
  "Homogeneous model" = "#6A3D9A"
)

n_bin_levels <- c("1", "2", "3–4", "5–8", "9–16", "17–32", "33+")

## ------------------------------------------------------------
## Unconditioned full-landscape calibration
## Fits p so that:
##   sum_ij [1 - (1 - p)^full_n_ij] = observed full regional links
## Then the same unconditioned formula is used after site removal:
##   expected retained links = sum_ij [1 - (1 - p)^retained_n_ij]
## ------------------------------------------------------------
fit_p_unconditioned <- function(pair_table){
  L_obs <- sum(pair_table$full_K > 0)
  n_vec <- pair_table$full_n
  
  expected_links <- function(p){
    sum(1 - (1 - p)^n_vec, na.rm = TRUE)
  }
  
  f <- function(p){
    expected_links(p) - L_obs
  }
  
  if(f(1e-12) >= 0) return(1e-12)
  if(f(1 - 1e-12) <= 0) return(1 - 1e-12)
  
  uniroot(f, c(1e-12, 1 - 1e-12), tol = 1e-12)$root
}

make_bins <- function(n){
  cut(
    n,
    breaks = c(0, 1, 2, 4, 8, 16, 32, Inf),
    labels = n_bin_levels,
    right = TRUE
  )
}

make_subsets <- function(sites){
  S <- length(sites)
  
  purrr::map_dfr(removal_levels, function(r){
    reps <- if(r == 0) 1 else n_site_reps
    m <- round(r * S)
    
    purrr::map_dfr(seq_len(reps), function(rep){
      removed <- if(m == 0) character(0) else sample(sites, m, replace = FALSE)
      tibble::tibble(
        removal_fraction = r,
        replicate = rep,
        site = setdiff(sites, removed)
      )
    })
  })
}

make_full_grid <- function(subset_index, pair_table){
  subset_index$.tmp_key <- 1L
  pair_table$.tmp_key <- 1L
  
  out <- dplyr::left_join(subset_index, pair_table, by = ".tmp_key") %>%
    dplyr::select(-.tmp_key)
  
  subset_index$.tmp_key <- NULL
  pair_table$.tmp_key <- NULL
  
  out
}

simulate_panel_a_support <- function(pair_table, dataset, p){
  dat <- pair_table %>%
    dplyr::mutate(
      n_bin = factor(make_bins(full_n), levels = n_bin_levels)
    )
  
  empirical <- dat %>%
    dplyr::group_by(dataset, n_bin) %>%
    dplyr::summarise(
      source = "Empirical data",
      mean_support = ifelse(
        sum(full_K > 0, na.rm = TRUE) >= 5,
        mean(full_K[full_K > 0], na.rm = TRUE),
        NA_real_
      ),
      q025 = NA_real_,
      q975 = NA_real_,
      empirical_n_links = sum(full_K > 0, na.rm = TRUE),
      n_pairs = dplyr::n(),
      .groups = "drop"
    )
  
  model_sims <- purrr::map_dfr(seq_len(n_model_sims_panel_a), function(sim_id){
    sim_K <- stats::rbinom(n = nrow(dat), size = dat$full_n, prob = p)
    
    dat %>%
      dplyr::mutate(sim_K = sim_K, sim_id = sim_id) %>%
      dplyr::group_by(dataset, n_bin, sim_id) %>%
      dplyr::summarise(
        model_mean_support = ifelse(
          sum(sim_K > 0, na.rm = TRUE) >= 5,
          mean(sim_K[sim_K > 0], na.rm = TRUE),
          NA_real_
        ),
        model_n_links = sum(sim_K > 0, na.rm = TRUE),
        n_pairs = dplyr::n(),
        .groups = "drop"
      )
  })
  
  model <- model_sims %>%
    dplyr::group_by(dataset, n_bin) %>%
    dplyr::summarise(
      source = "Homogeneous model",
      mean_support = mean(model_mean_support, na.rm = TRUE),
      q025 = stats::quantile(model_mean_support, 0.025, na.rm = TRUE),
      q975 = stats::quantile(model_mean_support, 0.975, na.rm = TRUE),
      empirical_n_links = NA_integer_,
      n_pairs = dplyr::first(n_pairs),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_support = ifelse(is.nan(mean_support), NA_real_, mean_support),
      q025 = ifelse(is.nan(q025), NA_real_, q025),
      q975 = ifelse(is.nan(q975), NA_real_, q975)
    )
  
  dplyr::bind_rows(empirical, model) %>%
    dplyr::mutate(
      n_bin = factor(n_bin, levels = n_bin_levels),
      source = factor(source, levels = c("Empirical data", "Homogeneous model"))
    )
}

run_one_dataset <- function(dataset){
  message("Figure 3 unconditioned baseline: ", dataset)
  
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
  
  ## Defensive: every observed interaction implies the pair is recorded together at that site.
  cooc <- dplyr::bind_rows(cooc, ints) %>%
    dplyr::distinct(site, consumer, resource)
  
  pair_table <- cooc %>%
    dplyr::count(consumer, resource, name = "full_n") %>%
    dplyr::left_join(
      ints %>% dplyr::count(consumer, resource, name = "full_K"),
      by = c("consumer", "resource")
    ) %>%
    dplyr::mutate(
      full_K = replace_na_num(full_K, 0L),
      dataset = dataset
    )
  
  if(any(pair_table$full_n < 1)) stop("Invalid full_n < 1 in ", dataset)
  if(any(pair_table$full_K > pair_table$full_n)) stop("Invalid full_K > full_n in ", dataset)
  
  p <- fit_p_unconditioned(pair_table)
  full_L <- sum(pair_table$full_K > 0)
  full_model_L <- sum(1 - (1 - p)^pair_table$full_n)
  
  message(
    "  p_unconditioned = ", signif(p, 5),
    "; empirical L0 = ", full_L,
    "; model expected L0 = ", round(full_model_L, 6),
    "; difference = ", signif(full_model_L - full_L, 5)
  )
  
  support_a <- simulate_panel_a_support(pair_table, dataset, p)
  
  sites <- sort(unique(cooc$site))
  subsets <- make_subsets(sites)
  
  ret_n <- subsets %>%
    dplyr::left_join(cooc %>% dplyr::mutate(cooc_here = 1L), by = "site") %>%
    dplyr::filter(!is.na(consumer), !is.na(resource)) %>%
    dplyr::group_by(removal_fraction, replicate, consumer, resource) %>%
    dplyr::summarise(retained_n = sum(cooc_here, na.rm = TRUE), .groups = "drop")
  
  ret_K <- subsets %>%
    dplyr::left_join(ints %>% dplyr::mutate(int_here = 1L), by = "site") %>%
    dplyr::filter(!is.na(consumer), !is.na(resource)) %>%
    dplyr::group_by(removal_fraction, replicate, consumer, resource) %>%
    dplyr::summarise(retained_K = sum(int_here, na.rm = TRUE), .groups = "drop")
  
  subset_index <- subsets %>% dplyr::distinct(removal_fraction, replicate)
  
  all_rep_pairs <- make_full_grid(
    subset_index,
    pair_table %>% dplyr::select(consumer, resource, full_n, full_K)
  )
  
  retention <- all_rep_pairs %>%
    dplyr::left_join(ret_n, by = c("removal_fraction", "replicate", "consumer", "resource")) %>%
    dplyr::left_join(ret_K, by = c("removal_fraction", "replicate", "consumer", "resource")) %>%
    dplyr::mutate(
      retained_n = replace_na_num(retained_n, 0L),
      retained_K = replace_na_num(retained_K, 0L),
      model_link_prob = 1 - (1 - p)^retained_n
    ) %>%
    dplyr::group_by(removal_fraction, replicate) %>%
    dplyr::summarise(
      dataset = dataset,
      empirical_links = sum(retained_K > 0, na.rm = TRUE),
      model_mean_links = sum(model_link_prob, na.rm = TRUE),
      model_sd_links = sqrt(sum(model_link_prob * (1 - model_link_prob), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      model_q025 = pmax(0, model_mean_links - 1.96 * model_sd_links),
      model_q975 = model_mean_links + 1.96 * model_sd_links,
      full_L = full_L,
      p_unconditioned = p,
      full_model_L = full_model_L
    )
  
  ## Check that the unconditioned full-landscape calibration matches at removal 0.
  zero_check <- retention %>% dplyr::filter(removal_fraction == 0)
  if(nrow(zero_check) != 1) stop("Unexpected number of zero-removal rows in ", dataset)
  if(abs(zero_check$model_mean_links - zero_check$empirical_links) > 1e-5){
    stop(
      "Unconditioned calibration failed at zero removal in ", dataset,
      ": empirical = ", zero_check$empirical_links,
      "; model = ", zero_check$model_mean_links
    )
  }
  
  list(support_a = support_a, retention = retention)
}

outs <- lapply(all_dataset_names, run_one_dataset)

support_a <- dplyr::bind_rows(lapply(outs, `[[`, "support_a"))
retention <- dplyr::bind_rows(lapply(outs, `[[`, "retention"))

position_panel_a <- ggplot2::position_dodge(width = 0.45)

panel_a <- ggplot2::ggplot(
  support_a,
  ggplot2::aes(x = n_bin, y = mean_support, colour = source, group = source)
) +
  ggplot2::geom_point(
    data = support_a %>% dplyr::filter(source == "Empirical data"),
    size = 1.8,
    position = position_panel_a,
    na.rm = TRUE
  ) +
  ggplot2::geom_errorbar(
    data = support_a %>% dplyr::filter(source == "Homogeneous model"),
    ggplot2::aes(ymin = q025, ymax = q975),
    width = 0.18,
    linewidth = 0.45,
    position = position_panel_a,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    data = support_a %>% dplyr::filter(source == "Homogeneous model"),
    size = 1.8,
    position = position_panel_a,
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~ dataset, ncol = 5, scales = "free") +
  ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
  ggplot2::scale_x_discrete(drop = TRUE) +
  theme_pub() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1)
  ) +
  ggplot2::labs(
    x = "Sites where the pair co-occurs",
    y = "Average number of sites supporting each regional interaction",
    colour = NULL,
    title = "A. Local support conditional on co-occurrence support"
  )

ret_emp <- retention %>%
  dplyr::group_by(dataset, removal_fraction) %>%
  dplyr::summarise(
    source = "Empirical data",
    mean_links = mean(empirical_links, na.rm = TRUE),
    q025 = stats::quantile(empirical_links, 0.025, na.rm = TRUE),
    q975 = stats::quantile(empirical_links, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ret_mod <- retention %>%
  dplyr::group_by(dataset, removal_fraction) %>%
  dplyr::summarise(
    source = "Homogeneous model",
    mean_links = mean(model_mean_links, na.rm = TRUE),
    ## Model ribbon: analytical binomial approximation around the fixed-p expectation.
    ## It is not variation in p and not a refitted model envelope.
    q025 = mean(model_q025, na.rm = TRUE),
    q975 = mean(model_q975, na.rm = TRUE),
    .groups = "drop"
  )

ret_plot <- dplyr::bind_rows(ret_emp, ret_mod)

panel_b <- ggplot2::ggplot(
  ret_plot,
  ggplot2::aes(
    x = removal_fraction,
    y = mean_links,
    colour = source,
    fill = source,
    group = source
  )
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q975), alpha = 0.14, colour = NA) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 1.3) +
  ggplot2::facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
  ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_pub() +
  ggplot2::labs(
    x = "Sites removed",
    y = "Number of regional interaction links remaining",
    colour = NULL,
    fill = NULL,
    title = "B. Regional interaction links remaining under site removal"
  )

title_grob <- grid::textGrob(
  "Repeated local realisation explains persistence beyond homogeneous co-occurrence opportunity",
  gp = grid::gpar(fontface = "bold", fontsize = 13)
)

fig_grid <- gridExtra::arrangeGrob(
  title_grob,
  panel_a,
  panel_b,
  ncol = 1,
  heights = c(0.08, 1, 1)
)

save_png(panel_a, "Figure3A_support_given_cooccurrence.png", 13, 7)
save_png(panel_b, "Figure3B_absolute_links_remaining.png", 13, 7)
save_png_grid(fig_grid, "Figure3_repeated_support_beyond_homogeneous_opportunity.png", 13, 13.5)

message("Saved Figure 3 PNG outputs in: ", out_dir)
message("Baseline used: unconditioned homogeneous per-site p fitted once to match the full-landscape number of regional links.")
