## ------------------------------------------------------------
## Quercus helper functions
## Run from the parent folder of the repository.
## Assumes "Quercus/" exists.
## ------------------------------------------------------------

packages <- c("dplyr", "ggplot2", "igraph", "bipartite", "rJava", "XLConnect", "tibble")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

read_quercus_metadata <- function(){
  metadata <- read.csv("Quercus/metadata.csv", stringsAsFactors = FALSE)
  metadata <- metadata[1:22, 1:4]
  metadata
}

read_quercus_webs <- function(){

  metadata <- read_quercus_metadata()
  webs <- list()

  for(t in metadata$Tree){

    file <- paste0("Quercus/raw-data/Web", t, " 2006 Kaartinen.xlsx")
    excel <- XLConnect::loadWorkbook(file)
    int_mat <- XLConnect::readWorksheet(object = excel, sheet = "Sheet1")

    row.names(int_mat) <- int_mat$Col1
    int_mat <- int_mat[-1]

    if(length(which(colnames(int_mat) %in% c("Unpar", "unpar"))) > 0){
      int_mat <- int_mat[, -which(colnames(int_mat) %in% c("Unpar", "unpar")), drop = FALSE]
    }

    int_mat[int_mat != 0] <- 1
    webs[[as.character(t)]] <- as.matrix(int_mat)
  }

  webs
}

clean_quercus_web <- function(web){

  web <- as.data.frame(web)

  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)){
    web <- web[, !bad_cols, drop = FALSE]
  }

  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)){
    web <- web[!bad_rows, , drop = FALSE]
  }

  web <- as.matrix(web)
  web[web != 0] <- 1

  if(nrow(web) > 0){
    r_remove <- which(rowSums(web, na.rm = TRUE) == 0)
    if(length(r_remove) > 0){
      web <- web[-r_remove, , drop = FALSE]
    }
  }

  if(ncol(web) > 0){
    c_remove <- which(colSums(web, na.rm = TRUE) == 0)
    if(length(c_remove) > 0){
      web <- web[, -c_remove, drop = FALSE]
    }
  }

  web
}

get_site_web <- function(webs, site_id){

  if(!is.null(names(webs))){
    if(as.character(site_id) %in% names(webs)){
      return(clean_quercus_web(webs[[as.character(site_id)]]))
    }
  }

  site_num <- suppressWarnings(as.integer(site_id))
  if(!is.na(site_num) && site_num >= 1 && site_num <= length(webs)){
    return(clean_quercus_web(webs[[site_num]]))
  }

  stop(paste("Could not find Quercus site:", site_id))
}

build_quercus_site_tables <- function(webs, all_sites){

  cooc_list <- list()
  interaction_list <- list()
  occupancy_list <- list()

  counter_cooc <- 1
  counter_int <- 1
  counter_occ <- 1

  for(site_id in all_sites){

    web <- get_site_web(webs, site_id)

    if(nrow(web) == 0 || ncol(web) == 0){
      next
    }

    consumers <- colnames(web)
    resources <- rownames(web)

    cur_occ <- rbind(
      data.frame(site = site_id,
                 species = consumers,
                 trophic_level = "consumer"),
      data.frame(site = site_id,
                 species = resources,
                 trophic_level = "resource")
    )

    occupancy_list[[counter_occ]] <- cur_occ
    counter_occ <- counter_occ + 1

    cur_cooc <- expand.grid(
      consumer = consumers,
      resource = resources,
      stringsAsFactors = FALSE
    )
    cur_cooc$site <- site_id

    cooc_list[[counter_cooc]] <- cur_cooc[, c("site", "consumer", "resource")]
    counter_cooc <- counter_cooc + 1

    positive_cells <- which(web > 0, arr.ind = TRUE)

    if(nrow(positive_cells) > 0){
      cur_int <- data.frame(
        site = site_id,
        consumer = colnames(web)[positive_cells[, 2]],
        resource = rownames(web)[positive_cells[, 1]]
      ) %>%
        filter(!is.na(consumer),
               !is.na(resource),
               consumer != "",
               resource != "") %>%
        distinct()

      interaction_list[[counter_int]] <- cur_int
      counter_int <- counter_int + 1
    }
  }

  list(
    occupancy = bind_rows(occupancy_list) %>%
      distinct(site, species, trophic_level),
    cooc_triples = bind_rows(cooc_list) %>%
      distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>%
      distinct(site, consumer, resource)
  )
}

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0){
    return(NA_real_)
  }
  (observed - expected) / expected
}

## ------------------------------------------------------------
## Script: quercus_model_envelope_site_removal.r
## ------------------------------------------------------------

set.seed(123)

p_fixed <- 0.0227
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100
n_model_reps <- 100

out_dir <- "Quercus/envelope"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

webs <- read_quercus_webs()
all_sites <- names(webs)

site_tables <- build_quercus_site_tables(webs, all_sites)
cooc_triples <- site_tables$cooc_triples
empirical_site_interactions <- site_tables$empirical_site_interactions
all_sites <- sort(unique(cooc_triples$site))

make_site_subsets <- function(all_sites, removal_levels, n_site_reps){

  site_subsets <- list()
  subset_index <- NULL
  counter <- 1

  for(removal in removal_levels){

    n_total <- length(all_sites)
    n_keep <- max(1, round(n_total * (1 - removal)))
    reps_here <- ifelse(removal == 0, 1, n_site_reps)

    for(r in seq_len(reps_here)){

      sites_keep <- if(removal == 0){
        all_sites
      } else {
        sample(all_sites, size = n_keep, replace = FALSE)
      }

      site_subsets[[counter]] <- sites_keep

      cur <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        site_rep = r,
        n_sites_kept = length(sites_keep)
      )

      subset_index <- rbind(subset_index, cur)
      counter <- counter + 1
    }
  }

  list(index = subset_index, subsets = site_subsets)
}

expected_links_from_subset <- function(cooc_triples, sites_keep, p_fixed){

  cur <- cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")

  if(nrow(cur) == 0){
    return(data.frame(
      n_cooc_links = 0,
      expected_links_raw = NA_real_,
      expected_links_conditioned = NA_real_
    ))
  }

  cur <- cur %>%
    group_by(consumer) %>%
    mutate(N_alpha = sum(n_cooc)) %>%
    ungroup() %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )

  data.frame(
    n_cooc_links = nrow(cur),
    expected_links_raw = sum(cur$prob_int, na.rm = TRUE),
    expected_links_conditioned = sum(cur$prob_expected, na.rm = TRUE)
  )
}

observed_links_from_subset <- function(site_interactions, sites_keep){
  site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource) %>%
    nrow()
}

simulate_model_site_interactions <- function(cooc_triples, p_fixed){
  cooc_triples %>%
    mutate(interaction = rbinom(n(), size = 1, prob = p_fixed)) %>%
    filter(interaction == 1) %>%
    select(site, consumer, resource) %>%
    distinct()
}

evaluate_network_against_subsets <- function(site_interactions,
                                             subset_object,
                                             subset_expected,
                                             source,
                                             model_rep = NA_integer_){

  out <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){

    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    exp_row <- subset_expected[subset_expected$subset_id == subset_id, ]

    observed_links <- observed_links_from_subset(site_interactions, sites_keep)

    out[[i]] <- data.frame(
      source = source,
      model_rep = model_rep,
      subset_id = subset_id,
      removal_fraction = subset_object$index$removal_fraction[i],
      site_rep = subset_object$index$site_rep[i],
      n_sites_kept = subset_object$index$n_sites_kept[i],
      n_observed_links = observed_links,
      n_cooc_links = exp_row$n_cooc_links,
      expected_links_raw = exp_row$expected_links_raw,
      expected_links_conditioned = exp_row$expected_links_conditioned,
      divergence_raw = observed_links - exp_row$expected_links_raw,
      relative_divergence_raw = safe_relative_difference(observed_links, exp_row$expected_links_raw),
      divergence_conditioned = observed_links - exp_row$expected_links_conditioned,
      relative_divergence_conditioned = safe_relative_difference(observed_links, exp_row$expected_links_conditioned)
    )
  }

  bind_rows(out)
}

subset_object <- make_site_subsets(
  all_sites = all_sites,
  removal_levels = removal_levels,
  n_site_reps = n_site_reps
)

subset_expected <- bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){

  subset_id <- subset_object$index$subset_id[i]
  sites_keep <- subset_object$subsets[[subset_id]]

  cur_exp <- expected_links_from_subset(
    cooc_triples = cooc_triples,
    sites_keep = sites_keep,
    p_fixed = p_fixed
  )

  cur_exp$subset_id <- subset_id
  cur_exp
}))

empirical_results <- evaluate_network_against_subsets(
  site_interactions = empirical_site_interactions,
  subset_object = subset_object,
  subset_expected = subset_expected,
  source = "empirical"
)

model_results_list <- vector("list", n_model_reps)

for(m in seq_len(n_model_reps)){

  message("Running Quercus model replicate ", m, " / ", n_model_reps)

  synthetic_site_interactions <- simulate_model_site_interactions(
    cooc_triples = cooc_triples,
    p_fixed = p_fixed
  )

  model_results_list[[m]] <- evaluate_network_against_subsets(
    site_interactions = synthetic_site_interactions,
    subset_object = subset_object,
    subset_expected = subset_expected,
    source = "model_generated",
    model_rep = m
  )
}

model_results <- bind_rows(model_results_list)

empirical_summary <- empirical_results %>%
  group_by(removal_fraction) %>%
  summarise(
    empirical_mean_relative_divergence_raw =
      mean(relative_divergence_raw, na.rm = TRUE),
    empirical_sd_relative_divergence_raw =
      sd(relative_divergence_raw, na.rm = TRUE),
    empirical_mean_divergence_raw =
      mean(divergence_raw, na.rm = TRUE),
    empirical_mean_observed_links =
      mean(n_observed_links, na.rm = TRUE),
    empirical_mean_expected_links_raw =
      mean(expected_links_raw, na.rm = TRUE),
    .groups = "drop"
  )

model_mean_curves <- model_results %>%
  group_by(model_rep, removal_fraction) %>%
  summarise(
    model_mean_relative_divergence_raw =
      mean(relative_divergence_raw, na.rm = TRUE),
    model_mean_divergence_raw =
      mean(divergence_raw, na.rm = TRUE),
    .groups = "drop"
  )

model_envelope <- model_mean_curves %>%
  group_by(removal_fraction) %>%
  summarise(
    model_q025_relative_divergence_raw =
      quantile(model_mean_relative_divergence_raw, 0.025, na.rm = TRUE),
    model_q500_relative_divergence_raw =
      quantile(model_mean_relative_divergence_raw, 0.5, na.rm = TRUE),
    model_q975_relative_divergence_raw =
      quantile(model_mean_relative_divergence_raw, 0.975, na.rm = TRUE),
    model_mean_relative_divergence_raw =
      mean(model_mean_relative_divergence_raw, na.rm = TRUE),
    .groups = "drop"
  )

comparison_summary <- left_join(empirical_summary,
                                model_envelope,
                                by = "removal_fraction") %>%
  mutate(
    empirical_above_model_975 =
      empirical_mean_relative_divergence_raw > model_q975_relative_divergence_raw,
    empirical_below_model_025 =
      empirical_mean_relative_divergence_raw < model_q025_relative_divergence_raw
  )

write.csv(empirical_results,
          file.path(out_dir, "quercus_empirical_divergence_by_subset.csv"),
          row.names = FALSE)

write.csv(model_results,
          file.path(out_dir, "quercus_model_generated_divergence_by_subset.csv"),
          row.names = FALSE)

write.csv(comparison_summary,
          file.path(out_dir, "quercus_empirical_vs_model_envelope_summary.csv"),
          row.names = FALSE)

p <- ggplot() +
  geom_ribbon(
    data = model_envelope,
    aes(x = removal_fraction,
        ymin = model_q025_relative_divergence_raw,
        ymax = model_q975_relative_divergence_raw),
    fill = "grey80"
  ) +
  geom_line(
    data = model_envelope,
    aes(x = removal_fraction,
        y = model_q500_relative_divergence_raw),
    linewidth = 1,
    linetype = 2
  ) +
  geom_line(
    data = empirical_summary,
    aes(x = removal_fraction,
        y = empirical_mean_relative_divergence_raw),
    linewidth = 1.2
  ) +
  geom_point(
    data = empirical_summary,
    aes(x = removal_fraction,
        y = empirical_mean_relative_divergence_raw),
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme_classic(base_size = 14) +
  xlab("Fraction of sites removed") +
  ylab("Mean relative divergence") +
  ggtitle(
    "Quercus: empirical divergence vs model-generated envelope",
    subtitle = "Ribbon = 95% envelope of model-generated mean curves; solid line = empirical"
  )

ggsave(file.path(out_dir, "quercus_empirical_vs_model_envelope.png"),
       p,
       width = 8,
       height = 5.5,
       dpi = 300)

print(comparison_summary)
