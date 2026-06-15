
## ------------------------------------------------------------
## Script: gottin_PP_degree_site_removal.R
## ------------------------------------------------------------
dataset_label <- "Gottin_PP"
out_dir <- "Gottin/degree_site_removal_PP"
p_fixed <- 0.044
load_dataset_webs <- function(){
  if(!requireNamespace("bipartite", quietly = TRUE)) install.packages("bipartite")
  library(bipartite)
  raw_data <- read.csv("Gottin/raw-data/plant_poll_all_interactions.csv", head = TRUE, sep = ";")
  webs <- frame2webs(raw_data, varnames = c("P.Genus.Species", "Genus.Species", "Site"))
  names(webs) <- as.character(names(webs))
  webs
}

## ---------------------------
## Shared helpers
## ---------------------------

packages <- c("dplyr", "tidyr", "ggplot2", "tibble")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 100
n_model_reps <- 100
selected_removal_levels <- c(0, 0.4, 0.8)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_web_matrix <- function(web){
  web <- as.data.frame(web)
  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)) web <- web[, !bad_cols, drop = FALSE]
  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)) web <- web[!bad_rows, , drop = FALSE]
  web <- as.matrix(web)
  suppressWarnings(storage.mode(web) <- "numeric")
  web[is.na(web)] <- 0
  web[web != 0] <- 1
  web
}

build_site_tables_from_webs <- function(webs){
  cooc_list <- list()
  interaction_list <- list()
  counter_cooc <- 1
  counter_int <- 1
  for(site_id in names(webs)){
    web <- clean_web_matrix(webs[[site_id]])
    if(nrow(web) == 0 || ncol(web) == 0) next
    consumers <- colnames(web)
    resources <- rownames(web)
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    if(length(consumers) == 0 || length(resources) == 0) next
    cur_cooc <- expand.grid(consumer = consumers, resource = resources, stringsAsFactors = FALSE)
    cur_cooc$site <- site_id
    cooc_list[[counter_cooc]] <- cur_cooc[, c("site", "consumer", "resource")]
    counter_cooc <- counter_cooc + 1
    positive_cells <- which(web > 0, arr.ind = TRUE)
    if(nrow(positive_cells) > 0){
      cur_int <- data.frame(
        site = site_id,
        consumer = colnames(web)[positive_cells[, 2]],
        resource = rownames(web)[positive_cells[, 1]],
        stringsAsFactors = FALSE
      ) %>%
        filter(!is.na(consumer), !is.na(resource), consumer != "", resource != "") %>%
        distinct(site, consumer, resource)
      interaction_list[[counter_int]] <- cur_int
      counter_int <- counter_int + 1
    }
  }
  list(
    cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>% distinct(site, consumer, resource)
  )
}

make_site_subsets <- function(all_sites, removal_levels, n_site_reps){
  site_subsets <- list()
  subset_index <- NULL
  counter <- 1
  for(removal in removal_levels){
    n_total <- length(all_sites)
    n_keep <- max(1, round(n_total * (1 - removal)))
    reps_here <- ifelse(removal == 0, 1, n_site_reps)
    for(r in seq_len(reps_here)){
      sites_keep <- if(removal == 0) all_sites else sample(all_sites, n_keep, replace = FALSE)
      site_subsets[[counter]] <- sites_keep
      subset_index <- rbind(subset_index, data.frame(
        subset_id = counter,
        removal_fraction = removal,
        site_rep = r,
        n_sites_kept = length(sites_keep)
      ))
      counter <- counter + 1
    }
  }
  list(index = subset_index, subsets = site_subsets)
}

simulate_model_site_interactions <- function(cooc_triples, p_fixed){
  cooc_triples %>%
    mutate(interaction = rbinom(n(), size = 1, prob = p_fixed)) %>%
    filter(interaction == 1) %>%
    select(site, consumer, resource) %>%
    distinct(site, consumer, resource)
}

gini_coefficient <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  if(all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_along(x) - n - 1) * x) / (n * sum(x))
}

degree_table_from_interactions <- function(site_interactions, sites_keep){
  cur <- site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource)
  consumer_degrees <- cur %>%
    group_by(species = consumer) %>%
    summarise(degree = n_distinct(resource), .groups = "drop") %>%
    mutate(trophic_level = "consumer")
  resource_degrees <- cur %>%
    group_by(species = resource) %>%
    summarise(degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(trophic_level = "resource")
  bind_rows(consumer_degrees, resource_degrees) %>%
    filter(!is.na(species), species != "", degree > 0) %>%
    select(trophic_level, species, degree)
}

summarise_degrees <- function(degrees){
  if(nrow(degrees) == 0){
    return(data.frame(
      trophic_level = c("consumer", "resource"),
      mean_degree = NA_real_, median_degree = NA_real_, variance_degree = NA_real_,
      maximum_degree = NA_real_, proportion_degree_1 = NA_real_,
      n_species_degree_gt0 = 0L, gini_degree = NA_real_
    ))
  }
  degrees %>%
    group_by(trophic_level) %>%
    summarise(
      mean_degree = mean(degree, na.rm = TRUE),
      median_degree = median(degree, na.rm = TRUE),
      variance_degree = ifelse(n() > 1, var(degree, na.rm = TRUE), 0),
      maximum_degree = max(degree, na.rm = TRUE),
      proportion_degree_1 = mean(degree == 1, na.rm = TRUE),
      n_species_degree_gt0 = n_distinct(species),
      gini_degree = gini_coefficient(degree),
      .groups = "drop"
    ) %>%
    right_join(data.frame(trophic_level = c("consumer", "resource")), by = "trophic_level")
}

frequency_degrees <- function(degrees){
  if(nrow(degrees) == 0){
    return(data.frame(trophic_level = character(), degree = integer(), n_species = integer()))
  }
  degrees %>% count(trophic_level, degree, name = "n_species") %>% arrange(trophic_level, degree)
}

evaluate_network_for_subsets <- function(site_interactions, subset_object, source, model_rep = NA_integer_){
  metric_list <- vector("list", nrow(subset_object$index))
  freq_list <- vector("list", nrow(subset_object$index))
  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    degrees <- degree_table_from_interactions(site_interactions, sites_keep)
    cur_metrics <- summarise_degrees(degrees) %>%
      mutate(source = source, model_rep = model_rep, subset_id = subset_id,
             removal_fraction = subset_object$index$removal_fraction[i],
             site_rep = subset_object$index$site_rep[i],
             n_sites_kept = subset_object$index$n_sites_kept[i]) %>%
      select(source, model_rep, subset_id, removal_fraction, site_rep, n_sites_kept, trophic_level, everything())
    cur_freq <- frequency_degrees(degrees) %>%
      mutate(source = source, model_rep = model_rep, subset_id = subset_id,
             removal_fraction = subset_object$index$removal_fraction[i],
             site_rep = subset_object$index$site_rep[i],
             n_sites_kept = subset_object$index$n_sites_kept[i]) %>%
      select(source, model_rep, subset_id, removal_fraction, site_rep, n_sites_kept, trophic_level, degree, n_species)
    metric_list[[i]] <- cur_metrics
    freq_list[[i]] <- cur_freq
  }
  list(metrics = bind_rows(metric_list), frequency = bind_rows(freq_list))
}

make_metric_envelope <- function(model_metrics, empirical_metrics){
  metric_cols <- c("mean_degree", "median_degree", "variance_degree", "maximum_degree",
                   "proportion_degree_1", "n_species_degree_gt0", "gini_degree")
  model_long <- model_metrics %>% pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value")
  model_mean_curves <- model_long %>%
    group_by(model_rep, trophic_level, removal_fraction, metric) %>%
    summarise(model_mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  model_envelope <- model_mean_curves %>%
    group_by(trophic_level, removal_fraction, metric) %>%
    summarise(model_q025 = quantile(model_mean_value, 0.025, na.rm = TRUE),
              model_q500 = quantile(model_mean_value, 0.500, na.rm = TRUE),
              model_q975 = quantile(model_mean_value, 0.975, na.rm = TRUE),
              model_mean = mean(model_mean_value, na.rm = TRUE), .groups = "drop")
  empirical_long <- empirical_metrics %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "empirical_value") %>%
    group_by(trophic_level, removal_fraction, metric) %>%
    summarise(empirical_mean = mean(empirical_value, na.rm = TRUE),
              empirical_sd = sd(empirical_value, na.rm = TRUE), .groups = "drop")
  empirical_long %>%
    left_join(model_envelope, by = c("trophic_level", "removal_fraction", "metric")) %>%
    mutate(empirical_above_model_975 = empirical_mean > model_q975,
           empirical_below_model_025 = empirical_mean < model_q025)
}

make_frequency_envelope <- function(model_frequency, empirical_frequency){
  model_summary <- model_frequency %>%
    group_by(model_rep, removal_fraction, trophic_level, degree) %>%
    summarise(n_species = mean(n_species, na.rm = TRUE), .groups = "drop")
  model_envelope <- model_summary %>%
    group_by(removal_fraction, trophic_level, degree) %>%
    summarise(model_q025 = quantile(n_species, 0.025, na.rm = TRUE),
              model_q500 = quantile(n_species, 0.500, na.rm = TRUE),
              model_q975 = quantile(n_species, 0.975, na.rm = TRUE), .groups = "drop")
  empirical_summary <- empirical_frequency %>%
    group_by(removal_fraction, trophic_level, degree) %>%
    summarise(empirical_mean = mean(n_species, na.rm = TRUE), .groups = "drop")
  degrees_all <- sort(unique(c(model_envelope$degree, empirical_summary$degree)))
  if(length(degrees_all) == 0) degrees_all <- 0
  full_grid <- expand.grid(removal_fraction = selected_removal_levels,
                           trophic_level = c("consumer", "resource"),
                           degree = degrees_all, stringsAsFactors = FALSE)
  full_grid %>%
    left_join(model_envelope, by = c("removal_fraction", "trophic_level", "degree")) %>%
    left_join(empirical_summary, by = c("removal_fraction", "trophic_level", "degree")) %>%
    mutate(across(c(model_q025, model_q500, model_q975, empirical_mean), ~replace_na(.x, 0)))
}

plot_metric_envelope <- function(comparison, metric_name, y_label, output_file){
  dat <- comparison %>% filter(metric == metric_name)
  p <- ggplot(dat, aes(x = removal_fraction)) +
    geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") +
    geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) +
    geom_line(aes(y = empirical_mean), linewidth = 1.2) +
    geom_point(aes(y = empirical_mean), size = 2.5) +
    facet_wrap(~ trophic_level, scales = "free_y") +
    theme_classic(base_size = 14) +
    xlab("Fraction of sites removed") + ylab(y_label) +
    ggtitle(paste0(dataset_label, ": ", y_label, " under random site removal"),
            subtitle = "Ribbon = 95% model-generated envelope; dashed = model median; solid = empirical")
  ggsave(output_file, p, width = 8.5, height = 5.5, dpi = 300)
}

plot_degree_frequency <- function(freq_comparison, output_file){
  dat <- freq_comparison %>% filter(removal_fraction %in% selected_removal_levels)
  p <- ggplot(dat, aes(x = degree)) +
    geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") +
    geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) +
    geom_line(aes(y = empirical_mean), linewidth = 1) +
    geom_point(aes(y = empirical_mean), size = 2) +
    facet_grid(trophic_level ~ removal_fraction, scales = "free_y") +
    theme_classic(base_size = 13) +
    xlab("Degree") + ylab("Number of species") +
    ggtitle(paste0(dataset_label, ": degree-frequency distributions"),
            subtitle = "Columns = removal fraction; ribbon = model envelope; solid = empirical")
  ggsave(output_file, p, width = 10, height = 6.5, dpi = 300)
}

run_degree_site_removal <- function(){
  message("Loading ", dataset_label)
  webs <- load_dataset_webs()
  site_tables <- build_site_tables_from_webs(webs)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions
  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  message(dataset_label, ": empirical degree metrics")
  empirical_out <- evaluate_network_for_subsets(empirical_site_interactions, subset_object, "empirical", NA_integer_)
  write.csv(empirical_out$metrics, file.path(out_dir, "degree_metrics_empirical.csv"), row.names = FALSE)
  write.csv(empirical_out$frequency, file.path(out_dir, "degree_frequency_empirical.csv"), row.names = FALSE)
  model_metrics_list <- vector("list", n_model_reps)
  model_frequency_list <- vector("list", n_model_reps)
  for(m in seq_len(n_model_reps)){
    message(dataset_label, ": model replicate ", m, " / ", n_model_reps)
    sim_interactions <- simulate_model_site_interactions(cooc_triples, p_fixed)
    model_out <- evaluate_network_for_subsets(sim_interactions, subset_object, "model_generated", m)
    model_metrics_list[[m]] <- model_out$metrics
    model_frequency_list[[m]] <- model_out$frequency
  }
  model_metrics <- bind_rows(model_metrics_list)
  model_frequency <- bind_rows(model_frequency_list)
  write.csv(model_metrics, file.path(out_dir, "degree_metrics_model_generated.csv"), row.names = FALSE)
  write.csv(model_frequency, file.path(out_dir, "degree_frequency_model_generated.csv"), row.names = FALSE)
  comparison <- make_metric_envelope(model_metrics, empirical_out$metrics)
  write.csv(comparison, file.path(out_dir, "degree_metrics_empirical_vs_model_envelope.csv"), row.names = FALSE)
  freq_comparison <- make_frequency_envelope(
    model_frequency %>% filter(removal_fraction %in% selected_removal_levels),
    empirical_out$frequency %>% filter(removal_fraction %in% selected_removal_levels)
  )
  write.csv(freq_comparison, file.path(out_dir, "degree_frequency_empirical_vs_model_envelope.csv"), row.names = FALSE)
  plot_metric_envelope(comparison, "mean_degree", "Mean degree", file.path(out_dir, "mean_degree_empirical_vs_model_envelope.png"))
  plot_metric_envelope(comparison, "maximum_degree", "Maximum degree", file.path(out_dir, "maximum_degree_empirical_vs_model_envelope.png"))
  plot_metric_envelope(comparison, "proportion_degree_1", "Proportion of degree-1 species", file.path(out_dir, "proportion_degree1_empirical_vs_model_envelope.png"))
  plot_metric_envelope(comparison, "gini_degree", "Degree Gini coefficient", file.path(out_dir, "degree_gini_empirical_vs_model_envelope.png"))
  plot_degree_frequency(freq_comparison, file.path(out_dir, "degree_frequency_selected_removals_empirical_vs_model_envelope.png"))
  message("Finished ", dataset_label, ". Outputs saved to ", out_dir)
}

run_degree_site_removal()
