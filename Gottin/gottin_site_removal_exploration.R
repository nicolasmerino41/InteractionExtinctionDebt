## ------------------------------------------------------------
## Script: gottin_site_removal_exploration.R
##
## Purpose:
## First exploratory analysis of model divergence under random
## removal of spatial units in the Gottin dataset.
##
## Can run either:
##   - Gottin-HP: host-parasitoid network
##   - Gottin-PP: plant-pollinator network
#### ------------------------------------------------------------
## Script: gottin_site_removal_exploration.R
##
## Purpose:
## Exploratory analysis of signed model divergence under random
## removal of spatial units in the Gottin dataset.
##
## Can run either:
##   - Gottin-HP: host-parasitoid network
##   - Gottin-PP: plant-pollinator network
##
## Run from the parent folder of the repository.
## Assumes the folder "Gottin/" exists.
## ------------------------------------------------------------


## ---------------------------
## 0. Choose network
## ---------------------------

## Choose either "HP" or "PP"
network_type <- "PP"


## ---------------------------
## 1. Load packages
## ---------------------------

packages <- c("dplyr", "ggplot2", "bipartite", "igraph")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


## ---------------------------
## 2. Read Gottin data
## ---------------------------

metadata <- read.csv("Gottin/metadata.csv", head = TRUE, sep = ",")

if(network_type == "HP"){
  
  raw_data <- read.csv("Gottin/raw-data/host_para_all_interactions.csv",
                       head = TRUE,
                       sep = ";")
  
  webs <- frame2webs(
    raw_data,
    varnames = c("Genus.Species",
                 "P1.Species.Genus",
                 "Site",
                 "P1.cells")
  )
  
  p_fixed <- 0.074
  output_prefix <- "gottin_HP"
  
} else if(network_type == "PP"){
  
  raw_data <- read.csv("Gottin/raw-data/plant_poll_all_interactions.csv",
                       head = TRUE,
                       sep = ";")
  
  webs <- frame2webs(
    raw_data,
    varnames = c("P.Genus.Species",
                 "Genus.Species",
                 "Site")
  )
  
  p_fixed <- 0.044
  output_prefix <- "gottin_PP"
  
} else {
  stop("network_type must be either 'HP' or 'PP'")
}


## ---------------------------
## 3. Helper functions
## ---------------------------

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0){
    return(NA_real_)
  } else {
    return((observed - expected) / expected)
  }
}


get_site_web <- function(webs, site_id){
  
  ## First try list names
  if(!is.null(names(webs))){
    if(as.character(site_id) %in% names(webs)){
      return(webs[[as.character(site_id)]])
    }
  }
  
  ## Then try numeric indexing, matching the original script style
  if(is.numeric(site_id) || suppressWarnings(!is.na(as.numeric(site_id)))){
    site_num <- as.numeric(site_id)
    if(site_num >= 1 && site_num <= length(webs)){
      return(webs[[site_num]])
    }
  }
  
  stop(paste("Could not find web for site:", site_id))
}


clean_web_matrix <- function(web){
  
  web <- as.data.frame(web)
  
  ## Remove columns with missing or empty names.
  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)){
    web <- web[, !bad_cols, drop = FALSE]
  }
  
  ## Remove rows with missing or empty names.
  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)){
    web <- web[!bad_rows, , drop = FALSE]
  }
  
  return(as.matrix(web))
}


extract_realised_pairs_from_web <- function(web){
  
  web <- clean_web_matrix(web)
  
  if(nrow(web) == 0 || ncol(web) == 0){
    return(data.frame(resource = character(),
                      consumer = character()))
  }
  
  positive_cells <- which(web > 0, arr.ind = TRUE)
  
  if(nrow(positive_cells) == 0){
    return(data.frame(resource = character(),
                      consumer = character()))
  }
  
  data.frame(
    resource = rownames(web)[positive_cells[, 1]],
    consumer = colnames(web)[positive_cells[, 2]]
  )
}


build_occupancy_from_webs <- function(webs, sites_keep){
  
  occupancy_list <- list()
  counter <- 1
  
  for(site_id in sites_keep){
    
    web <- get_site_web(webs, site_id)
    web <- clean_web_matrix(web)
    
    if(nrow(web) == 0 && ncol(web) == 0){
      next
    }
    
    consumers <- colnames(web)
    resources <- rownames(web)
    
    cur_consumers <- data.frame(
      site = site_id,
      species = consumers,
      trophic_level = "consumer"
    )
    
    cur_resources <- data.frame(
      site = site_id,
      species = resources,
      trophic_level = "resource"
    )
    
    occupancy_list[[counter]] <- rbind(cur_consumers, cur_resources)
    counter <- counter + 1
  }
  
  occupancy <- dplyr::bind_rows(occupancy_list)
  
  occupancy <- occupancy %>%
    filter(!is.na(species),
           species != "")
  
  return(unique(occupancy))
}


build_realised_pairs_from_webs <- function(webs, sites_keep){
  
  pair_list <- list()
  counter <- 1
  
  for(site_id in sites_keep){
    
    web <- get_site_web(webs, site_id)
    
    cur_pairs <- extract_realised_pairs_from_web(web)
    
    if(nrow(cur_pairs) > 0){
      cur_pairs$site <- site_id
      pair_list[[counter]] <- cur_pairs
      counter <- counter + 1
    }
  }
  
  realised_pairs <- dplyr::bind_rows(pair_list)
  
  if(nrow(realised_pairs) == 0){
    return(data.frame(resource = character(),
                      consumer = character(),
                      site = character()))
  }
  
  realised_pairs <- realised_pairs %>%
    filter(!is.na(resource),
           !is.na(consumer),
           resource != "",
           consumer != "")
  
  return(unique(realised_pairs))
}


## ---------------------------
## 4. Function for one site subset
## ---------------------------

run_gottin_subset <- function(webs,
                              sites_keep,
                              p_fixed,
                              removal_fraction = NA_real_,
                              rep_id = NA_integer_){
  
  ## Build occupancy from kept sites
  occupancy <- build_occupancy_from_webs(webs, sites_keep)
  
  if(nrow(occupancy) == 0){
    return(data.frame(
      removal_fraction = removal_fraction,
      rep = rep_id,
      n_sites_kept = length(unique(sites_keep)),
      n_consumers = 0,
      n_resources = 0,
      n_realised_links = 0,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_species_bias_raw = NA_real_,
      mean_species_bias_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }
  
  ## Build realised interactions from kept sites
  realised_site_pairs <- build_realised_pairs_from_webs(webs, sites_keep)
  
  realised_pairs <- realised_site_pairs %>%
    distinct(resource, consumer)
  
  n_realised_links <- nrow(realised_pairs)
  
  ## Build co-occurrences from occupancy
  consumer_sites <- occupancy %>%
    filter(trophic_level == "consumer") %>%
    select(site, consumer = species)
  
  resource_sites <- occupancy %>%
    filter(trophic_level == "resource") %>%
    select(site, resource = species)
  
  cooc <- inner_join(consumer_sites,
                     resource_sites,
                     by = "site") %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop") %>%
    transmute(
      species = consumer,
      potential_prey = resource,
      n_cooc = n_cooc
    )
  
  n_cooc_links <- nrow(cooc)
  
  if(n_cooc_links == 0){
    return(data.frame(
      removal_fraction = removal_fraction,
      rep = rep_id,
      n_sites_kept = length(unique(sites_keep)),
      n_consumers = length(unique(consumer_sites$consumer)),
      n_resources = length(unique(resource_sites$resource)),
      n_realised_links = n_realised_links,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_species_bias_raw = NA_real_,
      mean_species_bias_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }
  
  ## N_alpha = total co-occurrence frequency per consumer
  cooc <- cooc %>%
    group_by(species) %>%
    mutate(
      N_alpha = sum(n_cooc),
      total_potential_interactions = n_distinct(potential_prey)
    ) %>%
    ungroup()
  
  ## Interaction-rate model
  cooc <- cooc %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )
  
  ## Expected total links
  n_expected_links_raw <- sum(cooc$prob_int, na.rm = TRUE)
  n_expected_links_conditioned <- sum(cooc$prob_expected, na.rm = TRUE)
  
  ## Species-level expected degree
  expected_degree <- cooc %>%
    group_by(species) %>%
    summarise(
      expected_indegree_raw = sum(prob_int, na.rm = TRUE),
      expected_indegree_conditioned = sum(prob_expected, na.rm = TRUE),
      potential_degree = n_distinct(potential_prey),
      .groups = "drop"
    )
  
  observed_degree <- realised_pairs %>%
    group_by(species = consumer) %>%
    summarise(
      observed_indegree = n_distinct(resource),
      .groups = "drop"
    )
  
  degree_compare <- merge(expected_degree,
                          observed_degree,
                          by = "species",
                          all.x = TRUE)
  
  degree_compare$observed_indegree[is.na(degree_compare$observed_indegree)] <- 0
  
  ## Signed species-level bias:
  ## positive = observed > expected; negative = observed < expected.
  mean_species_bias_raw <- mean(
    degree_compare$observed_indegree -
      degree_compare$expected_indegree_raw,
    na.rm = TRUE
  )
  
  mean_species_bias_conditioned <- mean(
    degree_compare$observed_indegree -
      degree_compare$expected_indegree_conditioned,
    na.rm = TRUE
  )
  
  ## Absolute species-level error:
  ## magnitude of mismatch, ignoring direction.
  mean_abs_species_error_raw <- mean(
    abs(degree_compare$observed_indegree -
          degree_compare$expected_indegree_raw),
    na.rm = TRUE
  )
  
  mean_abs_species_error_conditioned <- mean(
    abs(degree_compare$observed_indegree -
          degree_compare$expected_indegree_conditioned),
    na.rm = TRUE
  )
  
  ## Signed network-level divergence:
  ## positive = more realised links than expected;
  ## negative = fewer realised links than expected.
  divergence_raw <- n_realised_links - n_expected_links_raw
  divergence_conditioned <- n_realised_links - n_expected_links_conditioned
  
  relative_divergence_raw <- safe_relative_difference(
    n_realised_links,
    n_expected_links_raw
  )
  
  relative_divergence_conditioned <- safe_relative_difference(
    n_realised_links,
    n_expected_links_conditioned
  )
  
  data.frame(
    removal_fraction = removal_fraction,
    rep = rep_id,
    n_sites_kept = length(unique(sites_keep)),
    n_consumers = length(unique(consumer_sites$consumer)),
    n_resources = length(unique(resource_sites$resource)),
    n_realised_links = n_realised_links,
    n_cooc_links = n_cooc_links,
    n_expected_links_raw = n_expected_links_raw,
    n_expected_links_conditioned = n_expected_links_conditioned,
    divergence_raw = divergence_raw,
    relative_divergence_raw = relative_divergence_raw,
    divergence_conditioned = divergence_conditioned,
    relative_divergence_conditioned = relative_divergence_conditioned,
    mean_species_bias_raw = mean_species_bias_raw,
    mean_species_bias_conditioned = mean_species_bias_conditioned,
    mean_abs_species_error_raw = mean_abs_species_error_raw,
    mean_abs_species_error_conditioned = mean_abs_species_error_conditioned
  )
}


## ---------------------------
## 5. Random site-removal experiment
## ---------------------------

set.seed(123)

all_sites <- sort(unique(metadata$site))

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 100

results_list <- list()
counter <- 1

for(removal in removal_levels){
  
  n_sites_total <- length(all_sites)
  n_sites_keep <- max(1, round(n_sites_total * (1 - removal)))
  
  reps_here <- ifelse(removal == 0, 1, n_replicates)
  
  for(r in seq_len(reps_here)){
    
    if(removal == 0){
      sites_keep <- all_sites
    } else {
      sites_keep <- sample(all_sites,
                           size = n_sites_keep,
                           replace = FALSE)
    }
    
    results_list[[counter]] <- run_gottin_subset(
      webs = webs,
      sites_keep = sites_keep,
      p_fixed = p_fixed,
      removal_fraction = removal,
      rep_id = r
    )
    
    counter <- counter + 1
  }
}

site_removal_results <- bind_rows(results_list)


## ---------------------------
## 6. Save results
## ---------------------------

write.csv(site_removal_results,
          file = paste0("Gottin/", output_prefix,
                        "_site_removal_divergence.csv"),
          row.names = FALSE)


## ---------------------------
## 7. Plot relative divergence: raw probabilities
## ---------------------------

plot_raw <- ggplot(site_removal_results,
                   aes(x = removal_fraction,
                       y = relative_divergence_raw)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.35, size = 2) +
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               linewidth = 1.2) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.02) +
  theme_classic(base_size = 14) +
  xlab("Fraction of sites removed") +
  ylab("Relative divergence: (observed - expected) / expected") +
  ggtitle(paste0("Gottin-", network_type,
                 ": signed divergence under random site removal"),
          subtitle = "Expected links calculated with raw pairwise interaction probabilities")

ggsave(paste0("Gottin/", output_prefix,
              "_site_removal_relative_divergence_raw.png"),
       plot_raw,
       width = 7,
       height = 5,
       dpi = 300)


## ---------------------------
## 8. Plot relative divergence: conditioned probabilities
## ---------------------------

plot_conditioned <- ggplot(site_removal_results,
                           aes(x = removal_fraction,
                               y = relative_divergence_conditioned)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.35, size = 2) +
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               linewidth = 1.2) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.02) +
  theme_classic(base_size = 14) +
  xlab("Fraction of sites removed") +
  ylab("Relative divergence: (observed - expected) / expected") +
  ggtitle(paste0("Gottin-", network_type,
                 ": signed divergence under random site removal"),
          subtitle = "Expected links calculated with conditioned probabilities")

ggsave(paste0("Gottin/", output_prefix,
              "_site_removal_relative_divergence_conditioned.png"),
       plot_conditioned,
       width = 7,
       height = 5,
       dpi = 300)


## ---------------------------
## 9. Plot signed absolute-scale divergence
## ---------------------------

plot_signed <- ggplot(site_removal_results,
                      aes(x = removal_fraction,
                          y = divergence_raw)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.35, size = 2) +
  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               linewidth = 1.2) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.02) +
  theme_classic(base_size = 14) +
  xlab("Fraction of sites removed") +
  ylab("Observed links - expected links") +
  ggtitle(paste0("Gottin-", network_type,
                 ": signed divergence under random site removal"))

ggsave(paste0("Gottin/", output_prefix,
              "_site_removal_signed_divergence_raw.png"),
       plot_signed,
       width = 7,
       height = 5,
       dpi = 300)


## ---------------------------
## 10. Console summary
## ---------------------------

summary_table <- site_removal_results %>%
  group_by(removal_fraction) %>%
  summarise(
    mean_realised_links = mean(n_realised_links, na.rm = TRUE),
    mean_expected_links_raw = mean(n_expected_links_raw, na.rm = TRUE),
    mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
    mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
    sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
    mean_species_bias_raw = mean(mean_species_bias_raw, na.rm = TRUE),
    mean_species_bias_conditioned = mean(mean_species_bias_conditioned, na.rm = TRUE),
    mean_abs_species_error_raw = mean(mean_abs_species_error_raw, na.rm = TRUE),
    mean_abs_species_error_conditioned = mean(mean_abs_species_error_conditioned, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

write.csv(summary_table,
          file = paste0("Gottin/", output_prefix,
                        "_site_removal_summary.csv"),
          row.names = FALSE)
