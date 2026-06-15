## ------------------------------------------------------------
## Nahuel helper functions
## Run from the parent folder of the repository.
## Assumes "Nahuel/" exists.
## ------------------------------------------------------------

packages <- c("dplyr", "ggplot2", "igraph", "bipartite")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

read_nahuel_webs <- function(){
  
  files <- c(
    "vaz_ag_matr_f.txt",
    "vaz_cl_matr_f.txt",
    "vaz_ll_matr_f.txt",
    "vaz_mh_matr_f.txt",
    "vaz_mnh_matr_f.txt",
    "vaz_qh_matr_f.txt",
    "vaz_qnh_matr_f.txt",
    "vaz_s_matr_f.txt"
  )
  
  webs <- lapply(files, function(f){
    t(read.csv(file.path("Nahuel/raw-data", f),
               sep = "\t",
               header = FALSE))
  })
  
  names(webs) <- seq_along(webs)
  webs
}

clean_nahuel_web <- function(web, reference_web){
  
  web <- as.matrix(web)
  
  colnames(web) <- paste0("Pol", seq_len(ncol(reference_web)))
  rownames(web) <- paste0("Plant", seq_len(nrow(reference_web)))
  
  r_remove <- which(rowSums(web, na.rm = TRUE) == 0)
  if(length(r_remove) > 0){
    web <- web[-r_remove, , drop = FALSE]
  }
  
  c_remove <- which(colSums(web, na.rm = TRUE) == 0)
  if(length(c_remove) > 0){
    web <- web[, -c_remove, drop = FALSE]
  }
  
  web
}

get_site_web <- function(webs, site_id){
  
  site_num <- suppressWarnings(as.integer(site_id))
  
  if(is.na(site_num) || site_num < 1 || site_num > length(webs)){
    stop(paste("Could not find Nahuel site:", site_id))
  }
  
  clean_nahuel_web(webs[[site_num]], webs[[1]])
}

build_nahuel_site_tables <- function(webs, all_sites){
  
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
## Script: data_processing_occurrence_nahuel_site_removal.r
##
## Purpose:
## Signed site-removal divergence exploration for the Nahuel dataset.
## ------------------------------------------------------------

set.seed(123)

p_fixed <- 0.065
output_prefix <- "nahuel"

webs <- read_nahuel_webs()
all_sites <- seq_along(webs)

run_nahuel_subset <- function(webs,
                              sites_keep,
                              p_fixed,
                              removal_fraction = NA_real_,
                              rep_id = NA_integer_){
  
  site_tables <- build_nahuel_site_tables(webs, sites_keep)
  
  occupancy <- site_tables$occupancy
  cooc_triples <- site_tables$cooc_triples
  realised_site_pairs <- site_tables$empirical_site_interactions
  
  if(nrow(occupancy) == 0 || nrow(cooc_triples) == 0){
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
  
  realised_pairs <- realised_site_pairs %>%
    distinct(resource, consumer)
  
  n_realised_links <- nrow(realised_pairs)
  
  cooc <- cooc_triples %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop") %>%
    transmute(
      species = consumer,
      potential_prey = resource,
      n_cooc = n_cooc
    )
  
  n_cooc_links <- nrow(cooc)
  
  cooc <- cooc %>%
    group_by(species) %>%
    mutate(
      N_alpha = sum(n_cooc),
      total_potential_interactions = n_distinct(potential_prey)
    ) %>%
    ungroup() %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )
  
  n_expected_links_raw <- sum(cooc$prob_int, na.rm = TRUE)
  n_expected_links_conditioned <- sum(cooc$prob_expected, na.rm = TRUE)
  
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
  
  divergence_raw <- n_realised_links - n_expected_links_raw
  divergence_conditioned <- n_realised_links - n_expected_links_conditioned
  
  data.frame(
    removal_fraction = removal_fraction,
    rep = rep_id,
    n_sites_kept = length(unique(sites_keep)),
    n_consumers = length(unique(occupancy$species[occupancy$trophic_level == "consumer"])),
    n_resources = length(unique(occupancy$species[occupancy$trophic_level == "resource"])),
    n_realised_links = n_realised_links,
    n_cooc_links = n_cooc_links,
    n_expected_links_raw = n_expected_links_raw,
    n_expected_links_conditioned = n_expected_links_conditioned,
    divergence_raw = divergence_raw,
    relative_divergence_raw = safe_relative_difference(n_realised_links, n_expected_links_raw),
    divergence_conditioned = divergence_conditioned,
    relative_divergence_conditioned = safe_relative_difference(n_realised_links, n_expected_links_conditioned),
    mean_species_bias_raw = mean_species_bias_raw,
    mean_species_bias_conditioned = mean_species_bias_conditioned,
    mean_abs_species_error_raw = mean_abs_species_error_raw,
    mean_abs_species_error_conditioned = mean_abs_species_error_conditioned
  )
}

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 100

results_list <- list()
counter <- 1

for(removal in removal_levels){
  
  n_sites_total <- length(all_sites)
  n_sites_keep <- max(1, round(n_sites_total * (1 - removal)))
  reps_here <- ifelse(removal == 0, 1, n_replicates)
  
  for(r in seq_len(reps_here)){
    
    sites_keep <- if(removal == 0){
      all_sites
    } else {
      sample(all_sites, size = n_sites_keep, replace = FALSE)
    }
    
    results_list[[counter]] <- run_nahuel_subset(
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

write.csv(site_removal_results,
          file = "Nahuel/nahuel_site_removal_divergence.csv",
          row.names = FALSE)

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
  ggtitle("Nahuel: signed divergence under random site removal",
          subtitle = "Expected links calculated with raw pairwise interaction probabilities")

ggsave("Nahuel/nahuel_site_removal_relative_divergence_raw.png",
       plot_raw,
       width = 7,
       height = 5,
       dpi = 300)

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
  ggtitle("Nahuel: signed divergence under random site removal")

ggsave("Nahuel/nahuel_site_removal_signed_divergence_raw.png",
       plot_signed,
       width = 7,
       height = 5,
       dpi = 300)

summary_table <- site_removal_results %>%
  group_by(removal_fraction) %>%
  summarise(
    mean_realised_links = mean(n_realised_links, na.rm = TRUE),
    mean_expected_links_raw = mean(n_expected_links_raw, na.rm = TRUE),
    mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
    mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
    sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
    mean_species_bias_raw = mean(mean_species_bias_raw, na.rm = TRUE),
    mean_abs_species_error_raw = mean(mean_abs_species_error_raw, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

write.csv(summary_table,
          file = "Nahuel/nahuel_site_removal_summary.csv",
          row.names = FALSE)
