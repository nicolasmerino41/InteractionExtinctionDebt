## ------------------------------------------------------------
## Script: galpar_site_removal_exploration.R
##
## Purpose:
## First exploratory analysis of model divergence under random
## removal of spatial units in the Galpar dataset.
##
## Run from the parent folder of the repository.
## Assumes the folder "Salix/" exists.
## ------------------------------------------------------------


## ---------------------------
## 0. Load packages
## ---------------------------

packages <- c("dplyr", "ggplot2", "magrittr", "reshape2", "igraph",
              "bipartite", "AICcmodavg", "broom", "data.table")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ---------------------------
## 1. Read Galpar / Salix data
## ---------------------------
## This follows the original Galpar processing script.

unlink("Salix/raw-data/csv", recursive = TRUE)
unlink("Salix/raw-data/rdata", recursive = TRUE)

source("Salix/lib/format4R.r")
get_formatData("Salix/raw-data/Salix_webs.csv")

df_site <- readRDS("Salix/raw-data/rdata/df_site.rds")
df_interact <- readRDS("Salix/raw-data/rdata/df_interact.rds")

df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR / df_interact$N_GALLS

site_interact <- merge(df_site, df_interact, by = "REARING_NUMBER")


## ---------------------------
## 2. Define helper function
## ---------------------------

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0){
    return(NA_real_)
  } else {
    return((observed - expected) / expected)
  }
}


## ---------------------------
## 3. Function for one site subset
## ---------------------------
## This rebuilds, for a subset of sites:
##   1. realised Galpar interactions
##   2. co-occurrence links
##   3. co-occurrence frequencies
##   4. expected links under the fixed p model
##   5. divergence metrics

run_galpar_subset <- function(site_interact,
                              sites_keep,
                              p_fixed = 0.121,
                              removal_fraction = NA_real_,
                              rep_id = NA_integer_){
  
  sub <- site_interact[site_interact$SITE %in% sites_keep, ]
  
  sub <- sub[!is.na(sub$SITE), ]
  sub <- sub[!is.na(sub$RGALLER), ]
  
  ## Realised galler-parasitoid interactions
  real_pairs <- sub %>%
    filter(!is.na(RPAR),
           RPAR != "none",
           RPAR != "") %>%
    distinct(species = RPAR,
             potential_prey = RGALLER)
  
  n_realised_links <- nrow(real_pairs)
  
  ## If no realised interactions remain, return empty summary
  if(n_realised_links == 0){
    return(data.frame(
      removal_fraction = removal_fraction,
      rep = rep_id,
      n_sites_kept = length(unique(sites_keep)),
      n_predators = 0,
      n_resources = length(unique(sub$RGALLER)),
      n_realised_links = 0,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }
  
  predators <- sort(unique(real_pairs$species))
  
  ## Build predator-resource co-occurrences.
  ## For each parasitoid, find sites where it occurs;
  ## then list all gallers present at those same sites.
  cooc_list <- list()
  
  for(pred in predators){
    
    pred_sites <- unique(sub$SITE[sub$RPAR == pred])
    
    if(length(pred_sites) == 0){
      next
    }
    
    cur <- sub[sub$SITE %in% pred_sites, c("SITE", "RGALLER")]
    cur <- unique(cur)
    
    cur_cooc <- cur %>%
      group_by(RGALLER) %>%
      summarise(n_cooc = n_distinct(SITE), .groups = "drop") %>%
      transmute(
        species = pred,
        potential_prey = RGALLER,
        n_cooc = n_cooc
      )
    
    cooc_list[[pred]] <- cur_cooc
  }
  
  cooc <- bind_rows(cooc_list)
  
  cooc <- cooc %>%
    filter(!is.na(potential_prey),
           potential_prey != "")
  
  n_cooc_links <- nrow(cooc)
  
  if(n_cooc_links == 0){
    return(data.frame(
      removal_fraction = removal_fraction,
      rep = rep_id,
      n_sites_kept = length(unique(sites_keep)),
      n_predators = length(predators),
      n_resources = length(unique(sub$RGALLER)),
      n_realised_links = n_realised_links,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }
  
  ## N_alpha = total co-occurrence frequency per predator
  cooc <- cooc %>%
    group_by(species) %>%
    mutate(
      N_alpha = sum(n_cooc),
      total_potential_interactions = n_distinct(potential_prey)
    ) %>%
    ungroup()
  
  ## Probability of interaction from the original model
  cooc <- cooc %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )
  
  ## Expected total links
  n_expected_links_raw <- sum(cooc$prob_int, na.rm = TRUE)
  n_expected_links_conditioned <- sum(cooc$prob_expected, na.rm = TRUE)
  
  ## Species-level expected indegree
  expected_degree <- cooc %>%
    group_by(species) %>%
    summarise(
      expected_indegree_raw = sum(prob_int, na.rm = TRUE),
      expected_indegree_conditioned = sum(prob_expected, na.rm = TRUE),
      potential_degree = n_distinct(potential_prey),
      .groups = "drop"
    )
  
  observed_degree <- real_pairs %>%
    group_by(species) %>%
    summarise(
      observed_indegree = n_distinct(potential_prey),
      .groups = "drop"
    )
  
  degree_compare <- merge(expected_degree,
                          observed_degree,
                          by = "species",
                          all.x = TRUE)
  
  degree_compare$observed_indegree[is.na(degree_compare$observed_indegree)] <- 0
  
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
  
  ## Divergence
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
  
  ## Return one-row summary
  data.frame(
    removal_fraction = removal_fraction,
    rep = rep_id,
    n_sites_kept = length(unique(sites_keep)),
    n_predators = length(unique(real_pairs$species)),
    n_resources = length(unique(real_pairs$potential_prey)),
    n_realised_links = n_realised_links,
    n_cooc_links = n_cooc_links,
    n_expected_links_raw = n_expected_links_raw,
    n_expected_links_conditioned = n_expected_links_conditioned,
    divergence_raw = divergence_raw,
    relative_divergence_raw = relative_divergence_raw,
    divergence_conditioned = divergence_conditioned,
    relative_divergence_conditioned = relative_divergence_conditioned,
    mean_abs_species_error_raw = mean_abs_species_error_raw,
    mean_abs_species_error_conditioned = mean_abs_species_error_conditioned
  )
}


## ---------------------------
## 4. Random site-removal experiment
## ---------------------------

set.seed(123)

all_sites <- sort(unique(site_interact$SITE))

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 100

p_fixed_galpar <- 0.121

results_list <- list()
counter <- 1

for(removal in removal_levels){
  
  n_sites_total <- length(all_sites)
  n_sites_keep <- max(1, round(n_sites_total * (1 - removal)))
  
  ## For removal = 0, one replicate is enough because all sites are kept.
  reps_here <- ifelse(removal == 0, 1, n_replicates)
  
  for(r in seq_len(reps_here)){
    
    if(removal == 0){
      sites_keep <- all_sites
    } else {
      sites_keep <- sample(all_sites,
                           size = n_sites_keep,
                           replace = FALSE)
    }
    
    results_list[[counter]] <- run_galpar_subset(
      site_interact = site_interact,
      sites_keep = sites_keep,
      p_fixed = p_fixed_galpar,
      removal_fraction = removal,
      rep_id = r
    )
    
    counter <- counter + 1
  }
}

site_removal_results <- bind_rows(results_list)


## ---------------------------
## 5. Save results
## ---------------------------

write.csv(site_removal_results,
          file = "Salix/galpar_site_removal_divergence.csv",
          row.names = FALSE)


## ---------------------------
## 6. Plot total-link divergence
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
  ylab("Relative divergence: observed - expected / expected") +
  ggtitle("Galpar: divergence under random site removal",
          subtitle = "Expected links calculated with raw pairwise interaction probabilities")

ggsave("Salix/galpar_site_removal_relative_divergence_raw.png",
       plot_raw,
       width = 7,
       height = 5,
       dpi = 300)


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
  ylab("Relative divergence: observed - expected / expected") +
  ggtitle("Galpar: divergence under random site removal",
          subtitle = "Expected links calculated with conditioned probabilities")

ggsave("Salix/galpar_site_removal_relative_divergence_conditioned.png",
       plot_conditioned,
       width = 7,
       height = 5,
       dpi = 300)


## ---------------------------
## 7. Plot absolute divergence
## ---------------------------

plot_abs <- ggplot(site_removal_results,
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
  ggtitle("Galpar: absolute divergence under random site removal")

ggsave("Salix/galpar_site_removal_absolute_divergence_raw.png",
       plot_abs,
       width = 7,
       height = 5,
       dpi = 300)


## ---------------------------
## 8. Quick console summary
## ---------------------------

summary_table <- site_removal_results %>%
  group_by(removal_fraction) %>%
  summarise(
    mean_realised_links = mean(n_realised_links, na.rm = TRUE),
    mean_expected_links_raw = mean(n_expected_links_raw, na.rm = TRUE),
    mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
    sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
    mean_abs_species_error_raw = mean(mean_abs_species_error_raw, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

write.csv(summary_table,
          file = "Salix/galpar_site_removal_summary.csv",
          row.names = FALSE)