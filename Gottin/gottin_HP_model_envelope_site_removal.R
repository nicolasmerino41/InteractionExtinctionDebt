## ------------------------------------------------------------
## Script: gottin_HP_model_envelope_site_removal.R
##
## Purpose:
## Compare empirical Gottin-HP site-removal divergence against
## divergence produced by model-generated site-level interaction
## networks.
##
## Run from the parent folder of the repository.
## Assumes "Gottin/" exists.
## ------------------------------------------------------------


## ---------------------------
## 0. Load packages
## ---------------------------

packages <- c("dplyr", "ggplot2", "bipartite", "igraph")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if(file.exists("utils.r")){
  source("utils.r")
}


## ---------------------------
## 1. Read Gottin-HP data
## ---------------------------

ant_raw <- read.csv("Gottin/raw-data/host_para_all_interactions.csv",
                    head = TRUE,
                    sep = ";")

metadata <- read.csv("Gottin/metadata.csv",
                     head = TRUE,
                     sep = ",")

ant_nets <- frame2webs(
  ant_raw,
  varnames = c("Genus.Species",
               "P1.Species.Genus",
               "Site",
               "P1.cells")
)


## ---------------------------
## 2. Parameters
## ---------------------------

set.seed(123)

p_fixed <- 0.074

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)

n_site_reps <- 100

n_model_reps <- 100

output_prefix <- "gottin_HP"


## ---------------------------
## 3. Helper functions
## ---------------------------

clean_web_matrix <- function(web){
  
  web <- as.data.frame(web)
  
  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)){
    web <- web[, !bad_cols, drop = FALSE]
  }
  
  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)){
    web <- web[!bad_rows, , drop = FALSE]
  }
  
  as.matrix(web)
}


get_site_web <- function(webs, site_id){
  
  if(!is.null(names(webs))){
    if(as.character(site_id) %in% names(webs)){
      return(clean_web_matrix(webs[[as.character(site_id)]]))
    }
  }
  
  site_num <- suppressWarnings(as.numeric(site_id))
  
  if(!is.na(site_num) &&
     site_num >= 1 &&
     site_num <= length(webs)){
    return(clean_web_matrix(webs[[site_num]]))
  }
  
  stop(paste("Could not find web for site:", site_id))
}


safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0){
    return(NA_real_)
  } else {
    return((observed - expected) / expected)
  }
}


build_gottin_HP_site_tables <- function(webs, all_sites){
  
  cooc_list <- list()
  interaction_list <- list()
  counter_cooc <- 1
  counter_int <- 1
  
  for(site_id in all_sites){
    
    web <- get_site_web(webs, site_id)
    
    if(nrow(web) == 0 || ncol(web) == 0){
      next
    }
    
    consumers <- colnames(web)
    resources <- rownames(web)
    
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    
    if(length(consumers) == 0 || length(resources) == 0){
      next
    }
    
    cur_cooc <- expand.grid(
      consumer = consumers,
      resource = resources,
      stringsAsFactors = FALSE
    )
    
    cur_cooc$site <- site_id
    
    cooc_list[[counter_cooc]] <- cur_cooc[, c("site",
                                              "consumer",
                                              "resource")]
    counter_cooc <- counter_cooc + 1
    
    positive_cells <- which(web > 0, arr.ind = TRUE)
    
    if(nrow(positive_cells) > 0){
      
      cur_int <- data.frame(
        site = site_id,
        consumer = colnames(web)[positive_cells[, 2]],
        resource = rownames(web)[positive_cells[, 1]]
      )
      
      cur_int <- cur_int %>%
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
    cooc_triples = bind_rows(cooc_list) %>%
      distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>%
      distinct(site, consumer, resource)
  )
}


make_site_subsets <- function(all_sites,
                              removal_levels,
                              n_site_reps){
  
  site_subsets <- list()
  subset_index <- NULL
  counter <- 1
  
  for(removal in removal_levels){
    
    n_total <- length(all_sites)
    n_keep <- max(1, round(n_total * (1 - removal)))
    
    reps_here <- ifelse(removal == 0, 1, n_site_reps)
    
    for(r in seq_len(reps_here)){
      
      if(removal == 0){
        sites_keep <- all_sites
      } else {
        sites_keep <- sample(all_sites,
                             size = n_keep,
                             replace = FALSE)
      }
      
      site_subsets[[counter]] <- sites_keep
      
      cur <- data.frame(
        subset_id = counter,
        removal_fraction = removal,
        site_rep = r,
        n_sites_kept = length(sites_keep)
      )
      
      if(is.null(subset_index)){
        subset_index <- cur
      } else {
        subset_index <- rbind(subset_index, cur)
      }
      
      counter <- counter + 1
    }
  }
  
  list(index = subset_index,
       subsets = site_subsets)
}


expected_links_from_subset <- function(cooc_triples,
                                       sites_keep,
                                       p_fixed){
  
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


observed_links_from_subset <- function(site_interactions,
                                       sites_keep){
  
  cur <- site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource)
  
  nrow(cur)
}


simulate_model_site_interactions <- function(cooc_triples,
                                             p_fixed){
  
  sim <- cooc_triples
  sim$interaction <- rbinom(n = nrow(sim),
                            size = 1,
                            prob = p_fixed)
  
  sim <- sim %>%
    filter(interaction == 1) %>%
    select(site, consumer, resource) %>%
    distinct()
  
  sim
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
    
    observed_links <- observed_links_from_subset(site_interactions,
                                                 sites_keep)
    
    divergence_raw <- observed_links - exp_row$expected_links_raw
    relative_divergence_raw <- safe_relative_difference(
      observed_links,
      exp_row$expected_links_raw
    )
    
    divergence_conditioned <- observed_links -
      exp_row$expected_links_conditioned
    
    relative_divergence_conditioned <- safe_relative_difference(
      observed_links,
      exp_row$expected_links_conditioned
    )
    
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
      divergence_raw = divergence_raw,
      relative_divergence_raw = relative_divergence_raw,
      divergence_conditioned = divergence_conditioned,
      relative_divergence_conditioned = relative_divergence_conditioned
    )
  }
  
  bind_rows(out)
}


## ---------------------------
## 4. Build tables
## ---------------------------

all_sites <- sort(unique(metadata$site))

site_tables <- build_gottin_HP_site_tables(
  webs = ant_nets,
  all_sites = all_sites
)

cooc_triples <- site_tables$cooc_triples

empirical_site_interactions <- site_tables$empirical_site_interactions

all_sites <- sort(unique(cooc_triples$site))


## ---------------------------
## 5. Create shared site-removal subsets
## ---------------------------

subset_object <- make_site_subsets(
  all_sites = all_sites,
  removal_levels = removal_levels,
  n_site_reps = n_site_reps
)


## ---------------------------
## 6. Compute model expectation for each subset
## ---------------------------

subset_expected_list <- vector("list", nrow(subset_object$index))

for(i in seq_len(nrow(subset_object$index))){
  
  subset_id <- subset_object$index$subset_id[i]
  sites_keep <- subset_object$subsets[[subset_id]]
  
  cur_exp <- expected_links_from_subset(
    cooc_triples = cooc_triples,
    sites_keep = sites_keep,
    p_fixed = p_fixed
  )
  
  cur_exp$subset_id <- subset_id
  
  subset_expected_list[[i]] <- cur_exp
}

subset_expected <- bind_rows(subset_expected_list)


## ---------------------------
## 7. Empirical divergence
## ---------------------------

empirical_results <- evaluate_network_against_subsets(
  site_interactions = empirical_site_interactions,
  subset_object = subset_object,
  subset_expected = subset_expected,
  source = "empirical",
  model_rep = NA_integer_
)


## ---------------------------
## 8. Model-generated divergence envelope
## ---------------------------

model_results_list <- vector("list", n_model_reps)

for(m in seq_len(n_model_reps)){
  
  message("Running model replicate ", m, " / ", n_model_reps)
  
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


## ---------------------------
## 9. Summaries
## ---------------------------

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
      empirical_mean_relative_divergence_raw >
      model_q975_relative_divergence_raw,
    empirical_below_model_025 =
      empirical_mean_relative_divergence_raw <
      model_q025_relative_divergence_raw
  )


## ---------------------------
## 10. Save outputs
## ---------------------------

write.csv(empirical_results,
          paste0("Gottin/envelope/",output_prefix, "_empirical_divergence_by_subset.csv"),
          row.names = FALSE)

write.csv(model_results,
          paste0("Gottin/envelope/",output_prefix, "_model_generated_divergence_by_subset.csv"),
          row.names = FALSE)

write.csv(comparison_summary,
          paste0("Gottin/envelope/",output_prefix, "_empirical_vs_model_envelope_summary.csv"),
          row.names = FALSE)


## ---------------------------
## 11. Plot
## ---------------------------

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
    "Gottin-HP: empirical divergence vs model-generated envelope",
    subtitle = "Ribbon = 95% envelope of model-generated mean curves; solid line = empirical"
  )

ggsave(paste0("Gottin/envelope/",output_prefix, "_empirical_vs_model_envelope.png"),
       p,
       width = 8,
       height = 5.5,
       dpi = 300)

print(comparison_summary)