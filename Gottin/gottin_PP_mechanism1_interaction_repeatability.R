## ------------------------------------------------------------
## Script: gottin_PP_mechanism1_interaction_repeatability.R
## ------------------------------------------------------------

packages <- c("dplyr", "ggplot2", "bipartite", "igraph")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "Gottin/mechanism1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mut_raw <- read.csv("Gottin/raw-data/plant_poll_all_interactions.csv",
                    head = TRUE, sep = ";")

metadata <- read.csv("Gottin/metadata.csv",
                     head = TRUE, sep = ",")

mut_nets <- frame2webs(
  mut_raw,
  varnames = c("P.Genus.Species",
               "Genus.Species",
               "Site")
)

set.seed(123)

p_fixed <- 0.044  # change if you have calibrated PP value
n_model_reps <- 1000

clean_web_matrix <- function(web){
  web <- as.data.frame(web)
  
  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)) web <- web[, !bad_cols, drop = FALSE]
  
  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)) web <- web[!bad_rows, , drop = FALSE]
  
  as.matrix(web)
}

get_site_web <- function(webs, site_id){
  if(!is.null(names(webs))){
    if(as.character(site_id) %in% names(webs)){
      return(clean_web_matrix(webs[[as.character(site_id)]]))
    }
  }
  
  site_num <- suppressWarnings(as.numeric(site_id))
  
  if(!is.na(site_num) && site_num >= 1 && site_num <= length(webs)){
    return(clean_web_matrix(webs[[site_num]]))
  }
  
  stop(paste("Could not find web for site:", site_id))
}

build_gottin_PP_site_tables <- function(webs, all_sites){
  
  cooc_list <- list()
  interaction_list <- list()
  counter_cooc <- 1
  counter_int <- 1
  
  for(site_id in all_sites){
    
    web <- get_site_web(webs, site_id)
    
    if(nrow(web) == 0 || ncol(web) == 0) next
    
    consumers <- colnames(web)
    resources <- rownames(web)
    
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    
    if(length(consumers) == 0 || length(resources) == 0) next
    
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
        filter(!is.na(consumer), !is.na(resource),
               consumer != "", resource != "") %>%
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

all_sites <- sort(unique(metadata$site))

site_tables <- build_gottin_PP_site_tables(
  webs = mut_nets,
  all_sites = all_sites
)

cooc_triples <- site_tables$cooc_triples
empirical_site_interactions <- site_tables$empirical_site_interactions

cooc_counts <- cooc_triples %>%
  group_by(consumer, resource) %>%
  summarise(n_cooccurring_sites = n_distinct(site), .groups = "drop")

empirical_repeatability <- empirical_site_interactions %>%
  group_by(consumer, resource) %>%
  summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
  left_join(cooc_counts, by = c("consumer", "resource")) %>%
  mutate(
    repeatability = n_interacting_sites / n_cooccurring_sites,
    source = "empirical",
    model_rep = NA_integer_
  )

write.csv(empirical_repeatability,
          file.path(out_dir, "gottin_PP_empirical_interaction_repeatability.csv"),
          row.names = FALSE)

simulate_model_repeatability <- function(cooc_triples, cooc_counts, p_fixed, model_rep){
  
  sim <- cooc_triples
  sim$interaction <- rbinom(nrow(sim), size = 1, prob = p_fixed)
  
  sim %>%
    filter(interaction == 1) %>%
    group_by(consumer, resource) %>%
    summarise(n_interacting_sites = n_distinct(site), .groups = "drop") %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      source = "model_generated",
      model_rep = model_rep
    )
}

model_repeatability_list <- vector("list", n_model_reps)

for(r in seq_len(n_model_reps)){
  message("Gottin-PP model replicate ", r, " / ", n_model_reps)
  
  model_repeatability_list[[r]] <- simulate_model_repeatability(
    cooc_triples = cooc_triples,
    cooc_counts = cooc_counts,
    p_fixed = p_fixed,
    model_rep = r
  )
}

model_repeatability <- bind_rows(model_repeatability_list)

write.csv(model_repeatability,
          file.path(out_dir, "gottin_PP_model_generated_interaction_repeatability.csv"),
          row.names = FALSE)

all_repeatability <- bind_rows(empirical_repeatability, model_repeatability)

summary_by_source <- all_repeatability %>%
  group_by(source) %>%
  summarise(
    n_links = n(),
    mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
    median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
    max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
    proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_by_source,
          file.path(out_dir, "gottin_PP_repeatability_summary_by_source.csv"),
          row.names = FALSE)

model_summary_by_rep <- model_repeatability %>%
  group_by(model_rep) %>%
  summarise(
    n_links = n(),
    mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
    median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
    max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
    proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE),
    .groups = "drop"
  )

empirical_summary <- empirical_repeatability %>%
  summarise(
    n_links = n(),
    mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
    median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
    max_interacting_sites = max(n_interacting_sites, na.rm = TRUE),
    proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE)
  )

comparison_to_model <- data.frame(
  metric = c("mean_interacting_sites", "median_interacting_sites",
             "q95_interacting_sites", "max_interacting_sites",
             "proportion_single_site_links", "mean_repeatability",
             "median_repeatability"),
  empirical_value = c(empirical_summary$mean_interacting_sites,
                      empirical_summary$median_interacting_sites,
                      empirical_summary$q95_interacting_sites,
                      empirical_summary$max_interacting_sites,
                      empirical_summary$proportion_single_site_links,
                      empirical_summary$mean_repeatability,
                      empirical_summary$median_repeatability),
  model_q025 = c(quantile(model_summary_by_rep$mean_interacting_sites, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_interacting_sites, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$q95_interacting_sites, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$max_interacting_sites, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$proportion_single_site_links, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$mean_repeatability, 0.025, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_repeatability, 0.025, na.rm = TRUE)),
  model_q500 = c(quantile(model_summary_by_rep$mean_interacting_sites, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_interacting_sites, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$q95_interacting_sites, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$max_interacting_sites, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$proportion_single_site_links, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$mean_repeatability, 0.5, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_repeatability, 0.5, na.rm = TRUE)),
  model_q975 = c(quantile(model_summary_by_rep$mean_interacting_sites, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_interacting_sites, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$q95_interacting_sites, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$max_interacting_sites, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$proportion_single_site_links, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$mean_repeatability, 0.975, na.rm = TRUE),
                 quantile(model_summary_by_rep$median_repeatability, 0.975, na.rm = TRUE))
) %>%
  mutate(
    empirical_above_model_975 = empirical_value > model_q975,
    empirical_below_model_025 = empirical_value < model_q025
  )

write.csv(model_summary_by_rep,
          file.path(out_dir, "gottin_PP_model_summary_by_rep.csv"),
          row.names = FALSE)

write.csv(comparison_to_model,
          file.path(out_dir, "gottin_PP_empirical_vs_model_repeatability_summary.csv"),
          row.names = FALSE)

plot_data <- all_repeatability %>%
  mutate(
    source = factor(source, levels = c("model_generated", "empirical")),
    n_interacting_sites_capped = pmin(n_interacting_sites, 20)
  )

p_hist <- ggplot(plot_data, aes(x = n_interacting_sites_capped, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 20) +
  theme_classic(base_size = 14) +
  xlab("Number of sites where interaction is observed, capped at 20") +
  ylab("Number of pairwise interactions") +
  ggtitle("Gottin-PP: spatial repetition of realised interactions",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "gottin_PP_hist_interacting_sites_empirical_vs_model.png"),
       p_hist, width = 8, height = 5.5, dpi = 300)

p_density <- ggplot(plot_data, aes(x = repeatability, color = source)) +
  geom_density(linewidth = 1.2, na.rm = TRUE) +
  theme_classic(base_size = 14) +
  xlab("Interaction repeatability: interacting sites / co-occurring sites") +
  ylab("Density") +
  ggtitle("Gottin-PP: interaction repeatability",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "gottin_PP_density_repeatability_empirical_vs_model.png"),
       p_density, width = 8, height = 5.5, dpi = 300)

p_box <- ggplot(plot_data, aes(x = source, y = n_interacting_sites)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_classic(base_size = 14) +
  xlab("") +
  ylab("Number of sites where interaction is observed") +
  ggtitle("Gottin-PP: interaction spatial repetition",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "gottin_PP_boxplot_interacting_sites_empirical_vs_model.png"),
       p_box, width = 7, height = 5.5, dpi = 300)

print(comparison_to_model)

