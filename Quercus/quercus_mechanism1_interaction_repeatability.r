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
## Script: quercus_mechanism1_interaction_repeatability.r
## ------------------------------------------------------------

set.seed(123)

p_fixed <- 0.0227
n_model_reps <- 1000

out_dir <- "Quercus/mechanism1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

webs <- read_quercus_webs()
all_sites <- names(webs)

site_tables <- build_quercus_site_tables(webs, all_sites)
cooc_triples <- site_tables$cooc_triples
empirical_site_interactions <- site_tables$empirical_site_interactions

cooc_counts <- cooc_triples %>%
  group_by(consumer, resource) %>%
  summarise(
    n_cooccurring_sites = n_distinct(site),
    .groups = "drop"
  )

empirical_repeatability <- empirical_site_interactions %>%
  group_by(consumer, resource) %>%
  summarise(
    n_interacting_sites = n_distinct(site),
    .groups = "drop"
  ) %>%
  left_join(cooc_counts, by = c("consumer", "resource")) %>%
  mutate(
    repeatability = n_interacting_sites / n_cooccurring_sites,
    source = "empirical",
    model_rep = NA_integer_
  )

write.csv(empirical_repeatability,
          file = file.path(out_dir, "quercus_empirical_interaction_repeatability.csv"),
          row.names = FALSE)

simulate_model_repeatability <- function(cooc_triples,
                                         cooc_counts,
                                         p_fixed,
                                         model_rep){

  sim <- cooc_triples

  sim$interaction <- rbinom(
    n = nrow(sim),
    size = 1,
    prob = p_fixed
  )

  sim %>%
    filter(interaction == 1) %>%
    group_by(consumer, resource) %>%
    summarise(
      n_interacting_sites = n_distinct(site),
      .groups = "drop"
    ) %>%
    left_join(cooc_counts, by = c("consumer", "resource")) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      source = "model_generated",
      model_rep = model_rep
    )
}

model_repeatability_list <- vector("list", n_model_reps)

for(r in seq_len(n_model_reps)){

  message("Quercus model replicate ", r, " / ", n_model_reps)

  model_repeatability_list[[r]] <- simulate_model_repeatability(
    cooc_triples = cooc_triples,
    cooc_counts = cooc_counts,
    p_fixed = p_fixed,
    model_rep = r
  )
}

model_repeatability <- bind_rows(model_repeatability_list)

write.csv(model_repeatability,
          file = file.path(out_dir, "quercus_model_generated_interaction_repeatability.csv"),
          row.names = FALSE)

all_repeatability <- bind_rows(
  empirical_repeatability,
  model_repeatability
)

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
          file = file.path(out_dir, "quercus_repeatability_summary_by_source.csv"),
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
  metric = c(
    "mean_interacting_sites",
    "median_interacting_sites",
    "q95_interacting_sites",
    "max_interacting_sites",
    "proportion_single_site_links",
    "mean_repeatability",
    "median_repeatability"
  ),
  empirical_value = c(
    empirical_summary$mean_interacting_sites,
    empirical_summary$median_interacting_sites,
    empirical_summary$q95_interacting_sites,
    empirical_summary$max_interacting_sites,
    empirical_summary$proportion_single_site_links,
    empirical_summary$mean_repeatability,
    empirical_summary$median_repeatability
  ),
  model_q025 = c(
    quantile(model_summary_by_rep$mean_interacting_sites, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$median_interacting_sites, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$q95_interacting_sites, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$max_interacting_sites, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$proportion_single_site_links, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$mean_repeatability, 0.025, na.rm = TRUE),
    quantile(model_summary_by_rep$median_repeatability, 0.025, na.rm = TRUE)
  ),
  model_q500 = c(
    quantile(model_summary_by_rep$mean_interacting_sites, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$median_interacting_sites, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$q95_interacting_sites, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$max_interacting_sites, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$proportion_single_site_links, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$mean_repeatability, 0.5, na.rm = TRUE),
    quantile(model_summary_by_rep$median_repeatability, 0.5, na.rm = TRUE)
  ),
  model_q975 = c(
    quantile(model_summary_by_rep$mean_interacting_sites, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$median_interacting_sites, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$q95_interacting_sites, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$max_interacting_sites, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$proportion_single_site_links, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$mean_repeatability, 0.975, na.rm = TRUE),
    quantile(model_summary_by_rep$median_repeatability, 0.975, na.rm = TRUE)
  )
) %>%
  mutate(
    empirical_above_model_975 = empirical_value > model_q975,
    empirical_below_model_025 = empirical_value < model_q025
  )

write.csv(model_summary_by_rep,
          file = file.path(out_dir, "quercus_model_summary_by_rep.csv"),
          row.names = FALSE)

write.csv(comparison_to_model,
          file = file.path(out_dir, "quercus_empirical_vs_model_repeatability_summary.csv"),
          row.names = FALSE)

plot_data <- all_repeatability %>%
  mutate(
    source = factor(source,
                    levels = c("model_generated", "empirical")),
    n_interacting_sites_capped = pmin(n_interacting_sites, 20)
  )

p_hist <- ggplot(plot_data,
                 aes(x = n_interacting_sites_capped,
                     fill = source)) +
  geom_histogram(position = "identity",
                 alpha = 0.45,
                 bins = 20) +
  theme_classic(base_size = 14) +
  xlab("Number of sites where interaction is observed, capped at 20") +
  ylab("Number of pairwise interactions") +
  ggtitle("Quercus: spatial repetition of realised interactions",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "quercus_hist_interacting_sites_empirical_vs_model.png"),
       p_hist,
       width = 8,
       height = 5.5,
       dpi = 300)

p_density <- ggplot(plot_data,
                    aes(x = repeatability,
                        color = source)) +
  geom_density(linewidth = 1.2, na.rm = TRUE) +
  theme_classic(base_size = 14) +
  xlab("Interaction repeatability: interacting sites / co-occurring sites") +
  ylab("Density") +
  ggtitle("Quercus: interaction repeatability",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "quercus_density_repeatability_empirical_vs_model.png"),
       p_density,
       width = 8,
       height = 5.5,
       dpi = 300)

p_box <- ggplot(plot_data,
                aes(x = source,
                    y = n_interacting_sites)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_classic(base_size = 14) +
  xlab("") +
  ylab("Number of sites where interaction is observed") +
  ggtitle("Quercus: interaction spatial repetition",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "quercus_boxplot_interacting_sites_empirical_vs_model.png"),
       p_box,
       width = 7,
       height = 5.5,
       dpi = 300)

print(comparison_to_model)
