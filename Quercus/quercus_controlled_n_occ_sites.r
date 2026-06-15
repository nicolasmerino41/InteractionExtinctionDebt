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
## Script: quercus_controlled_n_occ_sites.r
##
## Question:
## Is empirical repeatability still higher than expected after
## controlling for the number of co-occurring sites?
## ------------------------------------------------------------

set.seed(123)

p_fixed <- 0.0227
n_model_reps <- 1000
min_empirical_links <- 5

out_dir <- "Quercus/controlled_n_occ_sites"
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

empirical_pairs <- empirical_site_interactions %>%
  group_by(consumer, resource) %>%
  summarise(
    n_interacting_sites = n_distinct(site),
    .groups = "drop"
  ) %>%
  left_join(cooc_counts, by = c("consumer", "resource")) %>%
  mutate(
    repeatability = n_interacting_sites / n_cooccurring_sites
  )

write.csv(empirical_pairs,
          file.path(out_dir, "quercus_empirical_pair_repeatability.csv"),
          row.names = FALSE)

make_nocc_bin <- function(x){
  cut(
    x,
    breaks = c(0, 1, 2, 5, 10, 20, Inf),
    labels = c("1", "2", "3-5", "6-10", "11-20", "21+"),
    right = TRUE
  )
}

summarise_repeatability <- function(df, group_var){

  df %>%
    group_by({{ group_var }}) %>%
    summarise(
      n_links = n(),
      mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
      median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
      q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
      mean_repeatability = mean(repeatability, na.rm = TRUE),
      median_repeatability = median(repeatability, na.rm = TRUE),
      proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
      .groups = "drop"
    )
}

simulate_model_pairs <- function(cooc_counts, p_fixed, model_rep){

  simulated <- cooc_counts

  simulated$n_interacting_sites <- rbinom(
    n = nrow(simulated),
    size = simulated$n_cooccurring_sites,
    prob = p_fixed
  )

  simulated %>%
    filter(n_interacting_sites > 0) %>%
    mutate(
      repeatability = n_interacting_sites / n_cooccurring_sites,
      model_rep = model_rep,
      cooc_bin = make_nocc_bin(n_cooccurring_sites)
    )
}

empirical_pairs <- empirical_pairs %>%
  mutate(cooc_bin = make_nocc_bin(n_cooccurring_sites))

empirical_exact <- summarise_repeatability(
  empirical_pairs,
  n_cooccurring_sites
) %>%
  rename(empirical_n_links = n_links,
         empirical_mean_interacting_sites = mean_interacting_sites,
         empirical_median_interacting_sites = median_interacting_sites,
         empirical_q95_interacting_sites = q95_interacting_sites,
         empirical_mean_repeatability = mean_repeatability,
         empirical_median_repeatability = median_repeatability,
         empirical_proportion_single_site_links = proportion_single_site_links)

empirical_binned <- summarise_repeatability(
  empirical_pairs,
  cooc_bin
) %>%
  rename(empirical_n_links = n_links,
         empirical_mean_interacting_sites = mean_interacting_sites,
         empirical_median_interacting_sites = median_interacting_sites,
         empirical_q95_interacting_sites = q95_interacting_sites,
         empirical_mean_repeatability = mean_repeatability,
         empirical_median_repeatability = median_repeatability,
         empirical_proportion_single_site_links = proportion_single_site_links)

model_pairs_list <- vector("list", n_model_reps)

for(r in seq_len(n_model_reps)){

  message("Quercus controlled simulation ", r, " / ", n_model_reps)

  model_pairs_list[[r]] <- simulate_model_pairs(
    cooc_counts = cooc_counts,
    p_fixed = p_fixed,
    model_rep = r
  )
}

model_pairs <- bind_rows(model_pairs_list)

write.csv(model_pairs,
          file.path(out_dir, "quercus_model_generated_pair_repeatability.csv"),
          row.names = FALSE)

model_exact_by_rep <- model_pairs %>%
  group_by(model_rep, n_cooccurring_sites) %>%
  summarise(
    n_links = n(),
    mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
    median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE),
    proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
    .groups = "drop"
  )

model_binned_by_rep <- model_pairs %>%
  group_by(model_rep, cooc_bin) %>%
  summarise(
    n_links = n(),
    mean_interacting_sites = mean(n_interacting_sites, na.rm = TRUE),
    median_interacting_sites = median(n_interacting_sites, na.rm = TRUE),
    q95_interacting_sites = quantile(n_interacting_sites, 0.95, na.rm = TRUE),
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE),
    proportion_single_site_links = mean(n_interacting_sites == 1, na.rm = TRUE),
    .groups = "drop"
  )

make_envelope <- function(model_summary, group_var){

  model_summary %>%
    group_by({{ group_var }}) %>%
    summarise(
      model_n_links_median = median(n_links, na.rm = TRUE),

      model_mean_interacting_sites_q025 = quantile(mean_interacting_sites, 0.025, na.rm = TRUE),
      model_mean_interacting_sites_q500 = quantile(mean_interacting_sites, 0.500, na.rm = TRUE),
      model_mean_interacting_sites_q975 = quantile(mean_interacting_sites, 0.975, na.rm = TRUE),

      model_mean_repeatability_q025 = quantile(mean_repeatability, 0.025, na.rm = TRUE),
      model_mean_repeatability_q500 = quantile(mean_repeatability, 0.500, na.rm = TRUE),
      model_mean_repeatability_q975 = quantile(mean_repeatability, 0.975, na.rm = TRUE),

      model_prop_single_q025 = quantile(proportion_single_site_links, 0.025, na.rm = TRUE),
      model_prop_single_q500 = quantile(proportion_single_site_links, 0.500, na.rm = TRUE),
      model_prop_single_q975 = quantile(proportion_single_site_links, 0.975, na.rm = TRUE),

      .groups = "drop"
    )
}

model_exact_envelope <- make_envelope(model_exact_by_rep, n_cooccurring_sites)
model_binned_envelope <- make_envelope(model_binned_by_rep, cooc_bin)

comparison_exact <- empirical_exact %>%
  left_join(model_exact_envelope, by = "n_cooccurring_sites") %>%
  mutate(
    enough_empirical_links = empirical_n_links >= min_empirical_links,

    empirical_repeatability_above_model =
      empirical_mean_repeatability > model_mean_repeatability_q975,

    empirical_repeatability_below_model =
      empirical_mean_repeatability < model_mean_repeatability_q025,

    empirical_interacting_sites_above_model =
      empirical_mean_interacting_sites > model_mean_interacting_sites_q975,

    empirical_interacting_sites_below_model =
      empirical_mean_interacting_sites < model_mean_interacting_sites_q025,

    empirical_single_site_below_model =
      empirical_proportion_single_site_links < model_prop_single_q025,

    empirical_single_site_above_model =
      empirical_proportion_single_site_links > model_prop_single_q975
  )

comparison_binned <- empirical_binned %>%
  left_join(model_binned_envelope, by = "cooc_bin") %>%
  mutate(
    enough_empirical_links = empirical_n_links >= min_empirical_links,

    empirical_repeatability_above_model =
      empirical_mean_repeatability > model_mean_repeatability_q975,

    empirical_repeatability_below_model =
      empirical_mean_repeatability < model_mean_repeatability_q025,

    empirical_interacting_sites_above_model =
      empirical_mean_interacting_sites > model_mean_interacting_sites_q975,

    empirical_interacting_sites_below_model =
      empirical_mean_interacting_sites < model_mean_interacting_sites_q025,

    empirical_single_site_below_model =
      empirical_proportion_single_site_links < model_prop_single_q025,

    empirical_single_site_above_model =
      empirical_proportion_single_site_links > model_prop_single_q975
  )

write.csv(comparison_exact,
          file.path(out_dir, "quercus_controlled_exact_n_cooccurring_sites.csv"),
          row.names = FALSE)

write.csv(comparison_binned,
          file.path(out_dir, "quercus_controlled_binned_n_cooccurring_sites.csv"),
          row.names = FALSE)

plot_exact_data <- comparison_exact %>%
  filter(enough_empirical_links)

p_exact_repeatability <- ggplot(plot_exact_data,
                                aes(x = n_cooccurring_sites)) +
  geom_ribbon(aes(ymin = model_mean_repeatability_q025,
                  ymax = model_mean_repeatability_q975),
              fill = "grey80") +
  geom_line(aes(y = model_mean_repeatability_q500),
            linetype = 2,
            linewidth = 1) +
  geom_point(aes(y = empirical_mean_repeatability,
                 size = empirical_n_links)) +
  geom_line(aes(y = empirical_mean_repeatability),
            linewidth = 1.1) +
  theme_classic(base_size = 14) +
  xlab("Number of co-occurring sites") +
  ylab("Mean repeatability") +
  ggtitle("Quercus: repeatability controlled for co-occurring sites",
          subtitle = "Ribbon = 95% model envelope; points/solid line = empirical")

ggsave(file.path(out_dir, "quercus_exact_controlled_mean_repeatability.png"),
       p_exact_repeatability,
       width = 8,
       height = 5.5,
       dpi = 300)

p_exact_interacting <- ggplot(plot_exact_data,
                              aes(x = n_cooccurring_sites)) +
  geom_ribbon(aes(ymin = model_mean_interacting_sites_q025,
                  ymax = model_mean_interacting_sites_q975),
              fill = "grey80") +
  geom_line(aes(y = model_mean_interacting_sites_q500),
            linetype = 2,
            linewidth = 1) +
  geom_point(aes(y = empirical_mean_interacting_sites,
                 size = empirical_n_links)) +
  geom_line(aes(y = empirical_mean_interacting_sites),
            linewidth = 1.1) +
  theme_classic(base_size = 14) +
  xlab("Number of co-occurring sites") +
  ylab("Mean number of interacting sites") +
  ggtitle("Quercus: interacting sites controlled for co-occurring sites",
          subtitle = "Ribbon = 95% model envelope; points/solid line = empirical")

ggsave(file.path(out_dir, "quercus_exact_controlled_mean_interacting_sites.png"),
       p_exact_interacting,
       width = 8,
       height = 5.5,
       dpi = 300)

p_binned_repeatability <- ggplot(comparison_binned,
                                 aes(x = cooc_bin)) +
  geom_errorbar(aes(ymin = model_mean_repeatability_q025,
                    ymax = model_mean_repeatability_q975),
                width = 0.15,
                linewidth = 1) +
  geom_point(aes(y = model_mean_repeatability_q500),
             shape = 1,
             size = 3) +
  geom_point(aes(y = empirical_mean_repeatability,
                 size = empirical_n_links),
             size = 3) +
  theme_classic(base_size = 14) +
  xlab("Number of co-occurring sites, binned") +
  ylab("Mean repeatability") +
  ggtitle("Quercus: binned control for co-occurring sites",
          subtitle = "Error bars = 95% model envelope; filled points = empirical")

ggsave(file.path(out_dir, "quercus_binned_controlled_mean_repeatability.png"),
       p_binned_repeatability,
       width = 8,
       height = 5.5,
       dpi = 300)

print(comparison_binned)
