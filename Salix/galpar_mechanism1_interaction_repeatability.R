## ------------------------------------------------------------
## Script: galpar_mechanism1_interaction_repeatability.R
##
## Purpose:
## Test whether empirical realised interactions are more spatially
## repeated than model-generated interactions.
##
## Outputs saved in:
##   Salix/mechanism1/
##
## Run from the parent repository folder.
## ------------------------------------------------------------


## ---------------------------
## 0. Load packages
## ---------------------------

packages <- c("dplyr", "ggplot2", "magrittr", "reshape2",
              "igraph", "bipartite", "data.table")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


## ---------------------------
## 1. Output folder
## ---------------------------

out_dir <- "Salix/mechanism1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------
## 2. Read Galpar / Salix data
## ---------------------------

unlink("Salix/raw-data/csv", recursive = TRUE)
unlink("Salix/raw-data/rdata", recursive = TRUE)

source("Salix/lib/format4R.r")
get_formatData("Salix/raw-data/Salix_webs.csv")

df_site <- readRDS("Salix/raw-data/rdata/df_site.rds")
df_interact <- readRDS("Salix/raw-data/rdata/df_interact.rds")

df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR / df_interact$N_GALLS

site_interact <- merge(df_site, df_interact, by = "REARING_NUMBER")


## ---------------------------
## 3. Parameters
## ---------------------------

set.seed(123)

p_fixed <- 0.121
n_model_reps <- 1000


## ---------------------------
## 4. Build empirical co-occurrence and interaction tables
## ---------------------------

valid_rpar <- !is.na(site_interact$RPAR) &
  site_interact$RPAR != "" &
  site_interact$RPAR != "none"

valid_rgaller <- !is.na(site_interact$RGALLER) &
  site_interact$RGALLER != ""

consumer_sites <- site_interact[valid_rpar, ] %>%
  transmute(site = SITE,
            consumer = RPAR) %>%
  distinct()

resource_sites <- site_interact[valid_rgaller, ] %>%
  transmute(site = SITE,
            resource = RGALLER) %>%
  distinct()

cooc_triples <- inner_join(consumer_sites,
                           resource_sites,
                           by = "site") %>%
  distinct(site, consumer, resource)

empirical_site_interactions <- site_interact[valid_rpar & valid_rgaller, ] %>%
  transmute(site = SITE,
            consumer = RPAR,
            resource = RGALLER) %>%
  distinct()


## ---------------------------
## 5. Pair-level co-occurrence counts
## ---------------------------

cooc_counts <- cooc_triples %>%
  group_by(consumer, resource) %>%
  summarise(
    n_cooccurring_sites = n_distinct(site),
    .groups = "drop"
  )


## ---------------------------
## 6. Empirical interaction repeatability
## ---------------------------

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
          file = file.path(out_dir, "galpar_empirical_interaction_repeatability.csv"),
          row.names = FALSE)


## ---------------------------
## 7. Model-generated interaction repeatability
## ---------------------------

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
  
  sim_links <- sim %>%
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
  
  return(sim_links)
}


model_repeatability_list <- vector("list", n_model_reps)

for(r in seq_len(n_model_reps)){
  
  message("Galpar model replicate ", r, " / ", n_model_reps)
  
  model_repeatability_list[[r]] <- simulate_model_repeatability(
    cooc_triples = cooc_triples,
    cooc_counts = cooc_counts,
    p_fixed = p_fixed,
    model_rep = r
  )
}

model_repeatability <- bind_rows(model_repeatability_list)

write.csv(model_repeatability,
          file = file.path(out_dir, "galpar_model_generated_interaction_repeatability.csv"),
          row.names = FALSE)


## ---------------------------
## 8. Summary tables
## ---------------------------

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
          file = file.path(out_dir, "galpar_repeatability_summary_by_source.csv"),
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
)

comparison_to_model <- comparison_to_model %>%
  mutate(
    empirical_above_model_975 = empirical_value > model_q975,
    empirical_below_model_025 = empirical_value < model_q025
  )

write.csv(model_summary_by_rep,
          file = file.path(out_dir, "galpar_model_summary_by_rep.csv"),
          row.names = FALSE)

write.csv(comparison_to_model,
          file = file.path(out_dir, "galpar_empirical_vs_model_repeatability_summary.csv"),
          row.names = FALSE)


## ---------------------------
## 9. Plots
## ---------------------------

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
  ggtitle("Galpar: spatial repetition of realised interactions",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "galpar_hist_interacting_sites_empirical_vs_model.png"),
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
  ggtitle("Galpar: interaction repeatability",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "galpar_density_repeatability_empirical_vs_model.png"),
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
  ggtitle("Galpar: interaction spatial repetition",
          subtitle = "Empirical vs model-generated links")

ggsave(file.path(out_dir, "galpar_boxplot_interacting_sites_empirical_vs_model.png"),
       p_box,
       width = 7,
       height = 5.5,
       dpi = 300)


print(comparison_to_model)