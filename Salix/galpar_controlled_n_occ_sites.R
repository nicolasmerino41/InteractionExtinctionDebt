## ------------------------------------------------------------
## Script: galpar_controlled_n_occ_sites.R
##
## Question:
## Is empirical repeatability still higher than expected after
## controlling for the number of co-occurring sites?
##
## Outputs saved in:
##   Salix/controlled_n_occ_sites/
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

out_dir <- "Salix/controlled_n_occ_sites"
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
min_empirical_links <- 5


## ---------------------------
## 4. Build co-occurrence and empirical interaction tables
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
## 6. Empirical pair-level repeatability
## ---------------------------

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
          file.path(out_dir, "galpar_empirical_pair_repeatability.csv"),
          row.names = FALSE)


## ---------------------------
## 7. Helper functions
## ---------------------------

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


## ---------------------------
## 8. Empirical summaries, exact and binned
## ---------------------------

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


## ---------------------------
## 9. Model simulations
## ---------------------------

model_pairs_list <- vector("list", n_model_reps)

for(r in seq_len(n_model_reps)){
  
  message("Galpar controlled simulation ", r, " / ", n_model_reps)
  
  model_pairs_list[[r]] <- simulate_model_pairs(
    cooc_counts = cooc_counts,
    p_fixed = p_fixed,
    model_rep = r
  )
}

model_pairs <- bind_rows(model_pairs_list)

write.csv(model_pairs,
          file.path(out_dir, "galpar_model_generated_pair_repeatability.csv"),
          row.names = FALSE)


## ---------------------------
## 10. Model summaries, exact and binned
## ---------------------------

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


## ---------------------------
## 11. Empirical vs model comparison
## ---------------------------

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
          file.path(out_dir, "galpar_controlled_exact_n_cooccurring_sites.csv"),
          row.names = FALSE)

write.csv(comparison_binned,
          file.path(out_dir, "galpar_controlled_binned_n_cooccurring_sites.csv"),
          row.names = FALSE)


## ---------------------------
## 12. Plots
## ---------------------------

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
  ggtitle("Galpar: repeatability controlled for co-occurring sites",
          subtitle = "Ribbon = 95% model envelope; points/solid line = empirical")

ggsave(file.path(out_dir, "galpar_exact_controlled_mean_repeatability.png"),
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
  ggtitle("Galpar: interacting sites controlled for co-occurring sites",
          subtitle = "Ribbon = 95% model envelope; points/solid line = empirical")

ggsave(file.path(out_dir, "galpar_exact_controlled_mean_interacting_sites.png"),
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
  ggtitle("Galpar: binned control for co-occurring sites",
          subtitle = "Error bars = 95% model envelope; filled points = empirical")

ggsave(file.path(out_dir, "galpar_binned_controlled_mean_repeatability.png"),
       p_binned_repeatability,
       width = 8,
       height = 5.5,
       dpi = 300)


print(comparison_binned)