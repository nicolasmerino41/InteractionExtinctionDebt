## ------------------------------------------------------------
## Script: All/scripts/31_figure2_empirical_support_heterogeneity.R
## Empirical support heterogeneity.
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "patchwork", "scales")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, name, width = 12, height = 8){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

make_pair_support <- function(dataset){
  message("Loading pair support: ", dataset)
  site_tables <- get_dataset_site_tables(dataset)
  cooc <- site_tables$cooc_triples %>% distinct(site, consumer, resource)
  ints <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource)
  cooc_counts <- cooc %>% group_by(consumer, resource) %>% summarise(full_n = n_distinct(site), .groups = "drop")
  int_counts <- ints %>% group_by(consumer, resource) %>% summarise(full_K = n_distinct(site), .groups = "drop")
  cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(dataset = dataset, full_K = replace_na(full_K, 0L), full_B = full_K > 0, interaction_occupancy = full_K / full_n) %>%
    select(dataset, consumer, resource, full_n, full_K, full_B, interaction_occupancy)
}

all_pairs <- bind_rows(lapply(all_dataset_names, make_pair_support))
write.csv2(all_pairs, file.path(out_dir, "31_all_cooccurrence_pair_support.csv"), row.names = FALSE)

support_class <- function(K){
  case_when(
    K == 1 ~ "Observed in 1 site",
    K == 2 ~ "Observed in 2 sites",
    K >= 3 & K <= 4 ~ "Observed in 3–4 sites",
    K >= 5 & K <= 9 ~ "Observed in 5–9 sites",
    K >= 10 ~ "Observed in 10 or more sites",
    TRUE ~ NA_character_
  )
}

support_levels <- c("Observed in 1 site", "Observed in 2 sites", "Observed in 3–4 sites", "Observed in 5–9 sites", "Observed in 10 or more sites")
support_cols <- c("Observed in 1 site" = "#D55E00", "Observed in 2 sites" = "#E69F00", "Observed in 3–4 sites" = "#56B4E9", "Observed in 5–9 sites" = "#0072B2", "Observed in 10 or more sites" = "#004C6D")

links <- all_pairs %>% filter(full_K >= 1) %>% mutate(support_class = factor(support_class(full_K), levels = support_levels))

dataset_support <- links %>%
  count(dataset, support_class, name = "n_links") %>%
  group_by(dataset) %>% mutate(prop = n_links / sum(n_links)) %>% ungroup()

across_support <- dataset_support %>%
  group_by(support_class) %>% summarise(prop = mean(prop, na.rm = TRUE), n_links = NA_integer_, .groups = "drop") %>%
  mutate(dataset = "Across landscapes")

support_dist <- bind_rows(dataset_support, across_support)
dataset_order <- dataset_support %>% filter(support_class == "Observed in 1 site") %>% arrange(desc(prop)) %>% pull(dataset)
support_dist <- support_dist %>% mutate(dataset = factor(dataset, levels = rev(c("Across landscapes", dataset_order))))
write.csv2(support_dist, file.path(out_dir, "31_support_distribution_by_dataset.csv"), row.names = FALSE)

panel_a <- ggplot(support_dist, aes(x = dataset, y = prop, fill = support_class)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = support_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.title.y = element_blank()) +
  ylab("Share of regional interaction links") +
  ggtitle("A. Regional interaction links differ greatly in how many sites support them")

bin_n <- function(x){
  cut(x, breaks = c(0,1,2,4,8,16,32,Inf), labels = c("1", "2", "3–4", "5–8", "9–16", "17–32", "33+"), right = TRUE)
}
bin_K <- function(x){
  cut(x, breaks = c(-1,0,1,2,4,8,16,32,Inf), labels = c("0", "1", "2", "3–4", "5–8", "9–16", "17–32", "33+"), right = TRUE)
}

heat <- all_pairs %>%
  group_by(dataset) %>% mutate(dataset_weight = 1 / n()) %>% ungroup() %>%
  mutate(n_bin = bin_n(full_n), K_bin = bin_K(full_K)) %>%
  group_by(n_bin, K_bin) %>% summarise(weighted_count = sum(dataset_weight), .groups = "drop") %>%
  mutate(weighted_count = weighted_count / sum(weighted_count))

panel_b <- ggplot(heat, aes(x = n_bin, y = K_bin, fill = weighted_count)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", labels = percent_format(accuracy = 0.1)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom") +
  xlab("Sites where pair is recorded together") +
  ylab("Sites where interaction is observed") +
  ggtitle("B. Where species meet does not fully determine where interactions are sustained", subtitle = "All co-occurring pairs are included, including pairs never observed interacting.")

fig <- panel_a / panel_b + plot_layout(heights = c(1, 1.1)) + plot_annotation(title = "Empirical interaction links differ sharply in local support")
save_both(fig, "31_Figure2_empirical_support_heterogeneity", width = 12, height = 10)

checks <- data.frame(
  check = c("Figure 2B includes K=0 pairs", "all full_n >= 1", "all full_K <= full_n"),
  pass = c(any(all_pairs$full_K == 0), all(all_pairs$full_n >= 1), all(all_pairs$full_K <= all_pairs$full_n))
)
write.csv2(checks, file.path(out_dir, "31_Figure2_checks.csv"), row.names = FALSE)
message("Saved Figure 2 outputs in: ", out_dir)
