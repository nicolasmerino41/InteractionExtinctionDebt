## ------------------------------------------------------------
## Script: All/scripts/31_figure4_cooccurrence_baseline_and_support_persistence.R
## Homogeneous fixed-p baseline and support-mediated persistence.
## ------------------------------------------------------------
source("All/scripts/00_dataset_loaders_and_helpers_all.R")

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "patchwork", "future", "future.apply", "parallelly")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_site_reps <- 50
out_dir <- "All/outputs/31_hidden_support_main_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, name, width = 11, height = 7){
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
}
make_pair_id <- function(consumer, resource) paste(consumer, resource, sep = "___")
plink <- function(n, p) 1 - (1 - p)^n
fit_p_to_links <- function(n_vec, L_obs){
  if(L_obs <= 0) return(0)
  if(L_obs >= length(n_vec)) return(1)
  f <- function(p) sum(plink(n_vec, p), na.rm = TRUE) - L_obs
  uniroot(f, interval = c(0, 1), tol = 1e-10)$root
}
support_class <- function(n){
  case_when(n == 2 ~ "2 sites", n >= 3 & n <= 4 ~ "3–4 sites", n >= 5 & n <= 8 ~ "5–8 sites", n >= 9 & n <= 16 ~ "9–16 sites", n >= 17 ~ "17 or more sites", TRUE ~ NA_character_)
}
support_levels <- c("2 sites", "3–4 sites", "5–8 sites", "9–16 sites", "17 or more sites")

run_one_dataset <- function(dataset){
  message("Running Figure 4 baseline: ", dataset)
  site_tables <- get_dataset_site_tables(dataset)
  cooc <- site_tables$cooc_triples %>% distinct(site, consumer, resource) %>% mutate(pair_id = make_pair_id(consumer, resource))
  ints <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource) %>% mutate(pair_id = make_pair_id(consumer, resource))
  pair_table <- cooc %>%
    group_by(pair_id, consumer, resource) %>% summarise(full_n = n_distinct(site), .groups = "drop") %>%
    left_join(ints %>% group_by(pair_id) %>% summarise(full_K = n_distinct(site), .groups = "drop"), by = "pair_id") %>%
    mutate(full_K = replace_na(full_K, 0L), full_B = full_K > 0)
  L_obs <- sum(pair_table$full_B)
  p_full <- fit_p_to_links(pair_table$full_n, L_obs)
  expected_L <- sum(plink(pair_table$full_n, p_full))
  realised <- pair_table %>% filter(full_K >= 1)

  empirical_support <- realised %>%
    mutate(class = factor(support_class(full_n), levels = support_levels)) %>%
    filter(!is.na(class)) %>%
    group_by(dataset = dataset, class) %>%
    summarise(empirical_fraction_repeated = mean(full_K > 1), n_emp_links = n(), .groups = "drop")
  model_support <- pair_table %>%
    mutate(class = factor(support_class(full_n), levels = support_levels)) %>%
    filter(!is.na(class)) %>%
    group_by(dataset = dataset, class) %>%
    summarise(n_pairs = n(), model_fraction_repeated = sum(plink(full_n, p_full) - full_n * p_full * (1 - p_full)^(full_n - 1), na.rm = TRUE) / sum(plink(full_n, p_full), na.rm = TRUE), .groups = "drop")

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  removal_rows <- lapply(seq_len(nrow(subset_object$index)), function(i){
    idx <- subset_object$index[i,]
    sites_keep <- subset_object$subsets[[idx$subset_id]]
    ret_cooc <- cooc %>% filter(site %in% sites_keep) %>% group_by(pair_id) %>% summarise(retained_n = n_distinct(site), .groups = "drop")
    ret_int <- ints %>% filter(site %in% sites_keep) %>% group_by(pair_id) %>% summarise(retained_K = n_distinct(site), .groups = "drop")
    empirical_retained <- realised %>% left_join(ret_int, by = "pair_id") %>% mutate(retained_K = replace_na(retained_K, 0L)) %>% summarise(ret = mean(retained_K > 0)) %>% pull(ret)
    expected_retained_links <- ret_cooc %>% summarise(x = sum(plink(retained_n, p_full), na.rm = TRUE)) %>% pull(x)
    data.frame(dataset = dataset, removal_fraction = idx$removal_fraction, replicate = idx$site_rep, empirical_binary_retention = empirical_retained, model_expected_binary_retention = expected_retained_links / expected_L, empirical_minus_model = empirical_retained - expected_retained_links / expected_L)
  })

  list(calibration = data.frame(dataset = dataset, p_full = p_full, calibration_method = "unconditioned full-network link-count match: sum[1-(1-p)^n]=L_obs", empirical_full_regional_link_count = L_obs, expected_model_full_regional_link_count = expected_L, calibration_difference = L_obs - expected_L, notes = "p calibrated once on intact network; not refitted after site removal"), support = full_join(empirical_support, model_support, by = c("dataset", "class")), removal = bind_rows(removal_rows))
}

workers <- max(1, min(length(all_dataset_names), parallelly::availableCores() - 1))
future::plan(future::multisession, workers = workers)
outputs <- future.apply::future_lapply(all_dataset_names, run_one_dataset, future.seed = TRUE)
future::plan(future::sequential)

calibration <- bind_rows(lapply(outputs, `[[`, "calibration"))
support_summary <- bind_rows(lapply(outputs, `[[`, "support"))
retention <- bind_rows(lapply(outputs, `[[`, "removal"))
write.csv2(calibration, file.path(out_dir, "31_p_full_calibration.csv"), row.names = FALSE)
write.csv2(support_summary, file.path(out_dir, "31_homogeneous_support_by_cooccurrence_class.csv"), row.names = FALSE)
write.csv2(retention, file.path(out_dir, "31_empirical_minus_model_link_retention.csv"), row.names = FALSE)

support_long <- support_summary %>%
  select(dataset, class, empirical_fraction_repeated, model_fraction_repeated) %>%
  pivot_longer(cols = c(empirical_fraction_repeated, model_fraction_repeated), names_to = "source", values_to = "fraction") %>%
  mutate(source = recode(source, empirical_fraction_repeated = "Empirical data", model_fraction_repeated = "Homogeneous baseline"))
dataset_support <- support_long %>% group_by(dataset, class, source) %>% summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop")
cross_support <- dataset_support %>% group_by(class, source) %>% summarise(median_fraction = median(fraction, na.rm = TRUE), .groups = "drop")
cols <- c("Empirical data" = "#E69F00", "Homogeneous baseline" = "#6A3D9A")
panel_a <- ggplot() +
  geom_line(data = dataset_support, aes(x = class, y = fraction, group = interaction(dataset, source), colour = source), alpha = 0.20, linewidth = 0.45) +
  geom_line(data = cross_support, aes(x = class, y = median_fraction, group = source, colour = source), linewidth = 1.25) +
  geom_point(data = cross_support, aes(x = class, y = median_fraction, colour = source), size = 2) +
  scale_colour_manual(values = cols) + coord_cartesian(ylim = c(0,1)) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1)) +
  xlab("Sites where pair is recorded together") + ylab("Fraction of regional links observed in more than one site") +
  ggtitle("A. Given similar opportunities to meet, empirical interactions are more often repeatedly observed")

dataset_ret <- retention %>% group_by(dataset, removal_fraction) %>% summarise(empirical_minus_model = median(empirical_minus_model, na.rm = TRUE), .groups = "drop")
cross_ret <- dataset_ret %>% group_by(removal_fraction) %>% summarise(median_difference = median(empirical_minus_model, na.rm = TRUE), q25 = quantile(empirical_minus_model, 0.25, na.rm = TRUE), q75 = quantile(empirical_minus_model, 0.75, na.rm = TRUE), .groups = "drop")
panel_b <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_ribbon(data = cross_ret, aes(x = removal_fraction, ymin = q25, ymax = q75), fill = "#E69F00", alpha = 0.15) +
  geom_line(data = dataset_ret, aes(x = removal_fraction, y = empirical_minus_model, group = dataset), colour = "#E69F00", alpha = 0.25, linewidth = 0.45) +
  geom_line(data = cross_ret, aes(x = removal_fraction, y = median_difference), colour = "#D55E00", linewidth = 1.25) +
  geom_point(data = cross_ret, aes(x = removal_fraction, y = median_difference), colour = "#D55E00", size = 2) +
  theme_classic(base_size = 10) + xlab("Proportion of sites removed") + ylab("Extra regional links retained in empirical networks") +
  ggtitle("B. Excess binary-link retention after site removal", subtitle = "Positive values mean empirical regional links persist more than expected from a common per-site interaction probability.")
fig <- panel_a / panel_b + plot_annotation(title = "Co-occurrence can recover intact structure while missing support-mediated persistence")
save_both(fig, "31_Figure4_cooccurrence_baseline_and_support_persistence", width = 11, height = 8)

checks <- data.frame(check = c("p_calibrated_once", "zero_removal_difference_near_zero", "no_refit_after_removal"), pass = c(TRUE, all(abs(retention$empirical_minus_model[retention$removal_fraction == 0]) < 1e-6), TRUE))
write.csv2(checks, file.path(out_dir, "31_Figure4_checks.csv"), row.names = FALSE)
message("Saved Figure 4 outputs in: ", out_dir)
message("An interaction still present somewhere is not necessarily equally represented across the landscape. These analyses quantify the number of sampled sites still supporting each interaction and its persistence under random site loss. They do not directly measure the ecological effect of the interaction within a local site.")
