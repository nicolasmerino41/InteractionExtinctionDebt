
## ------------------------------------------------------------
## Script: nahuel_degree_shape_followup.R
## ------------------------------------------------------------

dataset_label <- "Nahuel"
p_fixed <- 0.065
previous_degree_dir <- "Nahuel/degree_site_removal"
out_dir <- "Nahuel/degree_shape_followup"

load_dataset_webs <- function(){
  files <- c("vaz_ag_matr_f.txt", "vaz_cl_matr_f.txt", "vaz_ll_matr_f.txt", "vaz_mh_matr_f.txt", "vaz_mnh_matr_f.txt", "vaz_qh_matr_f.txt", "vaz_qnh_matr_f.txt", "vaz_s_matr_f.txt")
  vaz <- lapply(files, function(f) t(read.csv(file.path("Nahuel/raw-data", f), sep = "\t", header = FALSE)))
  n_row <- dim(vaz[[1]])[1]; n_col <- dim(vaz[[1]])[2]
  webs <- list()
  for(i in seq_along(vaz)){
    n <- vaz[[i]]
    colnames(n) <- paste0("Pol", seq_len(n_col)); rownames(n) <- paste0("Plant", seq_len(n_row))
    n[n != 0] <- 1
    r_remove <- which(rowSums(n) == 0); if(length(r_remove) != 0) n <- n[-r_remove, , drop = FALSE]
    c_remove <- which(colSums(n) == 0); if(length(c_remove) != 0) n <- n[, -c_remove, drop = FALSE]
    webs[[as.character(i)]] <- n
  }
  webs
}

## ---------------------------
## Shared settings and helpers
## ---------------------------

packages <- c("dplyr", "tidyr", "ggplot2", "tibble")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
selected_removal_levels <- c(0, 0.4, 0.8)
n_site_reps <- 100
n_forced_model_reps <- 100

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_web_matrix <- function(web){
  web <- as.data.frame(web)
  bad_cols <- is.na(colnames(web)) | colnames(web) == ""
  if(any(bad_cols)) web <- web[, !bad_cols, drop = FALSE]
  bad_rows <- is.na(rownames(web)) | rownames(web) == ""
  if(any(bad_rows)) web <- web[!bad_rows, , drop = FALSE]
  web <- as.matrix(web)
  suppressWarnings(storage.mode(web) <- "numeric")
  web[is.na(web)] <- 0
  web[web != 0] <- 1
  web
}

build_site_tables_from_webs <- function(webs){
  cooc_list <- list(); interaction_list <- list(); counter_cooc <- 1; counter_int <- 1
  for(site_id in names(webs)){
    web <- clean_web_matrix(webs[[site_id]])
    if(nrow(web) == 0 || ncol(web) == 0) next
    consumers <- colnames(web); resources <- rownames(web)
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    if(length(consumers) == 0 || length(resources) == 0) next
    cur_cooc <- expand.grid(consumer = consumers, resource = resources, stringsAsFactors = FALSE)
    cur_cooc$site <- site_id
    cooc_list[[counter_cooc]] <- cur_cooc[, c("site", "consumer", "resource")]
    counter_cooc <- counter_cooc + 1
    positive_cells <- which(web > 0, arr.ind = TRUE)
    if(nrow(positive_cells) > 0){
      cur_int <- data.frame(site = site_id, consumer = colnames(web)[positive_cells[, 2]], resource = rownames(web)[positive_cells[, 1]], stringsAsFactors = FALSE) %>%
        filter(!is.na(consumer), !is.na(resource), consumer != "", resource != "") %>%
        distinct(site, consumer, resource)
      interaction_list[[counter_int]] <- cur_int
      counter_int <- counter_int + 1
    }
  }
  list(cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource), empirical_site_interactions = bind_rows(interaction_list) %>% distinct(site, consumer, resource))
}

make_site_subsets <- function(all_sites, removal_levels, n_site_reps){
  site_subsets <- list(); subset_index <- NULL; counter <- 1
  for(removal in removal_levels){
    n_total <- length(all_sites); n_keep <- max(1, round(n_total * (1 - removal))); reps_here <- ifelse(removal == 0, 1, n_site_reps)
    for(r in seq_len(reps_here)){
      sites_keep <- if(removal == 0) all_sites else sample(all_sites, n_keep, replace = FALSE)
      site_subsets[[counter]] <- sites_keep
      subset_index <- rbind(subset_index, data.frame(subset_id = counter, removal_fraction = removal, site_rep = r, n_sites_kept = length(sites_keep)))
      counter <- counter + 1
    }
  }
  list(index = subset_index, subsets = site_subsets)
}

gini_coefficient <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  if(all(x == 0)) return(0)
  x <- sort(x); n <- length(x)
  sum((2 * seq_along(x) - n - 1) * x) / (n * sum(x))
}

degree_table_from_pairs <- function(pairs){
  if(nrow(pairs) == 0) return(data.frame(trophic_level = character(), species = character(), degree = integer()))
  consumer_degrees <- pairs %>% group_by(species = consumer) %>% summarise(degree = n_distinct(resource), .groups = "drop") %>% mutate(trophic_level = "consumer")
  resource_degrees <- pairs %>% group_by(species = resource) %>% summarise(degree = n_distinct(consumer), .groups = "drop") %>% mutate(trophic_level = "resource")
  bind_rows(consumer_degrees, resource_degrees) %>% filter(!is.na(species), species != "", degree > 0) %>% select(trophic_level, species, degree)
}

summarise_degrees <- function(degrees){
  if(nrow(degrees) == 0){
    return(data.frame(trophic_level = c("consumer", "resource"), mean_degree = NA_real_, median_degree = NA_real_, variance_degree = NA_real_, maximum_degree = NA_real_, proportion_degree_1 = NA_real_, n_species_degree_gt0 = 0L, gini_degree = NA_real_))
  }
  degrees %>% group_by(trophic_level) %>% summarise(mean_degree = mean(degree, na.rm = TRUE), median_degree = median(degree, na.rm = TRUE), variance_degree = ifelse(n() > 1, var(degree, na.rm = TRUE), 0), maximum_degree = max(degree, na.rm = TRUE), proportion_degree_1 = mean(degree == 1, na.rm = TRUE), n_species_degree_gt0 = n_distinct(species), gini_degree = gini_coefficient(degree), .groups = "drop") %>% right_join(data.frame(trophic_level = c("consumer", "resource")), by = "trophic_level")
}

frequency_degrees <- function(degrees){
  if(nrow(degrees) == 0) return(data.frame(trophic_level = character(), degree = integer(), n_species = integer()))
  degrees %>% count(trophic_level, degree, name = "n_species") %>% arrange(trophic_level, degree)
}

get_model_prob_pairs_from_subset <- function(cooc_triples, sites_keep, p_fixed){
  cooc_triples %>% filter(site %in% sites_keep) %>% group_by(consumer, resource) %>% summarise(n_cooc = n_distinct(site), .groups = "drop") %>% mutate(prob_int = 1 - (1 - p_fixed)^n_cooc) %>% filter(prob_int > 0)
}

sample_forced_model_pairs <- function(prob_pairs, L_emp){
  if(nrow(prob_pairs) == 0 || L_emp <= 0) return(data.frame(consumer = character(), resource = character()))
  L_sample <- min(L_emp, nrow(prob_pairs))
  sampled_rows <- sample(seq_len(nrow(prob_pairs)), size = L_sample, replace = FALSE, prob = prob_pairs$prob_int)
  prob_pairs[sampled_rows, c("consumer", "resource"), drop = FALSE] %>% distinct(consumer, resource)
}

estimate_powerlaw_from_degrees <- function(degree_vector, xmin = 1L, min_n_reliable = 20L){
  x <- degree_vector[!is.na(degree_vector)]; x <- x[x >= xmin]; n <- length(x)
  attempted <- n >= 2 && length(unique(x)) >= 2
  reliable <- attempted && n >= min_n_reliable && max(x) >= 3
  if(!attempted){ alpha <- NA_real_ } else { denom <- sum(log(x / (xmin - 0.5))); alpha <- ifelse(denom > 0, 1 + n / denom, NA_real_) }
  data.frame(n_species_used = n, xmin_used = xmin, fit_attempted = attempted, fit_reliable = reliable, exponent_alpha = alpha, estimator = "continuous_MLE_approximation_xmin_1")
}

expand_degree_frequency <- function(freq_df){
  if(nrow(freq_df) == 0) return(integer())
  rep(freq_df$degree, times = freq_df$n_species)
}

estimate_powerlaw_from_frequency_file <- function(freq_df, source_label){
  out <- freq_df %>% group_by(source, model_rep, subset_id, removal_fraction, site_rep, n_sites_kept, trophic_level) %>% summarise(degree_vector = list(expand_degree_frequency(cur_data())), .groups = "drop")
  est_list <- vector("list", nrow(out))
  for(i in seq_len(nrow(out))){
    est <- estimate_powerlaw_from_degrees(out$degree_vector[[i]])
    est_list[[i]] <- bind_cols(out[i, -which(names(out) == "degree_vector")], est)
  }
  bind_rows(est_list) %>% mutate(source = source_label)
}

make_powerlaw_envelope <- function(empirical_pw, model_pw){
  model_envelope <- model_pw %>% group_by(trophic_level, removal_fraction) %>% summarise(model_q025 = quantile(exponent_alpha, 0.025, na.rm = TRUE), model_q500 = quantile(exponent_alpha, 0.500, na.rm = TRUE), model_q975 = quantile(exponent_alpha, 0.975, na.rm = TRUE), model_mean = mean(exponent_alpha, na.rm = TRUE), model_n_reliable = sum(fit_reliable, na.rm = TRUE), .groups = "drop")
  empirical_summary <- empirical_pw %>% group_by(trophic_level, removal_fraction) %>% summarise(empirical_mean = mean(exponent_alpha, na.rm = TRUE), empirical_sd = sd(exponent_alpha, na.rm = TRUE), empirical_n_reliable = sum(fit_reliable, na.rm = TRUE), .groups = "drop")
  empirical_summary %>% left_join(model_envelope, by = c("trophic_level", "removal_fraction")) %>% mutate(empirical_above_model_975 = empirical_mean > model_q975, empirical_below_model_025 = empirical_mean < model_q025)
}

plot_powerlaw <- function(comparison){
  p <- ggplot(comparison, aes(x = removal_fraction)) + geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") + geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) + geom_line(aes(y = empirical_mean), linewidth = 1.2) + geom_point(aes(y = empirical_mean), size = 2.5) + facet_wrap(~ trophic_level, scales = "free_y") + theme_classic(base_size = 14) + xlab("Fraction of sites removed") + ylab("Fitted power-law exponent") + ggtitle(paste0(dataset_label, ": power-law exponent under site removal"), subtitle = "Diagnostic continuous MLE approximation; ribbon = model envelope")
  ggsave(file.path(out_dir, "powerlaw_exponent_under_site_removal.png"), p, width = 8.5, height = 5.5, dpi = 300)
}

make_forced_metric_envelope <- function(model_metrics, empirical_metrics){
  metric_cols <- c("mean_degree", "median_degree", "variance_degree", "maximum_degree", "proportion_degree_1", "n_species_degree_gt0", "gini_degree")
  model_long <- model_metrics %>% pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value")
  model_envelope <- model_long %>% group_by(trophic_level, removal_fraction, metric) %>% summarise(model_q025 = quantile(value, 0.025, na.rm = TRUE), model_q500 = quantile(value, 0.500, na.rm = TRUE), model_q975 = quantile(value, 0.975, na.rm = TRUE), model_mean = mean(value, na.rm = TRUE), .groups = "drop")
  empirical_long <- empirical_metrics %>% pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "empirical_value") %>% group_by(trophic_level, removal_fraction, metric) %>% summarise(empirical_mean = mean(empirical_value, na.rm = TRUE), .groups = "drop")
  empirical_long %>% left_join(model_envelope, by = c("trophic_level", "removal_fraction", "metric")) %>% mutate(empirical_above_model_975 = empirical_mean > model_q975, empirical_below_model_025 = empirical_mean < model_q025)
}

make_forced_frequency_envelope <- function(model_frequency, empirical_frequency){
  model_envelope <- model_frequency %>% filter(removal_fraction %in% selected_removal_levels) %>% group_by(removal_fraction, trophic_level, degree) %>% summarise(model_q025 = quantile(n_species, 0.025, na.rm = TRUE), model_q500 = quantile(n_species, 0.500, na.rm = TRUE), model_q975 = quantile(n_species, 0.975, na.rm = TRUE), .groups = "drop")
  empirical_summary <- empirical_frequency %>% filter(removal_fraction %in% selected_removal_levels) %>% group_by(removal_fraction, trophic_level, degree) %>% summarise(empirical_mean = mean(n_species, na.rm = TRUE), .groups = "drop")
  all_degrees <- sort(unique(c(model_envelope$degree, empirical_summary$degree)))
  full_grid <- expand.grid(removal_fraction = selected_removal_levels, trophic_level = c("consumer", "resource"), degree = all_degrees, stringsAsFactors = FALSE)
  full_grid %>% left_join(model_envelope, by = c("removal_fraction", "trophic_level", "degree")) %>% left_join(empirical_summary, by = c("removal_fraction", "trophic_level", "degree")) %>% mutate(across(c(model_q025, model_q500, model_q975, empirical_mean), ~replace_na(.x, 0)))
}

plot_forced_metric <- function(envelope, metric_name, y_label, out_file){
  dat <- envelope %>% filter(metric == metric_name)
  p <- ggplot(dat, aes(x = removal_fraction)) + geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") + geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) + geom_line(aes(y = empirical_mean), linewidth = 1.2) + geom_point(aes(y = empirical_mean), size = 2.5) + facet_wrap(~ trophic_level, scales = "free_y") + theme_classic(base_size = 14) + xlab("Fraction of sites removed") + ylab(y_label) + ggtitle(paste0(dataset_label, ": forced-link-count model, ", y_label), subtitle = "Model has exactly the same number of links as empirical subset")
  ggsave(out_file, p, width = 8.5, height = 5.5, dpi = 300)
}

plot_forced_frequency <- function(freq_envelope){
  p <- ggplot(freq_envelope, aes(x = degree)) + geom_ribbon(aes(ymin = model_q025, ymax = model_q975), fill = "grey80") + geom_line(aes(y = model_q500), linetype = 2, linewidth = 1) + geom_line(aes(y = empirical_mean), linewidth = 1) + geom_point(aes(y = empirical_mean), size = 2) + facet_grid(trophic_level ~ removal_fraction, scales = "free_y") + theme_classic(base_size = 13) + xlab("Degree") + ylab("Number of species") + ggtitle(paste0(dataset_label, ": forced-link-count degree-frequency distributions"), subtitle = "Columns = removal fraction; ribbon = forced model envelope; solid = empirical")
  ggsave(file.path(out_dir, "forced_link_model_degree_frequency_selected_removals.png"), p, width = 10, height = 6.5, dpi = 300)
}

run_powerlaw_analysis <- function(){
  empirical_freq_file <- file.path(previous_degree_dir, "degree_frequency_empirical.csv")
  model_freq_file <- file.path(previous_degree_dir, "degree_frequency_model_generated.csv")
  if(!file.exists(empirical_freq_file) || !file.exists(model_freq_file)) stop("Previous degree-frequency files not found. Run the degree_site_removal script first.")
  empirical_freq <- read.csv(empirical_freq_file, stringsAsFactors = FALSE)
  model_freq <- read.csv(model_freq_file, stringsAsFactors = FALSE)
  empirical_pw <- estimate_powerlaw_from_frequency_file(empirical_freq, "empirical")
  model_pw <- estimate_powerlaw_from_frequency_file(model_freq, "model_generated")
  comparison_pw <- make_powerlaw_envelope(empirical_pw, model_pw)
  write.csv(empirical_pw, file.path(out_dir, "powerlaw_exponents_empirical.csv"), row.names = FALSE)
  write.csv(model_pw, file.path(out_dir, "powerlaw_exponents_model_generated.csv"), row.names = FALSE)
  write.csv(comparison_pw, file.path(out_dir, "powerlaw_exponents_empirical_vs_model_envelope.csv"), row.names = FALSE)
  plot_powerlaw(comparison_pw)
}

run_forced_link_analysis <- function(){
  message("Loading ", dataset_label, " for forced-link-count analysis")
  webs <- load_dataset_webs()
  site_tables <- build_site_tables_from_webs(webs)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions
  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_site_reps)
  emp_metric_list <- list(); emp_freq_list <- list(); model_metric_list <- list(); model_freq_list <- list()
  counter_emp <- 1; counter_model <- 1
  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    empirical_pairs <- empirical_site_interactions %>% filter(site %in% sites_keep) %>% distinct(consumer, resource)
    L_emp <- nrow(empirical_pairs)
    emp_degrees <- degree_table_from_pairs(empirical_pairs)
    emp_metric_list[[counter_emp]] <- summarise_degrees(emp_degrees) %>% mutate(source = "empirical", subset_id = subset_id, removal_fraction = subset_object$index$removal_fraction[i], site_rep = subset_object$index$site_rep[i], n_sites_kept = subset_object$index$n_sites_kept[i], L_emp = L_emp)
    emp_freq_list[[counter_emp]] <- frequency_degrees(emp_degrees) %>% mutate(source = "empirical", subset_id = subset_id, removal_fraction = subset_object$index$removal_fraction[i], site_rep = subset_object$index$site_rep[i], n_sites_kept = subset_object$index$n_sites_kept[i], L_emp = L_emp)
    counter_emp <- counter_emp + 1
    prob_pairs <- get_model_prob_pairs_from_subset(cooc_triples, sites_keep, p_fixed)
    for(m in seq_len(n_forced_model_reps)){
      forced_pairs <- sample_forced_model_pairs(prob_pairs, L_emp)
      model_degrees <- degree_table_from_pairs(forced_pairs)
      model_metric_list[[counter_model]] <- summarise_degrees(model_degrees) %>% mutate(source = "forced_link_model", forced_model_rep = m, subset_id = subset_id, removal_fraction = subset_object$index$removal_fraction[i], site_rep = subset_object$index$site_rep[i], n_sites_kept = subset_object$index$n_sites_kept[i], L_emp = L_emp, L_forced = nrow(forced_pairs))
      model_freq_list[[counter_model]] <- frequency_degrees(model_degrees) %>% mutate(source = "forced_link_model", forced_model_rep = m, subset_id = subset_id, removal_fraction = subset_object$index$removal_fraction[i], site_rep = subset_object$index$site_rep[i], n_sites_kept = subset_object$index$n_sites_kept[i], L_emp = L_emp, L_forced = nrow(forced_pairs))
      counter_model <- counter_model + 1
    }
    message(dataset_label, ": forced-link subset ", i, " / ", nrow(subset_object$index))
  }
  empirical_metrics <- bind_rows(emp_metric_list); empirical_frequency <- bind_rows(emp_freq_list); model_metrics <- bind_rows(model_metric_list); model_frequency <- bind_rows(model_freq_list)
  write.csv(empirical_metrics, file.path(out_dir, "degree_metrics_forced_link_model_empirical.csv"), row.names = FALSE)
  write.csv(model_metrics, file.path(out_dir, "degree_metrics_forced_link_model_generated.csv"), row.names = FALSE)
  write.csv(empirical_frequency, file.path(out_dir, "degree_frequency_forced_link_model_empirical.csv"), row.names = FALSE)
  write.csv(model_frequency, file.path(out_dir, "degree_frequency_forced_link_model_generated.csv"), row.names = FALSE)
  forced_envelope <- make_forced_metric_envelope(model_metrics, empirical_metrics)
  forced_freq_envelope <- make_forced_frequency_envelope(model_frequency, empirical_frequency)
  write.csv(forced_envelope, file.path(out_dir, "degree_metrics_forced_link_model_envelope.csv"), row.names = FALSE)
  write.csv(forced_freq_envelope, file.path(out_dir, "degree_frequency_forced_link_model_envelope.csv"), row.names = FALSE)
  plot_forced_metric(forced_envelope, "gini_degree", "Degree Gini coefficient", file.path(out_dir, "forced_link_model_degree_gini.png"))
  plot_forced_metric(forced_envelope, "maximum_degree", "Maximum degree", file.path(out_dir, "forced_link_model_maximum_degree.png"))
  plot_forced_metric(forced_envelope, "proportion_degree_1", "Proportion of degree-1 species", file.path(out_dir, "forced_link_model_proportion_degree1.png"))
  plot_forced_frequency(forced_freq_envelope)
}

message("Running power-law follow-up for ", dataset_label)
run_powerlaw_analysis()
message("Running forced-link-count follow-up for ", dataset_label)
run_forced_link_analysis()
message("Finished ", dataset_label, ". Outputs saved to ", out_dir)
