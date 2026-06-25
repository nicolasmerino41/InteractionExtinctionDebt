## ------------------------------------------------------------
## Script: All/scripts/20_degree_repeatability_and_partner_loss.R
##
## Purpose:
## Simple descriptive consumer-level analysis:
##   1. Do consumers with more realised resource partners have less
##      locally repeated realised links?
##   2. Do consumers with more realised resource partners lose a larger
##      fraction of their original partners under random site removal?
##   3. Is mean link repeatability associated with partner loss?
##
## Uses only the original 10 datasets.
## No probabilistic models, beta-binomial models, null models, or mediation models.
##
## Run from the parent repository folder.
## Outputs:
##   All/SeparatedResults/script20_degree_repeatability_partner_loss/<dataset>/
##   All/CombinedOutputs/
## ------------------------------------------------------------

packages <- c("dplyr", "ggplot2", "tidyr", "tibble", "bipartite")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

optional_packages <- c("rJava", "XLConnect", "magrittr")
for(pkg in optional_packages){
  suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))
}

set.seed(123)

## ---------------------------
## Settings
## ---------------------------

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 500

dataset_names <- c("Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
                   "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny")

base_out <- "All"
sep_out <- file.path(base_out, "SeparatedResults", "script20_degree_repeatability_partner_loss")
combined_out <- file.path(base_out, "CombinedOutputs")

dir.create(file.path(base_out, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(sep_out, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_out, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Loader helpers copied in the same operational style as the original
## site-removal scripts.
## ---------------------------

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
  if(nrow(web) > 0){
    r_remove <- which(rowSums(web) == 0)
    if(length(r_remove) > 0) web <- web[-r_remove, , drop = FALSE]
  }
  if(ncol(web) > 0){
    c_remove <- which(colSums(web) == 0)
    if(length(c_remove) > 0) web <- web[, -c_remove, drop = FALSE]
  }
  web
}

make_site_subsets <- function(all_sites, removal_levels, n_replicates){
  all_sites <- sort(unique(as.character(all_sites)))
  site_subsets <- list()
  subset_index <- NULL
  counter <- 1
  for(removal in removal_levels){
    n_total <- length(all_sites)
    n_keep <- max(1, round(n_total * (1 - removal)))
    reps_here <- ifelse(removal == 0, 1, n_replicates)
    for(r in seq_len(reps_here)){
      sites_keep <- if(removal == 0) all_sites else sample(all_sites, size = n_keep, replace = FALSE)
      site_subsets[[counter]] <- sites_keep
      subset_index <- rbind(subset_index, data.frame(
        subset_id = counter,
        removal_fraction = removal,
        replicate = r,
        n_sites_kept = length(sites_keep),
        stringsAsFactors = FALSE
      ))
      counter <- counter + 1
    }
  }
  list(index = subset_index, subsets = site_subsets)
}

build_site_tables_from_webs <- function(webs){
  cooc_list <- list()
  interaction_list <- list()
  counter_cooc <- 1
  counter_int <- 1
  for(site_id in names(webs)){
    web <- clean_web_matrix(webs[[site_id]])
    if(nrow(web) == 0 || ncol(web) == 0) next
    consumers <- colnames(web)
    resources <- rownames(web)
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    if(length(consumers) == 0 || length(resources) == 0) next

    cur_cooc <- expand.grid(consumer = consumers, resource = resources, stringsAsFactors = FALSE)
    cur_cooc$site <- as.character(site_id)
    cooc_list[[counter_cooc]] <- cur_cooc[, c("site", "consumer", "resource")]
    counter_cooc <- counter_cooc + 1

    positive_cells <- which(web > 0, arr.ind = TRUE)
    if(nrow(positive_cells) > 0){
      cur_int <- data.frame(
        site = as.character(site_id),
        consumer = colnames(web)[positive_cells[, 2]],
        resource = rownames(web)[positive_cells[, 1]],
        stringsAsFactors = FALSE
      ) %>%
        filter(!is.na(consumer), !is.na(resource), consumer != "", resource != "") %>%
        distinct(site, consumer, resource)
      interaction_list[[counter_int]] <- cur_int
      counter_int <- counter_int + 1
    }
  }
  cooc <- bind_rows(cooc_list) %>% distinct(site, consumer, resource)
  ints <- bind_rows(interaction_list) %>% distinct(site, consumer, resource)

  ## Defensive consistency: an observed interaction is also an operational
  ## pair-by-site co-occurrence record.
  cooc <- bind_rows(cooc, ints) %>% distinct(site, consumer, resource)

  list(cooc_triples = cooc, empirical_site_interactions = ints)
}

load_quercus_webs <- function(){
  if(!requireNamespace("XLConnect", quietly = TRUE)){
    stop("XLConnect is required for Quercus. Install binary if possible: install.packages('XLConnect', type = 'binary')")
  }
  metadata <- read.csv("Quercus/metadata.csv", stringsAsFactors = FALSE)
  metadata <- metadata[1:22, 1:4]
  webs <- list()
  for(t in metadata$Tree){
    excel <- XLConnect::loadWorkbook(paste0("Quercus/raw-data/Web", t, " 2006 Kaartinen.xlsx"))
    int_mat <- XLConnect::readWorksheet(object = excel, sheet = "Sheet1")
    row.names(int_mat) <- int_mat$Col1
    int_mat <- int_mat[-1]
    if(length(which(colnames(int_mat) %in% c("Unpar", "unpar"))) > 0){
      int_mat <- int_mat[-which(colnames(int_mat) %in% c("Unpar", "unpar"))]
    }
    webs[[as.character(t)]] <- clean_web_matrix(int_mat)
  }
  webs
}

load_nahuel_webs <- function(){
  files <- c("vaz_ag_matr_f.txt", "vaz_cl_matr_f.txt", "vaz_ll_matr_f.txt", "vaz_mh_matr_f.txt",
             "vaz_mnh_matr_f.txt", "vaz_qh_matr_f.txt", "vaz_qnh_matr_f.txt", "vaz_s_matr_f.txt")
  vaz <- lapply(files, function(f) t(read.csv(file.path("Nahuel/raw-data", f), sep = "\t", header = FALSE)))
  n_row <- dim(vaz[[1]])[1]
  n_col <- dim(vaz[[1]])[2]
  webs <- list()
  for(i in seq_along(vaz)){
    n <- vaz[[i]]
    colnames(n) <- paste0("Pol", seq_len(n_col))
    rownames(n) <- paste0("Plant", seq_len(n_row))
    webs[[as.character(i)]] <- clean_web_matrix(n)
  }
  webs
}

load_gottin_HP_webs <- function(){
  raw_data <- read.csv("Gottin/raw-data/host_para_all_interactions.csv", head = TRUE, sep = ";")
  webs <- bipartite::frame2webs(raw_data, varnames = c("Genus.Species", "P1.Species.Genus", "Site", "P1.cells"))
  names(webs) <- as.character(names(webs))
  lapply(webs, clean_web_matrix)
}

load_gottin_PP_webs <- function(){
  raw_data <- read.csv("Gottin/raw-data/plant_poll_all_interactions.csv", head = TRUE, sep = ";")
  webs <- bipartite::frame2webs(raw_data, varnames = c("P.Genus.Species", "Genus.Species", "Site"))
  names(webs) <- as.character(names(webs))
  lapply(webs, clean_web_matrix)
}

load_salix_galpar_site_interact <- function(){
  if(!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
  library(magrittr)

  dir.create("Salix/raw-data/csv", recursive = TRUE, showWarnings = FALSE)
  dir.create("Salix/raw-data/rdata", recursive = TRUE, showWarnings = FALSE)

  source("Salix/lib/format4R.r")
  get_formatData("Salix/raw-data/Salix_webs.csv")

  df_site <- readRDS("Salix/raw-data/rdata/df_site.rds")
  df_interact <- readRDS("Salix/raw-data/rdata/df_interact.rds")
  df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR / df_interact$N_GALLS

  merge(df_site, df_interact, by = "REARING_NUMBER")
}

build_salix_galpar_site_tables_original <- function(){
  site_interact <- load_salix_galpar_site_interact()
  site_interact <- site_interact %>%
    filter(!is.na(SITE), !is.na(RGALLER), RGALLER != "")

  ints <- site_interact %>%
    filter(!is.na(RPAR), RPAR != "none", RPAR != "") %>%
    transmute(site = as.character(SITE), consumer = as.character(RPAR), resource = as.character(RGALLER)) %>%
    distinct(site, consumer, resource)

  realised_pairs <- ints %>% distinct(consumer, resource)
  consumers <- sort(unique(realised_pairs$consumer))

  cooc_list <- list()
  for(pred in consumers){
    pred_sites <- unique(site_interact$SITE[site_interact$RPAR == pred &
                                            !is.na(site_interact$RPAR) &
                                            site_interact$RPAR != "none" &
                                            site_interact$RPAR != ""])
    if(length(pred_sites) == 0) next

    cur <- site_interact %>%
      filter(SITE %in% pred_sites, !is.na(RGALLER), RGALLER != "") %>%
      transmute(site = as.character(SITE), consumer = as.character(pred), resource = as.character(RGALLER)) %>%
      distinct(site, consumer, resource)

    cooc_list[[pred]] <- cur
  }

  cooc <- bind_rows(cooc_list) %>% distinct(site, consumer, resource)
  cooc <- bind_rows(cooc, ints) %>% distinct(site, consumer, resource)

  list(cooc_triples = cooc, empirical_site_interactions = ints)
}

load_garraf_montseny_olot_webs <- function(cur_dataset){
  if(!requireNamespace("XLConnect", quietly = TRUE)){
    stop("XLConnect is required for Garraf-Montseny-Olot. Install binary if possible: install.packages('XLConnect', type = 'binary')")
  }
  file_name <- switch(cur_dataset,
                      Garraf_HP = "garraf-hp",
                      Garraf_PP = "garraf-pp",
                      Garraf_PP2 = "garraf-pp-2",
                      Olot = "olot",
                      Montseny = "montseny")
  excel <- XLConnect::loadWorkbook(paste0("Garraf-Montseny-Olot/raw-data/", file_name, ".xlsx"))
  sheet_names <- XLConnect::getSheets(excel)
  names(sheet_names) <- sheet_names
  sheet_list <- lapply(sheet_names, function(.sheet) XLConnect::readWorksheet(object = excel, sheet = .sheet))
  site_sheets <- switch(cur_dataset,
                        Garraf_HP = unique(sheet_names)[-c(1, 27)],
                        Garraf_PP = unique(sheet_names)[-c(1, 42)],
                        Garraf_PP2 = unique(sheet_names),
                        Olot = unique(sheet_names)[-c(1, 16)],
                        Montseny = unique(sheet_names)[-c(1, 20)])
  webs <- list()
  for(p in site_sheets){
    n <- sheet_list[[p]]
    row.names(n) <- n[, 1]
    n <- n[-1]
    if(cur_dataset == "Garraf_PP2") n <- t(n)
    if(cur_dataset == "Montseny"){
      rownames(n) <- paste0("Plant-", rownames(n))
      colnames(n) <- paste0("Pol-", colnames(n))
    }
    webs[[as.character(p)]] <- clean_web_matrix(n)
  }
  webs
}

get_dataset_site_tables <- function(dataset){
  if(dataset == "Salix_Galpar") return(build_salix_galpar_site_tables_original())
  webs <- switch(dataset,
                 Quercus = load_quercus_webs(),
                 Nahuel = load_nahuel_webs(),
                 Gottin_HP = load_gottin_HP_webs(),
                 Gottin_PP = load_gottin_PP_webs(),
                 Garraf_HP = load_garraf_montseny_olot_webs("Garraf_HP"),
                 Garraf_PP = load_garraf_montseny_olot_webs("Garraf_PP"),
                 Garraf_PP2 = load_garraf_montseny_olot_webs("Garraf_PP2"),
                 Olot = load_garraf_montseny_olot_webs("Olot"),
                 Montseny = load_garraf_montseny_olot_webs("Montseny"),
                 stop("Unknown dataset: ", dataset))
  build_site_tables_from_webs(webs)
}

## ---------------------------
## Analysis helpers
## ---------------------------

make_pair_table <- function(dataset, cooc_triples, site_interactions){
  cooc_counts <- cooc_triples %>%
    distinct(site, consumer, resource) %>%
    group_by(consumer, resource) %>%
    summarise(n = n_distinct(site), .groups = "drop")

  int_counts <- site_interactions %>%
    distinct(site, consumer, resource) %>%
    group_by(consumer, resource) %>%
    summarise(K = n_distinct(site), .groups = "drop")

  cooc_counts %>%
    left_join(int_counts, by = c("consumer", "resource")) %>%
    mutate(
      dataset = dataset,
      K = replace_na(K, 0L),
      realised_link = K > 0,
      repeatability = ifelse(n > 0, K / n, NA_real_)
    ) %>%
    select(dataset, consumer, resource, n, K, realised_link, repeatability)
}

make_consumer_full_metrics <- function(pair_table){
  pair_table %>%
    filter(realised_link) %>%
    group_by(dataset, consumer) %>%
    summarise(
      full_degree = n_distinct(resource),
      mean_link_repeatability = mean(repeatability, na.rm = TRUE),
      median_link_repeatability = median(repeatability, na.rm = TRUE),
      number_of_links = n(),
      .groups = "drop"
    )
}

spearman_safe <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3) return(NA_real_)
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

make_summary_rows <- function(consumer_table){
  ds <- unique(consumer_table$dataset)

  data.frame(
    dataset = ds,
    test_name = c("Degree vs repeatability",
                  "Degree vs partner loss",
                  "Repeatability vs partner loss"),
    spearman_correlation_using_mean_repeatability = c(
      spearman_safe(consumer_table$full_degree, consumer_table$mean_link_repeatability),
      spearman_safe(consumer_table$full_degree, consumer_table$mean_degree_loss),
      spearman_safe(consumer_table$mean_link_repeatability, consumer_table$mean_degree_loss)
    ),
    spearman_correlation_using_median_repeatability = c(
      spearman_safe(consumer_table$full_degree, consumer_table$median_link_repeatability),
      NA_real_,
      spearman_safe(consumer_table$median_link_repeatability, consumer_table$mean_degree_loss)
    ),
    number_of_consumers = nrow(consumer_table),
    stringsAsFactors = FALSE
  )
}

## ---------------------------
## Dataset runner
## ---------------------------

run_one_dataset <- function(dataset){
  message("Running script 20 descriptive consumer analysis: ", dataset)

  dataset_dir <- file.path(sep_out, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  cooc_triples <- site_tables$cooc_triples %>% distinct(site, consumer, resource)
  site_interactions <- site_tables$empirical_site_interactions %>% distinct(site, consumer, resource)

  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_replicates)

  pair_table <- make_pair_table(dataset, cooc_triples, site_interactions)
  full_metrics <- make_consumer_full_metrics(pair_table)

  full_partner_table <- pair_table %>%
    filter(realised_link) %>%
    distinct(dataset, consumer, resource)

  full_degree_table <- full_partner_table %>%
    group_by(dataset, consumer) %>%
    summarise(full_degree = n_distinct(resource), .groups = "drop")

  retention_rows <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]

    retained_links <- site_interactions %>%
      filter(site %in% sites_keep) %>%
      distinct(consumer, resource) %>%
      mutate(retained = TRUE)

    cur <- full_partner_table %>%
      left_join(retained_links, by = c("consumer", "resource")) %>%
      mutate(retained = replace_na(retained, FALSE)) %>%
      group_by(dataset, consumer) %>%
      summarise(retained_degree = sum(retained), .groups = "drop") %>%
      left_join(full_degree_table, by = c("dataset", "consumer")) %>%
      mutate(
        removal_fraction = subset_object$index$removal_fraction[i],
        replicate = subset_object$index$replicate[i],
        subset_id = subset_id,
        degree_retention = retained_degree / full_degree,
        degree_loss = 1 - degree_retention
      ) %>%
      select(dataset, consumer, removal_fraction, replicate, subset_id,
             full_degree, retained_degree, degree_retention, degree_loss)

    retention_rows[[i]] <- cur
  }

  retention_by_subset <- bind_rows(retention_rows)

  mean_retention <- retention_by_subset %>%
    filter(removal_fraction > 0) %>%
    group_by(dataset, consumer) %>%
    summarise(
      mean_degree_retention = mean(degree_retention, na.rm = TRUE),
      mean_degree_loss = mean(degree_loss, na.rm = TRUE),
      .groups = "drop"
    )

  consumer_table <- full_metrics %>%
    left_join(mean_retention, by = c("dataset", "consumer")) %>%
    select(dataset, consumer, full_degree, mean_link_repeatability,
           median_link_repeatability, mean_degree_retention, mean_degree_loss)

  summary_table <- make_summary_rows(consumer_table)

  checks <- data.frame(
    dataset = dataset,
    all_consumers_full_degree_at_least_1 = all(consumer_table$full_degree >= 1),
    repeatability_values_between_0_and_1 = all(consumer_table$mean_link_repeatability >= 0 & consumer_table$mean_link_repeatability <= 1 &
                                                consumer_table$median_link_repeatability >= 0 & consumer_table$median_link_repeatability <= 1),
    degree_retention_and_loss_between_0_and_1 = all(consumer_table$mean_degree_retention >= 0 & consumer_table$mean_degree_retention <= 1 &
                                                     consumer_table$mean_degree_loss >= 0 & consumer_table$mean_degree_loss <= 1),
    retained_degree_never_exceeds_full_degree = all(retention_by_subset$retained_degree <= retention_by_subset$full_degree),
    degree_retention_equals_1_at_zero_removal = all(abs(retention_by_subset$degree_retention[retention_by_subset$removal_fraction == 0] - 1) < 1e-12),
    no_models_or_model_predictions_used = TRUE,
    n_consumers = nrow(consumer_table),
    n_site_removal_rows = nrow(retention_by_subset),
    stringsAsFactors = FALSE
  )

  write.csv2(consumer_table,
             file.path(dataset_dir, paste0(dataset, "_20_degree_repeatability_partner_loss_consumers.csv")),
             row.names = FALSE)
  write.csv2(summary_table,
             file.path(dataset_dir, paste0(dataset, "_20_degree_repeatability_partner_loss_summary.csv")),
             row.names = FALSE)
  write.csv2(checks,
             file.path(dataset_dir, paste0(dataset, "_20_degree_repeatability_partner_loss_checks.csv")),
             row.names = FALSE)

  message("Finished script 20 dataset: ", dataset)

  list(
    consumers = consumer_table,
    summary = summary_table,
    checks = checks,
    retention_by_subset = retention_by_subset
  )
}

## ---------------------------
## Run all datasets
## ---------------------------

all_outputs <- lapply(dataset_names, run_one_dataset)
names(all_outputs) <- dataset_names

consumer_all <- bind_rows(lapply(all_outputs, `[[`, "consumers"))
summary_all <- bind_rows(lapply(all_outputs, `[[`, "summary"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

consumer_all$dataset <- factor(consumer_all$dataset, levels = dataset_names)
summary_all$dataset <- factor(summary_all$dataset, levels = dataset_names)
checks_all$dataset <- factor(checks_all$dataset, levels = dataset_names)

write.csv2(
  consumer_all,
  file.path(combined_out, "20_degree_repeatability_partner_loss_consumers.csv"),
  row.names = FALSE
)

write.csv2(
  summary_all,
  file.path(combined_out, "20_degree_repeatability_partner_loss_summary.csv"),
  row.names = FALSE
)

write.csv2(
  checks_all,
  file.path(combined_out, "20_degree_repeatability_partner_loss_checks.csv"),
  row.names = FALSE
)

message("\nValidation checks for script 20:")
print(checks_all)

if(!all(checks_all$all_consumers_full_degree_at_least_1,
        checks_all$repeatability_values_between_0_and_1,
        checks_all$degree_retention_and_loss_between_0_and_1,
        checks_all$retained_degree_never_exceeds_full_degree,
        checks_all$degree_retention_equals_1_at_zero_removal,
        checks_all$no_models_or_model_predictions_used)){
  stop("Validation failed in script 20. Inspect 20_degree_repeatability_partner_loss_checks.csv.")
}

message("\nSpearman summary for script 20:")
print(summary_all)

## ---------------------------
## Figures
## ---------------------------

point_colour <- "grey20"
trend_colour <- "#0072B2"
test_colours <- c(
  "Degree vs repeatability" = "#0072B2",
  "Degree vs partner loss" = "#D55E00",
  "Repeatability vs partner loss" = "#009E73"
)

p_degree_repeat <- ggplot(
  consumer_all,
  aes(x = full_degree, y = mean_link_repeatability)
) +
  geom_point(aes(colour = "Consumers"), alpha = 0.8, size = 1.8) +
  geom_smooth(aes(colour = "Linear visual guide"), method = "lm", se = TRUE, linewidth = 0.9) +
  scale_colour_manual(values = c("Consumers" = point_colour,
                                 "Linear visual guide" = trend_colour)) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 10) +
  labs(x = "Number of realised resource partners",
       y = "Mean repeatability of realised links",
       colour = "")

ggsave(
  file.path(combined_out, "20_degree_repeatability.png"),
  p_degree_repeat,
  width = 14,
  height = 7,
  dpi = 300
)

p_degree_loss <- ggplot(
  consumer_all,
  aes(x = full_degree, y = mean_degree_loss)
) +
  geom_point(aes(colour = "Consumers"), alpha = 0.8, size = 1.8) +
  geom_smooth(aes(colour = "Linear visual guide"), method = "lm", se = TRUE, linewidth = 0.9) +
  scale_colour_manual(values = c("Consumers" = point_colour,
                                 "Linear visual guide" = trend_colour)) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 10) +
  labs(x = "Number of realised resource partners",
       y = "Mean fraction of original partners lost",
       colour = "")

ggsave(
  file.path(combined_out, "20_degree_partner_loss.png"),
  p_degree_loss,
  width = 14,
  height = 7,
  dpi = 300
)

p_repeat_loss <- ggplot(
  consumer_all,
  aes(x = mean_link_repeatability, y = mean_degree_loss)
) +
  geom_point(aes(colour = "Consumers"), alpha = 0.8, size = 1.8) +
  geom_smooth(aes(colour = "Linear visual guide"), method = "lm", se = TRUE, linewidth = 0.9) +
  scale_colour_manual(values = c("Consumers" = point_colour,
                                 "Linear visual guide" = trend_colour)) +
  facet_wrap(~ dataset, ncol = 5) +
  theme_classic(base_size = 10) +
  labs(x = "Mean repeatability of realised links",
       y = "Mean fraction of original partners lost",
       colour = "")

ggsave(
  file.path(combined_out, "20_repeatability_partner_loss.png"),
  p_repeat_loss,
  width = 14,
  height = 7,
  dpi = 300
)

summary_plot <- summary_all %>%
  mutate(
    dataset = factor(dataset, levels = dataset_names),
    test_name = factor(test_name,
                       levels = c("Degree vs repeatability",
                                  "Degree vs partner loss",
                                  "Repeatability vs partner loss"))
  )

p_summary <- ggplot(
  summary_plot,
  aes(x = dataset,
      y = spearman_correlation_using_mean_repeatability,
      colour = test_name)
) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(size = 2.4, position = position_dodge(width = 0.55), na.rm = TRUE) +
  scale_colour_manual(values = test_colours) +
  facet_wrap(~ test_name, ncol = 1) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Dataset",
       y = "Spearman correlation",
       colour = "Test")

ggsave(
  file.path(combined_out, "20_degree_repeatability_partner_loss_summary.png"),
  p_summary,
  width = 12,
  height = 9,
  dpi = 300
)

message("All done. Script 20 outputs saved in All/CombinedOutputs and ", sep_out, ".")
