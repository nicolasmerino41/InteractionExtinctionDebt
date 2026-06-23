
## ------------------------------------------------------------
## Script: All/scripts/01_site_removal_divergence_all.R
##
## Purpose:
## Unified site-removal divergence pipeline for 10 datasets:
## Quercus, Nahuel, Salix/Galpar, Gottin_HP, Gottin_PP,
## Garraf_HP, Garraf_PP, Garraf_PP2, Olot, Montseny.
##
## Run from the parent repository folder.
## Outputs:
##   All/SeparatedResults/site_removal_divergence/<dataset>/
##   All/CombinedOutputs/
## ------------------------------------------------------------

packages <- c("dplyr", "ggplot2", "tidyr", "tibble", "igraph", "bipartite")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

optional_packages <- c("rJava", "XLConnect", "magrittr", "reshape2", "data.table", "broom", "AICcmodavg")
for(pkg in optional_packages){
  suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 100

base_out <- "All"
sep_out <- file.path(base_out, "SeparatedResults", "site_removal_divergence")
combined_out <- file.path(base_out, "CombinedOutputs")

dir.create(file.path(base_out, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(sep_out, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_out, recursive = TRUE, showWarnings = FALSE)

dataset_config <- data.frame(
  dataset = c("Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
              "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"),
  p_fixed = c(0.0227, 0.065, 0.121, 0.074, 0.044,
              0.076, 0.100, 0.0962, 0.070, 0.099),
  stringsAsFactors = FALSE
)

datasets_to_run <- dataset_config$dataset

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0) NA_real_ else (observed - expected) / expected
}

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
        subset_id = counter, removal_fraction = removal,
        rep = r, n_sites_kept = length(sites_keep)
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
    cur_cooc$site <- site_id
    cooc_list[[counter_cooc]] <- cur_cooc[, c("site", "consumer", "resource")]
    counter_cooc <- counter_cooc + 1
    positive_cells <- which(web > 0, arr.ind = TRUE)
    if(nrow(positive_cells) > 0){
      cur_int <- data.frame(
        site = site_id,
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
  list(
    cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>% distinct(site, consumer, resource)
  )
}

expected_links_from_subset <- function(cooc_triples, sites_keep, p_fixed){
  cur <- cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")
  if(nrow(cur) == 0){
    return(data.frame(n_cooc_links = 0, n_expected_links_raw = NA_real_,
                      n_expected_links_conditioned = NA_real_))
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
    n_expected_links_raw = sum(cur$prob_int, na.rm = TRUE),
    n_expected_links_conditioned = sum(cur$prob_expected, na.rm = TRUE)
  )
}

observed_links_from_subset <- function(site_interactions, sites_keep){
  site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource) %>%
    nrow()
}

species_error_from_subset <- function(cooc_triples, site_interactions, sites_keep, p_fixed){
  cooc <- cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")
  if(nrow(cooc) == 0){
    return(data.frame(
      mean_species_bias_raw = NA_real_, mean_species_bias_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_, mean_abs_species_error_conditioned = NA_real_
    ))
  }
  cooc <- cooc %>%
    group_by(consumer) %>%
    mutate(N_alpha = sum(n_cooc)) %>%
    ungroup() %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )
  expected_degree <- cooc %>%
    group_by(species = consumer) %>%
    summarise(expected_degree_raw = sum(prob_int, na.rm = TRUE),
              expected_degree_conditioned = sum(prob_expected, na.rm = TRUE),
              .groups = "drop")
  observed_degree <- site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource) %>%
    group_by(species = consumer) %>%
    summarise(observed_degree = n_distinct(resource), .groups = "drop")
  degree_compare <- merge(expected_degree, observed_degree, by = "species", all.x = TRUE)
  degree_compare$observed_degree[is.na(degree_compare$observed_degree)] <- 0
  data.frame(
    mean_species_bias_raw = mean(degree_compare$observed_degree - degree_compare$expected_degree_raw, na.rm = TRUE),
    mean_species_bias_conditioned = mean(degree_compare$observed_degree - degree_compare$expected_degree_conditioned, na.rm = TRUE),
    mean_abs_species_error_raw = mean(abs(degree_compare$observed_degree - degree_compare$expected_degree_raw), na.rm = TRUE),
    mean_abs_species_error_conditioned = mean(abs(degree_compare$observed_degree - degree_compare$expected_degree_conditioned), na.rm = TRUE)
  )
}

## Original Salix/Galpar logic, kept separate because Galpar co-occurrence is not
## equivalent to a generic incidence-matrix loader.
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

run_galpar_subset_original <- function(site_interact,
                                       sites_keep,
                                       p_fixed = 0.121,
                                       removal_fraction = NA_real_,
                                       rep_id = NA_integer_,
                                       subset_id = NA_integer_){

  sub <- site_interact[site_interact$SITE %in% sites_keep, ]
  sub <- sub[!is.na(sub$SITE), ]
  sub <- sub[!is.na(sub$RGALLER), ]

  real_pairs <- sub %>%
    filter(!is.na(RPAR), RPAR != "none", RPAR != "") %>%
    distinct(species = RPAR, potential_prey = RGALLER)

  n_realised_links <- nrow(real_pairs)

  if(n_realised_links == 0){
    return(data.frame(
      dataset = "Salix_Galpar",
      removal_fraction = removal_fraction,
      rep = rep_id,
      subset_id = subset_id,
      n_sites_kept = length(unique(sites_keep)),
      n_realised_links = 0,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_species_bias_raw = NA_real_,
      mean_species_bias_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }

  predators <- sort(unique(real_pairs$species))
  cooc_list <- list()

  for(pred in predators){
    pred_sites <- unique(sub$SITE[sub$RPAR == pred])
    if(length(pred_sites) == 0) next

    cur <- sub[sub$SITE %in% pred_sites, c("SITE", "RGALLER")]
    cur <- unique(cur)

    cur_cooc <- cur %>%
      group_by(RGALLER) %>%
      summarise(n_cooc = n_distinct(SITE), .groups = "drop") %>%
      transmute(
        species = pred,
        potential_prey = RGALLER,
        n_cooc = n_cooc
      )

    cooc_list[[pred]] <- cur_cooc
  }

  cooc <- bind_rows(cooc_list) %>%
    filter(!is.na(potential_prey), potential_prey != "")

  n_cooc_links <- nrow(cooc)

  if(n_cooc_links == 0){
    return(data.frame(
      dataset = "Salix_Galpar",
      removal_fraction = removal_fraction,
      rep = rep_id,
      subset_id = subset_id,
      n_sites_kept = length(unique(sites_keep)),
      n_realised_links = n_realised_links,
      n_cooc_links = 0,
      n_expected_links_raw = NA_real_,
      n_expected_links_conditioned = NA_real_,
      divergence_raw = NA_real_,
      relative_divergence_raw = NA_real_,
      divergence_conditioned = NA_real_,
      relative_divergence_conditioned = NA_real_,
      mean_species_bias_raw = NA_real_,
      mean_species_bias_conditioned = NA_real_,
      mean_abs_species_error_raw = NA_real_,
      mean_abs_species_error_conditioned = NA_real_
    ))
  }

  cooc <- cooc %>%
    group_by(species) %>%
    mutate(
      N_alpha = sum(n_cooc),
      total_potential_interactions = n_distinct(potential_prey)
    ) %>%
    ungroup() %>%
    mutate(
      prob_int = 1 - (1 - p_fixed)^n_cooc,
      prob_expected = prob_int / (1 - (1 - p_fixed)^N_alpha)
    )

  n_expected_links_raw <- sum(cooc$prob_int, na.rm = TRUE)
  n_expected_links_conditioned <- sum(cooc$prob_expected, na.rm = TRUE)

  expected_degree <- cooc %>%
    group_by(species) %>%
    summarise(
      expected_degree_raw = sum(prob_int, na.rm = TRUE),
      expected_degree_conditioned = sum(prob_expected, na.rm = TRUE),
      .groups = "drop"
    )

  observed_degree <- real_pairs %>%
    group_by(species) %>%
    summarise(observed_degree = n_distinct(potential_prey), .groups = "drop")

  degree_compare <- merge(expected_degree, observed_degree, by = "species", all.x = TRUE)
  degree_compare$observed_degree[is.na(degree_compare$observed_degree)] <- 0

  divergence_raw <- n_realised_links - n_expected_links_raw
  divergence_conditioned <- n_realised_links - n_expected_links_conditioned

  data.frame(
    dataset = "Salix_Galpar",
    removal_fraction = removal_fraction,
    rep = rep_id,
    subset_id = subset_id,
    n_sites_kept = length(unique(sites_keep)),
    n_realised_links = n_realised_links,
    n_cooc_links = n_cooc_links,
    n_expected_links_raw = n_expected_links_raw,
    n_expected_links_conditioned = n_expected_links_conditioned,
    divergence_raw = divergence_raw,
    relative_divergence_raw = safe_relative_difference(n_realised_links, n_expected_links_raw),
    divergence_conditioned = divergence_conditioned,
    relative_divergence_conditioned = safe_relative_difference(n_realised_links, n_expected_links_conditioned),
    mean_species_bias_raw = mean(degree_compare$observed_degree - degree_compare$expected_degree_raw, na.rm = TRUE),
    mean_species_bias_conditioned = mean(degree_compare$observed_degree - degree_compare$expected_degree_conditioned, na.rm = TRUE),
    mean_abs_species_error_raw = mean(abs(degree_compare$observed_degree - degree_compare$expected_degree_raw), na.rm = TRUE),
    mean_abs_species_error_conditioned = mean(abs(degree_compare$observed_degree - degree_compare$expected_degree_conditioned), na.rm = TRUE)
  )
}

run_salix_galpar_original <- function(){
  dataset <- "Salix_Galpar"
  message("Running ", dataset)

  p_fixed <- dataset_config$p_fixed[dataset_config$dataset == dataset]
  dataset_dir <- file.path(sep_out, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  site_interact <- load_salix_galpar_site_interact()
  all_sites <- sort(unique(site_interact$SITE))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_replicates)

  out_list <- vector("list", nrow(subset_object$index))

  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]

    out_list[[i]] <- run_galpar_subset_original(
      site_interact = site_interact,
      sites_keep = sites_keep,
      p_fixed = p_fixed,
      removal_fraction = subset_object$index$removal_fraction[i],
      rep_id = subset_object$index$rep[i],
      subset_id = subset_id
    )
  }

  results <- bind_rows(out_list)

  summary_table <- results %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      mean_realised_links = mean(n_realised_links, na.rm = TRUE),
      mean_expected_links_raw = mean(n_expected_links_raw, na.rm = TRUE),
      mean_expected_links_conditioned = mean(n_expected_links_conditioned, na.rm = TRUE),
      mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
      mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
      sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
      mean_divergence_conditioned = mean(divergence_conditioned, na.rm = TRUE),
      mean_relative_divergence_conditioned = mean(relative_divergence_conditioned, na.rm = TRUE),
      mean_species_bias_raw = mean(mean_species_bias_raw, na.rm = TRUE),
      mean_abs_species_error_raw = mean(mean_abs_species_error_raw, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv2(results,
            file.path(dataset_dir, paste0(dataset, "_site_removal_divergence.csv")),
            row.names = FALSE)

  write.csv2(summary_table,
            file.path(dataset_dir, paste0(dataset, "_site_removal_summary.csv")),
            row.names = FALSE)

  message("Finished ", dataset)
  return(list(results = results, summary = summary_table))
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

load_salix_galpar_webs <- function(){
  if(!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
  library(magrittr)
  unlink("Salix/raw-data/csv", recursive = TRUE)
  unlink("Salix/raw-data/rdata", recursive = TRUE)
  source("Salix/lib/format4R.r")
  get_formatData("Salix/raw-data/Salix_webs.csv")
  df_site <- readRDS("Salix/raw-data/rdata/df_site.rds")
  df_interact <- readRDS("Salix/raw-data/rdata/df_interact.rds")
  df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR / df_interact$N_GALLS
  site_interact <- merge(df_site, df_interact, by = "REARING_NUMBER")
  site_interact <- site_interact %>%
    filter(!is.na(SITE), !is.na(RGALLER), !is.na(RPAR),
           RPAR != "none", RPAR != "", RGALLER != "")
  webs <- list()
  for(s in sort(unique(site_interact$SITE))){
    cur <- site_interact %>% filter(SITE == s)
    mat <- xtabs(NB_GALLS_PAR ~ RGALLER + RPAR, data = cur)
    webs[[as.character(s)]] <- clean_web_matrix(mat)
  }
  webs
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

load_dataset_webs <- function(dataset){
  if(dataset == "Quercus") return(load_quercus_webs())
  if(dataset == "Nahuel") return(load_nahuel_webs())
  if(dataset == "Salix_Galpar") return(load_salix_galpar_webs())
  if(dataset == "Gottin_HP") return(load_gottin_HP_webs())
  if(dataset == "Gottin_PP") return(load_gottin_PP_webs())
  if(dataset %in% c("Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny")){
    return(load_garraf_montseny_olot_webs(dataset))
  }
  stop("Unknown dataset: ", dataset)
}

run_one_dataset <- function(dataset){
  if(dataset == "Salix_Galpar"){
    return(run_salix_galpar_original())
  }
  message("Running ", dataset)
  p_fixed <- dataset_config$p_fixed[dataset_config$dataset == dataset]
  dataset_dir <- file.path(sep_out, dataset)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  webs <- load_dataset_webs(dataset)
  site_tables <- build_site_tables_from_webs(webs)
  cooc_triples <- site_tables$cooc_triples
  empirical_site_interactions <- site_tables$empirical_site_interactions
  all_sites <- sort(unique(cooc_triples$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_replicates)
  out_list <- vector("list", nrow(subset_object$index))
  for(i in seq_len(nrow(subset_object$index))){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    exp_row <- expected_links_from_subset(cooc_triples, sites_keep, p_fixed)
    n_realised_links <- observed_links_from_subset(empirical_site_interactions, sites_keep)
    species_error <- species_error_from_subset(cooc_triples, empirical_site_interactions, sites_keep, p_fixed)
    divergence_raw <- n_realised_links - exp_row$n_expected_links_raw
    divergence_conditioned <- n_realised_links - exp_row$n_expected_links_conditioned
    out_list[[i]] <- data.frame(
      dataset = dataset,
      removal_fraction = subset_object$index$removal_fraction[i],
      rep = subset_object$index$rep[i],
      subset_id = subset_id,
      n_sites_kept = subset_object$index$n_sites_kept[i],
      n_realised_links = n_realised_links,
      n_cooc_links = exp_row$n_cooc_links,
      n_expected_links_raw = exp_row$n_expected_links_raw,
      n_expected_links_conditioned = exp_row$n_expected_links_conditioned,
      divergence_raw = divergence_raw,
      relative_divergence_raw = safe_relative_difference(n_realised_links, exp_row$n_expected_links_raw),
      divergence_conditioned = divergence_conditioned,
      relative_divergence_conditioned = safe_relative_difference(n_realised_links, exp_row$n_expected_links_conditioned),
      species_error
    )
  }
  results <- bind_rows(out_list)
  summary_table <- results %>%
    group_by(dataset, removal_fraction) %>%
    summarise(
      mean_realised_links = mean(n_realised_links, na.rm = TRUE),
      mean_expected_links_raw = mean(n_expected_links_raw, na.rm = TRUE),
      mean_expected_links_conditioned = mean(n_expected_links_conditioned, na.rm = TRUE),
      mean_divergence_raw = mean(divergence_raw, na.rm = TRUE),
      mean_relative_divergence_raw = mean(relative_divergence_raw, na.rm = TRUE),
      sd_relative_divergence_raw = sd(relative_divergence_raw, na.rm = TRUE),
      mean_divergence_conditioned = mean(divergence_conditioned, na.rm = TRUE),
      mean_relative_divergence_conditioned = mean(relative_divergence_conditioned, na.rm = TRUE),
      mean_species_bias_raw = mean(mean_species_bias_raw, na.rm = TRUE),
      mean_abs_species_error_raw = mean(mean_abs_species_error_raw, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv2(results, file.path(dataset_dir, paste0(dataset, "_site_removal_divergence.csv")), row.names = FALSE)
  write.csv2(summary_table, file.path(dataset_dir, paste0(dataset, "_site_removal_summary.csv")), row.names = FALSE)
  return(list(results = results, summary = summary_table))
}

all_results <- list()
all_summaries <- list()

for(ds in datasets_to_run){
  cur <- run_one_dataset(ds)
  all_results[[ds]] <- cur$results
  all_summaries[[ds]] <- cur$summary
}

combined_results <- bind_rows(all_results)
combined_summary <- bind_rows(all_summaries)
combined_summary$dataset <- factor(combined_summary$dataset, levels = dataset_config$dataset)

write.csv2(combined_results, file.path(combined_out, "site_removal_divergence_all_datasets.csv"), row.names = FALSE)
write.csv2(combined_summary, file.path(combined_out, "site_removal_summary_all_datasets.csv"), row.names = FALSE)

p_combined_rel <- ggplot(combined_summary, aes(x = removal_fraction, y = mean_relative_divergence_raw)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 12) +
  xlab("Fraction of sites removed") +
  ylab("Mean relative divergence") +
  ggtitle("Site-removal relative divergence across all datasets",
          subtitle = "Raw expected links")

ggsave(file.path(combined_out, "combined_relative_divergence_raw.png"),
       p_combined_rel, width = 14, height = 7, dpi = 300)

p_combined_abs <- ggplot(combined_summary, aes(x = removal_fraction, y = mean_divergence_raw)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 12) +
  xlab("Fraction of sites removed") +
  ylab("Mean observed - expected links") +
  ggtitle("Site-removal signed absolute divergence across all datasets",
          subtitle = "Raw expected links; panels are free y-scale")

ggsave(file.path(combined_out, "combined_absolute_divergence_raw.png"),
       p_combined_abs, width = 14, height = 7, dpi = 300)

links_long <- combined_summary %>%
  select(dataset, removal_fraction, mean_realised_links, mean_expected_links_raw, mean_expected_links_conditioned) %>%
  pivot_longer(cols = c(mean_realised_links, mean_expected_links_raw, mean_expected_links_conditioned),
               names_to = "curve", values_to = "links")

p_combined_links <- ggplot(links_long, aes(x = removal_fraction, y = links, linetype = curve)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 12) +
  xlab("Fraction of sites removed") +
  ylab("Number of links") +
  ggtitle("Observed and expected links under site removal",
          subtitle = "Realised empirical links vs raw and conditioned model expectations")

ggsave(file.path(combined_out, "combined_observed_vs_expected_links.png"),
       p_combined_links, width = 14, height = 7, dpi = 300)

message("All done. Outputs saved in All/SeparatedResults/site_removal_divergence and All/CombinedOutputs.")
