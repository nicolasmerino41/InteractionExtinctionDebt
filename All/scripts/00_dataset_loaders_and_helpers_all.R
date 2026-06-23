
## ------------------------------------------------------------
## Script: All/scripts/00_dataset_loaders_and_helpers_all.R
##
## Purpose:
## Shared loaders and helper functions for the unified 10-dataset pipeline.
##
## Run analysis scripts from the parent repository folder.
## Inputs remain in original folders:
##   Quercus/, Nahuel/, Salix/, Gottin/, Garraf-Montseny-Olot/
## Outputs are written under All/.
## ------------------------------------------------------------

packages <- c("dplyr", "tidyr", "ggplot2", "tibble", "igraph", "bipartite")

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

dataset_config <- data.frame(
  dataset = c(
    "Quercus", "Nahuel", "Salix_Galpar",
    "Gottin_HP", "Gottin_PP",
    "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"
  ),
  p_fixed = c(
    0.0227, 0.065, 0.121,
    0.074, 0.044,
    0.076, 0.100, 0.0962, 0.070, 0.099
  ),
  stringsAsFactors = FALSE
)

all_dataset_names <- dataset_config$dataset

get_p_fixed <- function(dataset){
  p <- dataset_config$p_fixed[dataset_config$dataset == dataset]
  if(length(p) != 1 || is.na(p)) stop("Missing p_fixed for dataset: ", dataset)
  p
}

make_output_dirs <- function(result_type){
  base_out <- "All"
  sep_out <- file.path(base_out, "SeparatedResults", result_type)
  combined_out <- file.path(base_out, "CombinedOutputs")
  dir.create(file.path(base_out, "scripts"), recursive = TRUE, showWarnings = FALSE)
  dir.create(sep_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(combined_out, recursive = TRUE, showWarnings = FALSE)
  list(base = base_out, separated = sep_out, combined = combined_out)
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

safe_relative_difference <- function(observed, expected){
  if(is.na(expected) || expected == 0) return(NA_real_)
  (observed - expected) / expected
}

make_site_subsets <- function(all_sites, removal_levels, n_site_reps){
  site_subsets <- list()
  subset_index <- NULL
  counter <- 1

  for(removal in removal_levels){
    n_total <- length(all_sites)
    n_keep <- max(1, round(n_total * (1 - removal)))
    reps_here <- ifelse(removal == 0, 1, n_site_reps)

    for(r in seq_len(reps_here)){
      sites_keep <- if(removal == 0) all_sites else sample(all_sites, size = n_keep, replace = FALSE)
      site_subsets[[counter]] <- sites_keep
      subset_index <- rbind(
        subset_index,
        data.frame(
          subset_id = counter,
          removal_fraction = removal,
          site_rep = r,
          n_sites_kept = length(sites_keep)
        )
      )
      counter <- counter + 1
    }
  }

  list(index = subset_index, subsets = site_subsets)
}

build_site_tables_from_webs <- function(webs){
  cooc_list <- list()
  interaction_list <- list()
  occupancy_list <- list()
  counter_cooc <- 1
  counter_int <- 1
  counter_occ <- 1

  for(site_id in names(webs)){
    web <- clean_web_matrix(webs[[site_id]])
    if(nrow(web) == 0 || ncol(web) == 0) next

    consumers <- colnames(web)
    resources <- rownames(web)
    consumers <- consumers[!is.na(consumers) & consumers != ""]
    resources <- resources[!is.na(resources) & resources != ""]
    if(length(consumers) == 0 || length(resources) == 0) next

    occupancy_list[[counter_occ]] <- rbind(
      data.frame(site = site_id, species = consumers, trophic_level = "consumer"),
      data.frame(site = site_id, species = resources, trophic_level = "resource")
    )
    counter_occ <- counter_occ + 1

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
    occupancy = bind_rows(occupancy_list) %>% distinct(site, species, trophic_level),
    cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>% distinct(site, consumer, resource)
  )
}

## ---------------------------
## Dataset loaders
## ---------------------------

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
      int_mat <- int_mat[, -which(colnames(int_mat) %in% c("Unpar", "unpar")), drop = FALSE]
    }

    webs[[as.character(t)]] <- clean_web_matrix(int_mat)
  }

  webs
}

load_nahuel_webs <- function(){
  files <- c(
    "vaz_ag_matr_f.txt", "vaz_cl_matr_f.txt", "vaz_ll_matr_f.txt", "vaz_mh_matr_f.txt",
    "vaz_mnh_matr_f.txt", "vaz_qh_matr_f.txt", "vaz_qnh_matr_f.txt", "vaz_s_matr_f.txt"
  )

  vaz <- lapply(files, function(f){
    t(read.csv(file.path("Nahuel/raw-data", f), sep = "\t", header = FALSE))
  })

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
  webs <- bipartite::frame2webs(
    raw_data,
    varnames = c("Genus.Species", "P1.Species.Genus", "Site", "P1.cells")
  )
  names(webs) <- as.character(names(webs))
  lapply(webs, clean_web_matrix)
}

load_gottin_PP_webs <- function(){
  raw_data <- read.csv("Gottin/raw-data/plant_poll_all_interactions.csv", head = TRUE, sep = ";")
  webs <- bipartite::frame2webs(
    raw_data,
    varnames = c("P.Genus.Species", "Genus.Species", "Site")
  )
  names(webs) <- as.character(names(webs))
  lapply(webs, clean_web_matrix)
}

load_garraf_montseny_olot_webs <- function(dataset){
  if(!requireNamespace("XLConnect", quietly = TRUE)){
    stop("XLConnect is required for Garraf-Montseny-Olot. Install binary if possible: install.packages('XLConnect', type = 'binary')")
  }

  file_name <- switch(
    dataset,
    Garraf_HP = "garraf-hp",
    Garraf_PP = "garraf-pp",
    Garraf_PP2 = "garraf-pp-2",
    Olot = "olot",
    Montseny = "montseny"
  )

  excel <- XLConnect::loadWorkbook(paste0("Garraf-Montseny-Olot/raw-data/", file_name, ".xlsx"))
  sheet_names <- XLConnect::getSheets(excel)
  names(sheet_names) <- sheet_names
  sheet_list <- lapply(sheet_names, function(.sheet){
    XLConnect::readWorksheet(object = excel, sheet = .sheet)
  })

  site_sheets <- switch(
    dataset,
    Garraf_HP = unique(sheet_names)[-c(1, 27)],
    Garraf_PP = unique(sheet_names)[-c(1, 42)],
    Garraf_PP2 = unique(sheet_names),
    Olot = unique(sheet_names)[-c(1, 16)],
    Montseny = unique(sheet_names)[-c(1, 20)]
  )

  webs <- list()

  for(p in site_sheets){
    n <- sheet_list[[p]]
    row.names(n) <- n[, 1]
    n <- n[-1]

    if(dataset == "Garraf_PP2"){
      n <- t(n)
    }

    if(dataset == "Montseny"){
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
  if(dataset == "Gottin_HP") return(load_gottin_HP_webs())
  if(dataset == "Gottin_PP") return(load_gottin_PP_webs())
  if(dataset %in% c("Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny")){
    return(load_garraf_montseny_olot_webs(dataset))
  }
  stop("Dataset does not use generic web loader: ", dataset)
}

## ---------------------------
## Salix/Galpar original-style loader
## ---------------------------
load_salix_galpar_site_interact <- function(){
  
  df_site <- readRDS("Salix/raw-data/rdata/df_site.rds")
  df_interact <- readRDS("Salix/raw-data/rdata/df_interact.rds")
  
  df_interact$PAR_RATE <- df_interact$NB_GALLS_PAR /
    df_interact$N_GALLS
  
  merge(df_site,
        df_interact,
        by = "REARING_NUMBER")
}

build_salix_galpar_site_tables <- function(site_interact){
  site_interact <- site_interact %>%
    filter(!is.na(SITE), !is.na(RGALLER), RGALLER != "")

  empirical_site_interactions <- site_interact %>%
    filter(!is.na(RPAR), RPAR != "none", RPAR != "") %>%
    transmute(
      site = as.character(SITE),
      consumer = as.character(RPAR),
      resource = as.character(RGALLER)
    ) %>%
    distinct(site, consumer, resource)

  predators <- sort(unique(empirical_site_interactions$consumer))

  cooc_list <- list()

  for(pred in predators){
    pred_sites <- unique(site_interact$SITE[site_interact$RPAR == pred])

    if(length(pred_sites) == 0){
      next
    }

    cur <- site_interact[site_interact$SITE %in% pred_sites, c("SITE", "RGALLER")]
    cur <- unique(cur)

    cur_cooc <- cur %>%
      transmute(
        site = as.character(SITE),
        consumer = as.character(pred),
        resource = as.character(RGALLER)
      ) %>%
      filter(!is.na(resource), resource != "") %>%
      distinct(site, consumer, resource)

    cooc_list[[pred]] <- cur_cooc
  }

  cooc_triples <- bind_rows(cooc_list) %>%
    distinct(site, consumer, resource)

  occupancy <- bind_rows(
    empirical_site_interactions %>%
      transmute(site, species = consumer, trophic_level = "consumer"),
    cooc_triples %>%
      transmute(site, species = resource, trophic_level = "resource")
  ) %>%
    distinct(site, species, trophic_level)

  list(
    occupancy = occupancy,
    cooc_triples = cooc_triples,
    empirical_site_interactions = empirical_site_interactions
  )
}

get_dataset_site_tables <- function(dataset){
  if(dataset == "Salix_Galpar"){
    site_interact <- load_salix_galpar_site_interact()
    return(build_salix_galpar_site_tables(site_interact))
  }

  webs <- load_dataset_webs(dataset)
  build_site_tables_from_webs(webs)
}

## ---------------------------
## Analysis helpers
## ---------------------------

expected_links_from_subset <- function(cooc_triples, sites_keep, p_fixed){
  cur <- cooc_triples %>%
    filter(site %in% sites_keep) %>%
    group_by(consumer, resource) %>%
    summarise(n_cooc = n_distinct(site), .groups = "drop")

  if(nrow(cur) == 0){
    return(data.frame(
      n_cooc_links = 0,
      expected_links_raw = NA_real_,
      expected_links_conditioned = NA_real_
    ))
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
    expected_links_raw = sum(cur$prob_int, na.rm = TRUE),
    expected_links_conditioned = sum(cur$prob_expected, na.rm = TRUE)
  )
}

observed_links_from_subset <- function(site_interactions, sites_keep){
  site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource) %>%
    nrow()
}

simulate_model_site_interactions <- function(cooc_triples, p_fixed){
  cooc_triples %>%
    mutate(interaction = rbinom(n(), size = 1, prob = p_fixed)) %>%
    filter(interaction == 1) %>%
    select(site, consumer, resource) %>%
    distinct(site, consumer, resource)
}

gini_coefficient <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  if(all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_along(x) - n - 1) * x) / (n * sum(x))
}

degree_table_from_interactions <- function(site_interactions, sites_keep){
  cur <- site_interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(consumer, resource)

  if(nrow(cur) == 0){
    return(data.frame(trophic_level = character(), species = character(), degree = integer()))
  }

  consumer_degrees <- cur %>%
    group_by(species = consumer) %>%
    summarise(degree = n_distinct(resource), .groups = "drop") %>%
    mutate(trophic_level = "consumer")

  resource_degrees <- cur %>%
    group_by(species = resource) %>%
    summarise(degree = n_distinct(consumer), .groups = "drop") %>%
    mutate(trophic_level = "resource")

  bind_rows(consumer_degrees, resource_degrees) %>%
    filter(!is.na(species), species != "", degree > 0) %>%
    select(trophic_level, species, degree)
}

summarise_degrees <- function(degrees){
  if(nrow(degrees) == 0){
    return(data.frame(
      trophic_level = c("consumer", "resource"),
      mean_degree = NA_real_,
      median_degree = NA_real_,
      variance_degree = NA_real_,
      maximum_degree = NA_real_,
      proportion_degree_1 = NA_real_,
      n_species_degree_gt0 = 0L,
      gini_degree = NA_real_
    ))
  }

  degrees %>%
    group_by(trophic_level) %>%
    summarise(
      mean_degree = mean(degree, na.rm = TRUE),
      median_degree = median(degree, na.rm = TRUE),
      variance_degree = ifelse(n() > 1, var(degree, na.rm = TRUE), 0),
      maximum_degree = max(degree, na.rm = TRUE),
      proportion_degree_1 = mean(degree == 1, na.rm = TRUE),
      n_species_degree_gt0 = n_distinct(species),
      gini_degree = gini_coefficient(degree),
      .groups = "drop"
    ) %>%
    right_join(data.frame(trophic_level = c("consumer", "resource")), by = "trophic_level")
}

frequency_degrees <- function(degrees){
  if(nrow(degrees) == 0){
    return(data.frame(trophic_level = character(), degree = integer(), n_species = integer()))
  }

  degrees %>%
    count(trophic_level, degree, name = "n_species") %>%
    arrange(trophic_level, degree)
}

make_nocc_bin <- function(x){
  cut(
    x,
    breaks = c(0, 1, 2, 5, 10, 20, Inf),
    labels = c("1", "2", "3-5", "6-10", "11-20", "21+"),
    right = TRUE
  )
}
