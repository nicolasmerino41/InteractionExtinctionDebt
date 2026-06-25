## ------------------------------------------------------------
## Script: All/scripts/19_describe_layers_of_interaction_loss.R
##
## Purpose:
## Descriptive site-removal analysis of four layers of interaction loss.
## No fitted models, no beta-binomial model, no model predictions.
##
## Uses the same dataset-loading style as the original site-removal script:
##   - load site webs
##   - clean binary matrices
##   - build co-occurrence triples as site x consumer x resource records
##   - build empirical site interactions as positive matrix cells
##
## Outputs:
##   All/CombinedOutputs/19_*.csv
##   All/CombinedOutputs/19_*.png
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

## Same removal levels as earlier site-removal scripts.
removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
n_replicates <- 500

base_out <- "All"
combined_out <- file.path(base_out, "CombinedOutputs")
dir.create(combined_out, recursive = TRUE, showWarnings = FALSE)

dataset_config <- data.frame(
  dataset = c("Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
              "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"),
  stringsAsFactors = FALSE
)

datasets_to_run <- dataset_config$dataset

## ---------------------------
## Original-style helpers
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

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
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

    cur_cooc <- expand.grid(
      consumer = consumers,
      resource = resources,
      stringsAsFactors = FALSE
    )
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

  list(
    cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource),
    empirical_site_interactions = bind_rows(interaction_list) %>% distinct(site, consumer, resource)
  )
}

## ---------------------------
## Dataset loaders copied from original style
## ---------------------------

load_quercus_webs <- function(){
  if(!requireNamespace("XLConnect", quietly = TRUE)){
    stop("XLConnect is required for Quercus.")
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

## For Salix/Galpar, use the original special logic rather than the generic web loader.
## Consumer = parasitoid RPAR; resource = galler RGALLER.
## A consumer can only be operationally paired with resources at sites where that consumer was observed.
build_salix_site_tables_original <- function(){
  site_interact <- load_salix_galpar_site_interact()
  site_interact <- site_interact %>%
    filter(!is.na(SITE), !is.na(RGALLER), RGALLER != "")

  real_rows <- site_interact %>%
    filter(!is.na(RPAR), RPAR != "none", RPAR != "") %>%
    transmute(site = as.character(SITE), consumer = as.character(RPAR), resource = as.character(RGALLER)) %>%
    distinct(site, consumer, resource)

  consumers <- sort(unique(real_rows$consumer))
  cooc_list <- list()

  for(cons in consumers){
    cons_sites <- unique(real_rows$site[real_rows$consumer == cons])
    cur <- site_interact %>%
      filter(as.character(SITE) %in% cons_sites) %>%
      transmute(site = as.character(SITE), consumer = cons, resource = as.character(RGALLER)) %>%
      filter(!is.na(resource), resource != "") %>%
      distinct(site, consumer, resource)
    cooc_list[[cons]] <- cur
  }

  list(
    cooc_triples = bind_rows(cooc_list) %>% distinct(site, consumer, resource),
    empirical_site_interactions = real_rows
  )
}

load_garraf_montseny_olot_webs <- function(cur_dataset){
  if(!requireNamespace("XLConnect", quietly = TRUE)){
    stop("XLConnect is required for Garraf-Montseny-Olot.")
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
  if(dataset == "Gottin_HP") return(load_gottin_HP_webs())
  if(dataset == "Gottin_PP") return(load_gottin_PP_webs())
  if(dataset %in% c("Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny")){
    return(load_garraf_montseny_olot_webs(dataset))
  }
  stop("Unknown dataset: ", dataset)
}

get_site_tables_original_style <- function(dataset){
  if(dataset == "Salix_Galpar"){
    return(build_salix_site_tables_original())
  }
  build_site_tables_from_webs(load_dataset_webs(dataset))
}

## ---------------------------
## Core analysis
## ---------------------------

summarise_layers_one_subset <- function(dataset, cooc, ints, original_links, sites_keep, removal_fraction, replicate){
  cooc_ret <- cooc %>% filter(site %in% sites_keep)
  int_ret <- ints %>% filter(site %in% sites_keep)

  retained_cooc_records <- nrow(cooc_ret)
  retained_potential_pairs <- nrow(distinct(cooc_ret, pair_id))
  retained_interaction_records <- nrow(int_ret)
  retained_regional_links <- nrow(distinct(int_ret, pair_id))

  full_cooc_records <- attr(cooc, "full_cooc_records")
  full_potential_pairs <- attr(cooc, "full_potential_pairs")
  full_interaction_records <- attr(ints, "full_interaction_records")
  full_regional_links <- attr(ints, "full_regional_links")

  layer_rows <- data.frame(
    dataset = dataset,
    removal_fraction = removal_fraction,
    replicate = replicate,
    layer = c("Co-occurrence records", "Potential pairs", "Interaction records", "Regional interaction links"),
    full_count = c(full_cooc_records, full_potential_pairs, full_interaction_records, full_regional_links),
    retained_count = c(retained_cooc_records, retained_potential_pairs, retained_interaction_records, retained_regional_links),
    stringsAsFactors = FALSE
  ) %>%
    mutate(fraction_remaining = ifelse(full_count > 0, retained_count / full_count, NA_real_))

  retained_cooc_by_link <- cooc_ret %>%
    filter(pair_id %in% original_links$pair_id) %>%
    count(pair_id, name = "retained_n")

  retained_int_by_link <- int_ret %>%
    filter(pair_id %in% original_links$pair_id) %>%
    count(pair_id, name = "retained_K")

  link_level <- original_links %>%
    left_join(retained_cooc_by_link, by = "pair_id") %>%
    left_join(retained_int_by_link, by = "pair_id") %>%
    mutate(
      retained_n = ifelse(is.na(retained_n), 0L, retained_n),
      retained_K = ifelse(is.na(retained_K), 0L, retained_K),
      removal_fraction = removal_fraction,
      replicate = replicate,
      link_state = case_when(
        retained_K > 0 ~ "Interaction still observed",
        retained_n > 0 & retained_K == 0 ~ "Recorded together, interaction absent",
        retained_n == 0 ~ "Pair no longer recorded together"
      )
    ) %>%
    select(dataset, consumer, resource, pair_id, full_n, full_K, full_repeatability,
           removal_fraction, replicate, retained_n, retained_K, link_state, initial_support_group)

  state_rows <- link_level %>%
    count(dataset, removal_fraction, replicate, link_state, name = "number_in_state") %>%
    tidyr::complete(
      dataset,
      removal_fraction,
      replicate,
      link_state = c("Interaction still observed",
                     "Recorded together, interaction absent",
                     "Pair no longer recorded together"),
      fill = list(number_in_state = 0L)
    ) %>%
    mutate(
      number_of_original_regional_links = nrow(original_links),
      fraction_in_state = number_in_state / number_of_original_regional_links
    )

  survival_rows <- link_level %>%
    group_by(dataset, removal_fraction, replicate, initial_support_group) %>%
    summarise(
      number_original_links_in_group = n(),
      number_still_observed = sum(retained_K > 0),
      fraction_still_observed = number_still_observed / number_original_links_in_group,
      .groups = "drop"
    )

  list(layers = layer_rows, states = state_rows, survival = survival_rows, link_level = link_level)
}

run_one_dataset <- function(dataset){
  message("Running descriptive layers: ", dataset)

  site_tables <- get_site_tables_original_style(dataset)

  cooc <- site_tables$cooc_triples %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource)) %>%
    distinct(site, consumer, resource, pair_id)

  ints <- site_tables$empirical_site_interactions %>%
    mutate(site = as.character(site),
           consumer = as.character(consumer),
           resource = as.character(resource),
           pair_id = make_pair_id(consumer, resource)) %>%
    distinct(site, consumer, resource, pair_id)

  ## Defensive but original-compatible: an observed interaction cell is treated as a recorded-together cell.
  ## This should normally add zero rows for generic web-loaded datasets.
  cooc <- bind_rows(cooc, ints %>% select(site, consumer, resource, pair_id)) %>%
    distinct(site, consumer, resource, pair_id)

  attr(cooc, "full_cooc_records") <- nrow(cooc)
  attr(cooc, "full_potential_pairs") <- nrow(distinct(cooc, pair_id))
  attr(ints, "full_interaction_records") <- nrow(ints)
  attr(ints, "full_regional_links") <- nrow(distinct(ints, pair_id))

  full_n <- cooc %>%
    count(consumer, resource, pair_id, name = "full_n")

  full_K <- ints %>%
    count(consumer, resource, pair_id, name = "full_K")

  original_links <- full_K %>%
    left_join(full_n, by = c("consumer", "resource", "pair_id")) %>%
    mutate(
      dataset = dataset,
      full_repeatability = full_K / full_n,
      initial_support_group = case_when(
        full_K == 1 ~ "Observed in 1 site",
        full_K == 2 ~ "Observed in 2 sites",
        full_K %in% 3:4 ~ "Observed in 3–4 sites",
        full_K >= 5 ~ "Observed in 5 or more sites"
      )
    ) %>%
    select(dataset, consumer, resource, pair_id, full_n, full_K, full_repeatability, initial_support_group)

  all_sites <- sort(unique(cooc$site))
  subset_object <- make_site_subsets(all_sites, removal_levels, n_replicates)

  out <- lapply(seq_len(nrow(subset_object$index)), function(i){
    subset_id <- subset_object$index$subset_id[i]
    sites_keep <- subset_object$subsets[[subset_id]]
    summarise_layers_one_subset(
      dataset = dataset,
      cooc = cooc,
      ints = ints,
      original_links = original_links,
      sites_keep = sites_keep,
      removal_fraction = subset_object$index$removal_fraction[i],
      replicate = subset_object$index$replicate[i]
    )
  })

  layers <- bind_rows(lapply(out, `[[`, "layers"))
  states <- bind_rows(lapply(out, `[[`, "states"))
  survival <- bind_rows(lapply(out, `[[`, "survival"))
  link_level <- bind_rows(lapply(out, `[[`, "link_level"))

  ## Keep support groups only if the full dataset has at least 10 original links in that group.
  eligible_groups <- original_links %>%
    count(dataset, initial_support_group, name = "number_original_links_in_group_full") %>%
    filter(number_original_links_in_group_full >= 10)

  survival <- survival %>%
    semi_join(eligible_groups, by = c("dataset", "initial_support_group"))

  ## Validation summary for this dataset. No model is fitted or used.
  zero_layers_ok <- layers %>%
    filter(removal_fraction == 0) %>%
    summarise(ok = all(abs(fraction_remaining - 1) < 1e-12, na.rm = TRUE)) %>%
    pull(ok)

  layer_wide <- layers %>%
    select(dataset, removal_fraction, replicate, layer, retained_count, full_count) %>%
    pivot_wider(names_from = layer, values_from = c(retained_count, full_count))

  int_records_le_cooc <- all(layer_wide$`retained_count_Interaction records` <= layer_wide$`retained_count_Co-occurrence records`, na.rm = TRUE)
  links_le_pairs <- all(layer_wide$`retained_count_Regional interaction links` <= layer_wide$`retained_count_Potential pairs`, na.rm = TRUE)

  state_sum_ok <- states %>%
    group_by(dataset, removal_fraction, replicate) %>%
    summarise(total_fraction = sum(fraction_in_state), total_n = sum(number_in_state), .groups = "drop") %>%
    summarise(ok = all(abs(total_fraction - 1) < 1e-12) && all(total_n == nrow(original_links))) %>%
    pull(ok)

  counts_le_full <- layers %>%
    summarise(ok = all(retained_count <= full_count, na.rm = TRUE)) %>%
    pull(ok)

  link_state_exactly_one <- !any(is.na(link_level$link_state)) &&
    nrow(link_level) == nrow(original_links) * nrow(subset_object$index)

  checks <- data.frame(
    dataset = dataset,
    n_sites = length(all_sites),
    n_subsets = nrow(subset_object$index),
    n_original_regional_links = nrow(original_links),
    zero_removal_all_retention_equals_1 = zero_layers_ok,
    interaction_records_never_exceed_cooccurrence_records = int_records_le_cooc,
    regional_links_never_exceed_potential_pairs = links_le_pairs,
    each_original_link_assigned_exactly_one_state = link_state_exactly_one,
    three_state_proportions_sum_to_1 = state_sum_ok,
    all_retained_counts_le_full_counts = counts_le_full,
    no_models_fitted_or_used = TRUE,
    stringsAsFactors = FALSE
  )

  list(layers = layers, states = states, survival = survival, link_level = link_level, checks = checks)
}

## ---------------------------
## Run sequentially: this is descriptive and avoids future overhead.
## ---------------------------

all_outputs <- lapply(datasets_to_run, run_one_dataset)
names(all_outputs) <- datasets_to_run

layers_all <- bind_rows(lapply(all_outputs, `[[`, "layers"))
states_all <- bind_rows(lapply(all_outputs, `[[`, "states"))
survival_all <- bind_rows(lapply(all_outputs, `[[`, "survival"))
link_level_all <- bind_rows(lapply(all_outputs, `[[`, "link_level"))
checks_all <- bind_rows(lapply(all_outputs, `[[`, "checks"))

## ---------------------------
## Save tables
## ---------------------------

write.csv2(layers_all, file.path(combined_out, "19_layers_under_site_removal_summary.csv"), row.names = FALSE)
write.csv2(states_all, file.path(combined_out, "19_link_states_under_site_removal.csv"), row.names = FALSE)
write.csv2(survival_all, file.path(combined_out, "19_link_survival_by_initial_local_support.csv"), row.names = FALSE)
write.csv2(link_level_all %>% select(-pair_id, -initial_support_group),
           file.path(combined_out, "19_link_level_site_removal_data.csv"), row.names = FALSE)
write.csv2(checks_all, file.path(combined_out, "19_layers_under_site_removal_checks.csv"), row.names = FALSE)

cat("\nValidation checks for script 19:\n")
print(checks_all)

if(any(!checks_all$zero_removal_all_retention_equals_1) ||
   any(!checks_all$interaction_records_never_exceed_cooccurrence_records) ||
   any(!checks_all$regional_links_never_exceed_potential_pairs) ||
   any(!checks_all$each_original_link_assigned_exactly_one_state) ||
   any(!checks_all$three_state_proportions_sum_to_1) ||
   any(!checks_all$all_retained_counts_le_full_counts)){
  warning("At least one validation check failed. See 19_layers_under_site_removal_checks.csv")
}

## ---------------------------
## Plot summaries
## ---------------------------

layer_summary <- layers_all %>%
  group_by(dataset, removal_fraction, layer) %>%
  summarise(
    median_fraction = median(fraction_remaining, na.rm = TRUE),
    q025 = quantile(fraction_remaining, 0.025, na.rm = TRUE),
    q975 = quantile(fraction_remaining, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(dataset = factor(dataset, levels = datasets_to_run))

layer_colours <- c(
  "Co-occurrence records" = "#1b9e77",
  "Potential pairs" = "#d95f02",
  "Interaction records" = "#7570b3",
  "Regional interaction links" = "#e7298a"
)

p_layers <- ggplot(layer_summary,
                   aes(x = removal_fraction, y = median_fraction, colour = layer, fill = layer)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = layer_colours, name = "Layer") +
  scale_fill_manual(values = layer_colours, name = "Layer") +
  facet_wrap(~ dataset, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Fraction remaining")

ggsave(file.path(combined_out, "19_layers_under_site_removal.png"),
       p_layers, width = 14, height = 7, dpi = 300)

state_summary <- states_all %>%
  group_by(dataset, removal_fraction, link_state) %>%
  summarise(fraction_in_state = median(fraction_in_state, na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset, removal_fraction) %>%
  mutate(fraction_in_state = fraction_in_state / sum(fraction_in_state, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    dataset = factor(dataset, levels = datasets_to_run),
    link_state = factor(link_state,
                        levels = c("Pair no longer recorded together",
                                   "Recorded together, interaction absent",
                                   "Interaction still observed"))
  )

state_colours <- c(
  "Interaction still observed" = "#1b9e77",
  "Recorded together, interaction absent" = "#d95f02",
  "Pair no longer recorded together" = "#7570b3"
)

p_states <- ggplot(state_summary,
                   aes(x = removal_fraction, y = fraction_in_state, fill = link_state)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = state_colours, name = "Link state") +
  facet_wrap(~ dataset, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Fraction of original regional interaction links")

ggsave(file.path(combined_out, "19_link_states_under_site_removal.png"),
       p_states, width = 14, height = 7, dpi = 300)

survival_summary <- survival_all %>%
  group_by(dataset, removal_fraction, initial_support_group) %>%
  summarise(
    median_fraction = median(fraction_still_observed, na.rm = TRUE),
    q025 = quantile(fraction_still_observed, 0.025, na.rm = TRUE),
    q975 = quantile(fraction_still_observed, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dataset = factor(dataset, levels = datasets_to_run),
    initial_support_group = factor(initial_support_group,
                                   levels = c("Observed in 1 site",
                                              "Observed in 2 sites",
                                              "Observed in 3–4 sites",
                                              "Observed in 5 or more sites"))
  )

support_colours <- c(
  "Observed in 1 site" = "#1b9e77",
  "Observed in 2 sites" = "#d95f02",
  "Observed in 3–4 sites" = "#7570b3",
  "Observed in 5 or more sites" = "#e7298a"
)

p_survival <- ggplot(survival_summary,
                     aes(x = removal_fraction, y = median_fraction,
                         colour = initial_support_group, fill = initial_support_group)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = support_colours, name = "Initial support") +
  scale_fill_manual(values = support_colours, name = "Initial support") +
  facet_wrap(~ dataset, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  xlab("Proportion of sites removed") +
  ylab("Fraction of original links still observed")

ggsave(file.path(combined_out, "19_link_survival_by_initial_local_support.png"),
       p_survival, width = 14, height = 7, dpi = 300)

message("Finished script 19. Outputs saved in ", combined_out)
