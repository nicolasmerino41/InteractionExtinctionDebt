## ------------------------------------------------------------
## Script: Sabatino/scripts/01_sierra_level_basic_cooccurrence_analysis.R
##
## Purpose:
## First-pass Galiana-style co-occurrence/site-removal analysis
## for Sabatino plant-pollinator data.
##
## Input:
##   Sabatino/Network 2008_updated.xlsx
##   sheet: Matrix complete_2008_marzo2024
##
## Output:
##   Sabatino/outputs/01_sierra_level_basic_cooccurrence_analysis/
##
## Notes:
##   Site = Sierra.
##   Co-occurrence is operational: plant recorded in Sierra + valid insect recorded in Sierra.
##   SV / blank / missing / suspicious insect labels are excluded from primary networks.
##   Repeated records are collapsed to binary links for the primary analysis.
## ------------------------------------------------------------
packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "lubridate", "stringr", "purrr", "tibble")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

input_file <- "Sabatino/Network 2008_updated.xlsx"
input_sheet <- "Matrix complete_2008_marzo2024"

sv_handling <- "exclude"

min_sites_per_date <- 6
min_empirical_links_per_date <- 15

target_removal_fractions <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)

max_exact_subsets <- Inf
n_random_subsets <- 1000
min_sites_retained <- 2

n_model_reps <- 50

output_dir <- "Sabatino/outputs/01_sierra_level_basic_cooccurrence_analysis"
dir.create("Sabatino/scripts", recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_date <- function(x){
  out <- suppressWarnings(lubridate::ymd(x))
  if(all(is.na(out))) out <- suppressWarnings(as.Date(x))
  out
}

trim_to_na <- function(x){
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

safe_relative_error_observed_denominator <- function(observed, expected){
  if(is.na(observed) || observed == 0) return(NA_real_)
  (observed - expected) / observed
}

safe_relative_error_expected_denominator <- function(observed, expected){
  if(is.na(expected) || expected == 0) return(NA_real_)
  (observed - expected) / expected
}

is_suspicious_insect_label <- function(x){
  x_clean <- tolower(stringr::str_trim(as.character(x)))
  is.na(x_clean) |
    x_clean == "" |
    x_clean %in% c("sv", "s.v.", "s v", "na", "n/a", "none", "sin visita",
                   "no visitor", "no visitors", "unknown", "unidentified")
}

is_valid_insect_label <- function(x) !is_suspicious_insect_label(x)

make_pair_id <- function(plant, insect) paste(plant, insect, sep = "___")

expected_links_from_p <- function(cooc_counts, p){
  if(nrow(cooc_counts) == 0 || is.na(p)) return(NA_real_)
  sum(1 - (1 - p)^cooc_counts$n_cooccurring_sierras, na.rm = TRUE)
}

fit_p_calibrated <- function(cooc_counts, observed_links){
  if(nrow(cooc_counts) == 0 || is.na(observed_links)) return(NA_real_)
  if(observed_links <= 0) return(0)
  max_links <- nrow(cooc_counts)
  if(observed_links >= max_links) return(1)
  f <- function(p) expected_links_from_p(cooc_counts, p) - observed_links
  tryCatch(uniroot(f, interval = c(0, 1), tol = 1e-10)$root,
           error = function(e) NA_real_)
}

nearest_feasible_removal_levels <- function(n_sites, target_fractions, min_sites_retained = 2){
  n_removed_target <- round(n_sites * target_fractions)
  n_removed_target <- pmin(n_removed_target, n_sites - min_sites_retained)
  n_removed_target <- pmax(n_removed_target, 0)
  data.frame(
    target_removal_fraction = target_fractions,
    n_removed = n_removed_target,
    actual_removal_fraction = n_removed_target / n_sites,
    n_sites_retained = n_sites - n_removed_target
  ) %>% distinct(n_removed, .keep_all = TRUE)
}

generate_sierra_subsets <- function(sites, analysis_id){
  sites <- sort(unique(as.character(sites)))
  n_sites <- length(sites)
  feasible <- nearest_feasible_removal_levels(n_sites, target_removal_fractions, min_sites_retained)
  subset_list <- list()
  subset_index <- list()
  counter <- 1

  for(i in seq_len(nrow(feasible))){
    n_retained <- feasible$n_sites_retained[i]
    n_combinations <- choose(n_sites, n_retained)
    use_exact <- is.finite(n_combinations) && n_combinations <= max_exact_subsets

    if(use_exact){
      kept_sets <- combn(sites, n_retained, simplify = FALSE)
      subset_source <- "exact"
    } else {
      kept_sets <- replicate(n_random_subsets, sort(sample(sites, n_retained, replace = FALSE)), simplify = FALSE)
      kept_sets <- unique(kept_sets)
      subset_source <- "random"
    }

    for(j in seq_along(kept_sets)){
      kept <- sort(as.character(kept_sets[[j]]))
      subset_list[[counter]] <- kept
      subset_index[[counter]] <- data.frame(
        analysis_id = analysis_id,
        subset_id = counter,
        target_removal_fraction = feasible$target_removal_fraction[i],
        actual_removal_fraction = feasible$actual_removal_fraction[i],
        n_removed = n_sites - length(kept),
        n_sites_retained = length(kept),
        subset_source = subset_source
      )
      counter <- counter + 1
    }
  }
  list(subsets = subset_list, index = bind_rows(subset_index))
}

## ---------------------------
## 1. Read and audit raw data
## ---------------------------

raw <- readxl::read_excel(input_file, sheet = input_sheet) %>%
  select(!starts_with("..."))

required_cols <- c("Sierra", "Date", "Plant_species", "Number_of_flowers",
                   "Insect_species", "Insect_abundance", "Plant_ID")
missing_cols <- setdiff(required_cols, names(raw))
if(length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

raw_clean <- raw %>%
  mutate(
    row_id = row_number(),
    Sierra_original = Sierra,
    Date_original = Date,
    Plant_species_original = Plant_species,
    Insect_species_original = Insect_species,
    Sierra = trim_to_na(Sierra),
    Date_chr = trim_to_na(Date),
    Date_parsed = safe_date(Date_chr),
    Date_label = ifelse(is.na(Date_parsed), Date_chr, as.character(Date_parsed)),
    Plant_species = trim_to_na(Plant_species),
    Insect_species = trim_to_na(Insect_species),
    Plant_ID_chr = trim_to_na(Plant_ID),
    Number_of_flowers_numeric = safe_numeric(Number_of_flowers),
    Insect_abundance_numeric = safe_numeric(Insect_abundance),
    suspicious_insect_label = is_suspicious_insect_label(Insect_species),
    valid_insect_label = is_valid_insect_label(Insect_species),
    valid_interaction_record = !is.na(Sierra) &
      !is.na(Date_label) &
      !is.na(Plant_species) &
      !is.na(Insect_species) &
      valid_insect_label
  )

write.csv2(raw_clean, file.path(output_dir, "01_cleaned_raw_records_with_flags.csv"), row.names = FALSE)

dataset_overview <- data.frame(
  n_rows_raw = nrow(raw),
  n_rows_clean = nrow(raw_clean),
  n_sierras = n_distinct(raw_clean$Sierra, na.rm = TRUE),
  n_dates = n_distinct(raw_clean$Date_label, na.rm = TRUE),
  n_plant_species_raw = n_distinct(raw_clean$Plant_species, na.rm = TRUE),
  n_insect_labels_raw = n_distinct(raw_clean$Insect_species, na.rm = TRUE),
  n_valid_insect_species = n_distinct(raw_clean$Insect_species[raw_clean$valid_insect_label], na.rm = TRUE),
  n_unique_plant_ids = n_distinct(raw_clean$Plant_ID_chr, na.rm = TRUE),
  n_rows_suspicious_insect = sum(raw_clean$suspicious_insect_label, na.rm = TRUE),
  proportion_rows_suspicious_insect = mean(raw_clean$suspicious_insect_label, na.rm = TRUE),
  n_rows_missing_sierra = sum(is.na(raw_clean$Sierra)),
  n_rows_missing_date = sum(is.na(raw_clean$Date_label)),
  n_rows_missing_plant = sum(is.na(raw_clean$Plant_species)),
  n_rows_missing_insect = sum(is.na(raw_clean$Insect_species)),
  n_rows_non_numeric_flowers = sum(!is.na(raw_clean$Number_of_flowers) & is.na(raw_clean$Number_of_flowers_numeric)),
  n_rows_zero_or_negative_flowers = sum(!is.na(raw_clean$Number_of_flowers_numeric) & raw_clean$Number_of_flowers_numeric <= 0),
  n_rows_non_numeric_insect_abundance = sum(!is.na(raw_clean$Insect_abundance) & is.na(raw_clean$Insect_abundance_numeric)),
  n_rows_zero_or_negative_insect_abundance = sum(!is.na(raw_clean$Insect_abundance_numeric) & raw_clean$Insect_abundance_numeric <= 0),
  n_valid_interaction_records_primary = sum(raw_clean$valid_interaction_record, na.rm = TRUE)
)

write.csv2(dataset_overview, file.path(output_dir, "01_dataset_overview.csv"), row.names = FALSE)

sv_audit <- raw_clean %>%
  filter(suspicious_insect_label | is.na(Insect_species)) %>%
  group_by(Insect_species_original, Insect_species) %>%
  summarise(
    n_rows = n(),
    n_sierras = n_distinct(Sierra, na.rm = TRUE),
    n_dates = n_distinct(Date_label, na.rm = TRUE),
    insect_abundance_min = min(Insect_abundance_numeric, na.rm = TRUE),
    insect_abundance_mean = mean(Insect_abundance_numeric, na.rm = TRUE),
    insect_abundance_median = median(Insect_abundance_numeric, na.rm = TRUE),
    insect_abundance_max = max(Insect_abundance_numeric, na.rm = TRUE),
    n_valid_plant_species = sum(!is.na(Plant_species)),
    n_valid_plant_id = sum(!is.na(Plant_ID_chr)),
    n_valid_flower_counts = sum(!is.na(Number_of_flowers_numeric)),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA_real_, .x)))

write.csv2(sv_audit, file.path(output_dir, "01_SV_audit.csv"), row.names = FALSE)

repeated_audit <- raw_clean %>%
  group_by(Sierra, Date_label, Plant_species, Insect_species) %>%
  summarise(
    n_raw_records = n(),
    n_distinct_plant_ids = n_distinct(Plant_ID_chr, na.rm = TRUE),
    insect_abundance_sum = sum(Insect_abundance_numeric, na.rm = TRUE),
    insect_abundance_mean = mean(Insect_abundance_numeric, na.rm = TRUE),
    insect_abundance_max = max(Insect_abundance_numeric, na.rm = TRUE),
    n_distinct_flower_counts = n_distinct(Number_of_flowers_numeric, na.rm = TRUE),
    flower_count_consistent = n_distinct_flower_counts <= 1,
    .groups = "drop"
  ) %>%
  filter(n_raw_records > 1) %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA_real_, .x)))

write.csv2(repeated_audit, file.path(output_dir, "01_repeated_record_audit.csv"), row.names = FALSE)

message("Repeated plant-insect records within Sierra-Date: ", nrow(repeated_audit))
message("Repeated records with >1 Plant_ID: ", sum(repeated_audit$n_distinct_plant_ids > 1, na.rm = TRUE))
message("Repeated records with inconsistent flower counts: ", sum(!repeated_audit$flower_count_consistent, na.rm = TRUE))
message("Plant_ID patterns are reported only; plot identity is not inferred.")

## ---------------------------
## 2. Coverage and collapsed binary records
## ---------------------------

valid_records <- raw_clean %>% filter(valid_interaction_record)

sierra_date_coverage <- raw_clean %>%
  group_by(Sierra, Date_label) %>%
  summarise(
    raw_row_count = n(),
    valid_interaction_record_count = sum(valid_interaction_record, na.rm = TRUE),
    n_plant_species = n_distinct(Plant_species, na.rm = TRUE),
    n_valid_insect_species = n_distinct(Insect_species[valid_insect_label], na.rm = TRUE),
    n_unique_observed_pairs = n_distinct(make_pair_id(Plant_species[valid_interaction_record],
                                                      Insect_species[valid_interaction_record])),
    n_distinct_plant_ids = n_distinct(Plant_ID_chr, na.rm = TRUE),
    total_insect_abundance = sum(Insect_abundance_numeric[valid_interaction_record], na.rm = TRUE),
    .groups = "drop"
  )

flowers_dedup <- raw_clean %>%
  filter(!is.na(Sierra), !is.na(Date_label), !is.na(Plant_ID_chr)) %>%
  group_by(Sierra, Date_label, Plant_ID_chr, Plant_species) %>%
  summarise(flowers_per_plant_id = max(Number_of_flowers_numeric, na.rm = TRUE), .groups = "drop") %>%
  mutate(flowers_per_plant_id = ifelse(is.infinite(flowers_per_plant_id), NA_real_, flowers_per_plant_id)) %>%
  group_by(Sierra, Date_label) %>%
  summarise(total_flowers_deduplicated_by_plant_id = sum(flowers_per_plant_id, na.rm = TRUE), .groups = "drop")

sierra_date_coverage <- sierra_date_coverage %>%
  left_join(flowers_dedup, by = c("Sierra", "Date_label"))

write.csv2(sierra_date_coverage, file.path(output_dir, "01_sierra_date_coverage.csv"), row.names = FALSE)

date_coverage_summary <- sierra_date_coverage %>%
  group_by(Date_label) %>%
  summarise(
    n_sierras_sampled = n_distinct(Sierra, na.rm = TRUE),
    total_raw_rows = sum(raw_row_count, na.rm = TRUE),
    total_valid_interaction_records = sum(valid_interaction_record_count, na.rm = TRUE),
    mean_valid_interaction_records_per_sierra = mean(valid_interaction_record_count, na.rm = TRUE),
    .groups = "drop"
  )

write.csv2(date_coverage_summary, file.path(output_dir, "01_date_coverage_summary.csv"), row.names = FALSE)

cleaned_interaction_records <- valid_records %>%
  group_by(Sierra, Date_label, Plant_species, Insect_species) %>%
  summarise(
    interaction_observed = 1L,
    n_raw_records = n(),
    n_distinct_plant_ids = n_distinct(Plant_ID_chr, na.rm = TRUE),
    interaction_abundance_sum = sum(Insect_abundance_numeric, na.rm = TRUE),
    interaction_abundance_mean = mean(Insect_abundance_numeric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pair_id = make_pair_id(Plant_species, Insect_species))

write.csv2(cleaned_interaction_records, file.path(output_dir, "01_cleaned_interaction_records.csv"), row.names = FALSE)

## ---------------------------
## 3. Build pooled and date-stratified analysis units
## ---------------------------

build_analysis_unit <- function(records, analysis_id, analysis_type){
  records <- records %>% filter(!is.na(Sierra), !is.na(Plant_species), !is.na(Insect_species))

  plant_by_site <- records %>% distinct(site = Sierra, plant = Plant_species)
  insect_by_site <- records %>% distinct(site = Sierra, insect = Insect_species)

  cooc <- plant_by_site %>%
    inner_join(insect_by_site, by = "site") %>%
    transmute(analysis_id, analysis_type, site, plant, insect, pair_id = make_pair_id(plant, insect)) %>%
    distinct()

  interactions <- records %>%
    transmute(analysis_id, analysis_type, site = Sierra, plant = Plant_species,
              insect = Insect_species, pair_id = make_pair_id(plant, insect)) %>%
    distinct()

  cooc_counts <- cooc %>%
    group_by(plant, insect, pair_id) %>%
    summarise(n_cooccurring_sierras = n_distinct(site), .groups = "drop")

  interaction_counts <- interactions %>%
    group_by(plant, insect, pair_id) %>%
    summarise(n_interacting_sierras = n_distinct(site), .groups = "drop")

  list(cooc = cooc, interactions = interactions, cooc_counts = cooc_counts,
       interaction_counts = interaction_counts, sites = sort(unique(cooc$site)))
}

pooled_unit <- build_analysis_unit(cleaned_interaction_records, "all_dates_pooled", "pooled")

date_units <- list()
skipped_dates <- list()

for(d in sort(unique(cleaned_interaction_records$Date_label))){
  cur <- cleaned_interaction_records %>% filter(Date_label == d)
  tmp <- build_analysis_unit(cur, d, "date_stratified")
  n_sites <- length(tmp$sites)
  n_links <- n_distinct(tmp$interactions$pair_id)
  reasons <- c()
  if(n_sites < min_sites_per_date) reasons <- c(reasons, "too_few_sierras")
  if(n_links < min_empirical_links_per_date) reasons <- c(reasons, "too_few_empirical_links")

  if(length(reasons) > 0){
    skipped_dates[[d]] <- data.frame(analysis_id = d, analysis_type = "date_stratified",
                                     n_sierras = n_sites, n_empirical_links = n_links,
                                     skipped_reason = paste(reasons, collapse = ";"))
  } else {
    date_units[[d]] <- tmp
  }
}

analysis_units <- c(list(all_dates_pooled = pooled_unit), date_units)
skipped_date_table <- bind_rows(skipped_dates)
if(nrow(skipped_date_table) == 0){
  skipped_date_table <- data.frame(analysis_id = character(), analysis_type = character(),
                                   n_sierras = integer(), n_empirical_links = integer(),
                                   skipped_reason = character())
}
write.csv2(skipped_date_table, file.path(output_dir, "01_skipped_date_stratified_analyses.csv"), row.names = FALSE)

analysis_unit_summary <- bind_rows(lapply(names(analysis_units), function(id){
  u <- analysis_units[[id]]
  n_sites <- length(u$sites)
  n_plants <- n_distinct(u$cooc$plant)
  n_insects <- n_distinct(u$cooc$insect)
  n_cooc_pairs <- n_distinct(u$cooc$pair_id)
  n_links <- n_distinct(u$interactions$pair_id)
  data.frame(
    analysis_id = id,
    analysis_type = unique(u$cooc$analysis_type),
    n_sierras = n_sites,
    n_plants = n_plants,
    n_insects = n_insects,
    n_cooccurrence_pairs = n_cooc_pairs,
    n_realised_links = n_links,
    connectance = n_links / (n_plants * n_insects),
    proportion_cooccurrence_pairs_realised = n_links / n_cooc_pairs
  )
}))
write.csv2(analysis_unit_summary, file.path(output_dir, "01_analysis_unit_summary.csv"), row.names = FALSE)

## ---------------------------
## 4. Calibration
## ---------------------------

calibration_summary <- bind_rows(lapply(names(analysis_units), function(id){
  u <- analysis_units[[id]]
  observed_links <- n_distinct(u$interactions$pair_id)
  p_cal <- fit_p_calibrated(u$cooc_counts, observed_links)
  expected_links <- expected_links_from_p(u$cooc_counts, p_cal)
  data.frame(
    analysis_id = id,
    analysis_type = unique(u$cooc$analysis_type),
    p_calibrated = p_cal,
    observed_full_links = observed_links,
    expected_full_links = expected_links,
    calibration_error = observed_links - expected_links,
    n_cooccurrence_pairs = nrow(u$cooc_counts),
    mean_n_cooccurring_sierras = mean(u$cooc_counts$n_cooccurring_sierras, na.rm = TRUE),
    median_n_cooccurring_sierras = median(u$cooc_counts$n_cooccurring_sierras, na.rm = TRUE),
    min_n_cooccurring_sierras = min(u$cooc_counts$n_cooccurring_sierras, na.rm = TRUE),
    max_n_cooccurring_sierras = max(u$cooc_counts$n_cooccurring_sierras, na.rm = TRUE)
  )
}))
write.csv2(calibration_summary, file.path(output_dir, "01_calibration_summary.csv"), row.names = FALSE)

## ---------------------------
## 5. Sierra-removal analysis
## ---------------------------

analyse_subset <- function(u, sites_keep, subset_row, p_cal, observed_full_links, expected_full_links){
  cooc_counts_subset <- u$cooc %>%
    filter(site %in% sites_keep) %>%
    group_by(plant, insect, pair_id) %>%
    summarise(n_cooccurring_sierras = n_distinct(site), .groups = "drop")

  observed_retained_links <- u$interactions %>%
    filter(site %in% sites_keep) %>%
    distinct(pair_id) %>%
    nrow()

  expected_retained_links <- expected_links_from_p(cooc_counts_subset, p_cal)

  data.frame(
    analysis_id = subset_row$analysis_id,
    target_removal_fraction = subset_row$target_removal_fraction,
    actual_removal_fraction = subset_row$actual_removal_fraction,
    subset_id = subset_row$subset_id,
    n_removed = subset_row$n_removed,
    n_sites_retained = subset_row$n_sites_retained,
    subset_source = subset_row$subset_source,
    empirical_retained_links = observed_retained_links,
    expected_retained_links = expected_retained_links,
    link_error = observed_retained_links - expected_retained_links,
    relative_link_error_observed_denominator = safe_relative_error_observed_denominator(observed_retained_links, expected_retained_links),
    relative_link_error_expected_denominator = safe_relative_error_expected_denominator(observed_retained_links, expected_retained_links),
    empirical_retained_link_proportion = observed_retained_links / observed_full_links,
    expected_retained_link_proportion = expected_retained_links / expected_full_links
  )
}

subset_objects <- list()
all_removal_results <- list()

for(id in names(analysis_units)){
  u <- analysis_units[[id]]
  cal <- calibration_summary %>% filter(analysis_id == id)
  subsets <- generate_sierra_subsets(u$sites, id)
  subset_objects[[id]] <- subsets

  all_removal_results[[id]] <- bind_rows(lapply(seq_len(nrow(subsets$index)), function(i){
    analyse_subset(u, subsets$subsets[[i]], subsets$index[i, ], cal$p_calibrated,
                   cal$observed_full_links, cal$expected_full_links)
  })) %>%
    mutate(analysis_type = unique(u$cooc$analysis_type))
}

sierra_removal_by_subset <- bind_rows(all_removal_results)
sierra_removal_summary <- sierra_removal_by_subset %>%
  group_by(analysis_id, analysis_type, target_removal_fraction, actual_removal_fraction, n_removed, n_sites_retained) %>%
  summarise(
    mean_empirical_retained_links = mean(empirical_retained_links, na.rm = TRUE),
    sd_empirical_retained_links = sd(empirical_retained_links, na.rm = TRUE),
    mean_expected_retained_links = mean(expected_retained_links, na.rm = TRUE),
    sd_expected_retained_links = sd(expected_retained_links, na.rm = TRUE),
    mean_link_error = mean(link_error, na.rm = TRUE),
    sd_link_error = sd(link_error, na.rm = TRUE),
    mean_relative_link_error_observed_denominator = mean(relative_link_error_observed_denominator, na.rm = TRUE),
    mean_relative_link_error_expected_denominator = mean(relative_link_error_expected_denominator, na.rm = TRUE),
    mean_empirical_retained_link_proportion = mean(empirical_retained_link_proportion, na.rm = TRUE),
    mean_expected_retained_link_proportion = mean(expected_retained_link_proportion, na.rm = TRUE),
    n_subsets = n(),
    .groups = "drop"
  )

write.csv2(sierra_removal_by_subset, file.path(output_dir, "01_sierra_removal_by_subset.csv"), row.names = FALSE)
write.csv2(sierra_removal_summary, file.path(output_dir, "01_sierra_removal_summary.csv"), row.names = FALSE)

## ---------------------------
## 6. Model envelope and repeatability
## ---------------------------

simulate_site_interactions <- function(u, p_cal){
  u$cooc %>%
    mutate(interaction_observed = rbinom(n(), size = 1, prob = p_cal)) %>%
    filter(interaction_observed == 1) %>%
    select(analysis_id, analysis_type, site, plant, insect, pair_id)
}

retained_links_from_interactions <- function(sim_ints, sites_keep){
  sim_ints %>% filter(site %in% sites_keep) %>% distinct(pair_id) %>% nrow()
}

model_envelope_rows <- list()
model_pair_repeatability_all <- list()

for(id in names(analysis_units)){
  message("Model envelope/repeatability: ", id)
  u <- analysis_units[[id]]
  cal <- calibration_summary %>% filter(analysis_id == id)
  subsets <- subset_objects[[id]]

  for(r in seq_len(n_model_reps)){
    sim_ints <- simulate_site_interactions(u, cal$p_calibrated)
    sim_full_links <- n_distinct(sim_ints$pair_id)

    model_envelope_rows[[paste(id, r, sep = "_")]] <- bind_rows(lapply(seq_len(nrow(subsets$index)), function(i){
      sites_keep <- subsets$subsets[[i]]
      subset_row <- subsets$index[i, ]
      retained_links <- retained_links_from_interactions(sim_ints, sites_keep)
      data.frame(
        analysis_id = id,
        analysis_type = unique(u$cooc$analysis_type),
        model_rep = r,
        subset_id = subset_row$subset_id,
        target_removal_fraction = subset_row$target_removal_fraction,
        actual_removal_fraction = subset_row$actual_removal_fraction,
        n_removed = subset_row$n_removed,
        n_sites_retained = subset_row$n_sites_retained,
        model_retained_links = retained_links,
        model_retained_link_proportion = ifelse(sim_full_links > 0, retained_links / sim_full_links, NA_real_)
      )
    }))

    sim_counts <- sim_ints %>%
      group_by(pair_id) %>%
      summarise(n_interacting_sierras = n_distinct(site), .groups = "drop")

    model_pair_repeatability_all[[paste(id, r, sep = "_")]] <- u$cooc_counts %>%
      left_join(sim_counts, by = "pair_id") %>%
      mutate(
        n_interacting_sierras = replace_na(n_interacting_sierras, 0L),
        realised_link = n_interacting_sierras > 0,
        repeatability = n_interacting_sierras / n_cooccurring_sierras,
        analysis_id = id,
        analysis_type = unique(u$cooc$analysis_type),
        model_rep = r
      ) %>%
      filter(realised_link)
  }
}

model_envelope_by_subset <- bind_rows(model_envelope_rows)
model_envelope_summary <- model_envelope_by_subset %>%
  group_by(analysis_id, analysis_type, target_removal_fraction, actual_removal_fraction, n_removed, n_sites_retained) %>%
  summarise(
    model_mean_retained_links = mean(model_retained_links, na.rm = TRUE),
    model_q025_retained_links = quantile(model_retained_links, 0.025, na.rm = TRUE),
    model_q500_retained_links = quantile(model_retained_links, 0.500, na.rm = TRUE),
    model_q975_retained_links = quantile(model_retained_links, 0.975, na.rm = TRUE),
    model_mean_retained_link_proportion = mean(model_retained_link_proportion, na.rm = TRUE),
    model_q025_retained_link_proportion = quantile(model_retained_link_proportion, 0.025, na.rm = TRUE),
    model_q500_retained_link_proportion = quantile(model_retained_link_proportion, 0.500, na.rm = TRUE),
    model_q975_retained_link_proportion = quantile(model_retained_link_proportion, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
write.csv2(model_envelope_summary, file.path(output_dir, "01_model_envelope_summary.csv"), row.names = FALSE)

empirical_pair_repeatability <- bind_rows(lapply(names(analysis_units), function(id){
  u <- analysis_units[[id]]
  u$cooc_counts %>%
    left_join(u$interaction_counts, by = c("plant", "insect", "pair_id")) %>%
    mutate(
      n_interacting_sierras = replace_na(n_interacting_sierras, 0L),
      realised_link = n_interacting_sierras > 0,
      observed_repeatability = n_interacting_sierras / n_cooccurring_sierras,
      analysis_id = id,
      analysis_type = unique(u$cooc$analysis_type)
    ) %>%
    filter(realised_link) %>%
    select(analysis_id, analysis_type, plant, insect, pair_id,
           n_cooccurring_sierras, n_interacting_sierras, observed_repeatability)
}))
write.csv2(empirical_pair_repeatability, file.path(output_dir, "01_empirical_pair_repeatability.csv"), row.names = FALSE)

model_pair_repeatability <- bind_rows(model_pair_repeatability_all)
model_pair_repeatability_summary <- model_pair_repeatability %>%
  group_by(analysis_id, analysis_type, model_rep) %>%
  summarise(
    mean_repeatability = mean(repeatability, na.rm = TRUE),
    median_repeatability = median(repeatability, na.rm = TRUE),
    proportion_one_sierra_links = mean(n_interacting_sierras == 1, na.rm = TRUE),
    mean_interacting_sierras = mean(n_interacting_sierras, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(analysis_id, analysis_type) %>%
  summarise(
    model_mean_repeatability = mean(mean_repeatability, na.rm = TRUE),
    model_q025_mean_repeatability = quantile(mean_repeatability, 0.025, na.rm = TRUE),
    model_q500_mean_repeatability = quantile(mean_repeatability, 0.500, na.rm = TRUE),
    model_q975_mean_repeatability = quantile(mean_repeatability, 0.975, na.rm = TRUE),
    model_mean_proportion_one_sierra_links = mean(proportion_one_sierra_links, na.rm = TRUE),
    model_q025_proportion_one_sierra_links = quantile(proportion_one_sierra_links, 0.025, na.rm = TRUE),
    model_q500_proportion_one_sierra_links = quantile(proportion_one_sierra_links, 0.500, na.rm = TRUE),
    model_q975_proportion_one_sierra_links = quantile(proportion_one_sierra_links, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
write.csv2(model_pair_repeatability_summary, file.path(output_dir, "01_model_pair_repeatability_summary.csv"), row.names = FALSE)

emp_repeat_by_n <- empirical_pair_repeatability %>%
  group_by(analysis_id, analysis_type, n_cooccurring_sierras) %>%
  summarise(
    empirical_n_links = n(),
    empirical_mean_repeatability = mean(observed_repeatability, na.rm = TRUE),
    empirical_median_repeatability = median(observed_repeatability, na.rm = TRUE),
    empirical_proportion_one_sierra_links = mean(n_interacting_sierras == 1, na.rm = TRUE),
    .groups = "drop"
  )

model_repeat_by_n <- model_pair_repeatability %>%
  group_by(analysis_id, analysis_type, model_rep, n_cooccurring_sierras) %>%
  summarise(
    model_mean_repeatability_rep = mean(repeatability, na.rm = TRUE),
    model_proportion_one_sierra_links_rep = mean(n_interacting_sierras == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(analysis_id, analysis_type, n_cooccurring_sierras) %>%
  summarise(
    model_mean_repeatability = mean(model_mean_repeatability_rep, na.rm = TRUE),
    model_q025_repeatability = quantile(model_mean_repeatability_rep, 0.025, na.rm = TRUE),
    model_q500_repeatability = quantile(model_mean_repeatability_rep, 0.500, na.rm = TRUE),
    model_q975_repeatability = quantile(model_mean_repeatability_rep, 0.975, na.rm = TRUE),
    model_mean_proportion_one_sierra_links = mean(model_proportion_one_sierra_links_rep, na.rm = TRUE),
    .groups = "drop"
  )

repeatability_by_ncooccurrence <- emp_repeat_by_n %>%
  full_join(model_repeat_by_n, by = c("analysis_id", "analysis_type", "n_cooccurring_sierras"))
write.csv2(repeatability_by_ncooccurrence, file.path(output_dir, "01_repeatability_by_ncooccurrence.csv"), row.names = FALSE)

## ---------------------------
## 7. Figures
## ---------------------------

pooled_removal <- sierra_removal_summary %>%
  filter(analysis_id == "all_dates_pooled") %>%
  left_join(model_envelope_summary %>% filter(analysis_id == "all_dates_pooled"),
            by = c("analysis_id", "analysis_type", "target_removal_fraction",
                   "actual_removal_fraction", "n_removed", "n_sites_retained"))

p1 <- ggplot(pooled_removal, aes(x = actual_removal_fraction)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(y = mean_link_error), linewidth = 1.1) +
  geom_point(aes(y = mean_link_error), size = 2.5) +
  theme_classic(base_size = 13) +
  xlab("Actual fraction of Sierras removed") +
  ylab("Observed retained links - expected retained links") +
  ggtitle("Sabatino pooled analysis: observed minus expected retained links")
ggsave(file.path(output_dir, "01_pooled_observed_minus_expected_links.png"), p1, width = 7.5, height = 5, dpi = 300)

p2 <- ggplot(pooled_removal, aes(x = actual_removal_fraction)) +
  geom_ribbon(aes(ymin = model_q025_retained_link_proportion,
                  ymax = model_q975_retained_link_proportion), fill = "grey80") +
  geom_line(aes(y = model_q500_retained_link_proportion), linetype = 2, linewidth = 1) +
  geom_line(aes(y = mean_empirical_retained_link_proportion), linewidth = 1.1) +
  geom_point(aes(y = mean_empirical_retained_link_proportion), size = 2.5) +
  theme_classic(base_size = 13) +
  xlab("Actual fraction of Sierras removed") +
  ylab("Retained-link proportion") +
  ggtitle("Sabatino pooled analysis: empirical retained links vs model envelope")
ggsave(file.path(output_dir, "01_pooled_empirical_vs_model_retained_proportion.png"), p2, width = 7.5, height = 5, dpi = 300)

pooled_repeat <- empirical_pair_repeatability %>%
  filter(analysis_id == "all_dates_pooled") %>%
  mutate(source = "empirical") %>%
  select(source, repeatability = observed_repeatability)

pooled_model_repeat <- model_pair_repeatability %>%
  filter(analysis_id == "all_dates_pooled") %>%
  mutate(source = "model") %>%
  select(source, repeatability)

p3 <- ggplot(bind_rows(pooled_repeat, pooled_model_repeat), aes(x = repeatability, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 25) +
  theme_classic(base_size = 13) +
  xlab("Interaction repeatability") +
  ylab("Number of realised links") +
  ggtitle("Sabatino pooled analysis: empirical vs model repeatability")
ggsave(file.path(output_dir, "01_pooled_repeatability_empirical_vs_model.png"), p3, width = 7.5, height = 5, dpi = 300)

pooled_rep_by_n <- repeatability_by_ncooccurrence %>% filter(analysis_id == "all_dates_pooled")
p4 <- ggplot(pooled_rep_by_n, aes(x = n_cooccurring_sierras)) +
  geom_ribbon(aes(ymin = model_q025_repeatability, ymax = model_q975_repeatability), fill = "grey80") +
  geom_line(aes(y = model_q500_repeatability), linetype = 2, linewidth = 1) +
  geom_line(aes(y = empirical_mean_repeatability), linewidth = 1.1) +
  geom_point(aes(y = empirical_mean_repeatability, size = empirical_n_links)) +
  theme_classic(base_size = 13) +
  xlab("Number of co-occurring Sierras") +
  ylab("Mean repeatability") +
  ggtitle("Sabatino pooled analysis: repeatability controlled by co-occurrence frequency")
ggsave(file.path(output_dir, "01_pooled_repeatability_by_ncooccurrence.png"), p4, width = 7.5, height = 5, dpi = 300)

date_removal <- sierra_removal_summary %>% filter(analysis_type == "date_stratified")
if(nrow(date_removal) > 0){
  p5 <- ggplot(date_removal, aes(x = actual_removal_fraction, y = mean_link_error)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    facet_wrap(~ analysis_id, scales = "free_y") +
    theme_classic(base_size = 11) +
    xlab("Actual fraction of Sierras removed") +
    ylab("Observed - expected retained links") +
    ggtitle("Date-stratified Sabatino analyses: observed minus expected links")
  ggsave(file.path(output_dir, "01_date_stratified_observed_minus_expected_links.png"), p5, width = 12, height = 8, dpi = 300)

  date_repeat_summary <- empirical_pair_repeatability %>%
    filter(analysis_type == "date_stratified") %>%
    group_by(analysis_id) %>%
    summarise(empirical_mean_repeatability = mean(observed_repeatability, na.rm = TRUE), .groups = "drop") %>%
    left_join(model_pair_repeatability_summary %>% filter(analysis_type == "date_stratified"), by = "analysis_id")

  p6 <- ggplot(date_repeat_summary, aes(x = empirical_mean_repeatability, y = model_mean_repeatability, label = analysis_id)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(size = 2) +
    geom_text(check_overlap = TRUE, vjust = -0.5, size = 3) +
    theme_classic(base_size = 12) +
    xlab("Empirical mean repeatability") +
    ylab("Model mean repeatability") +
    ggtitle("Date-stratified repeatability summary")
  ggsave(file.path(output_dir, "01_date_stratified_repeatability_summary.png"), p6, width = 7.5, height = 5.5, dpi = 300)
}

coverage_plot_data <- sierra_date_coverage %>%
  mutate(Date_label = factor(Date_label, levels = sort(unique(Date_label))),
         Sierra = factor(Sierra, levels = sort(unique(Sierra))))

p_cov <- ggplot(coverage_plot_data, aes(x = Date_label, y = Sierra, fill = valid_interaction_record_count)) +
  geom_tile() +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Date") +
  ylab("Sierra") +
  ggtitle("Sabatino data coverage overview", subtitle = "Fill = valid interaction-record count per Sierra-Date")
ggsave(file.path(output_dir, "01_data_coverage_overview.png"), p_cov, width = 10, height = 6, dpi = 300)

## ---------------------------
## 8. Interpretation notes
## ---------------------------

notes <- c(
  "Sabatino first-pass Sierra-level Galiana-style analysis",
  "",
  "This is an operational replication of the Galiana-style local co-occurrence framework.",
  "The spatial unit used here is Sierra.",
  "",
  "Important interpretation constraints:",
  "- Insects were detected through observed flower visits, so this analysis does not estimate true interaction probability conditional on independently surveyed insect occurrence.",
  "- Operational co-occurrence means: plant recorded in Sierra + valid insect recorded in Sierra.",
  "- SV, blank, missing, and clearly non-taxon insect labels are excluded from the primary binary interaction-network analysis, but retained in the cleaned raw-data audit.",
  "- Repeated rows within Sierra-Date are collapsed to binary plant-insect links for the primary analysis because their meaning is not yet confirmed.",
  "- Interaction abundance and flower abundance are retained in descriptive outputs but are not used in the primary homogeneous co-occurrence model.",
  "- Date-stratified analyses are exploratory and are only run for dates with adequate Sierra replication and enough empirical links.",
  "- Random Sierra removal is a structural/sampling diagnostic, not a literal habitat-loss experiment.",
  "",
  "Central interpretation:",
  "At zero Sierra removal, observed and expected full metaweb link numbers match by calibration.",
  "Positive observed-minus-expected divergence after Sierra removal means empirical links are retained more strongly than expected under one homogeneous co-occurrence-to-interaction probability."
)

writeLines(notes, con = file.path(output_dir, "01_interpretation_notes.txt"))

message("Finished Sabatino Sierra-level basic co-occurrence analysis.")
message("Outputs saved in: ", output_dir)
