
## ------------------------------------------------------------
## Script: Sabatino/scripts/02_census_round_galiana_style_analysis.R
##
## Census-round Galiana-style analysis for Sabatino.
## Site = Sierra. Census round = chronological sampling order within Sierra.
## Uses write.csv2().
## ------------------------------------------------------------

packages <- c("readxl","dplyr","tidyr","ggplot2","lubridate","stringr","tibble","future","future.apply")
for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
set.seed(123)

parallel_workers <- suppressWarnings(as.integer(Sys.getenv("R_PARALLEL_WORKERS", unset=NA_character_)))
if(is.na(parallel_workers) || parallel_workers < 1){
  detected_cores <- parallel::detectCores(logical=FALSE)
  if(is.na(detected_cores) || detected_cores < 1) detected_cores <- 1
  parallel_workers <- max(1, detected_cores - 1)
}
future::plan(future::multisession, workers=parallel_workers)
on.exit(future::plan(future::sequential), add=TRUE)
message("Parallel workers: ", parallel_workers)

input_file <- "Sabatino/Network 2008_updated.xlsx"
input_sheet <- "Matrix complete_2008_marzo2024"

out_dir <- "Sabatino/outputs/02_census_round_galiana_style_analysis"
dir.create("Sabatino/scripts", recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_sierras_per_round <- 8
min_empirical_links_per_round <- 20
target_removal_fractions <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
max_exact_subsets <- 5000
n_random_subsets <- 1000
min_sites_retained <- 2
n_model_reps <- 50

trim_to_na <- function(x){x <- stringr::str_trim(as.character(x)); x[x==""] <- NA_character_; x}
safe_date <- function(x){out <- suppressWarnings(lubridate::ymd(x)); if(all(is.na(out))) out <- suppressWarnings(as.Date(x)); out}
safe_num <- function(x) suppressWarnings(as.numeric(x))
pair_id <- function(p,i) paste(p, i, sep="___")
valid_insect <- function(x){
  z <- tolower(stringr::str_trim(as.character(x)))
  !(is.na(z) | z=="" | z %in% c("sv","s.v.","s v","na","n/a","none","unknown","unidentified","no visitor","no visitors","sin visita"))
}
rel_exp <- function(obs, exp){ifelse(is.na(exp) | exp==0, NA_real_, (obs-exp)/exp)}
rel_obs <- function(obs, exp){ifelse(is.na(obs) | obs==0, NA_real_, (obs-exp)/obs)}
gini <- function(x){x <- x[is.finite(x)]; if(length(x)==0) return(NA_real_); if(all(x==0)) return(0); x <- sort(x); n <- length(x); sum((2*seq_along(x)-n-1)*x)/(n*sum(x))}
exp_links <- function(cooc_counts, p){if(nrow(cooc_counts)==0 || is.na(p)) return(NA_real_); sum(1-(1-p)^cooc_counts$n_cooccurring_sierras, na.rm=TRUE)}
fit_p <- function(cooc_counts, L){
  if(nrow(cooc_counts)==0 || is.na(L)) return(NA_real_)
  if(L <= 0) return(0)
  if(L >= nrow(cooc_counts)) return(1)
  f <- function(p) exp_links(cooc_counts, p) - L
  tryCatch(uniroot(f, c(0,1), tol=1e-10)$root, error=function(e) NA_real_)
}
auc <- function(x,y){ok <- is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]; if(length(x)<2) return(NA_real_); o<-order(x); x<-x[o]; y<-y[o]; sum(diff(x)*(head(y,-1)+tail(y,-1))/2)}

make_subsets <- function(sites, round_id){
  sites <- sort(unique(as.character(sites)))
  n <- length(sites)
  feasible <- data.frame(target_removal_fraction=target_removal_fractions) |>
    mutate(n_removed = pmax(0, pmin(round(n*target_removal_fraction), n-min_sites_retained)),
           actual_removal_fraction = n_removed/n,
           n_sites_retained = n-n_removed) |>
    distinct(n_removed, .keep_all=TRUE)
  sets <- list(); idx <- list(); k <- 1
  for(i in seq_len(nrow(feasible))){
    nkeep <- feasible$n_sites_retained[i]
    ncomb <- choose(n, nkeep)
    exact <- is.finite(ncomb) && ncomb <= max_exact_subsets
    kept <- if(exact) combn(sites, nkeep, simplify=FALSE) else unique(replicate(n_random_subsets, sort(sample(sites,nkeep)), simplify=FALSE))
    source <- if(exact) "exact" else "random"
    for(j in seq_along(kept)){
      sets[[k]] <- kept[[j]]
      idx[[k]] <- data.frame(census_round=round_id, subset_id=k,
                             target_removal_fraction=feasible$target_removal_fraction[i],
                             actual_removal_fraction=feasible$actual_removal_fraction[i],
                             n_removed=n-length(kept[[j]]), n_sites_retained=length(kept[[j]]),
                             subset_source=source)
      k <- k+1
    }
  }
  list(subsets=sets, index=bind_rows(idx))
}

build_unit <- function(rec, r){
  rec <- rec |> filter(census_round==r)
  plant_site <- rec |> distinct(site=Sierra, plant=Plant_species)
  insect_site <- rec |> distinct(site=Sierra, insect=Insect_species)
  cooc <- inner_join(plant_site, insect_site, by="site") |>
    transmute(census_round=r, site, plant, insect, pair_id=pair_id(plant,insect)) |> distinct()
  ints <- rec |> transmute(census_round=r, site=Sierra, plant=Plant_species, insect=Insect_species, pair_id=pair_id(plant,insect)) |> distinct()
  cooc_counts <- cooc |> group_by(plant,insect,pair_id) |> summarise(n_cooccurring_sierras=n_distinct(site), .groups="drop")
  int_counts <- ints |> group_by(plant,insect,pair_id) |> summarise(n_interacting_sierras=n_distinct(site), .groups="drop")
  list(cooc=cooc, interactions=ints, cooc_counts=cooc_counts, interaction_counts=int_counts, sites=sort(unique(cooc$site)))
}
sim_ints <- function(u,p) u$cooc |> mutate(hit=rbinom(n(),1,p)) |> filter(hit==1) |> select(census_round,site,plant,insect,pair_id)
deg_metrics <- function(pairs){
  if(nrow(pairs)==0) return(data.frame(trophic_side=c("plant","insect"), mean_degree=NA_real_, maximum_degree=NA_real_, degree_gini=NA_real_, proportion_degree_1=NA_real_))
  bind_rows(
    pairs |> count(species=plant, name="degree") |> mutate(trophic_side="plant"),
    pairs |> count(species=insect, name="degree") |> mutate(trophic_side="insect")
  ) |> group_by(trophic_side) |> summarise(mean_degree=mean(degree), maximum_degree=max(degree), degree_gini=gini(degree), proportion_degree_1=mean(degree==1), .groups="drop")
}
deg_freq <- function(pairs, round_id, source, rep=NA_integer_){
  if(nrow(pairs)==0) return(data.frame())
  bind_rows(
    pairs |> count(species=plant, name="degree") |> mutate(trophic_side="plant"),
    pairs |> count(species=insect, name="degree") |> mutate(trophic_side="insect")
  ) |> count(trophic_side, degree, name="frequency") |> mutate(census_round=round_id, source=source, model_rep=rep)
}

## Read and clean
raw <- readxl::read_excel(input_file, sheet=input_sheet) |> select(!starts_with("..."))
dat <- raw |>
  mutate(row_id=row_number(),
         Sierra=trim_to_na(Sierra),
         Date_chr=trim_to_na(Date),
         Date_parsed=safe_date(Date_chr),
         Date_label=ifelse(is.na(Date_parsed), Date_chr, as.character(Date_parsed)),
         Plant_species=trim_to_na(Plant_species),
         Insect_species=trim_to_na(Insect_species),
         Plant_ID_chr=trim_to_na(Plant_ID),
         Number_of_flowers_numeric=safe_num(Number_of_flowers),
         Insect_abundance_numeric=safe_num(Insect_abundance),
         valid_insect_label=valid_insect(Insect_species),
         valid_interaction_record=!is.na(Sierra)&!is.na(Date_label)&!is.na(Plant_species)&!is.na(Insect_species)&valid_insect_label)

date_lookup <- dat |> filter(!is.na(Sierra), !is.na(Date_label)) |> distinct(Sierra, Date_label, Date_parsed) |>
  arrange(Sierra, Date_parsed, Date_label) |> group_by(Sierra) |> mutate(census_round=row_number()) |> ungroup()
dat <- dat |> left_join(date_lookup |> select(Sierra,Date_label,census_round), by=c("Sierra","Date_label"))
write.csv2(dat, file.path(out_dir,"02_cleaned_raw_records_with_census_rounds.csv"), row.names=FALSE)

valid <- dat |> filter(valid_interaction_record, !is.na(census_round))
all_sierras <- sort(unique(dat$Sierra[!is.na(dat$Sierra)]))

mapping <- dat |> group_by(Sierra, Date_label, census_round) |>
  summarise(n_raw_records=n(),
            n_valid_interaction_records=sum(valid_interaction_record, na.rm=TRUE),
            n_plants=n_distinct(Plant_species, na.rm=TRUE),
            n_insects=n_distinct(Insect_species[valid_insect_label], na.rm=TRUE),
            n_observed_links=n_distinct(pair_id(Plant_species[valid_interaction_record], Insect_species[valid_interaction_record])),
            n_plant_ids=n_distinct(Plant_ID_chr, na.rm=TRUE),
            .groups="drop") |> arrange(census_round,Sierra)
write.csv2(mapping, file.path(out_dir,"02_census_round_date_mapping.csv"), row.names=FALSE)

coverage <- mapping |> group_by(census_round) |>
  summarise(n_sierras_represented=n_distinct(Sierra),
            represented_sierras=paste(sort(unique(Sierra)), collapse="; "),
            missing_sierras=paste(setdiff(all_sierras, sort(unique(Sierra))), collapse="; "),
            total_valid_interaction_records=sum(n_valid_interaction_records, na.rm=TRUE),
            total_plants=n_distinct(valid$Plant_species[valid$census_round==first(census_round)]),
            total_insects=n_distinct(valid$Insect_species[valid$census_round==first(census_round)]),
            total_observed_links=n_distinct(pair_id(valid$Plant_species[valid$census_round==first(census_round)], valid$Insect_species[valid$census_round==first(census_round)])),
            .groups="drop") |>
  mutate(missing_sierras=ifelse(missing_sierras=="", NA_character_, missing_sierras))
write.csv2(coverage, file.path(out_dir,"02_census_round_coverage.csv"), row.names=FALSE)

skipped <- coverage |> mutate(skipped_reason=paste(ifelse(n_sierras_represented<min_sierras_per_round,"too_few_sierras",NA),
                                                   ifelse(total_observed_links<min_empirical_links_per_round,"too_few_empirical_links",NA),
                                                   sep=";"),
                              skipped_reason=gsub("NA;|;NA|NA","",skipped_reason)) |>
  filter(skipped_reason!="") |> select(census_round,n_sierras_represented,total_observed_links,skipped_reason)
write.csv2(skipped, file.path(out_dir,"02_skipped_census_rounds.csv"), row.names=FALSE)

eligible <- coverage |> filter(n_sierras_represented>=min_sierras_per_round, total_observed_links>=min_empirical_links_per_round) |> pull(census_round)
message("Eligible rounds: ", paste(eligible, collapse=", "))

collapsed <- valid |> group_by(census_round,Sierra,Plant_species,Insect_species) |>
  summarise(interaction_observed=1L, n_raw_records=n(), n_distinct_plant_ids=n_distinct(Plant_ID_chr, na.rm=TRUE),
            interaction_abundance_sum=sum(Insect_abundance_numeric, na.rm=TRUE),
            interaction_abundance_mean=mean(Insect_abundance_numeric, na.rm=TRUE), .groups="drop") |>
  mutate(pair_id=pair_id(Plant_species, Insect_species))
write.csv2(collapsed, file.path(out_dir,"02_cleaned_census_round_interaction_records.csv"), row.names=FALSE)

units <- future.apply::future_lapply(eligible, function(r) build_unit(collapsed,r), future.seed=TRUE); names(units) <- as.character(eligible)

cal <- bind_rows(future.apply::future_lapply(seq_along(units), function(i){
  r <- eligible[i]; u <- units[[i]]
  L <- n_distinct(u$interactions$pair_id)
  p <- fit_p(u$cooc_counts, L)
  E <- exp_links(u$cooc_counts, p)
  np <- n_distinct(u$cooc$plant); ni <- n_distinct(u$cooc$insect); nc <- nrow(u$cooc_counts)
  data.frame(census_round=r, n_sierras=length(u$sites), n_plants=np, n_insects=ni, n_cooccurrence_pairs=nc,
             observed_full_links=L, expected_full_links=E, calibration_error=L-E, p_calibrated=p,
             mean_n_cooccurring_sierras=mean(u$cooc_counts$n_cooccurring_sierras),
             median_n_cooccurring_sierras=median(u$cooc_counts$n_cooccurring_sierras),
             max_n_cooccurring_sierras=max(u$cooc_counts$n_cooccurring_sierras),
             connectance=L/(np*ni), conversion_rate=L/nc)
}, future.seed=TRUE))
write.csv2(cal, file.path(out_dir,"02_calibration_summary.csv"), row.names=FALSE)

## Site removal
subset_objects <- lapply(seq_along(units), function(i){
  r <- eligible[i]; u <- units[[i]]
  make_subsets(u$sites, r)
})
names(subset_objects) <- as.character(eligible)

rem_list <- future.apply::future_lapply(seq_along(units), function(i){
  r <- eligible[i]; u <- units[[i]]; c0 <- cal |> filter(census_round==r)
  subsets <- subset_objects[[as.character(r)]]
  bind_rows(lapply(seq_len(nrow(subsets$index)), function(j){
    keep <- subsets$subsets[[j]]; sr <- subsets$index[j,]
    cc <- u$cooc |> filter(site %in% keep) |> group_by(plant,insect,pair_id) |> summarise(n_cooccurring_sierras=n_distinct(site), .groups="drop")
    obs <- u$interactions |> filter(site %in% keep) |> distinct(pair_id) |> nrow()
    exp <- exp_links(cc, c0$p_calibrated)
    data.frame(census_round=r, subset_id=sr$subset_id, target_removal_fraction=sr$target_removal_fraction,
               actual_removal_fraction=sr$actual_removal_fraction, n_removed=sr$n_removed,
               n_sites_retained=sr$n_sites_retained, subset_source=sr$subset_source,
               empirical_retained_links=obs, expected_retained_links=exp, link_error=obs-exp,
               relative_divergence_expected_denominator=rel_exp(obs,exp),
               relative_divergence_observed_denominator=rel_obs(obs,exp),
               empirical_retained_link_proportion=obs/c0$observed_full_links,
               expected_retained_link_proportion=exp/c0$expected_full_links)
  }))
}, future.seed=TRUE)
names(rem_list) <- as.character(eligible)
rem_by_subset <- bind_rows(rem_list)
rem_summary <- rem_by_subset |> group_by(census_round,target_removal_fraction,actual_removal_fraction,n_removed,n_sites_retained) |>
  summarise(mean_empirical_retained_links=mean(empirical_retained_links), sd_empirical_retained_links=sd(empirical_retained_links),
            mean_expected_retained_links=mean(expected_retained_links), sd_expected_retained_links=sd(expected_retained_links),
            mean_link_error=mean(link_error), sd_link_error=sd(link_error),
            mean_relative_divergence_expected_denominator=mean(relative_divergence_expected_denominator, na.rm=TRUE),
            mean_relative_divergence_observed_denominator=mean(relative_divergence_observed_denominator, na.rm=TRUE),
            mean_empirical_retained_link_proportion=mean(empirical_retained_link_proportion),
            mean_expected_retained_link_proportion=mean(expected_retained_link_proportion), n_subsets=n(), .groups="drop")
write.csv2(rem_by_subset, file.path(out_dir,"02_sierra_removal_by_subset.csv"), row.names=FALSE)
write.csv2(rem_summary, file.path(out_dir,"02_sierra_removal_summary.csv"), row.names=FALSE)

## Model simulations, envelope, repeatability, degree
sim_results <- future.apply::future_lapply(seq_along(units), function(i){
  r <- eligible[i]; u <- units[[i]]; c0 <- cal |> filter(census_round==r); subsets <- subset_objects[[as.character(r)]]
  message("Simulating census round ", r)
  emp_pairs <- u$interactions |> distinct(plant,insect,pair_id)
  deg_emp_i <- deg_metrics(emp_pairs) |> mutate(census_round=r, source="empirical", model_rep=NA_integer_)
  freq_emp_i <- deg_freq(emp_pairs, r, "empirical")

  reps <- lapply(seq_len(n_model_reps), function(m){
    sim <- sim_ints(u, c0$p_calibrated)
    sim_full <- n_distinct(sim$pair_id)
    env_i <- bind_rows(lapply(seq_len(nrow(subsets$index)), function(j){
      sr <- subsets$index[j,]; retained <- sim |> filter(site %in% subsets$subsets[[j]]) |> distinct(pair_id) |> nrow()
      data.frame(census_round=r, model_rep=m, subset_id=sr$subset_id, target_removal_fraction=sr$target_removal_fraction,
                 actual_removal_fraction=sr$actual_removal_fraction, n_removed=sr$n_removed, n_sites_retained=sr$n_sites_retained,
                 model_retained_links=retained, model_retained_link_proportion=ifelse(sim_full>0, retained/sim_full, NA_real_))
    }))
    sim_counts <- sim |> group_by(pair_id) |> summarise(n_interacting_sierras=n_distinct(site), .groups="drop")
    mod_rep_i <- u$cooc_counts |> left_join(sim_counts, by="pair_id") |>
      mutate(n_interacting_sierras=replace_na(n_interacting_sierras,0L), realised_link=n_interacting_sierras>0,
             repeatability=n_interacting_sierras/n_cooccurring_sierras, census_round=r, model_rep=m) |> filter(realised_link)
    sim_pairs <- sim |> distinct(plant,insect,pair_id)
    deg_mod_i <- deg_metrics(sim_pairs) |> mutate(census_round=r, source="model", model_rep=m)
    freq_mod_i <- deg_freq(sim_pairs, r, "model", m)
    list(env=env_i, mod_rep=mod_rep_i, deg_mod=deg_mod_i, freq_mod=freq_mod_i)
  })

  list(deg_emp=deg_emp_i, freq_emp=freq_emp_i,
       env=bind_rows(lapply(reps, `[[`, "env")),
       mod_rep=bind_rows(lapply(reps, `[[`, "mod_rep")),
       deg_mod=bind_rows(lapply(reps, `[[`, "deg_mod")),
       freq_mod=bind_rows(lapply(reps, `[[`, "freq_mod")))
}, future.seed=TRUE)

env <- bind_rows(lapply(sim_results, `[[`, "env"))
mod_rep_rows <- lapply(sim_results, `[[`, "mod_rep")
deg_mod <- lapply(sim_results, `[[`, "deg_mod")
freq_mod <- lapply(sim_results, `[[`, "freq_mod")
deg_emp <- lapply(sim_results, `[[`, "deg_emp")
freq_emp <- lapply(sim_results, `[[`, "freq_emp")
env_summary <- env |> group_by(census_round,target_removal_fraction,actual_removal_fraction,n_removed,n_sites_retained) |>
  summarise(model_mean_retained_links=mean(model_retained_links, na.rm=TRUE),
            model_q025_retained_links=quantile(model_retained_links,0.025,na.rm=TRUE),
            model_q500_retained_links=quantile(model_retained_links,0.5,na.rm=TRUE),
            model_q975_retained_links=quantile(model_retained_links,0.975,na.rm=TRUE),
            model_mean_retained_link_proportion=mean(model_retained_link_proportion,na.rm=TRUE),
            model_q025_retained_link_proportion=quantile(model_retained_link_proportion,0.025,na.rm=TRUE),
            model_q500_retained_link_proportion=quantile(model_retained_link_proportion,0.5,na.rm=TRUE),
            model_q975_retained_link_proportion=quantile(model_retained_link_proportion,0.975,na.rm=TRUE), .groups="drop")
write.csv2(env_summary, file.path(out_dir,"02_model_envelope_summary.csv"), row.names=FALSE)

emp_rep <- bind_rows(lapply(seq_along(units), function(i){
  r <- eligible[i]; u <- units[[i]]
  u$cooc_counts |> left_join(u$interaction_counts, by=c("plant","insect","pair_id")) |>
    mutate(n_interacting_sierras=replace_na(n_interacting_sierras,0L), realised_link=n_interacting_sierras>0,
           observed_repeatability=n_interacting_sierras/n_cooccurring_sierras, census_round=r) |> filter(realised_link)
}))
write.csv2(emp_rep, file.path(out_dir,"02_empirical_pair_repeatability.csv"), row.names=FALSE)

mod_rep <- bind_rows(mod_rep_rows)
mod_rep_summary <- mod_rep |> group_by(census_round,model_rep) |>
  summarise(model_mean_repeatability_rep=mean(repeatability,na.rm=TRUE), model_median_repeatability_rep=median(repeatability,na.rm=TRUE),
            model_proportion_one_sierra_links_rep=mean(n_interacting_sierras==1,na.rm=TRUE),
            model_mean_interacting_sierras_rep=mean(n_interacting_sierras,na.rm=TRUE), .groups="drop") |>
  group_by(census_round) |> summarise(model_mean_repeatability=mean(model_mean_repeatability_rep,na.rm=TRUE),
                                      model_q025_mean_repeatability=quantile(model_mean_repeatability_rep,0.025,na.rm=TRUE),
                                      model_q500_mean_repeatability=quantile(model_mean_repeatability_rep,0.5,na.rm=TRUE),
                                      model_q975_mean_repeatability=quantile(model_mean_repeatability_rep,0.975,na.rm=TRUE),
                                      model_median_repeatability=mean(model_median_repeatability_rep,na.rm=TRUE),
                                      model_proportion_one_sierra_links=mean(model_proportion_one_sierra_links_rep,na.rm=TRUE),
                                      model_mean_interacting_sierras=mean(model_mean_interacting_sierras_rep,na.rm=TRUE), .groups="drop")
write.csv2(mod_rep_summary, file.path(out_dir,"02_model_pair_repeatability_summary.csv"), row.names=FALSE)

emp_by_n <- emp_rep |> group_by(census_round,n_cooccurring_sierras) |>
  summarise(empirical_n_links=n(), empirical_mean_repeatability=mean(observed_repeatability),
            empirical_median_repeatability=median(observed_repeatability),
            empirical_proportion_one_sierra_links=mean(n_interacting_sierras==1), .groups="drop")
mod_by_n <- mod_rep |> group_by(census_round,model_rep,n_cooccurring_sierras) |>
  summarise(model_mean_repeatability_rep=mean(repeatability), model_proportion_one_sierra_links_rep=mean(n_interacting_sierras==1), .groups="drop") |>
  group_by(census_round,n_cooccurring_sierras) |> summarise(model_mean_repeatability=mean(model_mean_repeatability_rep,na.rm=TRUE),
                                                           model_q025_repeatability=quantile(model_mean_repeatability_rep,0.025,na.rm=TRUE),
                                                           model_q500_repeatability=quantile(model_mean_repeatability_rep,0.5,na.rm=TRUE),
                                                           model_q975_repeatability=quantile(model_mean_repeatability_rep,0.975,na.rm=TRUE),
                                                           model_mean_proportion_one_sierra_links=mean(model_proportion_one_sierra_links_rep,na.rm=TRUE), .groups="drop")
rep_by_n <- full_join(emp_by_n, mod_by_n, by=c("census_round","n_cooccurring_sierras"))
write.csv2(rep_by_n, file.path(out_dir,"02_repeatability_by_ncooccurrence.csv"), row.names=FALSE)

deg_emp_df <- bind_rows(deg_emp)
deg_mod_df <- bind_rows(deg_mod)
deg_mod_summary <- deg_mod_df |> group_by(census_round,trophic_side) |>
  summarise(across(c(mean_degree, maximum_degree, degree_gini, proportion_degree_1),
                   list(model_mean=~mean(.x,na.rm=TRUE), model_q025=~quantile(.x,0.025,na.rm=TRUE),
                        model_q500=~quantile(.x,0.5,na.rm=TRUE), model_q975=~quantile(.x,0.975,na.rm=TRUE)),
                   .names="{.col}_{.fn}"), .groups="drop")
deg_out <- deg_emp_df |> rename(mean_degree_empirical=mean_degree, maximum_degree_empirical=maximum_degree,
                                degree_gini_empirical=degree_gini, proportion_degree_1_empirical=proportion_degree_1) |>
  select(census_round,trophic_side,mean_degree_empirical,maximum_degree_empirical,degree_gini_empirical,proportion_degree_1_empirical) |>
  left_join(deg_mod_summary, by=c("census_round","trophic_side"))
write.csv2(deg_out, file.path(out_dir,"02_zero_removal_degree_metrics.csv"), row.names=FALSE)

freq_emp_df <- bind_rows(freq_emp)
freq_mod_df <- bind_rows(freq_mod)
freq_mod_summary <- freq_mod_df |> group_by(census_round,trophic_side,degree,model_rep) |> summarise(frequency=sum(frequency),.groups="drop") |>
  group_by(census_round,trophic_side,degree) |> summarise(model_frequency_mean=mean(frequency), model_frequency_q025=quantile(frequency,0.025),
                                                          model_frequency_q500=quantile(frequency,0.5), model_frequency_q975=quantile(frequency,0.975), .groups="drop")
freq_emp_summary <- freq_emp_df |> group_by(census_round,trophic_side,degree) |> summarise(empirical_frequency=sum(frequency),.groups="drop")
freq_out <- full_join(freq_emp_summary, freq_mod_summary, by=c("census_round","trophic_side","degree")) |>
  mutate(across(c(empirical_frequency,model_frequency_mean,model_frequency_q025,model_frequency_q500,model_frequency_q975), ~replace_na(.x,0)))
write.csv2(freq_out, file.path(out_dir,"02_zero_removal_degree_frequency.csv"), row.names=FALSE)

## Summary metrics
emp_rep_summary <- emp_rep |> group_by(census_round) |> summarise(empirical_mean_repeatability=mean(observed_repeatability),
                                                                  empirical_median_repeatability=median(observed_repeatability),
                                                                  empirical_proportion_one_sierra_links=mean(n_interacting_sierras==1),
                                                                  empirical_mean_interacting_sierras=mean(n_interacting_sierras), .groups="drop")
div_overall <- rem_summary |> left_join(env_summary, by=c("census_round","target_removal_fraction","actual_removal_fraction","n_removed","n_sites_retained")) |>
  group_by(census_round) |> summarise(maximum_absolute_divergence=max(abs(mean_link_error),na.rm=TRUE),
                                     divergence_at_closest_0_4=mean_link_error[which.min(abs(actual_removal_fraction-0.4))],
                                     divergence_at_closest_0_6=mean_link_error[which.min(abs(actual_removal_fraction-0.6))],
                                     divergence_at_closest_0_8=mean_link_error[which.min(abs(actual_removal_fraction-0.8))],
                                     divergence_auc=auc(actual_removal_fraction, mean_link_error),
                                     divergence_slope=ifelse(n_distinct(actual_removal_fraction)>=2, coef(lm(mean_link_error~actual_removal_fraction))[2], NA_real_),
                                     fraction_levels_above_model_975=mean(mean_empirical_retained_link_proportion > model_q975_retained_link_proportion, na.rm=TRUE),
                                     largest_exceedance_above_model_envelope=max(mean_empirical_retained_link_proportion-model_q975_retained_link_proportion, na.rm=TRUE), .groups="drop")
summary_metrics <- cal |> left_join(emp_rep_summary, by="census_round") |> left_join(mod_rep_summary, by="census_round") |>
  mutate(repeatability_excess=empirical_mean_repeatability-model_mean_repeatability) |> left_join(div_overall, by="census_round")
write.csv2(summary_metrics, file.path(out_dir,"02_census_round_summary_metrics.csv"), row.names=FALSE)

## Figures
p_base <- summary_metrics |> select(census_round,p_calibrated,observed_full_links,connectance,empirical_mean_repeatability,repeatability_excess) |>
  pivot_longer(-census_round, names_to="metric", values_to="value") |>
  ggplot(aes(census_round,value)) + geom_line() + geom_point() + facet_wrap(~metric, scales="free_y") +
  theme_classic() + xlab("Census round") + ylab("Value") + ggtitle("Census-round baseline metrics")
ggsave(file.path(out_dir,"02_census_round_p_and_baseline_metrics.png"), p_base, width=10, height=7, dpi=300)

p_div <- ggplot(rem_summary, aes(actual_removal_fraction, mean_link_error)) + geom_hline(yintercept=0, linetype=2) +
  geom_line() + geom_point() + facet_wrap(~census_round, scales="free_y") + theme_classic() +
  xlab("Actual fraction of Sierras removed") + ylab("Observed - expected retained links") +
  ggtitle("Observed-minus-expected divergence across census rounds")
ggsave(file.path(out_dir,"02_census_round_divergence_curves.png"), p_div, width=12, height=8, dpi=300)

ret_plot <- rem_summary |> left_join(env_summary, by=c("census_round","target_removal_fraction","actual_removal_fraction","n_removed","n_sites_retained"))
p_ret <- ggplot(ret_plot, aes(actual_removal_fraction)) + geom_ribbon(aes(ymin=model_q025_retained_link_proportion,ymax=model_q975_retained_link_proportion), fill="grey80") +
  geom_line(aes(y=model_q500_retained_link_proportion), linetype=2) + geom_line(aes(y=mean_empirical_retained_link_proportion)) +
  geom_point(aes(y=mean_empirical_retained_link_proportion)) + facet_wrap(~census_round) + theme_classic() +
  xlab("Actual fraction of Sierras removed") + ylab("Retained-link proportion") +
  ggtitle("Empirical vs homogeneous-model retained-link proportion")
ggsave(file.path(out_dir,"02_census_round_empirical_vs_model_retention.png"), p_ret, width=12, height=8, dpi=300)

p_rep <- ggplot(rep_by_n, aes(n_cooccurring_sierras)) + geom_ribbon(aes(ymin=model_q025_repeatability,ymax=model_q975_repeatability), fill="grey80") +
  geom_line(aes(y=model_q500_repeatability), linetype=2) + geom_line(aes(y=empirical_mean_repeatability)) +
  geom_point(aes(y=empirical_mean_repeatability, size=empirical_n_links)) + facet_wrap(~census_round) +
  coord_cartesian(ylim=c(0,1)) + theme_classic() + xlab("Number of co-occurring Sierras") + ylab("Mean repeatability") +
  ggtitle("Repeatability controlled by co-occurrence frequency")
ggsave(file.path(out_dir,"02_census_round_repeatability_by_ncooccurrence.png"), p_rep, width=12, height=8, dpi=300)

p_pred <- ggplot(summary_metrics, aes(p_calibrated, divergence_auc, label=census_round)) + geom_point() + geom_smooth(method="lm", se=FALSE, linetype=2) +
  geom_text(vjust=-0.6, check_overlap=TRUE) + theme_classic() + xlab("Calibrated p") + ylab("Divergence AUC") +
  ggtitle("Cross-round divergence predictors")
ggsave(file.path(out_dir,"02_cross_round_divergence_predictors.png"), p_pred, width=7, height=5, dpi=300)

p_freq <- ggplot(freq_out, aes(degree)) + geom_ribbon(aes(ymin=model_frequency_q025,ymax=model_frequency_q975), fill="grey80") +
  geom_line(aes(y=model_frequency_q500), linetype=2) + geom_line(aes(y=empirical_frequency)) + geom_point(aes(y=empirical_frequency), size=1) +
  facet_grid(trophic_side~census_round, scales="free_y") + theme_classic() + xlab("Degree") + ylab("Frequency") +
  ggtitle("Removal-0 degree distributions")
ggsave(file.path(out_dir,"02_zero_removal_degree_distributions_empirical_vs_model.png"), p_freq, width=14, height=7, dpi=300)

deg_long <- deg_out |> select(census_round,trophic_side,
                              mean_degree_empirical, maximum_degree_empirical, degree_gini_empirical, proportion_degree_1_empirical,
                              mean_degree_model_q025, mean_degree_model_q500, mean_degree_model_q975,
                              maximum_degree_model_q025, maximum_degree_model_q500, maximum_degree_model_q975,
                              degree_gini_model_q025, degree_gini_model_q500, degree_gini_model_q975,
                              proportion_degree_1_model_q025, proportion_degree_1_model_q500, proportion_degree_1_model_q975) |>
  pivot_longer(-c(census_round,trophic_side), names_to="raw", values_to="value") |>
  mutate(metric=case_when(grepl("mean_degree",raw)~"mean_degree", grepl("maximum_degree",raw)~"maximum_degree", grepl("degree_gini",raw)~"degree_gini", TRUE~"proportion_degree_1"),
         stat=case_when(grepl("empirical",raw)~"empirical", grepl("q025",raw)~"q025", grepl("q500",raw)~"q500", grepl("q975",raw)~"q975")) |>
  select(census_round,trophic_side,metric,stat,value) |> pivot_wider(names_from=stat, values_from=value)
p_deg <- ggplot(deg_long, aes(census_round)) + geom_errorbar(aes(ymin=q025,ymax=q975), width=0.1) + geom_point(aes(y=q500), shape=1) + geom_point(aes(y=empirical)) +
  facet_grid(metric~trophic_side, scales="free_y") + theme_classic() + xlab("Census round") + ylab("Metric") +
  ggtitle("Removal-0 degree summaries")
ggsave(file.path(out_dir,"02_zero_removal_degree_summary_empirical_vs_model.png"), p_deg, width=11, height=9, dpi=300)

notes <- c(
  "Sabatino census-round Galiana-style analysis",
  "",
  "1. Census rounds align each Sierra's first, second, third, etc. sampling occasion; they do not represent identical calendar dates or identical phenological conditions across Sierras.",
  "2. Differences among rounds may reflect seasonal turnover, phenology, sampling timing, variable effort, or variation in local community composition. They are not automatically independent replicates in a strict statistical sense.",
  "3. Cross-round regressions are descriptive and should not be overinterpreted.",
  "4. The primary value is to test whether pooled Sabatino results are stable across sampling rounds.",
  "5. Interaction counts and flower abundance are preserved but not used in the homogeneous model.",
  "6. Sierra removal is a structural/sampling diagnostic, not a habitat-loss experiment.",
  "7. Co-occurrence is operational: plant recorded in a Sierra during a census round + valid insect recorded in that Sierra during that census round."
)
writeLines(notes, file.path(out_dir,"02_interpretation_notes.txt"))

message("Finished census-round analysis. Outputs saved in: ", out_dir)
