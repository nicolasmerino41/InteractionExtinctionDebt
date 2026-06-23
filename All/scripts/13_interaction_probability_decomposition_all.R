
## ------------------------------------------------------------
## Script: All/scripts/13_interaction_probability_decomposition_all.R
##
## Purpose:
## Decompose site-level interaction probability:
## Why are realised interactions arranged across sites in a way that
## makes metaweb links more persistent than expected from one global p?
##
## Outputs:
##   All/SeparatedResults/script13_interaction_probability_decomposition/
##   All/CombinedOutputs/
##
## Rules:
##   - all datasets in all_dataset_names
##   - write.csv2() for all tables
##   - combined plots only
##   - robust model failure handling
##   - Windows-safe parallelisation across datasets
## ------------------------------------------------------------

source("All/scripts/00_dataset_loaders_and_helpers_all.R")

parallel_packages <- c("future", "future.apply", "parallelly")
for(pkg in parallel_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(123)

removal_levels <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8)
run_cross_validation <- FALSE
cv_pair_models <- FALSE
n_site_reps <- 20

cv_folds <- 5
max_cv_rows <- Inf
prob_eps <- 1e-6

dirs <- make_output_dirs("script13_interaction_probability_decomposition")
sep_out <- dirs$separated
combined_out <- dirs$combined

n_workers <- max(1, parallelly::availableCores() - 1)
future::plan(future::multisession, workers = n_workers)
message("Using ", n_workers, " parallel workers across datasets.")

model_specs <- data.frame(
  model_type = c("M0_global", "M1_consumer", "M2_resource", "M3_site", "M6_pair"),
  formula_text = c(
    "interaction_observed ~ 1",
    "interaction_observed ~ consumer",
    "interaction_observed ~ resource",
    "interaction_observed ~ site",
    "interaction_observed ~ pair_id"
  ),
  is_pair_model = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

clamp_prob <- function(p){
  pmin(pmax(p, prob_eps), 1 - prob_eps)
}

log_loss_binary <- function(y, p){
  p <- clamp_prob(p)
  -mean(y * log(p) + (1 - y) * log(1 - p), na.rm = TRUE)
}

brier_binary <- function(y, p){
  p <- clamp_prob(p)
  mean((y - p)^2, na.rm = TRUE)
}

safe_relative_error <- function(observed, expected){
  if(is.na(expected) || expected == 0) return(NA_real_)
  (observed - expected) / expected
}

make_pair_id <- function(consumer, resource){
  paste(consumer, resource, sep = "___")
}

make_site_level_table <- function(dataset, site_tables){

  cooc <- site_tables$cooc_triples %>%
    distinct(site, consumer, resource) %>%
    mutate(
      dataset = dataset,
      pair_id = make_pair_id(consumer, resource)
    )

  ints <- site_tables$empirical_site_interactions %>%
    distinct(site, consumer, resource) %>%
    mutate(interaction_observed = 1L)

  cooc %>%
    left_join(ints, by = c("site", "consumer", "resource")) %>%
    mutate(
      interaction_observed = replace_na(interaction_observed, 0L),
      site = factor(site),
      consumer = factor(consumer),
      resource = factor(resource),
      pair_id = factor(pair_id)
    ) %>%
    select(dataset, site, consumer, resource, pair_id, interaction_observed)
}

safe_fit_glm <- function(dat, formula_text){

  form <- as.formula(formula_text)

  fit <- tryCatch(
    suppressWarnings(glm(form, data = dat, family = binomial())),
    error = function(e) e
  )

  if(inherits(fit, "error")){
    return(list(success = FALSE, fit = NULL, message = fit$message))
  }

  list(success = TRUE, fit = fit, message = "")
}

safe_predict <- function(fit, newdata, fallback_p){

  if(is.null(fit)){
    return(rep(fallback_p, nrow(newdata)))
  }

  p <- tryCatch(
    predict(fit, newdata = newdata, type = "response"),
    error = function(e) rep(fallback_p, nrow(newdata))
  )

  clamp_prob(as.numeric(p))
}

fit_one_model_summary <- function(dataset, dat, model_type, formula_text){

  n_opportunities <- nrow(dat)
  n_interactions <- sum(dat$interaction_observed)
  site_level_conversion <- n_interactions / n_opportunities
  fallback_p <- clamp_prob(site_level_conversion)

  fit_obj <- safe_fit_glm(dat, formula_text)

  if(!fit_obj$success){
    return(list(
      fit = NULL,
      summary = data.frame(
        dataset = dataset,
        model_type = model_type,
        n_opportunities = n_opportunities,
        n_interactions = n_interactions,
        site_level_conversion = site_level_conversion,
        AIC = NA_real_,
        residual_deviance = NA_real_,
        null_deviance = NA_real_,
        deviance_explained = NA_real_,
        log_loss = NA_real_,
        brier_score = NA_real_,
        fit_success = FALSE,
        fit_message = fit_obj$message
      ),
      predicted = rep(fallback_p, nrow(dat))
    ))
  }

  fit <- fit_obj$fit
  predicted <- safe_predict(fit, dat, fallback_p)

  residual_deviance <- tryCatch(fit$deviance, error = function(e) NA_real_)
  null_deviance <- tryCatch(fit$null.deviance, error = function(e) NA_real_)

  deviance_explained <- ifelse(
    is.na(null_deviance) || null_deviance == 0,
    NA_real_,
    1 - residual_deviance / null_deviance
  )

  summary <- data.frame(
    dataset = dataset,
    model_type = model_type,
    n_opportunities = n_opportunities,
    n_interactions = n_interactions,
    site_level_conversion = site_level_conversion,
    AIC = tryCatch(AIC(fit), error = function(e) NA_real_),
    residual_deviance = residual_deviance,
    null_deviance = null_deviance,
    deviance_explained = deviance_explained,
    log_loss = log_loss_binary(dat$interaction_observed, predicted),
    brier_score = brier_binary(dat$interaction_observed, predicted),
    fit_success = TRUE,
    fit_message = fit_obj$message
  )

  list(fit = fit, summary = summary, predicted = predicted)
}

make_cv_folds <- function(n, k){
  sample(rep(seq_len(k), length.out = n))
}

cv_one_model <- function(dataset, dat, model_type, formula_text, k = 5){

  if(nrow(dat) < k * 2 || length(unique(dat$interaction_observed)) < 2){
    return(data.frame(
      dataset = dataset,
      model_type = model_type,
      cv_log_loss = NA_real_,
      cv_brier_score = NA_real_,
      cv_success = FALSE,
      cv_message = "too few rows or only one response class"
    ))
  }

  dat_cv <- dat

  if(is.finite(max_cv_rows) && nrow(dat_cv) > max_cv_rows){
    dat_cv <- dat_cv[sample(seq_len(nrow(dat_cv)), max_cv_rows), , drop = FALSE]
  }

  folds <- make_cv_folds(nrow(dat_cv), k)
  fold_metrics <- list()

  for(fold in seq_len(k)){

    train <- droplevels(dat_cv[folds != fold, , drop = FALSE])
    test <- droplevels(dat_cv[folds == fold, , drop = FALSE])
    fallback_p <- clamp_prob(mean(train$interaction_observed))

    fit_obj <- safe_fit_glm(train, formula_text)

    if(!fit_obj$success){
      fold_metrics[[fold]] <- data.frame(
        fold = fold,
        log_loss = NA_real_,
        brier_score = NA_real_,
        success = FALSE,
        message = fit_obj$message
      )
      next
    }

    p <- safe_predict(fit_obj$fit, test, fallback_p)

    fold_metrics[[fold]] <- data.frame(
      fold = fold,
      log_loss = log_loss_binary(test$interaction_observed, p),
      brier_score = brier_binary(test$interaction_observed, p),
      success = TRUE,
      message = ""
    )
  }

  folds_out <- bind_rows(fold_metrics)

  data.frame(
    dataset = dataset,
    model_type = model_type,
    cv_log_loss = mean(folds_out$log_loss, na.rm = TRUE),
    cv_brier_score = mean(folds_out$brier_score, na.rm = TRUE),
    cv_success = any(folds_out$success),
    cv_message = paste(unique(folds_out$message[folds_out$message != ""]), collapse = " | ")
  ) %>%
    mutate(
      cv_log_loss = ifelse(is.nan(cv_log_loss), NA_real_, cv_log_loss),
      cv_brier_score = ifelse(is.nan(cv_brier_score), NA_real_, cv_brier_score)
    )
}

run_cv_for_dataset <- function(dataset, dat){

  if(!run_cross_validation){
    return(data.frame(
      dataset = dataset,
      model_type = model_specs$model_type,
      cv_log_loss = NA_real_,
      cv_brier_score = NA_real_,
      cv_success = FALSE,
      cv_message = "cross-validation disabled"
    ))
  }

  bind_rows(lapply(seq_len(nrow(model_specs)), function(i){

    spec <- model_specs[i, ]

    if(spec$is_pair_model && !cv_pair_models){
      return(data.frame(
        dataset = dataset,
        model_type = spec$model_type,
        cv_log_loss = NA_real_,
        cv_brier_score = NA_real_,
        cv_success = FALSE,
        cv_message = "pair-model CV disabled"
      ))
    }

    message(dataset, ": CV ", spec$model_type)

    cv_one_model(
      dataset = dataset,
      dat = dat,
      model_type = spec$model_type,
      formula_text = spec$formula_text,
      k = cv_folds
    )
  }))
}

make_site_subsets_for_dataset <- function(sites){
  make_site_subsets(as.character(sites), removal_levels, n_site_reps)
}

empirical_retained_links <- function(dat, sites_keep){
  dat %>%
    filter(site %in% sites_keep, interaction_observed == 1) %>%
    distinct(pair_id) %>%
    nrow()
}

expected_retained_links_from_predicted_p <- function(dat, sites_keep, predicted_p){

  cur <- dat %>%
    mutate(predicted_p = predicted_p) %>%
    filter(site %in% sites_keep)

  if(nrow(cur) == 0){
    return(NA_real_)
  }

  cur %>%
    group_by(pair_id) %>%
    summarise(
      p_pair_retained = 1 - prod(1 - predicted_p),
      .groups = "drop"
    ) %>%
    summarise(expected_retained_links = sum(p_pair_retained, na.rm = TRUE)) %>%
    pull(expected_retained_links)
}

retained_link_prediction_for_model <- function(dataset, dat, predicted_p, model_type,
                                               subset_object,
                                               empirical_full_links,
                                               expected_full_links){

  bind_rows(lapply(seq_len(nrow(subset_object$index)), function(i){

    sites_keep <- subset_object$subsets[[i]]
    subset_row <- subset_object$index[i, ]

    empirical_links <- empirical_retained_links(dat, sites_keep)

    expected_links <- expected_retained_links_from_predicted_p(
      dat = dat,
      sites_keep = sites_keep,
      predicted_p = predicted_p
    )

    data.frame(
      dataset = dataset,
      model_type = model_type,
      subset_id = subset_row$subset_id,
      removal_fraction = subset_row$removal_fraction,
      site_rep = subset_row$site_rep,
      n_sites_kept = subset_row$n_sites_kept,
      empirical_retained_links = empirical_links,
      expected_retained_links = expected_links,
      absolute_error = empirical_links - expected_links,
      relative_error = safe_relative_error(empirical_links, expected_links),
      empirical_retained_link_proportion = empirical_links / empirical_full_links,
      expected_retained_link_proportion_relative_to_expected_full = expected_links / expected_full_links,
      expected_retained_link_proportion_relative_to_empirical_full = expected_links / empirical_full_links
    )
  }))
}

run_one_dataset <- function(dataset){

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  message("Running script 13: ", dataset)

  out_dir <- file.path(sep_out, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  site_tables <- get_dataset_site_tables(dataset)
  dat <- make_site_level_table(dataset, site_tables)

  empirical_full_links <- dat %>%
    filter(interaction_observed == 1) %>%
    distinct(pair_id) %>%
    nrow()

  all_sites <- sort(unique(as.character(dat$site)))
  subset_object <- make_site_subsets_for_dataset(all_sites)

  fit_summaries <- list()
  predicted_list <- list()

  for(i in seq_len(nrow(model_specs))){
    spec <- model_specs[i, ]
    message(dataset, ": fitting ", spec$model_type)

    fit_res <- fit_one_model_summary(
      dataset = dataset,
      dat = dat,
      model_type = spec$model_type,
      formula_text = spec$formula_text
    )

    fit_summaries[[spec$model_type]] <- fit_res$summary
    predicted_list[[spec$model_type]] <- fit_res$predicted
  }

  fit_summary <- bind_rows(fit_summaries)

  m0_log <- fit_summary$log_loss[fit_summary$model_type == "M0_global"][1]
  m0_brier <- fit_summary$brier_score[fit_summary$model_type == "M0_global"][1]
  m0_dev <- fit_summary$deviance_explained[fit_summary$model_type == "M0_global"][1]

  fit_summary <- fit_summary %>%
    mutate(
      delta_log_loss_vs_M0 = m0_log - log_loss,
      delta_brier_vs_M0 = m0_brier - brier_score,
      delta_deviance_explained_vs_M0 = deviance_explained - m0_dev
    )

  cv_summary <- run_cv_for_dataset(dataset, dat)

  m0_cv_log <- cv_summary$cv_log_loss[cv_summary$model_type == "M0_global"][1]
  m0_cv_brier <- cv_summary$cv_brier_score[cv_summary$model_type == "M0_global"][1]

  cv_summary <- cv_summary %>%
    mutate(
      delta_cv_log_loss_vs_M0 = m0_cv_log - cv_log_loss,
      delta_cv_brier_vs_M0 = m0_cv_brier - cv_brier_score
    )

  retained_predictions <- bind_rows(lapply(names(predicted_list), function(model_type){

    predicted_p <- predicted_list[[model_type]]

    expected_full <- expected_retained_links_from_predicted_p(
      dat = dat,
      sites_keep = all_sites,
      predicted_p = predicted_p
    )

    retained_link_prediction_for_model(
      dataset = dataset,
      dat = dat,
      predicted_p = predicted_p,
      model_type = model_type,
      subset_object = subset_object,
      empirical_full_links = empirical_full_links,
      expected_full_links = expected_full
    )
  }))

  retained_summary <- retained_predictions %>%
    group_by(dataset, model_type, removal_fraction) %>%
    summarise(
      mean_empirical_retained_links = mean(empirical_retained_links, na.rm = TRUE),
      sd_empirical_retained_links = sd(empirical_retained_links, na.rm = TRUE),
      mean_expected_retained_links = mean(expected_retained_links, na.rm = TRUE),
      sd_expected_retained_links = sd(expected_retained_links, na.rm = TRUE),
      mean_absolute_error = mean(absolute_error, na.rm = TRUE),
      sd_absolute_error = sd(absolute_error, na.rm = TRUE),
      mean_relative_error = mean(relative_error, na.rm = TRUE),
      sd_relative_error = sd(relative_error, na.rm = TRUE),
      RMSE = sqrt(mean(absolute_error^2, na.rm = TRUE)),
      MAE = mean(abs(absolute_error), na.rm = TRUE),
      bias = mean(absolute_error, na.rm = TRUE),
      mean_empirical_retained_link_proportion = mean(empirical_retained_link_proportion, na.rm = TRUE),
      mean_expected_retained_link_proportion_relative_to_expected_full =
        mean(expected_retained_link_proportion_relative_to_expected_full, na.rm = TRUE),
      mean_expected_retained_link_proportion_relative_to_empirical_full =
        mean(expected_retained_link_proportion_relative_to_empirical_full, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv2(
    fit_summary,
    file.path(out_dir, paste0(dataset, "_script13_site_level_model_fit_summary.csv")),
    row.names = FALSE
  )

  write.csv2(
    cv_summary,
    file.path(out_dir, paste0(dataset, "_script13_site_level_model_cv_summary.csv")),
    row.names = FALSE
  )

  write.csv2(
    retained_summary,
    file.path(out_dir, paste0(dataset, "_script13_retained_link_prediction_summary.csv")),
    row.names = FALSE
  )

  list(
    fit_summary = fit_summary,
    cv_summary = cv_summary,
    retained_predictions = retained_predictions,
    retained_summary = retained_summary
  )
}

all_outputs <- future.apply::future_lapply(
  all_dataset_names,
  run_one_dataset,
  future.seed = TRUE
)

names(all_outputs) <- all_dataset_names

fit_summary_all <- bind_rows(lapply(all_outputs, `[[`, "fit_summary"))
cv_summary_all <- bind_rows(lapply(all_outputs, `[[`, "cv_summary"))
retained_prediction_all <- bind_rows(lapply(all_outputs, `[[`, "retained_predictions"))
retained_summary_all <- bind_rows(lapply(all_outputs, `[[`, "retained_summary"))

write.csv2(
  fit_summary_all,
  file.path(combined_out, "script13_site_level_model_fit_summary_combined.csv"),
  row.names = FALSE
)

write.csv2(
  cv_summary_all,
  file.path(combined_out, "script13_site_level_model_cv_summary_combined.csv"),
  row.names = FALSE
)

write.csv2(
  retained_prediction_all,
  file.path(combined_out, "script13_retained_link_prediction_by_subset_combined.csv"),
  row.names = FALSE
)

write.csv2(
  retained_summary_all,
  file.path(combined_out, "script13_retained_link_prediction_summary_combined.csv"),
  row.names = FALSE
)

fit_summary_all$dataset <- factor(fit_summary_all$dataset, levels = all_dataset_names)
cv_summary_all$dataset <- factor(cv_summary_all$dataset, levels = all_dataset_names)
retained_summary_all$dataset <- factor(retained_summary_all$dataset, levels = all_dataset_names)

plot_fit <- cv_summary_all %>%
  select(dataset, model_type, delta_cv_log_loss_vs_M0, delta_cv_brier_vs_M0) %>%
  left_join(
    fit_summary_all %>%
      select(dataset, model_type, delta_log_loss_vs_M0, delta_brier_vs_M0),
    by = c("dataset", "model_type")
  ) %>%
  mutate(
    improvement_metric = ifelse(!is.na(delta_cv_log_loss_vs_M0),
                                delta_cv_log_loss_vs_M0,
                                delta_log_loss_vs_M0),
    metric_source = ifelse(!is.na(delta_cv_log_loss_vs_M0),
                           "CV log loss",
                           "in-sample log loss")
  )

p_fit <- ggplot(plot_fit,
                aes(x = model_type,
                    y = improvement_metric)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_col() +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Improvement over M0_global") +
  ggtitle("Site-level model fit improvement relative to M0",
          subtitle = "Uses CV log loss when available; otherwise in-sample log loss")

ggsave(
  file.path(combined_out, "13_model_fit_improvement_relative_to_M0.png"),
  p_fit,
  width = 14,
  height = 7,
  dpi = 300
)

p_error <- ggplot(retained_summary_all,
                  aes(x = removal_fraction,
                      y = mean_absolute_error,
                      linetype = model_type)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  xlab("Fraction of sites removed") +
  ylab("Empirical retained links - expected retained links") +
  ggtitle("Retained-link prediction error by site-level model")

ggsave(
  file.path(combined_out, "13_retained_link_prediction_error_by_model.png"),
  p_error,
  width = 14,
  height = 7,
  dpi = 300
)

p_obs_exp <- ggplot(retained_summary_all,
                    aes(x = mean_empirical_retained_links,
                        y = mean_expected_retained_links,
                        shape = model_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free") +
  theme_classic(base_size = 10) +
  xlab("Empirical retained links") +
  ylab("Expected retained links") +
  ggtitle("Observed vs expected retained links by site-level model")

ggsave(
  file.path(combined_out, "13_observed_vs_expected_retained_links_by_model.png"),
  p_obs_exp,
  width = 14,
  height = 7,
  dpi = 300
)

component_importance <- plot_fit %>%
  mutate(
    component = case_when(
      model_type == "M1_consumer" ~ "consumer",
      model_type == "M2_resource" ~ "resource",
      model_type == "M3_site" ~ "site",
      model_type == "M4_consumer_resource" ~ "consumer+resource",
      model_type == "M5_consumer_resource_site" ~ "consumer+resource+site",
      model_type == "M6_pair" ~ "pair",
      model_type == "M7_pair_site" ~ "pair+site",
      TRUE ~ "global"
    )
  ) %>%
  filter(model_type != "M0_global")

write.csv2(
  component_importance,
  file.path(combined_out, "script13_component_importance_summary.csv"),
  row.names = FALSE
)

p_component <- ggplot(component_importance,
                      aes(x = component,
                          y = improvement_metric)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2) +
  facet_wrap(~ dataset, ncol = 5, scales = "free_y") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model component") +
  ylab("Improvement over M0_global") +
  ggtitle("Which component explains site-level interaction occurrence?",
          subtitle = "Higher values mean better performance than M0")

ggsave(
  file.path(combined_out, "13_component_importance_summary.png"),
  p_component,
  width = 14,
  height = 7,
  dpi = 300
)

future::plan(future::sequential)

message("Finished script 13 interaction probability decomposition.")
