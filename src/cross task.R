# Diagnostic nested cross-validation: all-50 feature set, fixed folds and seeds.
# The out-of-fold score is then applied unchanged to mortality in `cross_task`.
# Analysis-only all-50 ridge models; no composite-score comparator or graphics.
suppressPackageStartupMessages(library(glmnet))
seed_master<-20260909L;n_outer_folds<-5L;n_inner_folds<-5L;n_repeats<-20L;n_bootstrap<-2000L
paths<-list(tables=OUT)
write_csv<-function(x,path)fwrite(x,path,bom=TRUE,na='')
log_line<-function(...) {cat(...,'\n');flush.console()}
d236_all<-readinput('cross_task_GSE236713');d630_diag_all<-readinput('cross_task_GSE63042_diagnosis');d630_mort_all<-readinput('cross_task_GSE63042_mortality')
hallmark50<-grep('^HALLMARK_',names(d236_all),value=TRUE)
stopifnot(length(hallmark50)==50,setequal(hallmark50,grep('^HALLMARK_',names(d630_diag_all),value=TRUE)),setequal(hallmark50,grep('^HALLMARK_',names(d630_mort_all),value=TRUE)))
diag_sets<-list(GSE236713=copy(d236_all[disease=='Sepsis'|(disease=='SIRS'&disease_group=='OOHCA')]),GSE63042=copy(d630_diag_all[diagnostic_group%in%c('Sepsis','SIRS')]))
diag_sets$GSE236713[,diagnosis_y:=as.integer(disease=='Sepsis')];diag_sets$GSE63042[,diagnosis_y:=as.integer(diagnostic_group=='Sepsis')]
mort_sets<-list(GSE236713=copy(d236_all[disease=='Sepsis'&outcome%in%c('Died','Survived')]),GSE63042=copy(d630_mort_all[diagnostic_group=='Sepsis'&mortality_28day%in%c(0,1)]))
mort_sets$GSE236713[,mortality_y:=as.integer(outcome=='Died')];mort_sets$GSE63042[,mortality_y:=as.integer(mortality_28day)]
for(cohort in names(diag_sets)){
 for(nm in c('patient_id','sample_id')){set(diag_sets[[cohort]],j=nm,value=as.character(diag_sets[[cohort]][[nm]]));set(mort_sets[[cohort]],j=nm,value=as.character(mort_sets[[cohort]][[nm]]))}
 stopifnot(!anyDuplicated(diag_sets[[cohort]]$patient_id),!anyDuplicated(mort_sets[[cohort]]$patient_id),!anyNA(diag_sets[[cohort]][,..hallmark50]),!anyNA(mort_sets[[cohort]][,..hallmark50]),all(mort_sets[[cohort]]$patient_id%in%diag_sets[[cohort]]$patient_id))
}
stopifnot(nrow(diag_sets$GSE236713)==162,nrow(diag_sets$GSE63042)==129,nrow(mort_sets$GSE236713)==124,nrow(mort_sets$GSE63042)==106)
make_stratified_folds <- function(y, k, seed) {
  if (min(table(y)) < k) stop("A class has fewer observations than folds")
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in sort(unique(y))) {
    idx <- which(y == cl)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

auc_empirical <- function(score, y) {
  pos <- score[y == 1L]; neg <- score[y == 0L]
  ranks <- rank(c(pos, neg), ties.method = "average")
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

auc_boot_ci <- function(score, y, b = n_bootstrap, seed = seed_master) {
  pos <- score[y == 1L]; neg <- score[y == 0L]
  set.seed(seed)
  boots <- replicate(b, {
    bp <- sample(pos, length(pos), replace = TRUE)
    bn <- sample(neg, length(neg), replace = TRUE)
    ranks <- rank(c(bp, bn), ties.method = "average")
    (sum(ranks[seq_along(bp)]) - length(bp) * (length(bp) + 1) / 2) / (length(bp) * length(bn))
  })
  as.numeric(quantile(boots, c(.025, .975), names = FALSE, type = 6))
}

performance_row <- function(cohort, analysis, score, y, seed) {
  ci <- auc_boot_ci(score, y, seed = seed)
  data.table(
    cohort = cohort, analysis = analysis, n = length(y), n_class0 = sum(y == 0L), n_class1 = sum(y == 1L),
    AUC = auc_empirical(score, y), AUC_CI_low = ci[1], AUC_CI_high = ci[2],
    auc_direction = "Higher score predicts class 1", bootstrap_resamples = n_bootstrap, bootstrap_unit = "patient"
  )
}

nested_ridge_oof <- function(d, y_col, cohort, task, seed_base) {
  x <- as.matrix(d[, ..hallmark50]); y <- as.integer(d[[y_col]])
  ids <- as.character(d$patient_id); sample_ids <- as.character(d$sample_id)
  prediction_records <- list(); fold_records <- list(); parameter_records <- list()
  p_idx <- 0L; f_idx <- 0L; s_idx <- 0L
  for (repeat_id in seq_len(n_repeats)) {
    outer_seed <- seed_base + repeat_id * 1000L
    outer_fold <- make_stratified_folds(y, n_outer_folds, outer_seed)
    f_idx <- f_idx + 1L
    fold_records[[f_idx]] <- data.table(
      cohort = cohort, task = task, repeat_id = repeat_id, outer_fold = outer_fold, record_type = "outer_assignment",
      inner_for_outer_fold = NA_integer_, inner_fold = NA_integer_, patient_id = ids, sample_id = sample_ids, label = y,
      used_for_training = FALSE, used_for_outer_test = TRUE, seed = outer_seed
    )
    for (outer_id in seq_len(n_outer_folds)) {
      train_idx <- which(outer_fold != outer_id); test_idx <- which(outer_fold == outer_id)
      mu <- colMeans(x[train_idx, , drop = FALSE])
      sigma <- apply(x[train_idx, , drop = FALSE], 2, sd)
      sigma[!is.finite(sigma) | sigma == 0] <- 1
      x_train <- sweep(sweep(x[train_idx, , drop = FALSE], 2, mu, "-"), 2, sigma, "/")
      x_test <- sweep(sweep(x[test_idx, , drop = FALSE], 2, mu, "-"), 2, sigma, "/")
      inner_seed <- seed_base + repeat_id * 1000L + outer_id * 10L + 1L
      inner_fold <- make_stratified_folds(y[train_idx], n_inner_folds, inner_seed)
      cvfit <- cv.glmnet(
        x_train, y[train_idx], family = "binomial", alpha = 0, foldid = inner_fold,
        nfolds = n_inner_folds, type.measure = "deviance", standardize = FALSE, intercept = TRUE,
        parallel = FALSE, keep = FALSE
      )
      lambda_selected <- cvfit$lambda.min
      pred <- as.numeric(predict(cvfit, newx = x_test, s = "lambda.min", type = "response"))
      if (any(!is.finite(pred))) stop(sprintf("Non-finite OOF prediction: %s %s repeat %d fold %d", cohort, task, repeat_id, outer_id))
      p_idx <- p_idx + 1L
      prediction_records[[p_idx]] <- data.table(
        cohort = cohort, task = task, repeat_id = repeat_id, outer_fold = outer_id,
        patient_id = ids[test_idx], sample_id = sample_ids[test_idx], label = y[test_idx], oof_probability = pred
      )
      f_idx <- f_idx + 1L
      fold_records[[f_idx]] <- data.table(
        cohort = cohort, task = task, repeat_id = repeat_id, outer_fold = outer_fold[train_idx], record_type = "inner_assignment",
        inner_for_outer_fold = outer_id, inner_fold = inner_fold, patient_id = ids[train_idx], sample_id = sample_ids[train_idx], label = y[train_idx],
        used_for_training = TRUE, used_for_outer_test = FALSE, seed = inner_seed
      )
      s_idx <- s_idx + 1L
      parameter_records[[s_idx]] <- data.table(
        cohort = cohort, task = task, repeat_id = repeat_id, outer_fold = outer_id,
        pathway = hallmark50, outer_training_mean = as.numeric(mu), outer_training_sd = as.numeric(sigma),
        lambda_selected = lambda_selected, alpha = 0, outer_folds = n_outer_folds, inner_folds = n_inner_folds,
        standardize_in_glmnet = FALSE, preprocessing_fit_scope = "outer-training patients only", inner_type_measure = "deviance"
      )
    }
  }
  preds <- rbindlist(prediction_records)
  if (preds[, .N, by = .(repeat_id, patient_id)][, any(N != 1L)]) stop("Each patient must have exactly one OOF prediction per repeat")
  preds[, mean_oof_probability := mean(oof_probability), by = patient_id]
  list(
    predictions = preds,
    folds = rbindlist(fold_records, fill = TRUE),
    parameters = rbindlist(parameter_records)
  )
}

model_runs <- list()
for (i in seq_along(diag_sets)) {
  cohort <- names(diag_sets)[i]
  log_line("Starting nested CV diagnosis model: ", cohort)
  model_runs[[paste0(cohort, "_diagnosis")]] <- nested_ridge_oof(diag_sets[[cohort]], "diagnosis_y", cohort, "Diagnosis", seed_master + i * 100000L)
}

all_oof <- rbindlist(lapply(model_runs, `[[`, "predictions"))
all_folds <- rbindlist(lapply(model_runs, `[[`, "folds"), fill = TRUE)
all_parameters <- rbindlist(lapply(model_runs, `[[`, "parameters"))
write_csv(all_oof, file.path(paths$tables, "08_patient_level_oof_predictions.csv"))
write_csv(all_folds, file.path(paths$tables, "10_cross_task_fold_assignments.csv"))
write_csv(all_parameters, file.path(paths$tables, "11_cross_task_outer_scaling_and_lambda.csv"))

get_final_oof <- function(cohort_name, task_name) unique(all_oof[all_oof$cohort == cohort_name & all_oof$task == task_name, .(patient_id, sample_id, label, mean_oof_probability)])

diagnosis_perf <- rbindlist(lapply(seq_along(diag_sets), function(i) {
  cohort <- names(diag_sets)[i]; p <- get_final_oof(cohort, "Diagnosis")
  performance_row(cohort, "Diagnostic model OOF -> diagnosis", p$mean_oof_probability, p$label, seed_master + 700000L + i)
}))
write_csv(diagnosis_perf, file.path(paths$tables, "02_within_cohort_diagnosis_oof_performance.csv"))

transfer_mortality <- rbindlist(lapply(seq_along(mort_sets), function(i) {
  cohort <- names(mort_sets)[i]
  p <- get_final_oof(cohort, "Diagnosis")
  m <- mort_sets[[cohort]][, .(patient_id, mortality_y)]
  z <- merge(m, p[, .(patient_id, OOF_Diagnostic_Score = mean_oof_probability)], by = "patient_id", all.x = TRUE, sort = FALSE)
  if (anyNA(z$OOF_Diagnostic_Score)) stop(sprintf("%s is missing strict OOF diagnostic scores for mortality patients", cohort))
  performance_row(cohort, "Same diagnostic OOF score -> mortality", z$OOF_Diagnostic_Score, z$mortality_y, seed_master + 710000L + i)
}))
write_csv(transfer_mortality, file.path(paths$tables, "03_diagnostic_score_to_mortality_transfer.csv"))


for(i in seq_len(nrow(rbindlist(list(diagnosis_perf,transfer_mortality))))){
 d<-rbindlist(list(diagnosis_perf,transfer_mortality));q<-d[i];e<-expected('cross_task_diagnosis')[cohort==q$cohort&analysis==q$analysis];stopifnot(nrow(e)==1)
 for(nm in c('n','n_class0','n_class1','AUC','AUC_CI_low','AUC_CI_high','bootstrap_resamples'))check('FIGURE1_NESTED_CV',paste(q$cohort,q$analysis),nm,q[[nm]],e[[nm]],inputs=c('cross_task_GSE236713','cross_task_GSE63042_diagnosis','cross_task_GSE63042_mortality'),n=q$n,estimand='Diagnostic-only nested OOF; unchanged diagnostic predictions evaluated for mortality')
}
