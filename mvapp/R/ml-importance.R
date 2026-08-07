#' Cross-method phenotype importance  [Roadmap Phase 4a]
#'
#' Answers "which phenotypes predict/explain a target phenotype, and how much".
#' TODO (Phase 4a): fit random forest (ranger), LASSO/elastic net (glmnet), and
#' gradient boosting (xgboost) under a common tidymodels resampling scheme;
#' return a tidy importance table on a shared footing plus out-of-sample metrics.
#'
#' @param data A data frame (or `mvapp_dataset`).
#' @param target Name of the phenotype to explain.
#' @param predictors Character vector of predictor phenotypes.
#' @param methods Which learners to run.
#' @return A tidy data frame of per-method importance (to be defined in Phase 4).
#' @export
rf_importance <- function(data, target, predictors,
                          methods = c("rf", "lasso", "gbm")) {
  stop("rf_importance(): not yet implemented -- see Roadmap Phase 4a", call. = FALSE)
}
