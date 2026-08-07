#' Structural-equation / path model  [Roadmap Phase 4b]
#'
#' Tests a *directed* hypothesis among traits (complementary to the predictive
#' importance in [rf_importance()]). TODO (Phase 4b): wrap lavaan; return the
#' fitted model, path coefficients, fit indices (CFI/RMSEA/SRMR), and a diagram.
#'
#' @param data A data frame (or `mvapp_dataset`).
#' @param model A lavaan model specification (trait ~ trait ...).
#' @return A fitted path model (to be defined in Phase 4).
#' @export
path_model <- function(data, model) {
  stop("path_model(): not yet implemented -- see Roadmap Phase 4b", call. = FALSE)
}
