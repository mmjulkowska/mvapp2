#' Within-group Stress Tolerance Index (STI)
#'
#' Each stressed observation divided by the mean of *its own group's* controls
#' (the regime the lab uses when comparing genotypes in a Control-vs-Drought
#' design). Controls return `NA` (the index is defined on stressed plants).
#'
#' @param value Numeric phenotype (e.g. digital biomass).
#' @param group Grouping factor identifying the comparison unit (e.g. genotype).
#' @param condition Treatment factor (control vs stress).
#' @param control_level The baseline level of `condition`.
#' @return Numeric vector of STI values (NA for control rows).
#' @export
within_group_sti <- function(value, group, condition, control_level = "Control") {
  n <- length(value)
  stopifnot(length(group) == n, length(condition) == n)
  is_ctrl <- condition == control_level
  ctrl_means <- tapply(value[is_ctrl], group[is_ctrl], mean, na.rm = TRUE)
  denom <- ctrl_means[as.character(group)]
  as.numeric(ifelse(is_ctrl, NA_real_, value / denom))
}

#' Classic genotype-ranking tolerance indices
#'
#' @param drought,control Numeric vectors of the trait under drought and control.
#' @return A list with `STI1 = D/C` and `STI2 = D/sqrt(C)`.
#' @export
sti_indices <- function(drought, control) {
  list(STI1 = drought / control,
       STI2 = drought / sqrt(control))
}
