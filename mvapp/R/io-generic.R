#' Read a long-format phenotype CSV into an mvapp_dataset
#'
#' The generic MVApp input: one row per observation, with identifier and factor
#' columns followed by measured phenotypes (e.g. LINE, ACCESSION, TREATMENT,
#' DAY, then trait columns).
#'
#' @param path Path to a CSV file.
#' @param id,group,time,phenotypes Column-role assignments (see [mvapp_dataset()]).
#'   If `phenotypes` is NULL, every column not named as id/group/time is treated
#'   as a phenotype.
#' @param source Provenance tag stored on the dataset.
#' @param ... Passed to [utils::read.csv()].
#' @return An `mvapp_dataset`.
#' @export
read_long_csv <- function(path, id, group, time = NULL, phenotypes = NULL,
                          source = "generic", ...) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
  if (is.null(phenotypes)) {
    phenotypes <- setdiff(names(df), c(id, group, time))
  }
  mvapp_dataset(df, id = id, group = group, time = time,
                phenotypes = phenotypes, source = source)
}
