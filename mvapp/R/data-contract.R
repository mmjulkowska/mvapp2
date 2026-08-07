#' The canonical MVApp dataset
#'
#' Every analysis in MVApp consumes one structure: a tidy long-format data
#' frame plus a *role map* that names which columns are identifiers, grouping
#' factors, an optional time column, and the measured phenotypes. Both the
#' generic CSV importer and the PhenoSight importer produce this object, so all
#' downstream code (stats, plots, ML) is agnostic to the data source.
#'
#' @param data A data frame in long format.
#' @param id Character vector of identifier column names (e.g. plant/line ID).
#' @param group Character vector of grouping/factor columns (e.g. genotype, treatment).
#' @param time Optional single column name holding time (day or TOE hours).
#' @param phenotypes Character vector of measured phenotype columns.
#' @param source Provenance tag, e.g. "generic", "phenosight_ps2".
#' @return An object of class `mvapp_dataset`.
#' @export
new_mvapp_dataset <- function(data, id = character(), group = character(),
                              time = NULL, phenotypes = character(),
                              source = "generic") {
  stopifnot(is.data.frame(data))
  roles <- list(id = id, group = group, time = time, phenotypes = phenotypes)
  structure(list(data = as.data.frame(data, check.names = FALSE),
                 roles = roles, source = source),
            class = "mvapp_dataset")
}

#' @rdname new_mvapp_dataset
#' @export
validate_mvapp_dataset <- function(x) {
  stopifnot(inherits(x, "mvapp_dataset"))
  cols <- names(x$data)
  named <- unlist(x$roles[c("id", "group", "time", "phenotypes")],
                  use.names = FALSE)
  missing <- setdiff(named, cols)
  if (length(missing)) {
    stop("mvapp_dataset: columns not found in data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(x$roles$phenotypes) == 0) {
    stop("mvapp_dataset: at least one phenotype column is required",
         call. = FALSE)
  }
  invisible(x)
}

#' Construct and validate an mvapp_dataset
#' @inheritParams new_mvapp_dataset
#' @export
mvapp_dataset <- function(data, id = character(), group = character(),
                          time = NULL, phenotypes = character(),
                          source = "generic") {
  validate_mvapp_dataset(
    new_mvapp_dataset(data, id, group, time, phenotypes, source)
  )
}

#' @export
print.mvapp_dataset <- function(x, ...) {
  cat("<mvapp_dataset>  source:", x$source, "\n")
  cat("  rows:", nrow(x$data), " cols:", ncol(x$data), "\n")
  cat("  id:        ", paste(x$roles$id, collapse = ", "), "\n")
  cat("  group:     ", paste(x$roles$group, collapse = ", "), "\n")
  cat("  time:      ", if (is.null(x$roles$time)) "<none>" else x$roles$time, "\n")
  cat("  phenotypes:", paste(x$roles$phenotypes, collapse = ", "), "\n")
  invisible(x)
}
