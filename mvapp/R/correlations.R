# Correlation backend for MVApp 2.0 -- ported from the legacy Correlations tab.
# All statistics live here (tested); the Shiny module only renders.

#' Coerce a data frame or mvapp_dataset to a numeric phenotype matrix
#' @keywords internal
.cor_matrix_input <- function(data, phenotypes) {
  if (inherits(data, "mvapp_dataset")) {
    if (missing(phenotypes) || is.null(phenotypes)) {
      phenotypes <- data$roles$phenotypes
    }
    data <- data$data
  }
  stopifnot(is.data.frame(data))
  miss <- setdiff(phenotypes, names(data))
  if (length(miss)) {
    stop("correlations: phenotypes not found: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  X <- data[, phenotypes, drop = FALSE]
  num <- vapply(X, is.numeric, logical(1))
  if (!all(num)) {
    stop("correlations: non-numeric phenotype columns: ",
         paste(phenotypes[!num], collapse = ", "), call. = FALSE)
  }
  if (length(phenotypes) < 2L) {
    stop("correlations: need at least two phenotypes", call. = FALSE)
  }
  as.matrix(X)
}

#' Correlation matrix of phenotype columns
#'
#' @param data A data frame or [mvapp_dataset].
#' @param phenotypes Numeric columns to correlate (defaults to the dataset's
#'   phenotype roles when `data` is an mvapp_dataset).
#' @param method "pearson" (default) or "spearman".
#' @param use Passed to [stats::cor()]; default "pairwise.complete.obs".
#' @return A square numeric correlation matrix.
#' @export
cor_matrix <- function(data, phenotypes = NULL, method = c("pearson", "spearman"),
                       use = "pairwise.complete.obs") {
  method <- match.arg(method)
  X <- .cor_matrix_input(data, phenotypes)
  stats::cor(X, method = method, use = use)
}

#' Matrix of correlation p-values (pairwise)
#'
#' Replaces `corrplot::cor.mtest` with a tested base-R implementation, used to
#' cross out non-significant correlations. Diagonal is 0.
#'
#' @inheritParams cor_matrix
#' @return A square matrix of two-sided p-values with the same dim/names as
#'   [cor_matrix()].
#' @export
cor_pmatrix <- function(data, phenotypes = NULL, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  X <- .cor_matrix_input(data, phenotypes)
  p <- ncol(X)
  P <- matrix(0, p, p, dimnames = list(colnames(X), colnames(X)))
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      pv <- tryCatch(
        stats::cor.test(X[, i], X[, j], method = method)$p.value,
        error = function(e) NA_real_)
      P[i, j] <- P[j, i] <- pv
    }
  }
  P
}

#' Number of complete-case observations across the chosen phenotypes
#' @inheritParams cor_matrix
#' @return Integer count of rows with no missing values in `phenotypes`.
#' @export
cor_n <- function(data, phenotypes = NULL) {
  X <- .cor_matrix_input(data, phenotypes)
  sum(stats::complete.cases(X))
}

#' Correlation matrices split by a factor (for an upper/lower-triangle overlay)
#'
#' Computes two correlation matrices over the *same* phenotypes -- one per level
#' of a splitting factor -- so a single correlation plot can show one level in
#' its upper triangle and the other in its lower triangle (a compact way to
#' compare, e.g., Control vs Drought trait correlations side by side).
#'
#' @param data A data frame or [mvapp_dataset].
#' @param phenotypes Numeric columns to correlate.
#' @param by Name of the splitting factor column.
#' @param upper_level,lower_level The two levels of `by` for the upper and lower
#'   triangles. Default to the first two levels present in the data.
#' @param method "pearson" (default) or "spearman".
#' @return A list with `upper`/`lower` correlation matrices (identical row/col
#'   order), the chosen `upper_level`/`lower_level`, and `n_upper`/`n_lower`.
#' @export
cor_split <- function(data, phenotypes = NULL, by,
                      upper_level = NULL, lower_level = NULL,
                      method = c("pearson", "spearman")) {
  method <- match.arg(method)
  if (inherits(data, "mvapp_dataset")) {
    if (is.null(phenotypes)) phenotypes <- data$roles$phenotypes
    data <- data$data
  }
  stopifnot(is.data.frame(data))
  if (!by %in% names(data)) {
    stop("cor_split(): split factor not found: ", by, call. = FALSE)
  }
  levs <- unique(as.character(data[[by]]))
  if (is.null(upper_level)) upper_level <- levs[1]
  if (is.null(lower_level)) {
    lower_level <- levs[if (length(levs) >= 2L) 2L else 1L]
  }
  if (identical(as.character(upper_level), as.character(lower_level))) {
    stop("cor_split(): upper and lower levels must differ", call. = FALSE)
  }
  for (lv in c(upper_level, lower_level)) {
    if (!as.character(lv) %in% levs) {
      stop("cor_split(): level '", lv, "' not present in '", by, "'",
           call. = FALSE)
    }
  }
  du <- data[as.character(data[[by]]) == as.character(upper_level), , drop = FALSE]
  dl <- data[as.character(data[[by]]) == as.character(lower_level), , drop = FALSE]
  list(
    upper = cor_matrix(du, phenotypes, method = method),
    lower = cor_matrix(dl, phenotypes, method = method),
    p_upper = cor_pmatrix(du, phenotypes, method = method),
    p_lower = cor_pmatrix(dl, phenotypes, method = method),
    upper_level = upper_level, lower_level = lower_level,
    n_upper = cor_n(du, phenotypes), n_lower = cor_n(dl, phenotypes)
  )
}

#' Per-group regression stats for a scatter of two traits
#'
#' For the Scatterplots view: reports n, Pearson r, R-squared (the coefficient
#' of determination of the least-squares line shown on the plot, = r^2), and the
#' p-value -- overall, or one row per level of `group` when the scatter is
#' coloured by a factor (so every regression line displayed gets its own R^2).
#'
#' @param data A data frame or [mvapp_dataset].
#' @param x,y Trait column names (x = predictor, y = response).
#' @param group Optional factor column; one row of stats per level.
#' @return A data frame with columns `group`, `n`, `r`, `R2`, `p`.
#' @export
scatter_stats <- function(data, x, y, group = NULL) {
  if (inherits(data, "mvapp_dataset")) data <- data$data
  stopifnot(is.data.frame(data))
  for (v in c(x, y)) {
    if (!v %in% names(data)) {
      stop("scatter_stats(): column not found: ", v, call. = FALSE)
    }
  }
  one <- function(df) {
    df <- df[stats::complete.cases(df[c(x, y)]), , drop = FALSE]
    n <- nrow(df)
    if (n < 3 || stats::sd(df[[x]]) == 0 || stats::sd(df[[y]]) == 0) {
      return(data.frame(n = n, r = NA_real_, R2 = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(stats::cor.test(df[[x]], df[[y]], method = "pearson"))
    r <- unname(ct$estimate)
    data.frame(n = n, r = r, R2 = r^2, p = ct$p.value)
  }
  if (is.null(group) || !group %in% names(data)) {
    out <- cbind(group = "all", one(data), stringsAsFactors = FALSE)
  } else {
    parts <- split(data, as.character(data[[group]]))
    out <- do.call(rbind, lapply(names(parts), function(g)
      cbind(group = g, one(parts[[g]]), stringsAsFactors = FALSE)))
  }
  rownames(out) <- NULL
  out
}
