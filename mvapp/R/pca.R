#' Principal Component Analysis (MVApp backend)
#'
#' Ports the legacy MVApp PCA tab. It reproduces the outputs the old app
#' computed with FactoMineR (eigenvalues, individual coordinates, variable
#' coordinates, and percentage contributions) using base [stats::prcomp()], so
#' the backend carries no heavy dependency and is fully unit-tested. The results
#' are equivalent to `FactoMineR::PCA()` on the same scaled matrix.
#'
#' Rows with missing phenotype values are dropped by default -- the legacy app
#' handled NAs the same way, via its "missing values removed" data option.
#'
#' Deviation from the legacy code (intentional): the old tab could `scale()` the
#' matrix *and then* call `PCA()` with its default `scale.unit = TRUE`, scaling
#' twice. Here `scale` maps once onto `prcomp(scale. = )`, which is the correct
#' single scaling.
#'
#' @param data A data frame or an [mvapp_dataset].
#' @param phenotypes Character vector of numeric phenotype columns to analyse.
#' @param id_cols Optional columns pasted together to label individuals (e.g.
#'   genotype, treatment, time, ID). If `data` is an mvapp_dataset and this is
#'   NULL, its id/group/time roles are used. Falls back to row numbers.
#' @param scale Scale variables to unit variance (default TRUE; recommended when
#'   traits are on different measurement scales).
#' @param na_action "omit" (drop incomplete rows, default) or "fail".
#' @return An object of class `mvapp_pca`: a list with data frames `eig`, `ind`,
#'   `var_coord`, `contrib`, plus `n`, `scaled`, and the underlying `prcomp` fit.
#' @export
run_pca <- function(data, phenotypes, id_cols = NULL, scale = TRUE,
                    na_action = c("omit", "fail")) {
  na_action <- match.arg(na_action)

  if (inherits(data, "mvapp_dataset")) {
    if (is.null(id_cols)) {
      id_cols <- c(data$roles$id, data$roles$group, data$roles$time)
    }
    data <- data$data
  }
  stopifnot(is.data.frame(data))

  miss <- setdiff(phenotypes, names(data))
  if (length(miss)) {
    stop("run_pca(): phenotypes not found: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (length(phenotypes) < 2L) {
    stop("run_pca(): need at least two phenotypes", call. = FALSE)
  }

  X <- data[, phenotypes, drop = FALSE]
  num <- vapply(X, is.numeric, logical(1))
  if (!all(num)) {
    stop("run_pca(): non-numeric phenotype columns: ",
         paste(phenotypes[!num], collapse = ", "), call. = FALSE)
  }
  X <- as.matrix(X)

  labels <- if (is.null(id_cols) || length(id_cols) == 0L) {
    as.character(seq_len(nrow(X)))
  } else {
    do.call(paste, c(data[id_cols], sep = "_"))
  }
  rownames(X) <- make.unique(labels)

  ok <- stats::complete.cases(X)
  if (any(!ok)) {
    if (na_action == "fail") {
      stop("run_pca(): ", sum(!ok), " rows contain NAs", call. = FALSE)
    }
    X <- X[ok, , drop = FALSE]
  }
  if (nrow(X) < 3L) {
    stop("run_pca(): too few complete rows for PCA (need >= 3)", call. = FALSE)
  }

  fit  <- stats::prcomp(X, center = TRUE, scale. = scale)
  sdev <- fit$sdev
  eig  <- sdev^2
  k    <- length(eig)
  pc   <- paste0("Dim.", seq_len(k))

  eig_df <- data.frame(
    component          = seq_len(k),
    eigenvalue         = eig,
    variance_percent   = 100 * eig / sum(eig),
    cumulative_percent = cumsum(100 * eig / sum(eig)),
    row.names          = NULL
  )

  ind <- as.data.frame(fit$x)
  names(ind) <- pc
  ind <- cbind(id = rownames(fit$x), ind)
  rownames(ind) <- NULL

  rot <- fit$rotation                      # variables x PCs (unit-norm columns)
  var_coord <- sweep(rot, 2, sdev, `*`)    # variable-PC correlation
  var_coord_df <- data.frame(trait = rownames(rot), var_coord, row.names = NULL)
  names(var_coord_df) <- c("trait", pc)

  contrib_df <- data.frame(trait = rownames(rot), 100 * rot^2, row.names = NULL)
  names(contrib_df) <- c("trait", pc)       # contributions sum to 100 per PC

  structure(
    list(eig = eig_df, ind = ind, var_coord = var_coord_df,
         contrib = contrib_df, n = nrow(X), scaled = scale, prcomp = fit),
    class = "mvapp_pca"
  )
}

#' @export
print.mvapp_pca <- function(x, ...) {
  cat("<mvapp_pca>  n =", x$n, " variables =", nrow(x$contrib),
      " scaled =", x$scaled, "\n")
  cat("  variance explained (first dims): ",
      paste0(sprintf("%.1f%%", utils::head(x$eig$variance_percent, 5)),
             collapse = ", "), "\n")
  invisible(x)
}
