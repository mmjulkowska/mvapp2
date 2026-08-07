# Multidimensional scaling (classical / metric MDS) -- ports the MVApp MDS tab.
# Base R (dist + cmdscale); no heavy dependency; fully testable.

# shared: build a labelled, complete-case (optionally scaled) trait matrix
.trait_matrix <- function(data, phenotypes, id_cols = NULL, scale = TRUE) {
  if (inherits(data, "mvapp_dataset")) {
    if (is.null(id_cols)) {
      id_cols <- c(data$roles$id, data$roles$group, data$roles$time)
    }
    data <- data$data
  }
  stopifnot(is.data.frame(data))
  miss <- setdiff(phenotypes, names(data))
  if (length(miss)) {
    stop("phenotypes not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (length(phenotypes) < 2L) {
    stop("need at least two phenotypes", call. = FALSE)
  }
  X <- data[, phenotypes, drop = FALSE]
  num <- vapply(X, is.numeric, logical(1))
  if (!all(num)) {
    stop("non-numeric phenotype columns: ",
         paste(phenotypes[!num], collapse = ", "), call. = FALSE)
  }
  X <- as.matrix(X)
  labels <- if (is.null(id_cols) || length(id_cols) == 0L) {
    as.character(seq_len(nrow(X)))
  } else {
    do.call(paste, c(data[id_cols], sep = "_"))
  }
  rownames(X) <- make.unique(labels)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < 3L) stop("too few complete rows (need >= 3)", call. = FALSE)
  if (isTRUE(scale)) X <- scale(X)
  X
}

#' Classical (metric) multidimensional scaling
#'
#' Ports the MVApp MDS tab: distances between samples across the chosen traits,
#' embedded into `k` dimensions with [stats::cmdscale()].
#'
#' @param data A data frame or [mvapp_dataset].
#' @param phenotypes Numeric trait columns.
#' @param id_cols Columns pasted to label samples (from dataset roles if NULL).
#' @param scale Scale traits to unit variance first (default TRUE).
#' @param k Number of MDS dimensions (default 2).
#' @param dist_method Distance metric for [stats::dist()] (default "euclidean").
#' @return An object of class `mvapp_mds`: `points` (id + MDS dims), `eig`,
#'   goodness-of-fit `gof`, plus `n`, `k`, `scaled`.
#' @export
run_mds <- function(data, phenotypes, id_cols = NULL, scale = TRUE, k = 2,
                    dist_method = "euclidean") {
  X <- .trait_matrix(data, phenotypes, id_cols, scale)
  k <- min(k, nrow(X) - 1L)
  D <- stats::dist(X, method = dist_method)
  fit <- stats::cmdscale(D, k = k, eig = TRUE)
  pts <- as.data.frame(fit$points)
  names(pts) <- paste0("MDS", seq_len(ncol(pts)))
  pts <- cbind(id = rownames(X), pts)
  rownames(pts) <- NULL
  structure(list(points = pts, eig = fit$eig, gof = fit$GOF,
                 n = nrow(X), k = ncol(fit$points), scaled = scale,
                 dist_method = dist_method),
            class = "mvapp_mds")
}

#' @export
print.mvapp_mds <- function(x, ...) {
  cat("<mvapp_mds>  n =", x$n, " dims =", x$k, " scaled =", x$scaled,
      " GOF =", paste(round(x$gof, 3), collapse = "/"), "\n")
  invisible(x)
}
