# Clustering -- ports the MVApp hierarchical + k-means tabs. Base R
# (dist/hclust/cutree/kmeans); the trait-matrix prep is shared with MDS via
# .trait_matrix() (see mds.R).

#' Hierarchical clustering of samples
#'
#' @param data A data frame or [mvapp_dataset].
#' @param phenotypes Numeric trait columns.
#' @param id_cols Columns pasted to label samples (from dataset roles if NULL).
#' @param scale Scale traits to unit variance first (default TRUE).
#' @param dist_method Distance metric for [stats::dist()] (default "euclidean").
#' @param link_method Linkage for [stats::hclust()] (default "complete").
#' @param k Number of clusters to cut the tree into (default 3).
#' @return An object of class `mvapp_hclust`: the `hclust` fit, a `clusters`
#'   data frame (id + cluster), the `dist` object, plus `k`, `n`, `scaled`.
#' @export
run_hclust <- function(data, phenotypes, id_cols = NULL, scale = TRUE,
                       dist_method = "euclidean", link_method = "complete",
                       k = 3) {
  X <- .trait_matrix(data, phenotypes, id_cols, scale)
  k <- max(1L, min(as.integer(k), nrow(X) - 1L))
  D <- stats::dist(X, method = dist_method)
  hc <- stats::hclust(D, method = link_method)
  cl <- stats::cutree(hc, k = k)
  clusters <- data.frame(id = rownames(X), cluster = unname(cl),
                         row.names = NULL)
  structure(list(hclust = hc, clusters = clusters, dist = D, x = X, k = k,
                 n = nrow(X), scaled = scale, dist_method = dist_method,
                 link_method = link_method),
            class = "mvapp_hclust")
}

#' K-means clustering of samples
#'
#' @inheritParams run_hclust
#' @param k Number of clusters (default 3).
#' @param nstart Random restarts for [stats::kmeans()] (default 25).
#' @param max_k Largest k evaluated for the elbow (within-SS) curve (default 10).
#' @param seed Optional RNG seed for reproducibility (default 1; set NULL to skip).
#' @return An object of class `mvapp_kmeans`: the `kmeans` fit, a `clusters`
#'   data frame (id + cluster), an elbow `wss` data frame (k + total within-SS),
#'   plus `k`, `n`, `scaled`.
#' @export
run_kmeans <- function(data, phenotypes, id_cols = NULL, scale = TRUE, k = 3,
                       nstart = 25, max_k = 10, seed = 1) {
  X <- .trait_matrix(data, phenotypes, id_cols, scale)
  k <- max(1L, min(as.integer(k), nrow(X) - 1L))
  if (!is.null(seed)) set.seed(seed)
  km <- stats::kmeans(X, centers = k, nstart = nstart, iter.max = 50)
  clusters <- data.frame(id = rownames(X), cluster = unname(km$cluster),
                         row.names = NULL)
  max_k <- min(as.integer(max_k), nrow(X) - 1L)
  wss <- vapply(seq_len(max_k), function(kk) {
    if (!is.null(seed)) set.seed(seed)
    stats::kmeans(X, centers = kk, nstart = nstart, iter.max = 50)$tot.withinss
  }, numeric(1))
  structure(list(kmeans = km, clusters = clusters,
                 wss = data.frame(k = seq_len(max_k), tot_withinss = wss),
                 k = k, n = nrow(X), scaled = scale),
            class = "mvapp_kmeans")
}

#' @export
print.mvapp_hclust <- function(x, ...) {
  cat("<mvapp_hclust>  n =", x$n, " k =", x$k, " linkage =", x$link_method, "\n")
  invisible(x)
}

#' @export
print.mvapp_kmeans <- function(x, ...) {
  cat("<mvapp_kmeans>  n =", x$n, " k =", x$k,
      " between/total SS =", round(x$kmeans$betweenss / x$kmeans$totss, 3), "\n")
  invisible(x)
}
