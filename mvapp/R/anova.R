# ANOVA + Tukey HSD with compact-letter display, and cluster diagnostics.
# This is the Julkowska-lab signature test: one-way ANOVA of a trait across
# groups (e.g. clusters), Tukey HSD post-hoc, and significance letters ordered
# by mean -- ready to drop onto a boxplot. Reused by both cluster tabs (and,
# later, the Data Exploration ANOVA tab).

#' One-way ANOVA + Tukey HSD with compact-letter display
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric response column.
#' @param group Grouping factor column (e.g. "cluster", "TREATMENT").
#' @return An object of class `mvapp_anova_cld`: `summary` (one row per group:
#'   group, n, mean, sd, se, ymax, letter -- letters ordered so "a" is the
#'   highest mean), the overall `anova_p`, and the `tukey` object.
#' @export
anova_tukey_letters <- function(data, trait, group) {
  if (inherits(data, "mvapp_dataset")) data <- data$data
  stopifnot(is.data.frame(data))
  for (v in c(trait, group)) {
    if (!v %in% names(data)) {
      stop("anova_tukey_letters(): column not found: ", v, call. = FALSE)
    }
  }
  d <- data[stats::complete.cases(data[c(trait, group)]), c(trait, group)]
  names(d) <- c("value", "grp")
  d$grp <- factor(d$grp)
  # keep only replicated groups (Tukey needs >= 2 obs per group)
  keep <- names(which(table(d$grp) >= 2))
  d <- d[d$grp %in% keep, ]
  d$grp <- droplevels(d$grp)
  if (nlevels(d$grp) < 2) {
    stop("anova_tukey_letters(): need >= 2 groups with >= 2 observations each",
         call. = FALSE)
  }

  fit <- stats::aov(value ~ grp, data = d)
  anova_p <- summary(fit)[[1]][["Pr(>F)"]][1]
  tuk <- stats::TukeyHSD(fit)
  cld <- multcompView::multcompLetters4(fit, tuk)
  letters <- cld$grp$Letters   # named by group, ordered by decreasing mean

  agg <- do.call(rbind, lapply(split(d$value, d$grp), function(v) {
    data.frame(n = length(v), mean = mean(v), sd = stats::sd(v),
               se = stats::sd(v) / sqrt(length(v)), ymax = max(v))
  }))
  agg$group <- rownames(agg); rownames(agg) <- NULL
  agg$letter <- unname(letters[agg$group])
  agg <- agg[order(-agg$mean), c("group", "n", "mean", "sd", "se", "ymax", "letter")]

  structure(list(summary = agg, anova_p = anova_p, tukey = tuk),
            class = "mvapp_anova_cld")
}

#' @export
print.mvapp_anova_cld <- function(x, ...) {
  cat("<mvapp_anova_cld>  ANOVA p =", signif(x$anova_p, 3), "\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Mean between- and within-cluster distances
#'
#' @param hres An `mvapp_hclust` object (from [run_hclust()]).
#' @return A symmetric cluster-by-cluster matrix: off-diagonal = mean distance
#'   between members of the two clusters; diagonal = mean within-cluster distance.
#' @export
cluster_distances <- function(hres) {
  stopifnot(inherits(hres, "mvapp_hclust"))
  D <- as.matrix(hres$dist)
  cl <- hres$clusters$cluster
  ks <- sort(unique(cl))
  M <- matrix(NA_real_, length(ks), length(ks),
              dimnames = list(paste0("C", ks), paste0("C", ks)))
  for (i in seq_along(ks)) {
    for (j in seq_along(ks)) {
      sub <- D[cl == ks[i], cl == ks[j], drop = FALSE]
      if (ks[i] == ks[j]) {
        M[i, j] <- if (sum(cl == ks[i]) > 1) mean(sub[upper.tri(sub)]) else 0
      } else {
        M[i, j] <- mean(sub)
      }
    }
  }
  M
}

#' K-means diagnostics: elbow (WSS), average silhouette, and gap statistic
#'
#' Ports the MVApp "optimal number of clusters" tab -- three ways to pick k.
#'
#' @param data A data frame or [mvapp_dataset].
#' @param phenotypes Numeric trait columns.
#' @param id_cols,scale See [run_kmeans()].
#' @param max_k Largest k to evaluate (default 10).
#' @param B Bootstrap samples for the gap statistic (default 25).
#' @param gap Whether to compute the (expensive) gap statistic. Elbow and
#'   silhouette are always computed; set `FALSE` to skip only the gap column.
#' @param nstart Random restarts per k (default 5 -- diagnostics don't need the
#'   heavier restart count used for the final fit).
#' @param seed RNG seed (default 1).
#' @return A list with `table` (k, wss, silhouette, gap, gap_SE) and `best`
#'   (suggested k by silhouette and by the gap firstSEmax rule).
#' @export
cluster_diagnostics <- function(data, phenotypes, id_cols = NULL, scale = TRUE,
                                max_k = 10, B = 25, gap = TRUE, nstart = 5,
                                seed = 1) {
  X <- .trait_matrix(data, phenotypes, id_cols, scale)
  max_k <- max(2L, min(as.integer(max_k), nrow(X) - 1L))
  D <- stats::dist(X)

  wss <- numeric(max_k); sil <- rep(NA_real_, max_k)
  for (k in seq_len(max_k)) {
    if (!is.null(seed)) set.seed(seed)
    km <- stats::kmeans(X, centers = k, nstart = nstart, iter.max = 50)
    wss[k] <- km$tot.withinss
    if (k >= 2) sil[k] <- mean(cluster::silhouette(km$cluster, D)[, 3])
  }

  gap_col <- rep(NA_real_, max_k); gap_se <- rep(NA_real_, max_k)
  best_gap <- NA_integer_
  if (isTRUE(gap)) {
    if (!is.null(seed)) set.seed(seed)
    g <- cluster::clusGap(X, FUN = stats::kmeans, K.max = max_k, B = B,
                          nstart = nstart, iter.max = 50)
    gap_col <- g$Tab[, "gap"]; gap_se <- g$Tab[, "SE.sim"]
    best_gap <- cluster::maxSE(gap_col, gap_se, method = "firstSEmax")
  }
  tab <- data.frame(k = seq_len(max_k), wss = wss, silhouette = sil,
                    gap = gap_col, gap_SE = gap_se)
  list(table = tab,
       best = list(silhouette = which.max(sil), gap = best_gap))
}
