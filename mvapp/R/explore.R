# Data Exploration backends: normality, homogeneity of variance, t-tests, and
# two-way ANOVA. One-way ANOVA + Tukey letters is provided by
# anova_tukey_letters() (see anova.R) and reused by the module.

.explore_df <- function(data) {
  if (inherits(data, "mvapp_dataset")) data <- data$data
  stopifnot(is.data.frame(data))
  data
}

#' Shapiro-Wilk normality test, overall or per group
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric column to test.
#' @param group Optional factor; a test per level when supplied.
#' @return A data frame: group, n, W (Shapiro statistic), p. Groups with n < 3
#'   (or n > 5000, randomly subsampled to 5000, Shapiro's limit) are handled.
#' @export
shapiro_by_group <- function(data, trait, group = NULL) {
  data <- .explore_df(data)
  if (!trait %in% names(data)) {
    stop("shapiro_by_group(): trait not found: ", trait, call. = FALSE)
  }
  one <- function(v) {
    v <- v[is.finite(v)]
    n <- length(v)
    if (n < 3) return(data.frame(n = n, W = NA_real_, p = NA_real_))
    vv <- if (n > 5000) sample(v, 5000) else v
    s <- stats::shapiro.test(vv)
    data.frame(n = n, W = unname(s$statistic), p = s$p.value)
  }
  if (is.null(group) || !group %in% names(data)) {
    out <- cbind(group = "all", one(data[[trait]]), stringsAsFactors = FALSE)
  } else {
    parts <- split(data[[trait]], as.character(data[[group]]))
    out <- do.call(rbind, lapply(names(parts), function(g)
      cbind(group = g, one(parts[[g]]), stringsAsFactors = FALSE)))
  }
  rownames(out) <- NULL
  out
}

#' Test homogeneity of variance across groups
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric response.
#' @param group Grouping factor.
#' @param method "levene" (Brown-Forsythe, median-centred; base-R implementation)
#'   or "bartlett".
#' @return A list: method, statistic, df, p.
#' @export
variance_test <- function(data, trait, group, method = c("levene", "bartlett")) {
  method <- match.arg(method)
  data <- .explore_df(data)
  for (v in c(trait, group)) {
    if (!v %in% names(data)) stop("variance_test(): column not found: ", v, call. = FALSE)
  }
  d <- data[stats::complete.cases(data[c(trait, group)]), c(trait, group)]
  names(d) <- c("value", "grp"); d$grp <- factor(d$grp)
  if (nlevels(d$grp) < 2) stop("variance_test(): need >= 2 groups", call. = FALSE)

  if (method == "bartlett") {
    b <- stats::bartlett.test(value ~ grp, data = d)
    return(list(method = "Bartlett", statistic = unname(b$statistic),
                df = unname(b$parameter), p = b$p.value))
  }
  # Levene's test, Brown-Forsythe variant: ANOVA on |value - group median|
  med <- tapply(d$value, d$grp, stats::median)
  d$z <- abs(d$value - med[as.character(d$grp)])
  tab <- summary(stats::aov(z ~ grp, data = d))[[1]]
  list(method = "Levene (Brown-Forsythe)", statistic = tab[["F value"]][1],
       df = paste(tab[["Df"]], collapse = ", "), p = tab[["Pr(>F)"]][1])
}

#' Two-sample t-test of a trait between two groups
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric response.
#' @param group Factor with exactly two levels.
#' @param var_equal Assume equal variances (Student) vs Welch (default FALSE).
#' @param paired Paired t-test (default FALSE).
#' @param pair_by When paired, the categorical column that defines the matched
#'   pairs (e.g. "ACCESSION"). Values are averaged per `pair_by` x `group` cell,
#'   then the two group means are paired across `pair_by` levels (only levels
#'   present in both groups are used). Required for a paired test.
#' @return A list with group levels, per-group means, the mean difference, the
#'   t statistic, df, p, 95% CI, the test label, and (paired) `n_pairs`.
#' @export
t_test_two <- function(data, trait, group, var_equal = FALSE, paired = FALSE,
                       pair_by = NULL) {
  data <- .explore_df(data)
  cols <- c(trait, group, if (!is.null(pair_by)) pair_by)
  for (v in cols) {
    if (!v %in% names(data)) stop("t_test_two(): column not found: ", v, call. = FALSE)
  }
  d <- data[stats::complete.cases(data[cols]), cols]
  names(d)[1:2] <- c("value", "grp"); d$grp <- factor(d$grp)
  if (nlevels(d$grp) != 2) {
    stop("t_test_two(): need exactly 2 groups, got ", nlevels(d$grp), call. = FALSE)
  }
  lv <- levels(d$grp)

  if (isTRUE(paired)) {
    if (is.null(pair_by)) {
      stop("t_test_two(): a paired test needs 'pair_by' (the pairing variable, ",
           "e.g. genotype).", call. = FALSE)
    }
    names(d)[3] <- "pairid"
    agg <- stats::aggregate(value ~ pairid + grp, d, mean)   # average per cell
    w <- stats::reshape(agg, idvar = "pairid", timevar = "grp", direction = "wide")
    x <- w[[paste0("value.", lv[1])]]; y <- w[[paste0("value.", lv[2])]]
    ok <- is.finite(x) & is.finite(y)                        # complete pairs only
    x <- x[ok]; y <- y[ok]
    if (length(x) < 2) stop("t_test_two(): too few complete pairs", call. = FALSE)
    tt <- stats::t.test(x, y, paired = TRUE)
    return(list(groups = lv, mean1 = mean(x), mean2 = mean(y),
      diff = mean(x) - mean(y), statistic = unname(tt$statistic),
      df = unname(tt$parameter), p = tt$p.value,
      conf_int = as.numeric(tt$conf.int), n_pairs = length(x),
      method = sprintf("Paired t-test (paired by %s, %d pairs)", pair_by, length(x))))
  }

  tt <- stats::t.test(value ~ grp, data = d, var.equal = var_equal)
  m <- tapply(d$value, d$grp, mean)
  list(groups = lv, mean1 = unname(m[1]), mean2 = unname(m[2]),
       diff = unname(m[1] - m[2]), statistic = unname(tt$statistic),
       df = unname(tt$parameter), p = tt$p.value,
       conf_int = as.numeric(tt$conf.int), method = tt$method)
}

#' Kruskal-Wallis test with pairwise-Wilcoxon compact-letter display
#'
#' The non-parametric analogue of one-way ANOVA + Tukey: it ranks the data (no
#' normality or equal-variance assumption), tests whether groups differ
#' (Kruskal-Wallis), and derives significance letters from pairwise Wilcoxon
#' (Mann-Whitney) tests with multiplicity adjustment. Returns the same shape as
#' [anova_tukey_letters()] so it plugs into the same plot.
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric response.
#' @param group Grouping factor (2 groups -> Wilcoxon/Mann-Whitney; more ->
#'   Kruskal-Wallis).
#' @param p_adjust Multiple-testing adjustment for the pairwise tests (default "BH").
#' @return An object of class `mvapp_anova_cld`: `summary` (group, n, mean, sd,
#'   se, median, ymax, letter -- letters ordered so "a" is the highest median),
#'   `anova_p` (the Kruskal-Wallis p), and `test`.
#' @export
kruskal_letters <- function(data, trait, group, p_adjust = "BH") {
  data <- .explore_df(data)
  for (v in c(trait, group)) {
    if (!v %in% names(data)) stop("kruskal_letters(): column not found: ", v, call. = FALSE)
  }
  d <- data[stats::complete.cases(data[c(trait, group)]), c(trait, group)]
  names(d) <- c("value", "grp"); d$grp <- factor(d$grp)
  keep <- names(which(table(d$grp) >= 2))
  d <- d[d$grp %in% keep, ]; d$grp <- droplevels(d$grp)
  if (nlevels(d$grp) < 2) {
    stop("kruskal_letters(): need >= 2 groups with >= 2 observations", call. = FALSE)
  }
  kw <- stats::kruskal.test(value ~ grp, data = d)
  # order groups by descending median so letter "a" marks the highest
  meds <- tapply(d$value, d$grp, stats::median)
  d$grp <- factor(d$grp, levels = names(sort(meds, decreasing = TRUE)))
  levs <- levels(d$grp)
  pw <- suppressWarnings(
    stats::pairwise.wilcox.test(d$value, d$grp, p.adjust.method = p_adjust))
  combos <- utils::combn(levs, 2)
  pvec <- apply(combos, 2, function(cc) {
    val <- suppressWarnings(pw$p.value[cc[2], cc[1]])
    if (is.null(val) || is.na(val)) 1 else val
  })
  names(pvec) <- apply(combos, 2, paste, collapse = "-")
  letters <- multcompView::multcompLetters(pvec)$Letters

  agg <- do.call(rbind, lapply(split(d$value, d$grp), function(v) {
    data.frame(n = length(v), mean = mean(v), sd = stats::sd(v),
               se = stats::sd(v) / sqrt(length(v)), median = stats::median(v),
               ymax = max(v))
  }))
  agg$group <- rownames(agg); rownames(agg) <- NULL
  agg$letter <- unname(letters[agg$group])
  agg <- agg[match(levs, agg$group), ]
  structure(list(
    summary = agg[c("group", "n", "mean", "sd", "se", "median", "ymax", "letter")],
    anova_p = kw$p.value, test = "Kruskal-Wallis + pairwise Wilcoxon"),
    class = "mvapp_anova_cld")
}

#' Screen many traits with a Kruskal-Wallis test
#'
#' Non-parametric analogue of [anova_screen()]: `kruskal.test(trait ~ group)`
#' for every trait, ranked by p-value.
#'
#' @param data A data frame or [mvapp_dataset].
#' @param traits Character vector of numeric traits.
#' @param group Grouping factor.
#' @return A data frame (trait, n_groups, kruskal_p) ordered by increasing p.
#' @export
kruskal_screen <- function(data, traits, group) {
  data <- .explore_df(data)
  if (!group %in% names(data)) stop("kruskal_screen(): group not found: ", group, call. = FALSE)
  rows <- lapply(traits, function(tr) {
    if (!tr %in% names(data) || !is.numeric(data[[tr]])) return(NULL)
    d <- data[stats::complete.cases(data[c(tr, group)]), c(tr, group)]
    names(d) <- c("value", "grp"); d$grp <- factor(d$grp)
    keep <- names(which(table(d$grp) >= 2))
    d <- d[d$grp %in% keep, ]; d$grp <- droplevels(d$grp)
    p <- if (nlevels(d$grp) < 2) NA_real_ else
      tryCatch(stats::kruskal.test(value ~ grp, d)$p.value,
               error = function(e) NA_real_)
    data.frame(trait = tr, n_groups = nlevels(d$grp), kruskal_p = p)
  })
  out <- do.call(rbind, rows)
  out[order(out$kruskal_p), , drop = FALSE]
}

#' Kruskal-Wallis trait screen, split by a facet variable
#'
#' @inheritParams anova_screen_faceted
#' @return A data frame (trait + one p column per facet level).
#' @export
kruskal_screen_faceted <- function(data, traits, group, facet) {
  data <- .explore_df(data)
  if (!facet %in% names(data)) stop("kruskal_screen_faceted(): facet not found: ",
                                    facet, call. = FALSE)
  levs <- sort(unique(as.character(data[[facet]])))
  merged <- data.frame(trait = traits, stringsAsFactors = FALSE)
  for (lv in levs) {
    sub <- data[as.character(data[[facet]]) == lv, , drop = FALSE]
    sc <- kruskal_screen(sub, traits, group)
    m <- sc[c("trait", "kruskal_p")]; names(m)[2] <- lv
    merged <- merge(merged, m, by = "trait", all.x = TRUE, sort = FALSE)
  }
  pcols <- setdiff(names(merged), "trait")
  key <- suppressWarnings(apply(merged[pcols], 1,
    function(p) min(p, na.rm = TRUE)))
  merged[order(key), , drop = FALSE]
}

#' Screen many traits with a two-sample t-test
#'
#' @param data A data frame or [mvapp_dataset].
#' @param traits Character vector of numeric traits.
#' @param group Grouping factor with two levels.
#' @return A data frame (trait, p) ordered by increasing p.
#' @export
ttest_screen <- function(data, traits, group) {
  data <- .explore_df(data)
  if (!group %in% names(data)) stop("ttest_screen(): group not found: ", group, call. = FALSE)
  rows <- lapply(traits, function(tr) {
    if (!tr %in% names(data) || !is.numeric(data[[tr]])) return(NULL)
    d <- data[stats::complete.cases(data[c(tr, group)]), c(tr, group)]
    names(d) <- c("value", "grp"); d$grp <- factor(d$grp)
    p <- if (nlevels(d$grp) != 2) NA_real_ else
      tryCatch(stats::t.test(value ~ grp, d)$p.value, error = function(e) NA_real_)
    data.frame(trait = tr, p = p)
  })
  out <- do.call(rbind, rows)
  out[order(out$p), , drop = FALSE]
}

#' Screen many traits with a one-way ANOVA
#'
#' Runs `trait ~ group` for every trait and returns them ranked by p-value, so
#' you can see up front which traits differ across the groups.
#'
#' @param data A data frame or [mvapp_dataset].
#' @param traits Character vector of numeric traits to screen.
#' @param group Grouping factor.
#' @return A data frame (trait, n_groups, anova_p) ordered by increasing p.
#' @export
anova_screen <- function(data, traits, group) {
  data <- .explore_df(data)
  if (!group %in% names(data)) stop("anova_screen(): group not found: ", group, call. = FALSE)
  rows <- lapply(traits, function(tr) {
    if (!tr %in% names(data) || !is.numeric(data[[tr]])) return(NULL)
    d <- data[stats::complete.cases(data[c(tr, group)]), c(tr, group)]
    names(d) <- c("value", "grp"); d$grp <- factor(d$grp)
    keep <- names(which(table(d$grp) >= 2))
    d <- d[d$grp %in% keep, ]; d$grp <- droplevels(d$grp)
    p <- if (nlevels(d$grp) < 2) NA_real_ else
      summary(stats::aov(value ~ grp, d))[[1]][["Pr(>F)"]][1]
    data.frame(trait = tr, n_groups = nlevels(d$grp), anova_p = p)
  })
  out <- do.call(rbind, rows)
  out[order(out$anova_p), , drop = FALSE]
}

#' One-way ANOVA trait screen, split by a facet variable
#'
#' Runs [anova_screen()] separately within each level of `facet` and returns a
#' wide table: one row per trait, one p-value column per facet level.
#'
#' @param data A data frame or [mvapp_dataset].
#' @param traits Character vector of numeric traits.
#' @param group Grouping factor.
#' @param facet Column giving the facet subsets.
#' @return A data frame (trait + one p column per facet level), rows ordered by
#'   the smallest p across facets.
#' @export
anova_screen_faceted <- function(data, traits, group, facet) {
  data <- .explore_df(data)
  if (!facet %in% names(data)) stop("anova_screen_faceted(): facet not found: ",
                                    facet, call. = FALSE)
  levs <- sort(unique(as.character(data[[facet]])))
  merged <- data.frame(trait = traits, stringsAsFactors = FALSE)
  for (lv in levs) {
    sub <- data[as.character(data[[facet]]) == lv, , drop = FALSE]
    sc <- anova_screen(sub, traits, group)
    m <- sc[c("trait", "anova_p")]; names(m)[2] <- lv
    merged <- merge(merged, m, by = "trait", all.x = TRUE, sort = FALSE)
  }
  pcols <- setdiff(names(merged), "trait")
  key <- suppressWarnings(apply(merged[pcols], 1,
    function(p) min(p, na.rm = TRUE)))
  merged[order(key), , drop = FALSE]
}

#' Two-sample t-test trait screen, split by a facet variable
#'
#' @inheritParams anova_screen_faceted
#' @return A data frame (trait + one p column per facet level).
#' @export
ttest_screen_faceted <- function(data, traits, group, facet) {
  data <- .explore_df(data)
  if (!facet %in% names(data)) stop("ttest_screen_faceted(): facet not found: ",
                                    facet, call. = FALSE)
  levs <- sort(unique(as.character(data[[facet]])))
  merged <- data.frame(trait = traits, stringsAsFactors = FALSE)
  for (lv in levs) {
    sub <- data[as.character(data[[facet]]) == lv, , drop = FALSE]
    sc <- ttest_screen(sub, traits, group)
    names(sc)[2] <- lv
    merged <- merge(merged, sc, by = "trait", all.x = TRUE, sort = FALSE)
  }
  pcols <- setdiff(names(merged), "trait")
  key <- suppressWarnings(apply(merged[pcols], 1,
    function(p) min(p, na.rm = TRUE)))
  merged[order(key), , drop = FALSE]
}

#' Screen many traits with a two-way ANOVA
#'
#' @param data A data frame or [mvapp_dataset].
#' @param traits Character vector of numeric traits to screen.
#' @param factor_a,factor_b The two grouping factors.
#' @return A data frame (trait, p_A, p_B, p_interaction) ordered by the smallest
#'   of the three p-values.
#' @export
anova2_screen <- function(data, traits, factor_a, factor_b) {
  data <- .explore_df(data)
  rows <- lapply(traits, function(tr) {
    if (!tr %in% names(data) || !is.numeric(data[[tr]])) return(NULL)
    r <- tryCatch(anova_two_way(data, tr, factor_a, factor_b),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    t <- r$table
    getp <- function(term) { p <- t$p[t$term == term]; if (length(p)) p else NA_real_ }
    data.frame(trait = tr, p_A = getp(factor_a), p_B = getp(factor_b),
               p_interaction = getp(paste0(factor_a, ":", factor_b)))
  })
  out <- do.call(rbind, rows)
  out[order(pmin(out$p_A, out$p_B, out$p_interaction, na.rm = TRUE)), , drop = FALSE]
}

#' Two-way ANOVA of a trait across two factors (with interaction)
#'
#' @param data A data frame or [mvapp_dataset].
#' @param trait Numeric response.
#' @param factor_a,factor_b The two grouping factors.
#' @return A list: `table` (tidy ANOVA table, terms named after the factors),
#'   the fitted `aov`, and `means` (cell means for an interaction plot).
#' @export
anova_two_way <- function(data, trait, factor_a, factor_b) {
  data <- .explore_df(data)
  for (v in c(trait, factor_a, factor_b)) {
    if (!v %in% names(data)) stop("anova_two_way(): column not found: ", v, call. = FALSE)
  }
  if (factor_a == factor_b) stop("anova_two_way(): choose two different factors", call. = FALSE)
  d <- data[stats::complete.cases(data[c(trait, factor_a, factor_b)]),
            c(trait, factor_a, factor_b)]
  names(d) <- c("value", "A", "B"); d$A <- factor(d$A); d$B <- factor(d$B)
  fit <- stats::aov(value ~ A * B, data = d)
  tab <- as.data.frame(summary(fit)[[1]])
  tab$term <- trimws(rownames(tab)); rownames(tab) <- NULL
  tab$term <- c(factor_a, factor_b, paste0(factor_a, ":", factor_b),
                "Residuals")[match(tab$term, c("A", "B", "A:B", "Residuals"))]
  names(tab) <- c("Df", "Sum_Sq", "Mean_Sq", "F_value", "p", "term")
  tab <- tab[c("term", "Df", "Sum_Sq", "Mean_Sq", "F_value", "p")]
  means <- stats::aggregate(value ~ A + B, d, mean)
  names(means) <- c(factor_a, factor_b, "mean")
  list(table = tab, fit = fit, factor_a = factor_a, factor_b = factor_b,
       means = means)
}
