# Reproducibility helpers: generate runnable R code and default figure legends
# for each analysis. Pure + tested, so the Shiny "Show me the code" / "Show
# figure legend" features (and the log-book) all draw from one source of truth.

.rlines <- function(x) paste(x[!vapply(x, is.null, logical(1))], collapse = "\n")
.rvec   <- function(x) paste0("c(", paste0('"', x, '"', collapse = ", "), ")")

# palette line for generated corrplot code (Brewer or ggsci)
.code_palette_line <- function(palette, levels) {
  if (startsWith(palette, "ggsci:")) {
    key <- trimws(sub("ggsci:", "", palette))
    if (identical(key, "GSEA")) {
      'pal <- colorRampPalette(ggsci::pal_gsea()(12))(200)'
    } else {
      sprintf('pal <- colorRampPalette(ggsci::pal_%s()(10))(200)',
              tolower(key))
    }
  } else {
    sprintf('pal <- colorRampPalette(RColorBrewer::brewer.pal(%d, "%s"))(200)',
            as.integer(levels), palette)
  }
}

#' Generate R code that reproduces the correlation plot
#'
#' Returns a runnable, copy-pasteable R snippet (base R + corrplot) matching the
#' current settings, for the "Show me the code" educational feature.
#' @return A single character string of R code.
#' @export
code_correlation <- function(phenotypes, method = "pearson", plot_style = "circle",
                             type = "full", order = "original",
                             palette = "Spectral", levels = 10,
                             sig = FALSE, sig_level = 0.05,
                             split = FALSE, by = NULL,
                             upper_level = NULL, lower_level = NULL) {
  head <- c(
    "# --- Correlation plot (reproduces the MVApp Correlations tab) ---",
    "library(corrplot)",
    'data <- read.csv("your_data.csv", check.names = FALSE)  # your long-format data',
    paste0("traits <- ", .rvec(phenotypes)),
    .code_palette_line(palette, levels))

  if (isTRUE(split)) {
    body <- c(
      "",
      sprintf("# Split by %s: %s (upper triangle) vs %s (lower triangle)",
              by, upper_level, lower_level),
      sprintf('du <- data[data$%s == "%s", traits]; du <- du[complete.cases(du), ]',
              by, upper_level),
      sprintf('dl <- data[data$%s == "%s", traits]; dl <- dl[complete.cases(dl), ]',
              by, lower_level),
      sprintf('Mu <- cor(du, method = "%s"); Ml <- cor(dl, method = "%s")',
              method, method),
      sprintf('corrplot(Mu, type = "upper", method = "%s", tl.pos = "lt", col = pal)',
              plot_style),
      sprintf('corrplot(Ml, type = "lower", method = "%s", add = TRUE, tl.pos = "n", cl.pos = "n", col = pal)',
              plot_style))
  } else {
    call <- sprintf(
      'corrplot(M, method = "%s", type = "%s", order = "%s", col = pal%s)',
      plot_style, type, order,
      if (isTRUE(sig))
        sprintf(', p.mat = corrplot::cor.mtest(d)$p, sig.level = %s, insig = "pch"',
                sig_level) else "")
    body <- c(
      "",
      "d <- data[, traits]",
      "d <- d[complete.cases(d), ]  # drop rows with missing values",
      sprintf('M <- cor(d, method = "%s")', method),
      if (isTRUE(sig)) "# non-significant correlations are crossed out:" else NULL,
      call)
  }
  .rlines(c(head, body))
}

#' Generate R code that reproduces the scatterplot
#' @return A single character string of R code.
#' @export
code_scatter <- function(x, y, colour = NULL, lm = TRUE, ci = TRUE, rug = FALSE) {
  lines <- c(
    "# --- Scatterplot (reproduces the MVApp Scatterplots tab) ---",
    "library(ggplot2)",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    "",
    sprintf("p <- ggplot(data, aes(x = `%s`, y = `%s`%s)) +", x, y,
            if (!is.null(colour)) sprintf(", colour = `%s`", colour) else ""),
    "  geom_point(size = 2, alpha = 0.85) +",
    if (isTRUE(lm))
      sprintf('  geom_smooth(method = "lm", se = %s, formula = y ~ x) +',
              if (isTRUE(ci)) "TRUE" else "FALSE") else NULL,
    if (isTRUE(rug)) "  geom_rug(alpha = 0.4) +" else NULL,
    if (!is.null(colour)) "  ggsci::scale_colour_d3() +" else NULL,
    "  theme_minimal()",
    "p",
    "",
    if (!is.null(colour))
      sprintf('by(data, data$`%s`, function(df) summary(lm(`%s` ~ `%s`, df))$r.squared)  # R-squared per group',
              colour, y, x)
    else
      sprintf('summary(lm(`%s` ~ `%s`, data))$r.squared  # R-squared', y, x))
  .rlines(lines)
}

#' Generate R code that reproduces the PCA
#' @return A single character string of R code.
#' @export
code_pca <- function(phenotypes, scale = TRUE) {
  .rlines(c(
    "# --- PCA (reproduces the MVApp PCA tab) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    paste0("traits <- ", .rvec(phenotypes)),
    "X <- data[, traits]",
    "X <- X[complete.cases(X), ]  # drop rows with missing values",
    sprintf("pca <- prcomp(X, center = TRUE, scale. = %s)",
            if (isTRUE(scale)) "TRUE" else "FALSE"),
    "",
    "# Scree plot (variance explained per component)",
    "varexp <- 100 * pca$sdev^2 / sum(pca$sdev^2)",
    'barplot(varexp, names.arg = paste0("PC", seq_along(varexp)),',
    '        ylab = "% variance explained", main = "Scree plot")',
    "",
    "# Scores (samples) and loadings (variable contributions)",
    "scores   <- as.data.frame(pca$x)",
    "loadings <- as.data.frame(pca$rotation)"))
}

# ---------------------------------------------------------------- legends -----

#' Default figure legend for the correlation plot
#' @return A single character string.
#' @export
legend_correlation <- function(phenotypes, method = "pearson", n = NA,
                               split = FALSE, upper_level = NULL,
                               lower_level = NULL, sig = FALSE, sig_level = 0.05) {
  s <- sprintf("Correlation matrix of the traits %s, computed using %s's method%s.",
               paste(phenotypes, collapse = ", "),
               tools::toTitleCase(method),
               if (!is.na(n)) sprintf(" (n = %d complete observations)", n) else "")
  s <- paste0(s, " Each cell's colour and size encode the correlation ",
    "coefficient, which ranges from -1 to +1: one end of the diverging colour ",
    "scale marks strong negative correlations (traits moving in opposite ",
    "directions), the pale middle marks no correlation (near 0), and the other ",
    "end marks strong positive correlations (traits moving together), as shown ",
    "in the colour key.")
  if (isTRUE(split)) {
    s <- paste0(s, sprintf(
      " The upper triangle shows %s and the lower triangle shows %s.",
      upper_level, lower_level))
  }
  if (isTRUE(sig)) {
    s <- paste0(s, sprintf(
      " Correlations not significant at p < %s are crossed out.", sig_level))
  }
  s
}

#' Default figure legend for the scatterplot
#'
#' @param group_r2 Optional named numeric vector of R-squared per colour group
#'   (names = group levels). When supplied (and `colour` is set) the legend
#'   reports one R-squared per group instead of a single overall value.
#' @param group_p Optional named numeric vector of correlation p-values per
#'   colour group; when supplied alongside `group_r2` each group reports both.
#' @param p Optional overall correlation p-value (used when not grouped).
#' @return A single character string.
#' @export
legend_scatter <- function(x, y, colour = NULL, lm = TRUE, ci = TRUE, r2 = NULL,
                           group_r2 = NULL, group_p = NULL, p = NULL) {
  s <- sprintf("Scatterplot of %s (y-axis) against %s (x-axis)%s.", y, x,
               if (!is.null(colour)) sprintf(", coloured by %s", colour) else "")
  if (isTRUE(lm)) {
    s <- paste0(s, sprintf(" A linear regression line is shown%s.",
      if (isTRUE(ci)) " with its 95% confidence interval" else ""))
  }
  if (!is.null(colour) && !is.null(group_r2) && length(group_r2)) {
    grp <- names(group_r2)
    parts <- vapply(grp, function(g) {
      if (!is.null(group_p) && g %in% names(group_p)) {
        sprintf("%s = %.3f (p = %.3g)", g, group_r2[[g]], group_p[[g]])
      } else {
        sprintf("%s = %.3f", g, group_r2[[g]])
      }
    }, character(1))
    s <- paste0(s, sprintf(" Pearson R-squared (p-value) per %s: %s.",
                           colour, paste(parts, collapse = ", ")))
  } else if (!is.null(r2) && !is.na(r2)) {
    s <- paste0(s, sprintf(" Pearson R-squared = %.3f%s.", r2,
      if (!is.null(p) && !is.na(p)) sprintf(" (p = %.3g)", p) else ""))
  }
  s
}

#' Default figure legend for the PCA
#' @return A single character string.
#' @export
legend_pca <- function(phenotypes, n, scale = TRUE, var1 = NA, var2 = NA,
                       what = "overview", pcx = 1, pcy = 2, varx = NA, vary = NA,
                       colour = NULL) {
  base <- sprintf(
    "Principal component analysis of %d traits (%s) across %d samples, %s.",
    length(phenotypes), paste(phenotypes, collapse = ", "), n,
    if (isTRUE(scale)) "scaled to unit variance" else "on their original scale")
  switch(what,
    individuals = sprintf(paste0("%s Samples are projected onto PC%d (%.1f%%) ",
      "and PC%d (%.1f%%)%s; each point is one sample and the percentages are the ",
      "variance each component explains."),
      base, pcx, varx, pcy, vary,
      if (!is.null(colour)) sprintf(", coloured by %s", colour) else ""),
    scree = sprintf(paste0("%s Scree plot: bars and line show the %% of total ",
      "variance explained by each component in decreasing order; use it to decide ",
      "how many components to retain."), base),
    contributions = sprintf(paste0("%s Correlation circle: arrows show each ",
      "trait's loading on PC%d (%.1f%%) and PC%d (%.1f%%); traits with long arrows ",
      "pointing the same way co-vary, and arrow colour encodes each trait's ",
      "combined contribution to the two components."),
      base, pcx, varx, pcy, vary),
    eigenvalues = sprintf(paste0("%s Table of eigenvalues with the %% and ",
      "cumulative %% of total variance explained by each principal component."),
      base),
    sprintf("%s PC1 and PC2 explain %.1f%% and %.1f%% of the total variance, respectively.",
            base, var1, var2))
}

# ---------------------------------------------------- MDS & clustering --------

.scale_line <- function(scale) if (isTRUE(scale)) "X <- scale(X)" else NULL
.fmtp <- function(p) if (is.na(p)) "NA" else sprintf("%.3g", p)

# ------------------------------------------------ Curve fitting ---------------

#' @rdname code_pca
#' @export
code_curvefit <- function(time, response, samples, model = "linear") {
  fmla <- switch(model,
    linear      = sprintf("`%s` ~ `%s`", response, time),
    quadratic   = sprintf("`%s` ~ `%s` + I(`%s`^2)", response, time, time),
    exponential = sprintf("log(`%s`) ~ `%s`", response, time),
    square_root = sprintf("`%s` ~ sqrt(`%s`)", response, time),
    logarithmic = sprintf("`%s` ~ log(`%s`)", response, time),
    sprintf("`%s` ~ `%s`", response, time))
  .rlines(c(
    "# --- Curve fitting per sample/plant (reproduces the MVApp fitting tab) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    sprintf('data$sample <- interaction(%s, sep = " | ", drop = TRUE)',
      paste0("data$`", samples, "`", collapse = ", ")),
    "fits <- lapply(split(data, data$sample), function(s) {",
    sprintf("  m <- lm(%s, data = s)", fmla),
    sprintf("  pred <- %s",
      if (model == "exponential") "exp(fitted(m))" else "fitted(m)"),
    sprintf("  r2 <- cor(s$`%s`, pred)^2   # goodness-of-fit on the original scale", response),
    "  data.frame(sample = s$sample[1], INTERCEPT = coef(m)[1],",
    "             DELTA = coef(m)[2], r2 = r2)",
    "})",
    "result <- do.call(rbind, fits)"))
}

#' @rdname legend_pca
#' @export
legend_curvefit <- function(response, time, model, n_samples = NA,
                            median_r2 = NA, r2_cut = 0.7, n_below = NA) {
  mlab <- c(linear = "linear", quadratic = "quadratic",
    exponential = "exponential", square_root = "square-root",
    logarithmic = "logarithmic")[model]
  sprintf(paste0("A %s model of %s versus %s was fitted separately to each ",
    "individual sample / plant%s. Each fit yields an intercept and the primary ",
    "rate coefficient (DELTA), and the goodness-of-fit R-squared (the squared ",
    "correlation between observed and fitted values, on the original scale) ",
    "shows how well the model describes that plant%s. This table has one row per ",
    "plant. Plants with R-squared below %.2f%s are flagged as poorly fit and can ",
    "be excluded from downstream analysis."),
    mlab, response, time,
    if (!is.na(n_samples)) sprintf(" (%d samples)", n_samples) else "",
    if (!is.na(median_r2)) sprintf(" (median R-squared = %.2f)", median_r2) else "",
    r2_cut,
    if (!is.na(n_below)) sprintf(" - %d here", n_below) else "")
}

# ------------------------------------------------ Data Exploration ------------

#' @rdname code_pca
#' @export
code_normality <- function(trait, group = NULL) {
  .rlines(c(
    "# --- Normality: Shapiro-Wilk test + Q-Q plot ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    if (!is.null(group))
      sprintf('by(data$`%s`, data$`%s`, shapiro.test)  # test per group', trait, group)
    else sprintf('shapiro.test(data$`%s`)', trait),
    sprintf('qqnorm(data$`%s`); qqline(data$`%s`, col = "red")', trait, trait)))
}

#' @rdname code_pca
#' @export
code_variance <- function(trait, group, method = "levene") {
  test <- if (identical(method, "bartlett")) {
    sprintf('bartlett.test(`%s` ~ factor(`%s`), data = data)', trait, group)
  } else {
    c("# Levene's test (Brown-Forsythe, median-centred):",
      sprintf('z <- abs(data$`%s` - ave(data$`%s`, data$`%s`, FUN = median))', trait, trait, group),
      sprintf('summary(aov(z ~ factor(data$`%s`)))', group))
  }
  .rlines(c("# --- Homogeneity of variance ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)', test))
}

#' @rdname code_pca
#' @export
code_ttest <- function(trait, group, var_equal = FALSE, paired = FALSE,
                       pair_by = NULL) {
  if (isTRUE(paired) && !is.null(pair_by)) {
    return(.rlines(c(
      "# --- Paired t-test (average per pairing level, then pair across groups) ---",
      'data <- read.csv("your_data.csv", check.names = FALSE)',
      sprintf('m <- tapply(data$`%s`, list(data$`%s`, data$`%s`), mean)',
        trait, pair_by, group),
      "t.test(m[, 1], m[, 2], paired = TRUE)   # the two columns are the groups")))
  }
  .rlines(c("# --- Two-sample t-test ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    sprintf('t.test(`%s` ~ `%s`, data = data, var.equal = %s)',
      trait, group, if (var_equal) "TRUE" else "FALSE")))
}

#' @rdname code_pca
#' @export
code_anova2way <- function(trait, factor_a, factor_b) {
  .rlines(c("# --- Two-way ANOVA (with interaction) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    sprintf('fit <- aov(`%s` ~ `%s` * `%s`, data = data)', trait, factor_a, factor_b),
    "summary(fit)   # main effects + interaction",
    sprintf('interaction.plot(data$`%s`, data$`%s`, data$`%s`)',
      factor_a, factor_b, trait)))
}

#' @rdname code_pca
#' @export
code_kruskal <- function(trait, group) {
  .rlines(c(
    "# --- Non-parametric: Kruskal-Wallis + pairwise Wilcoxon ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    sprintf('kruskal.test(`%s` ~ `%s`, data = data)   # overall rank-based test', trait, group),
    sprintf('pairwise.wilcox.test(data$`%s`, data$`%s`, p.adjust.method = "BH")', trait, group),
    "# Compact-letter groups come from the pairwise p-value matrix",
    "# (e.g. multcompView::multcompLetters on the significant comparisons)."))
}

#' @rdname legend_pca
#' @export
legend_nonparam <- function(trait, group, kruskal_p = NA, n_groups = NA,
                            geom = "box") {
  gd <- switch(geom,
    box = "Box plots show the median, interquartile range and range",
    box_jitter = "Box plots (median, IQR, range) with the individual points overlaid",
    violin = "Violin plots show the distribution density",
    violin_jitter = "Violin plots with the individual points overlaid",
    bar_se = "Bars show the group mean +/- 1 standard error",
    bar_sd = "Bars show the group mean +/- 1 standard deviation",
    points_se = "Points show the individual measurements with the mean +/- 1 standard error",
    "Group distributions")
  sprintf(paste0("%s of %s across %s%s. Because the data need not be normally ",
    "distributed or have equal variances, groups were compared with the ",
    "rank-based Kruskal-Wallis test%s; post-hoc pairwise Wilcoxon (Mann-Whitney) ",
    "tests with BH adjustment give the letters - groups sharing a letter do not ",
    "differ significantly (p > 0.05). This is the non-parametric alternative to ",
    "one-way ANOVA + Tukey."),
    gd, trait, group, if (!is.na(n_groups)) sprintf(" (%d groups)", n_groups) else "",
    if (!is.na(kruskal_p)) sprintf(" (p = %.3g)", kruskal_p) else "")
}

#' @rdname legend_pca
#' @export
legend_normality <- function(trait, group = NULL) {
  sprintf(paste0("Distribution of %s%s. The histogram/density shows its shape ",
    "and the Q-Q plot compares the data quantiles against a normal distribution ",
    "(points on the line indicate normality). The Shapiro-Wilk test is reported%s; ",
    "a small p-value (< 0.05) indicates a significant departure from normality."),
    trait, if (!is.null(group)) sprintf(" per %s", group) else "",
    if (!is.null(group)) " for each group" else "")
}

#' @rdname legend_pca
#' @export
legend_variance <- function(trait, group, method = "Levene", p = NA) {
  sprintf(paste0("Homogeneity of variance of %s across %s, tested with %s's ",
    "test%s. A small p-value (< 0.05) indicates the group variances differ, ",
    "violating the equal-variance assumption of ANOVA and Student's t-test ",
    "(prefer Welch's t-test or a non-parametric test in that case)."),
    trait, group, method, if (!is.na(p)) sprintf(" (p = %.3g)", p) else "")
}

#' @rdname legend_pca
#' @export
legend_ttest <- function(trait, group, p = NA, method = "Welch two-sample t-test") {
  sprintf(paste0("Comparison of %s between the two levels of %s using a %s%s. A ",
    "small p-value (< 0.05) indicates the group means differ significantly."),
    trait, group, method, if (!is.na(p)) sprintf(" (p = %.3g)", p) else "")
}

#' @rdname legend_pca
#' @export
legend_anova2way <- function(trait, factor_a, factor_b, p_a = NA, p_b = NA,
                             p_int = NA) {
  sprintf(paste0("Two-way ANOVA of %s testing the main effects of %s and %s and ",
    "their interaction. The interaction plot shows the mean %s for each ",
    "combination; roughly parallel lines suggest no interaction, crossing or ",
    "diverging lines suggest one. Main effects: %s p = %s, %s p = %s; ",
    "interaction p = %s."),
    trait, factor_a, factor_b, trait, factor_a, .fmtp(p_a), factor_b,
    .fmtp(p_b), .fmtp(p_int))
}

#' @rdname code_pca
#' @export
code_mds <- function(phenotypes, scale = TRUE, dist_method = "euclidean") {
  .rlines(c(
    "# --- Classical MDS (reproduces the MVApp MDS tab) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    paste0("traits <- ", .rvec(phenotypes)),
    "X <- data[, traits]; X <- X[complete.cases(X), ]",
    .scale_line(scale),
    sprintf('D <- dist(X, method = "%s")', dist_method),
    "fit <- cmdscale(D, k = 2, eig = TRUE)",
    "mds <- as.data.frame(fit$points); names(mds) <- c(\"MDS1\", \"MDS2\")",
    "plot(mds$MDS1, mds$MDS2, xlab = \"MDS1\", ylab = \"MDS2\", pch = 19)"))
}

#' @rdname code_pca
#' @export
code_hclust <- function(phenotypes, scale = TRUE, dist_method = "euclidean",
                        link_method = "complete", k = 3) {
  .rlines(c(
    "# --- Hierarchical clustering (reproduces the MVApp clustering tab) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    paste0("traits <- ", .rvec(phenotypes)),
    "X <- data[, traits]; X <- X[complete.cases(X), ]",
    .scale_line(scale),
    sprintf('hc <- hclust(dist(X, method = "%s"), method = "%s")',
            dist_method, link_method),
    "plot(as.dendrogram(hc), horiz = TRUE)",
    sprintf("clusters <- cutree(hc, k = %d)", as.integer(k))))
}

#' @rdname code_pca
#' @export
code_kmeans <- function(phenotypes, scale = TRUE, k = 3, nstart = 25) {
  .rlines(c(
    "# --- K-means clustering (reproduces the MVApp clustering tab) ---",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    paste0("traits <- ", .rvec(phenotypes)),
    "X <- data[, traits]; X <- X[complete.cases(X), ]",
    .scale_line(scale),
    "set.seed(1)  # for reproducible cluster assignment",
    sprintf("km <- kmeans(X, centers = %d, nstart = %d)",
            as.integer(k), as.integer(nstart)),
    "clusters <- km$cluster",
    "# Elbow curve:",
    "wss <- sapply(1:10, function(kk) kmeans(X, kk, nstart = 25)$tot.withinss)",
    'plot(1:10, wss, type = "b", xlab = "k", ylab = "within-cluster SS")'))
}

#' @rdname legend_pca
#' @export
legend_mds <- function(phenotypes, n, scale = TRUE, dist_method = "euclidean",
                       gof = NA) {
  sprintf(paste0("Classical multidimensional scaling of %d samples based on %d ",
    "traits (%s), using %s distances%s. Points that are close together are ",
    "similar across the traits.%s"),
    n, length(phenotypes), paste(phenotypes, collapse = ", "), dist_method,
    if (isTRUE(scale)) " on scaled traits" else "",
    if (!is.na(gof)) sprintf(paste0(" Goodness-of-fit = %.3f: the proportion of ",
      "the original (full-dimensional) distance information preserved in this ",
      "2-D map - higher is better, with values around 0.8 or above indicating ",
      "the map faithfully represents the true distances between samples."),
      gof) else "")
}

#' Figure legend explaining a k-selection method (elbow / silhouette / gap)
#'
#' @param method One of "wss" (elbow), "silhouette", "gap".
#' @param best Optional suggested k to mention.
#' @param B Bootstrap samples (gap only), for the description.
#' @return A single character string.
#' @export
legend_optimal_k <- function(method, best = NA, B = 25) {
  sugg <- if (!is.na(best))
    sprintf(" The dashed line marks the suggested k = %s.", best) else ""
  switch(method,
    wss = paste0("Elbow plot: the total within-cluster sum of squares (a measure ",
      "of how tight the clusters are) as the number of clusters k increases. It ",
      "always decreases with more clusters, so choose the k at the 'elbow' - the ",
      "point where adding another cluster gives only a small further reduction. ",
      "The dashed line marks your currently selected k."),
    silhouette = paste0("Average silhouette width for each k. For every sample ",
      "the silhouette measures how well it fits its own cluster versus the ",
      "nearest neighbouring cluster (range -1 to 1; higher is better). The k with ",
      "the highest average silhouette gives the best-separated clustering.", sugg),
    gap = sprintf(paste0("Gap statistic: compares the within-cluster dispersion ",
      "to that expected under a null reference (no real clusters), estimated from ",
      "%d bootstrap samples. Larger gaps indicate more convincing clustering; the ",
      "suggested k uses the 'firstSEmax' rule - the smallest k whose gap is within ",
      "one standard error of the first peak.%s"), B, sugg),
    "Diagnostic for choosing the number of clusters.")
}

#' @rdname legend_pca
#' @export
legend_hclust <- function(phenotypes, n, k, link_method = "complete",
                          dist_method = "euclidean", scale = TRUE) {
  sprintf(paste0("Hierarchical clustering of %d samples based on %d traits (%s), ",
    "using %s distances and %s linkage%s. The tree was cut into %d clusters."),
    n, length(phenotypes), paste(phenotypes, collapse = ", "),
    dist_method, link_method,
    if (isTRUE(scale)) " on scaled traits" else "", as.integer(k))
}

#' @rdname code_pca
#' @export
code_anova_tukey <- function(trait, group) {
  .rlines(c(
    "# --- ANOVA + Tukey HSD with significance letters (Julkowska-lab style) ---",
    "library(ggplot2); library(multcompView)",
    'data <- read.csv("your_data.csv", check.names = FALSE)',
    sprintf('data$%s <- factor(data$%s)', group, group),
    "",
    sprintf("fit <- aov(`%s` ~ `%s`, data = data)", trait, group),
    "print(summary(fit))                       # overall ANOVA",
    "tuk <- TukeyHSD(fit)                       # pairwise post-hoc",
    "cld <- multcompLetters4(fit, tuk)          # letters ordered by mean",
    "",
    "# per-group summary + letter placement",
    sprintf("lab <- aggregate(`%s` ~ `%s`, data, function(v) max(v))", trait, group),
    sprintf('lab$letter <- cld$`%s`$Letters[as.character(lab$`%s`)]', group, group),
    "",
    sprintf("ggplot(data, aes(`%s`, `%s`, fill = `%s`)) +", group, trait, group),
    "  geom_boxplot(alpha = 0.8) +",
    sprintf("  geom_text(data = lab, aes(y = `%s`, label = letter), vjust = -0.5) +", trait),
    "  ggsci::scale_fill_d3() + theme_minimal()"))
}

#' @rdname legend_pca
#' @param geom The plot geom shown: one of "box", "box_jitter", "violin",
#'   "violin_jitter", "bar_se", "bar_sd", "points_se".
#' @param points_are What an individual plotted point represents (e.g.
#'   "individual measurements", or "the mean per genotype").
#' @export
legend_anova_tukey <- function(trait, group, anova_p = NA, n_groups = NA,
                               geom = "box", points_are = "individual measurements") {
  desc <- switch(geom,
    box = "Box plots show the median, interquartile range and full range",
    box_jitter = sprintf(paste0("Box plots (median, interquartile range, range) ",
      "with the underlying %s overlaid as points"), points_are),
    violin = "Violin plots show the distribution density (mirrored)",
    violin_jitter = sprintf(paste0("Violin plots (distribution density) with the ",
      "underlying %s overlaid as points"), points_are),
    bar_se = "Bars show the group mean with error bars of +/- 1 standard error",
    bar_sd = "Bars show the group mean with error bars of +/- 1 standard deviation",
    points_se = sprintf(paste0("Points show the %s, with the group mean +/- 1 ",
      "standard error overlaid in black"), points_are),
    "Group distributions")
  sprintf(paste0("%s of %s across %s%s. Groups were compared by one-way ANOVA%s ",
    "followed by Tukey's HSD; groups sharing a letter are not significantly ",
    "different (p > 0.05)."),
    desc, trait, group,
    if (!is.na(n_groups)) sprintf(" (%d groups)", n_groups) else "",
    if (!is.na(anova_p)) sprintf(" (p = %.3g)", anova_p) else "")
}

#' @rdname legend_pca
#' @export
legend_kmeans <- function(phenotypes, n, k, scale = TRUE, betweenss_ratio = NA) {
  sprintf(paste0("K-means clustering of %d samples into %d clusters based on %d ",
    "traits (%s)%s.%s"),
    n, as.integer(k), length(phenotypes), paste(phenotypes, collapse = ", "),
    if (isTRUE(scale)) " (scaled to unit variance)" else "",
    if (!is.na(betweenss_ratio))
      sprintf(" Between-cluster SS accounts for %.1f%% of the total.",
              100 * betweenss_ratio) else "")
}
