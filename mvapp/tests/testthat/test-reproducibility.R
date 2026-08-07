# The strongest test for generated code is that it PARSES as valid R.

test_that("code_correlation generates parseable R for all modes", {
  ph <- c("AREA", "BIOMASS", "ROUNDNESS")
  plain <- code_correlation(ph, method = "spearman", plot_style = "ellipse")
  expect_silent(parse(text = plain))
  expect_match(plain, "corrplot")
  expect_match(plain, "AREA")

  sigc <- code_correlation(ph, sig = TRUE, sig_level = 0.01)
  expect_silent(parse(text = sigc))
  expect_match(sigc, "cor.mtest", fixed = TRUE)

  splitc <- code_correlation(ph, split = TRUE, by = "TREATMENT",
                             upper_level = "Control", lower_level = "Drought")
  expect_silent(parse(text = splitc))
  expect_match(splitc, "Control", fixed = TRUE)
  expect_match(splitc, "type = \"lower\"", fixed = TRUE)

  ggc <- code_correlation(ph, palette = "ggsci: GSEA")
  expect_silent(parse(text = ggc))
  expect_match(ggc, "pal_gsea", fixed = TRUE)
})

test_that("code_scatter generates parseable R with and without colour", {
  a <- code_scatter("AREA", "BIOMASS", colour = "TREATMENT", lm = TRUE, ci = TRUE)
  expect_silent(parse(text = a))
  expect_match(a, "geom_smooth", fixed = TRUE)
  expect_match(a, "r.squared", fixed = TRUE)

  b <- code_scatter("AREA", "BIOMASS", colour = NULL, lm = FALSE, rug = TRUE)
  expect_silent(parse(text = b))
  expect_false(grepl("geom_smooth", b))
  expect_match(b, "geom_rug", fixed = TRUE)
})

test_that("code_pca generates parseable R", {
  p <- code_pca(c("AREA", "BIOMASS", "FvFm"), scale = TRUE)
  expect_silent(parse(text = p))
  expect_match(p, "prcomp", fixed = TRUE)
  expect_match(p, "scale. = TRUE", fixed = TRUE)
})

test_that("legends mention the key elements", {
  lc <- legend_correlation(c("AREA", "BIOMASS"), method = "pearson", n = 80,
                           split = TRUE, upper_level = "Control",
                           lower_level = "Drought", sig = TRUE, sig_level = 0.05)
  expect_match(lc, "Pearson")
  expect_match(lc, "n = 80", fixed = TRUE)
  expect_match(lc, "upper triangle shows Control")
  expect_match(lc, "crossed out")
  expect_match(lc, "strong negative correlations")
  expect_match(lc, "strong positive correlations")

  ls <- legend_scatter("AREA", "BIOMASS", colour = "TREATMENT", r2 = 0.973)
  expect_match(ls, "coloured by TREATMENT")
  expect_match(ls, "R-squared = 0.973", fixed = TRUE)

  # per-group R-squared (with p-values) overrides the single overall value
  lg <- legend_scatter("AREA", "BIOMASS", colour = "TREATMENT", r2 = 0.9,
                       group_r2 = c(Control = 0.955, Drought = 0.975),
                       group_p = c(Control = 1e-27, Drought = 4e-32))
  expect_match(lg, "per TREATMENT")
  expect_match(lg, "Control = 0.955 (p = 1e-27)", fixed = TRUE)
  expect_match(lg, "Drought = 0.975 (p = 4e-32)", fixed = TRUE)
  expect_false(grepl("R-squared = 0.900", lg))

  lp <- legend_pca(c("AREA", "BIOMASS", "FvFm"), n = 80, scale = TRUE,
                   var1 = 80.9, var2 = 12.3)
  expect_match(lp, "3 traits")
  expect_match(lp, "80.9%", fixed = TRUE)
  expect_match(lp, "scaled to unit variance")
})

test_that("legend_pca is specific per sub-plot and PC-aware", {
  ph <- c("AREA", "BIOMASS", "FvFm")
  ind <- legend_pca(ph, n = 80, what = "individuals", pcx = 2, pcy = 3,
                    varx = 12.3, vary = 4.1, colour = "TREATMENT")
  expect_match(ind, "PC2 (12.3%)", fixed = TRUE)
  expect_match(ind, "PC3 (4.1%)", fixed = TRUE)
  expect_match(ind, "coloured by TREATMENT")
  expect_match(legend_pca(ph, n = 80, what = "scree"), "Scree")
  expect_match(legend_pca(ph, n = 80, what = "contributions", pcx = 1, pcy = 2,
                          varx = 80, vary = 12), "Correlation circle")
  expect_match(legend_pca(ph, n = 80, what = "eigenvalues"), "eigenvalues")
})

test_that("legend_mds explains goodness-of-fit", {
  lm <- legend_mds(c("a", "b", "c"), n = 40, dist_method = "euclidean",
                   gof = 0.82)
  expect_match(lm, "Goodness-of-fit = 0.820", fixed = TRUE)
  expect_match(lm, "distance information preserved")
  expect_match(lm, "close together are similar")
})

test_that("legend_optimal_k explains each k-selection method", {
  expect_match(legend_optimal_k("wss"), "elbow")
  expect_match(legend_optimal_k("silhouette", best = 2), "silhouette")
  expect_match(legend_optimal_k("silhouette", best = 2), "k = 2", fixed = TRUE)
  g <- legend_optimal_k("gap", best = 4, B = 25)
  expect_match(g, "Gap statistic")
  expect_match(g, "25 bootstrap", fixed = TRUE)
  expect_match(g, "firstSEmax")
})

test_that("legend_anova_tukey describes the chosen geom", {
  expect_match(legend_anova_tukey("AREA", "cluster", geom = "bar_se"),
               "standard error")
  expect_match(legend_anova_tukey("AREA", "cluster", geom = "bar_sd"),
               "standard deviation")
  expect_match(legend_anova_tukey("AREA", "cluster", geom = "violin"), "Violin")
  expect_match(legend_anova_tukey("AREA", "cluster", geom = "box_jitter",
               points_are = "the mean per genotype"), "mean per genotype")
})
