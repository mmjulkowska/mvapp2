test_that("shapiro_by_group flags non-normal and normal data", {
  set.seed(1)
  df <- data.frame(
    g = rep(c("norm", "skew"), each = 60),
    v = c(rnorm(60), rexp(60)))
  res <- shapiro_by_group(df, "v", "g")
  expect_setequal(res$group, c("norm", "skew"))
  pn <- res$p[res$group == "norm"]; ps <- res$p[res$group == "skew"]
  expect_gt(pn, 0.05)     # normal -> not rejected
  expect_lt(ps, 0.05)     # exponential -> rejected
  # overall (no group)
  o <- shapiro_by_group(df, "v")
  expect_equal(nrow(o), 1); expect_equal(o$group, "all")
})

test_that("variance_test detects unequal variances (Levene and Bartlett)", {
  set.seed(2)
  df <- data.frame(g = rep(c("a", "b"), each = 50),
                   v = c(rnorm(50, 0, 1), rnorm(50, 0, 4)))
  lev <- variance_test(df, "v", "g", method = "levene")
  bar <- variance_test(df, "v", "g", method = "bartlett")
  expect_match(lev$method, "Levene")
  expect_lt(lev$p, 0.05)
  expect_lt(bar$p, 0.05)
  # equal variances -> not rejected
  df2 <- data.frame(g = rep(c("a", "b"), each = 50), v = rnorm(100))
  expect_gt(variance_test(df2, "v", "g")$p, 0.05)
})

test_that("t_test_two matches stats::t.test and validates group count", {
  set.seed(3)
  df <- data.frame(g = rep(c("ctrl", "trt"), each = 30),
                   v = c(rnorm(30, 5), rnorm(30, 7)))
  res <- t_test_two(df, "v", "g")
  ref <- t.test(v ~ g, data = df)
  expect_equal(res$p, ref$p.value)
  expect_equal(res$statistic, unname(ref$statistic))
  expect_lt(res$p, 0.05)
  three <- data.frame(g = rep(c("a", "b", "c"), each = 5), v = rnorm(15))
  expect_error(t_test_two(three, "v", "g"), "exactly 2 groups")
})

test_that("t_test_two paired-by-variable averages per cell and pairs", {
  set.seed(9)
  geno <- paste0("g", 1:12)
  d <- rbind(
    data.frame(ACC = rep(geno, 3), TRT = "Control", v = rep(1:12, 3) + rnorm(36, 0, 0.2)),
    data.frame(ACC = rep(geno, 3), TRT = "Salt",    v = rep(1:12, 3) + 2 + rnorm(36, 0, 0.2)))
  res <- t_test_two(d, "v", "TRT", paired = TRUE, pair_by = "ACC")
  expect_equal(res$n_pairs, 12)                 # one pair per genotype
  expect_match(res$method, "Paired")
  expect_lt(res$p, 1e-6)                         # consistent +2 shift per genotype
  # matches a manual paired test on the per-genotype means
  m <- tapply(d$v, list(d$ACC, d$TRT), mean)
  expect_equal(res$p, t.test(m[, "Control"], m[, "Salt"], paired = TRUE)$p.value)
  expect_error(t_test_two(d, "v", "TRT", paired = TRUE), "pair_by")
})

test_that("anova_two_way reports main effects and interaction", {
  set.seed(4)
  d <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), rep = 1:20)
  d$v <- 2 * (d$A == "a2") + 3 * (d$B == "b2") + rnorm(nrow(d), 0, 0.5)
  res <- anova_two_way(d, "v", "A", "B")
  expect_true(all(c("A", "B", "A:B", "Residuals") %in% res$table$term))
  pA <- res$table$p[res$table$term == "A"]
  pB <- res$table$p[res$table$term == "B"]
  expect_lt(pA, 0.05); expect_lt(pB, 0.05)          # both main effects real
  expect_equal(nrow(res$means), 4)                   # 2x2 cell means
  expect_error(anova_two_way(d, "v", "A", "A"), "different factors")
})

test_that("anova_screen ranks traits by one-way ANOVA p-value", {
  set.seed(5)
  df <- data.frame(g = rep(c("a", "b", "c"), each = 20),
    strong = c(rnorm(20, 0), rnorm(20, 5), rnorm(20, 10)),  # big group effect
    none = rnorm(60))                                        # no group effect
  sc <- anova_screen(df, c("strong", "none"), "g")
  expect_equal(sc$trait[1], "strong")          # most significant first
  expect_lt(sc$anova_p[sc$trait == "strong"], 1e-6)
  expect_gt(sc$anova_p[sc$trait == "none"], 0.05)
})

test_that("kruskal_letters is the non-parametric analogue with letters", {
  set.seed(8)
  df <- data.frame(
    g = rep(c("lo", "mid", "hi"), each = 25),
    v = c(rexp(25, 1), rexp(25, 1) + 3, rexp(25, 1) + 8))   # skewed, clearly ordered
  res <- kruskal_letters(df, "v", "g")
  expect_s3_class(res, "mvapp_anova_cld")
  expect_lt(res$anova_p, 1e-3)
  s <- res$summary
  expect_true("median" %in% names(s))
  expect_equal(s$group[1], "hi")            # ordered by descending median
  expect_match(s$letter[1], "a")
  expect_false(s$letter[s$group == "hi"] == s$letter[s$group == "lo"])
})

test_that("kruskal_screen ranks traits and faceted gives per-facet columns", {
  set.seed(13)
  d <- expand.grid(g = c("a", "b", "c"), day = c("d1", "d2"), r = 1:15)
  d$hit <- ifelse(d$day == "d2", 4 * as.integer(d$g), 0) + rexp(nrow(d))
  d$flat <- rexp(nrow(d))
  sc <- kruskal_screen(d, c("hit", "flat"), "g")
  expect_equal(sc$trait[1], "hit")
  fc <- kruskal_screen_faceted(d, c("hit", "flat"), "g", "day")
  expect_true(all(c("trait", "d1", "d2") %in% names(fc)))
  expect_gt(fc$d1[fc$trait == "hit"], 0.05)
  expect_lt(fc$d2[fc$trait == "hit"], 0.01)
})

test_that("ttest_screen ranks traits by two-group t-test p-value", {
  set.seed(10)
  df <- data.frame(g = rep(c("a", "b"), each = 30),
    big = c(rnorm(30, 0), rnorm(30, 4)), flat = rnorm(60))
  sc <- ttest_screen(df, c("big", "flat"), "g")
  expect_equal(sc$trait[1], "big")
  expect_lt(sc$p[sc$trait == "big"], 1e-6)
  expect_gt(sc$p[sc$trait == "flat"], 0.05)
})

test_that("anova_screen_faceted gives one p-column per facet level", {
  set.seed(11)
  d <- expand.grid(g = c("a", "b", "c"), day = c("d1", "d2"), r = 1:15)
  # trait differs by group only on d2
  d$tr <- ifelse(d$day == "d2", 5 * as.integer(d$g), 0) + rnorm(nrow(d), 0, 0.3)
  d$flat <- rnorm(nrow(d))
  sc <- anova_screen_faceted(d, c("tr", "flat"), "g", "day")
  expect_true(all(c("trait", "d1", "d2") %in% names(sc)))
  expect_gt(sc$d1[sc$trait == "tr"], 0.05)   # no group effect on day 1
  expect_lt(sc$d2[sc$trait == "tr"], 1e-6)   # strong effect on day 2
})

test_that("ttest_screen_faceted gives one p-column per facet level", {
  set.seed(12)
  d <- expand.grid(g = c("ctrl", "trt"), day = c("d1", "d2"), r = 1:20)
  d$tr <- ifelse(d$g == "trt" & d$day == "d2", 3, 0) + rnorm(nrow(d), 0, 0.4)
  sc <- ttest_screen_faceted(d, "tr", "g", "day")
  expect_true(all(c("trait", "d1", "d2") %in% names(sc)))
  expect_gt(sc$d1[1], 0.05)
  expect_lt(sc$d2[1], 1e-6)
})

test_that("anova2_screen reports per-term p-values across traits", {
  set.seed(6)
  d <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), r = 1:15)
  d$hit <- 3 * (d$A == "a2") + rnorm(nrow(d), 0, 0.5)   # A effect only
  d$flat <- rnorm(nrow(d))
  sc <- anova2_screen(d, c("hit", "flat"), "A", "B")
  expect_true(all(c("p_A", "p_B", "p_interaction") %in% names(sc)))
  expect_lt(sc$p_A[sc$trait == "hit"], 0.05)
  expect_gt(sc$p_A[sc$trait == "flat"], 0.05)
})
