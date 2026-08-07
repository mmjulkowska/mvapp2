test_that("anova_tukey_letters flags real differences and shares letters", {
  set.seed(1)
  df <- data.frame(
    value = c(rnorm(12, 5), rnorm(12, 5.2), rnorm(12, 10)),
    cluster = rep(c("A", "B", "C"), each = 12))
  res <- anova_tukey_letters(df, "value", "cluster")
  expect_s3_class(res, "mvapp_anova_cld")
  expect_lt(res$anova_p, 1e-3)
  s <- res$summary
  expect_setequal(s$group, c("A", "B", "C"))
  # C is far higher -> unique letter; A and B overlap -> share a letter
  lC <- s$letter[s$group == "C"]
  lA <- s$letter[s$group == "A"]
  lB <- s$letter[s$group == "B"]
  expect_false(lC == lA)
  expect_equal(lA, lB)
  # summary ordered by decreasing mean, "a" on the top group
  expect_equal(s$group[1], "C")
  expect_match(s$letter[1], "a")
})

test_that("anova_tukey_letters validates inputs", {
  df <- data.frame(v = 1:6, g = rep(c("x", "y"), 3))
  expect_error(anova_tukey_letters(df, "nope", "g"), "not found")
  one <- data.frame(v = 1:4, g = rep("only", 4))
  expect_error(anova_tukey_letters(one, "v", "g"), "2 groups")
})

test_that("cluster_distances is symmetric with sensible diagonal", {
  df <- data.frame(id = paste0("s", 1:20), grp = "x",
                   a = c(rnorm(10, 0), rnorm(10, 10)),
                   b = c(rnorm(10, 0), rnorm(10, 10)),
                   c = c(rnorm(10, 0), rnorm(10, 10)))
  hres <- run_hclust(df, c("a", "b", "c"), id_cols = "id", k = 2)
  M <- cluster_distances(hres)
  expect_equal(dim(M), c(2, 2))
  expect_true(isSymmetric(M))
  # two well-separated clusters: between-distance > each within-distance
  expect_gt(M[1, 2], M[1, 1])
  expect_gt(M[1, 2], M[2, 2])
})

test_that("cluster_diagnostics returns elbow/silhouette/gap and suggestions", {
  set.seed(2)
  df <- data.frame(id = paste0("s", 1:40),
                   a = c(rnorm(20, 0), rnorm(20, 8)),
                   b = c(rnorm(20, 0), rnorm(20, 8)),
                   c = c(rnorm(20, 0), rnorm(20, 8)),
                   d = rnorm(40))
  diag <- cluster_diagnostics(df, c("a", "b", "c", "d"), id_cols = "id",
                              max_k = 6, B = 20)
  tab <- diag$table
  expect_true(all(c("k", "wss", "silhouette", "gap") %in% names(tab)))
  expect_true(all(diff(tab$wss) <= 1e-6))               # WSS non-increasing
  expect_true(all(tab$silhouette[-1] >= -1 & tab$silhouette[-1] <= 1))
  expect_equal(diag$best$silhouette, 2)                 # 2 true clusters
})
