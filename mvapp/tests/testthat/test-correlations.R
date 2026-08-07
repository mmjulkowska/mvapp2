test_that("cor_matrix reproduces stats::cor on the phenotype columns", {
  set.seed(3)
  df <- data.frame(g = rep(c("A", "B"), 10),
                   a = rnorm(20), b = rnorm(20), c = rnorm(20))
  M <- cor_matrix(df, c("a", "b", "c"), method = "pearson")
  expect_equal(unname(M), unname(cor(df[c("a", "b", "c")], method = "pearson")))
  expect_equal(diag(M), c(a = 1, b = 1, c = 1))
  expect_equal(dim(M), c(3, 3))
})

test_that("cor_matrix supports spearman and mvapp_dataset input", {
  df <- data.frame(LINE = paste0("l", 1:6), TREATMENT = rep(c("C", "D"), 3),
                   x = c(1, 2, 3, 4, 5, 6), y = c(2, 1, 4, 3, 6, 5),
                   z = c(6, 5, 4, 3, 2, 1))
  ds <- mvapp_dataset(df, id = "LINE", group = "TREATMENT",
                      phenotypes = c("x", "y", "z"))
  M <- cor_matrix(ds, method = "spearman")           # phenotypes from roles
  expect_equal(dim(M), c(3, 3))
  expect_equal(M["x", "z"], -1)                       # perfectly anti-monotonic
})

test_that("cor_pmatrix is symmetric, zero-diagonal, and flags strong pairs", {
  set.seed(7)
  x <- rnorm(40); df <- data.frame(x = x, y = x + rnorm(40, 0, 0.05),
                                   z = rnorm(40))
  P <- cor_pmatrix(df, c("x", "y", "z"))
  expect_true(isSymmetric(P))
  expect_equal(diag(P), c(x = 0, y = 0, z = 0))
  expect_lt(P["x", "y"], 0.001)                       # x,y near-identical
  # p-value matches a direct cor.test
  expect_equal(P["x", "z"], cor.test(df$x, df$z)$p.value)
})

test_that("cor_n counts complete cases", {
  df <- data.frame(a = c(1, 2, NA, 4, 5), b = c(1, 2, 3, NA, 5),
                   c = c(1, 2, 3, 4, 5))
  expect_equal(cor_n(df, c("a", "b", "c")), 3)
})

test_that("correlation inputs are validated", {
  df <- data.frame(a = 1:5, b = letters[1:5], c = 6:10)
  expect_error(cor_matrix(df, c("a", "zzz")), "not found")
  expect_error(cor_matrix(df, c("a", "b")), "non-numeric")
  expect_error(cor_matrix(df, "a"), "at least two")
})
