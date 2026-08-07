test_that("run_pca returns coherent eigen/contribution structure", {
  set.seed(1)
  df <- data.frame(
    id  = paste0("s", 1:30),
    grp = rep(c("A", "B"), each = 15),
    t1  = rnorm(30),
    t2  = rnorm(30),
    t3  = rnorm(30),
    t4  = rnorm(30)
  )
  res <- run_pca(df, phenotypes = c("t1", "t2", "t3", "t4"),
                 id_cols = c("id", "grp"))

  expect_s3_class(res, "mvapp_pca")
  # variance percentages sum to 100
  expect_equal(sum(res$eig$variance_percent), 100)
  # scaled PCA: eigenvalues sum to the number of variables
  expect_equal(sum(res$eig$eigenvalue), 4)
  # contributions sum to 100 within each principal component
  contrib_cols <- res$contrib[, -1]
  expect_equal(unname(colSums(contrib_cols)), rep(100, ncol(contrib_cols)))
  # one row of scores per (complete) individual
  expect_equal(nrow(res$ind), 30)
  expect_equal(res$ind$id[1], "s1_A")
})

test_that("run_pca drops rows with NAs by default and can fail instead", {
  df <- data.frame(a = c(1, 2, 3, NA, 5, 6),
                   b = c(2, 1, 4, 3, 6, 5),
                   c = c(1, 3, 2, 4, 5, 7))
  res <- run_pca(df, phenotypes = c("a", "b", "c"))
  expect_equal(res$n, 5)
  expect_error(run_pca(df, phenotypes = c("a", "b", "c"), na_action = "fail"),
               "contain NAs")
})

test_that("run_pca validates inputs", {
  df <- data.frame(a = 1:5, b = letters[1:5], c = 6:10)
  expect_error(run_pca(df, phenotypes = c("a", "z")), "not found")
  expect_error(run_pca(df, phenotypes = "a"), "at least two")
  expect_error(run_pca(df, phenotypes = c("a", "b")), "non-numeric")
})

test_that("run_pca scores match a direct prcomp on the same matrix", {
  df <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8),
                   y = c(2, 1, 4, 3, 6, 5, 8, 7),
                   z = c(1, 4, 9, 16, 25, 36, 49, 64))
  res <- run_pca(df, phenotypes = c("x", "y", "z"), scale = TRUE)
  ref <- stats::prcomp(as.matrix(df), center = TRUE, scale. = TRUE)
  expect_equal(unname(as.matrix(res$ind[, -1])), unname(ref$x))
})

test_that("run_pca accepts an mvapp_dataset and uses its roles for labels", {
  df <- data.frame(LINE = paste0("L", 1:6),
                   TREATMENT = rep(c("Control", "Drought"), 3),
                   AREA = c(1, 2, 3, 4, 5, 6),
                   PERIMETER = c(6, 5, 4, 3, 2, 1),
                   ROUNDNESS = c(1, 2, 1, 2, 1, 2))
  ds <- mvapp_dataset(df, id = "LINE", group = "TREATMENT",
                      phenotypes = c("AREA", "PERIMETER", "ROUNDNESS"))
  res <- run_pca(ds, phenotypes = c("AREA", "PERIMETER", "ROUNDNESS"))
  expect_s3_class(res, "mvapp_pca")
  expect_true(grepl("^L1_Control", res$ind$id[1]))
})
