make_toy <- function(n = 30) {
  set.seed(11)
  data.frame(
    id = paste0("s", seq_len(n)),
    grp = rep(c("A", "B"), length.out = n),
    t1 = rnorm(n), t2 = rnorm(n), t3 = rnorm(n), t4 = rnorm(n))
}

test_that("run_mds embeds samples and reports GOF", {
  df <- make_toy()
  res <- run_mds(df, c("t1", "t2", "t3", "t4"), id_cols = c("id", "grp"), k = 2)
  expect_s3_class(res, "mvapp_mds")
  expect_equal(nrow(res$points), 30)
  expect_true(all(c("MDS1", "MDS2") %in% names(res$points)))
  expect_equal(res$points$id[1], "s1_A")
  expect_true(res$gof[1] >= 0 && res$gof[1] <= 1)
})

test_that("run_mds on scaled Euclidean matches PCA scores up to sign", {
  df <- make_toy()
  ph <- c("t1", "t2", "t3", "t4")
  mds <- run_mds(df, ph, scale = TRUE, k = 2)
  pca <- run_pca(df, ph, scale = TRUE)
  # classical MDS of Euclidean distances == PCA; compare absolute coordinates
  expect_equal(unname(abs(as.matrix(mds$points[, c("MDS1", "MDS2")]))),
               unname(abs(as.matrix(pca$ind[, c("Dim.1", "Dim.2")]))),
               tolerance = 1e-6)
})

test_that("run_hclust cuts into k clusters", {
  df <- make_toy()
  res <- run_hclust(df, c("t1", "t2", "t3", "t4"), id_cols = c("id"),
                    k = 3, link_method = "ward.D2")
  expect_s3_class(res, "mvapp_hclust")
  expect_s3_class(res$hclust, "hclust")
  expect_equal(length(unique(res$clusters$cluster)), 3)
  expect_equal(nrow(res$clusters), 30)
})

test_that("run_kmeans is reproducible with a seed and returns an elbow curve", {
  df <- make_toy()
  ph <- c("t1", "t2", "t3", "t4")
  a <- run_kmeans(df, ph, k = 3, seed = 42)
  b <- run_kmeans(df, ph, k = 3, seed = 42)
  expect_equal(a$clusters$cluster, b$clusters$cluster)     # deterministic
  expect_equal(sort(unique(a$clusters$cluster)), 1:3)
  expect_equal(nrow(a$wss), min(10, nrow(df) - 1))
  expect_true(all(diff(a$wss$tot_withinss) <= 1e-6))        # WSS non-increasing
})

test_that("ordination inputs are validated", {
  df <- make_toy()
  expect_error(run_mds(df, c("t1", "zzz")), "not found")
  expect_error(run_hclust(df, "t1"), "at least two")
})
