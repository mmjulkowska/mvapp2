test_that("cor_split returns aligned per-level correlation matrices", {
  df <- data.frame(
    trt = rep(c("Control", "Drought"), each = 20),
    a = c(1:20, 20:1) + rnorm(40, 0, 0.001),
    b = c(1:20, 1:20) + rnorm(40, 0, 0.001),
    c = rnorm(40)
  )
  sp <- cor_split(df, phenotypes = c("a", "b", "c"), by = "trt",
                  upper_level = "Control", lower_level = "Drought")
  expect_equal(dim(sp$upper), c(3, 3))
  expect_equal(dim(sp$lower), c(3, 3))
  expect_identical(rownames(sp$upper), rownames(sp$lower))  # aligned
  expect_equal(sp$upper_level, "Control")
  expect_equal(sp$n_upper, 20)
  # upper matrix equals cor_matrix on the Control subset
  ctrl <- df[df$trt == "Control", ]
  expect_equal(sp$upper, cor_matrix(ctrl, c("a", "b", "c")))
  # a,b move together under Control, oppositely under Drought
  expect_gt(sp$upper["a", "b"], 0.9)
  expect_lt(sp$lower["a", "b"], -0.9)
  # per-level p-value matrices are returned, aligned, and significant for a~b
  expect_equal(dim(sp$p_upper), c(3, 3))
  expect_identical(rownames(sp$p_upper), rownames(sp$upper))
  expect_true(isSymmetric(sp$p_lower))
  expect_lt(sp$p_upper["a", "b"], 0.001)
  expect_lt(sp$p_lower["a", "b"], 0.001)
})

test_that("cor_split defaults to the first two levels and validates", {
  df <- data.frame(trt = rep(c("C", "D"), each = 5),
                   x = rnorm(10), y = rnorm(10), z = rnorm(10))
  sp <- cor_split(df, c("x", "y", "z"), by = "trt")
  expect_setequal(c(sp$upper_level, sp$lower_level), c("C", "D"))
  expect_error(cor_split(df, c("x", "y", "z"), by = "nope"), "not found")
  expect_error(cor_split(df, c("x", "y", "z"), by = "trt",
                         upper_level = "C", lower_level = "C"), "must differ")
  expect_error(cor_split(df, c("x", "y", "z"), by = "trt",
                         upper_level = "Z", lower_level = "D"), "not present")
})
