test_that("scatter_stats computes R2 = r^2 overall", {
  set.seed(1)
  x <- rnorm(50); y <- 2 * x + rnorm(50, 0, 0.5)
  df <- data.frame(x = x, y = y)
  st <- scatter_stats(df, "x", "y")
  expect_equal(nrow(st), 1)
  expect_equal(st$group, "all")
  expect_equal(st$n, 50)
  expect_equal(st$R2, cor(x, y)^2)
  expect_equal(st$p, cor.test(x, y)$p.value)
})

test_that("scatter_stats returns one row per group", {
  df <- data.frame(
    g = rep(c("Control", "Drought"), each = 15),
    x = c(1:15, 1:15),
    y = c(1:15 + rnorm(15, 0, 0.01), 15:1 + rnorm(15, 0, 0.01))
  )
  st <- scatter_stats(df, "x", "y", group = "g")
  expect_equal(nrow(st), 2)
  expect_setequal(st$group, c("Control", "Drought"))
  ctrl <- st[st$group == "Control", ]
  drt  <- st[st$group == "Drought", ]
  expect_gt(ctrl$r, 0.99)     # positive slope
  expect_lt(drt$r, -0.99)     # negative slope
  expect_gt(ctrl$R2, 0.98)    # both near-perfect fits
  expect_gt(drt$R2, 0.98)
})

test_that("scatter_stats degrades gracefully and validates", {
  df <- data.frame(x = c(1, 2, NA), y = c(2, NA, 6), z = c(1, 1, 1))
  st <- scatter_stats(df, "x", "y")          # <3 complete pairs
  expect_true(is.na(st$R2))
  st2 <- scatter_stats(data.frame(x = 1:5, z = rep(1, 5)), "x", "z")  # zero variance
  expect_true(is.na(st2$R2))
  expect_error(scatter_stats(df, "x", "nope"), "not found")
})
