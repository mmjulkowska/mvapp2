make_growth <- function() {
  set.seed(1)
  # two samples: A grows linearly, B grows exponentially
  a <- do.call(rbind, lapply(0:7, function(d)
    data.frame(geno = "A", day = d, area = 10 + 5 * d + rnorm(3, 0, 0.5))))
  b <- do.call(rbind, lapply(0:7, function(d)
    data.frame(geno = "B", day = d, area = 10 * exp(0.3 * d) + rnorm(3, 0, 0.5))))
  rbind(a, b)
}

test_that("fit_curves returns per-sample intercept/delta/r2", {
  d <- make_growth()
  res <- fit_curves(d, "day", "area", "geno", model = "linear")
  expect_s3_class(res, "mvapp_curvefit")
  expect_equal(nrow(res$table), 2)
  expect_true(all(c("geno", "n", "INTERCEPT", "DELTA", "r2") %in% names(res$table)))
  # sample A is truly linear -> near-perfect linear fit and slope ~ 5
  a <- res$table[res$table$geno == "A", ]
  expect_gt(a$r2, 0.99)
  expect_equal(a$DELTA, 5, tolerance = 0.3)
  expect_equal(a$INTERCEPT, 10, tolerance = 1)
})

test_that("estimate_models picks the right shape per sample", {
  d <- make_growth()
  em <- estimate_models(d, "day", "area", "geno")
  expect_true(all(c("linear", "exponential", "best_model") %in% names(em)))
  expect_equal(em$best_model[em$sample == "A"], "linear")
  expect_equal(em$best_model[em$sample == "B"], "exponential")
})

test_that("curve_fit_line returns a fitted curve over the x range", {
  d <- make_growth()
  b <- d[d$geno == "B", ]
  ln <- curve_fit_line(b$day, b$area, "exponential", n = 50)
  expect_equal(nrow(ln), 50)
  expect_true(all(diff(ln$y) > 0))       # exponential growth is increasing
  expect_null(curve_fit_line(c(1, 1), c(2, 3), "linear"))  # too few distinct x
})

test_that("fit_curves validates and handles multiple grouping columns", {
  d <- make_growth(); d$trt <- rep(c("c", "s"), length.out = nrow(d))
  res <- fit_curves(d, "day", "area", c("geno", "trt"))
  expect_true(all(c("geno", "trt") %in% names(res$table)))
  expect_error(fit_curves(d, "nope", "area", "geno"), "not found")
})

test_that("fit_curves fits one curve per plant and carries group metadata", {
  set.seed(2)
  # two genotypes, three plants each, linear growth; plant is the sample ID
  d <- do.call(rbind, lapply(c("A", "B"), function(g) {
    slope <- if (g == "A") 5 else 2
    do.call(rbind, lapply(1:3, function(p) do.call(rbind, lapply(0:7, function(dy)
      data.frame(geno = g, plant = paste0(g, p), day = dy,
                 area = 10 + slope * dy + rnorm(1, 0, 0.4))))))
  }))
  res <- fit_curves(d, "day", "area", samples = "plant", groups = "geno")
  # one row per plant (6 plants), genotype carried along as metadata
  expect_equal(nrow(res$table), 6)
  expect_true(all(c("plant", "geno") %in% names(res$table)))
  expect_setequal(res$samples, "plant")
  expect_setequal(res$groups, "geno")
  # each plant fits near-perfectly and A plants grow ~5, B plants ~2
  expect_true(all(res$table$r2 > 0.99))
  expect_equal(mean(res$table$DELTA[res$table$geno == "A"]), 5, tolerance = 0.3)
  expect_equal(mean(res$table$DELTA[res$table$geno == "B"]), 2, tolerance = 0.3)
})
