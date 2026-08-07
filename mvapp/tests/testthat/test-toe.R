test_that("compute_toe measures hours from the earliest timestamp", {
  ts <- as.POSIXct(c("2024-01-01 00:00", "2024-01-01 12:00",
                     "2024-01-02 00:00"), tz = "UTC")
  expect_equal(compute_toe(ts), c(0, 12, 24))
})

test_that("compute_toe accepts character timestamps and a custom origin", {
  x <- compute_toe(c("2024-01-01 06:00", "2024-01-01 18:00"),
                   origin = as.POSIXct("2024-01-01 00:00", tz = "UTC"))
  expect_equal(x, c(6, 18))
})
