test_that("mvapp_dataset builds and validates", {
  df <- data.frame(LINE = c("a", "b"), TREATMENT = c("Control", "Drought"),
                   AREA = c(1, 2))
  ds <- mvapp_dataset(df, id = "LINE", group = "TREATMENT",
                      phenotypes = "AREA")
  expect_s3_class(ds, "mvapp_dataset")
  expect_equal(ds$roles$phenotypes, "AREA")
})

test_that("validation catches missing columns and empty phenotypes", {
  df <- data.frame(LINE = "a", AREA = 1)
  expect_error(mvapp_dataset(df, id = "NOPE", phenotypes = "AREA"),
               "columns not found")
  expect_error(mvapp_dataset(df, id = "LINE", phenotypes = character()),
               "at least one phenotype")
})
