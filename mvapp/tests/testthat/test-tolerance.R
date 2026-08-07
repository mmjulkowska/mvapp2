test_that("within_group_sti divides by each group's own control mean", {
  value     <- c(10, 12, 5, 8, 20, 24, 10, 16)
  group     <- c("g1", "g1", "g1", "g1", "g2", "g2", "g2", "g2")
  condition <- c("Control", "Control", "Drought", "Drought",
                 "Control", "Control", "Drought", "Drought")
  sti <- within_group_sti(value, group, condition)
  # g1 control mean = 11 -> drought 5,8 => 5/11, 8/11
  # g2 control mean = 22 -> drought 10,16 => 10/22, 16/22
  expect_true(all(is.na(sti[condition == "Control"])))
  expect_equal(sti[3:4], c(5/11, 8/11))
  expect_equal(sti[7:8], c(10/22, 16/22))
})

test_that("sti_indices returns STI1 and STI2", {
  out <- sti_indices(drought = c(4, 9), control = c(16, 9))
  expect_equal(out$STI1, c(0.25, 1))
  expect_equal(out$STI2, c(1, 3))
})
