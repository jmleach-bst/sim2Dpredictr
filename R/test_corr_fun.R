test_that("Set corr < corr.min to 0", {
  expect_equal(
    corr_fun(corr.structure = "ar1",
             i = 1, j = 1, k = 5, v = 5,
             rho = 0.5, corr.min = 0.02),
    0
  )
})

test_that("Leaves corr > corr.min alone ", {
  expect_gt(
    corr_fun(corr.structure = "ar1",
             i = 1, j = 1, k = 2, v = 2,
             rho = 0.5, corr.min = 0.02),
    0
  )
})
