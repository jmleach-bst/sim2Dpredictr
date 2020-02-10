test_that("adjacency matrix is binary", {
  expect_equal(length(table(proximity_builder(im.res = c(5, 5), type = "full",
                                              weight = "binary"))),
               2)
})

test_that("proximity matrix is not binary", {
  expect_gt(length(table(proximity_builder(im.res = c(5, 5),
                                           neighborhood = "round", r = 3,
                                           type = "full",
                                           weight = "distance"))),
            2)
})

test_that("proximity matrix has max of phi", {
  expect_equal(max(proximity_builder(im.res = c(5, 5),
                                     neighborhood = "round", r = 3,
                                     phi = 2,
                                     type = "full",
                                     weight = "distance")),
               2)
})

test_that("alpha = 0 gives independent variables", {
  expect_equal(diag(diag(precision_builder(im.res = c(2, 2), alpha = 0))),
               precision_builder(im.res = c(2, 2), alpha = 0))
})

test_that("tau > 0", {
  expect_error(precision_builder(im.res = c(2, 2), tau = 0))
})

test_that("phi > 0", {
  expect_error(precision_builder(im.res = c(2, 2), phi = -0.01))
})

test_that("alpha >= 0", {
  expect_error(precision_builder(im.res = c(2, 2), alpha = -0.01))
})

test_that("alpha < 1", {
  expect_error(precision_builder(im.res = c(2, 2), alpha = 1.01))
})

test_that("Corr matrix is positive definite", {
  expect_error(correlation_builder(im.res = c(5, 5),
                                   neighborhood = "round",
                                   r = 1, rho = 0.95))
})
