test_that("Error if # Categories misaligned with # of parameter sets.", {
  expect_error(generate_multinom_probs(V = 3, 
                                       B = list(b1 = 1:4,
                                                b2 = 1:4,
                                                b3 = 1:4),
                                       X = matrix(rnorm(10*3),
                                                  nrow = 10,
                                                  ncol = 3)))
})

test_that("Error if inconsistent # of parameters.", {
  expect_error(generate_multinom_probs(V = 4, 
                                       B = list(b1 = 1:4,
                                                b2 = 1:4,
                                                b3 = 1:6),
                                       X = matrix(rnorm(10*3),
                                                  nrow = 10,
                                                  ncol = 3)))
})

test_that("Valid probabilities are generated.", {
  expect_equal(rep(1, 10), 
               rowSums(generate_multinom_probs(V = 3, X = matrix(rnorm(10 * 2),
                                                                 ncol = 2, nrow = 10), 
                                               B = list(b1 = c(1, 0.25, -0.25),
                                                        b2 = c(-0.5, 0.15, 0.15)))))
})

test_that("Error for X.incl.X0 is TRUE", {
  expect_error(generate_multinom_probs(V = 3, X = matrix(rnorm(10 * 2),
                                                           ncol = 2, nrow = 10), 
                                         B = list(b1 = c(1, 0.25, -0.25),
                                                  b2 = c(-0.5, 0.15, 0.15)),
                                         X.incl.X0 = TRUE),
                 "1st column of X does not contain all 1's. Did you intend X.incl.X0 = FALSE?")
})

test_that("Warning for X.incl.X0 is FALSE", {
  expect_warning(generate_multinom_probs(V = 3, X = cbind(1, matrix(rnorm(10 * 2),
                                                           ncol = 2, nrow = 10)), 
                                         B = list(b1 = c(1, 0.25, -0.25, 1),
                                                  b2 = c(-0.5, 0.15, 0.15, 1)),
                                         X.incl.X0 = FALSE),
                 "1st column of X contains all 1's. Did you intend X.incl.X0 = TRUE?")
})