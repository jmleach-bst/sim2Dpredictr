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