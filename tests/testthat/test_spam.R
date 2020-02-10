test_that("create spam object S", {
  expect_true(spam::is.spam(chol_s2Dp(im.res = c(5, 5), use.spam = TRUE, rho = 0.5)$S))
})

test_that("create spam object Q", {
  expect_true(spam::is.spam(chol_s2Dp(im.res = c(5, 5), use.spam = TRUE,
                                      matrix.type = "prec")$Q))
})

test_that("correct class for chol. factor", {
  expect_true(class(chol_s2Dp(im.res = c(5, 5), use.spam = TRUE, rho = 0.5)$R) == 'spam.chol.NgPeyton')
})
