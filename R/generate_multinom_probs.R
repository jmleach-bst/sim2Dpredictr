#' Generate Probabiltiies for Multinomial Draws
#'
#' Obtain probabilities for each category of a multinomial distribution
#' based on covariate and parameter values based on the logit models
#' for the multinomial distribution.
#' 
#' @param V A numeric value stating the number of categories desired.
#' @param B A list, each element of which contains a parameter vector. The list should
#' have length \code{V - 1}. 
#' @param X A matrix, each row of which contains subject covariate/predictor values. 
#' @param X.incl.X0 Logical. When \code{TRUE}, \code{X} contains a column of 1's for
#' the intercept. Otherwise, a column of 1's is generated internally. Default is \code{FALSE}.
#' @return A matrix containing subject-specific probabilities for each category of the
#' multinomial distribution. The number of rows equals \code{nrow(X)} and the number of 
#' columns equals \code{V}. 
#' @examples 
#' ## number of categories
#' vt <- 3
#' 
#' ## covariate values
#' xt <- matrix(rnorm(10 * 2), ncol = 2, nrow = 10)
#' 
#' ## list of parameter vectors
#' bt <- list(b1 = c(1, 0.25, -0.25),
#'            b2 = c(-0.5, 0.15, 0.15))
#'
#' ## subject specific probabilities for each category
#' generate_multinom_probs(V = vt, X = xt, B = bt)
#' 
#' @export 
generate_multinom_probs <- function(V = NULL, B = NULL, X = NULL,
                                    X.incl.X0 = FALSE) {
  
  if (length(B) != V - 1) {
    stop("Length of B must equal V - 1")
  }
  
  for (b in 1:length(B)) {
    if (ncol(X) != length(B[[b]]) - 1) {
      stop("X and B should have the same number of predictors.")
    }
  }
  
  if (X.incl.X0 == FALSE) {
    X <- cbind(1, X)
  }
  
  
  # store exp(XB) for each of v = 1, ..., V
  exp.v.old <- NULL
  for (v in 1:(V - 1)) {
    exp.vv <- exp(X %*% B[[v]])
    exp.v <- cbind(exp.v.old, exp.vv)
    exp.v.old <- exp.vv
  }
  
  # reference group is last level
  exp.v <- cbind(exp.v, 1)
  
  # store probabilities
  p <- matrix(NA, nrow = nrow(exp.v), ncol = V)
  for (n in 1:nrow(exp.v)) {
    for (v in 1:V) {
      p[n, v] <- exp.v[n, v] / sum(exp.v[n, ])
    }
  }
  
  return(p)
}












