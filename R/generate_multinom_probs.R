#' Generate Probabiltiies for Multinomial Draws
#'
#' Obtain probabilities for each category of a multinomial distribution
#' based on covariate and parameter values based on the logit models
#' for the multinomial distribution.
#' 
#' @param V A numeric value stating the number of categories desired.
#' @param B A list, each element of which contains a parameter vector. The list should
#' have length \code{V - 1}. 
#' @param X A matrix, each row of which contains subject covariate/predictor values. \code{X}
#' should not contain a column for the intercept, as this will be added internally.
#' @return A matrix containing subject-specific probabilities for each category of the
#' multinomial distribution. The number of rows equals \code{nrow(X)} and the number of 
#' columns equals \code{V}. 
#' 
#' @export 
generate_multinom_probs <- function(V = NULL, B = NULL, X = NULL) {
  
  if (length(B) != V - 1) {
    stop("Length of B must equal V - 1")
  }
  
  for (b in 1:length(B)) {
    if (ncol(X) != length(B[[b]]) - 1) {
      stop("X and B should have the same number of predictors.")
    }
  }
  
  X <- cbind(1, X)
  
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












