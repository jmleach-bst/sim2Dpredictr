#' Generate Probabilities for Multinomial Draws
#'
#' Obtain probabilities for each category of a multinomial distribution
#' based on covariate and parameter values based on the logit models
#' for the multinomial distribution.
#' 
#' @param V A numeric value stating the number of categories desired.
#' @param B A list, each element of which contains a parameter vector. The 
#' list should have length \code{V - 1}, i.e., should contain parameter values
#' associated with all categories except the reference category, following 
#' Agresti (2007). Alternatively, \code{B} may be a list of length \code{V} 
#' if one desires to specify parameters for every category, i.e., the
#' over-parameterized model used in Friedman (2010). 
#' @param X A matrix, each row of which contains subject covariate/predictor 
#' values. 
#' @param X.incl.X0 Logical. When \code{TRUE}, \code{X} should contain column
#' of 1's for the intercept. Otherwise, a column of 1's is generated 
#' internally. Default is \code{FALSE}.
#' @return A matrix containing subject-specific probabilities for each 
#' category of the multinomial distribution. The number of rows equals 
#' \code{nrow(X)} and the number of columns equals \code{V}. 
#' @references
#' \insertRef{Agresti:2007}{sim2Dpredictr}
#' 
#' \insertRef{Friedman:2010}{sim2Dpredictr}
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
#' ## list of parameter vectors (over-parameterized model)
#' bu <- list(b1 = c(1, 0.25, -0.25),
#'            b2 = c(-0.5, 0.15, 0.15),
#'            b3 = c(-1, 0.1, -0.20))
#'
#' ## subject specific probabilities for each category
#' generate_multinom_probs(V = vt, X = xt, B = bt)
#' 
#' ## subject specific probabilities for each category 
#' ## (over-parameterized model)
#' generate_multinom_probs(V = vt, X = xt, B = bu)
#' 
#' @export 
generate_multinom_probs <- function(V = NULL, B = NULL, X = NULL,
                                    X.incl.X0 = FALSE) {
  
  if (!(length(B) %in% c(V - 1, V))) {
    stop("Length of B must equal V - 1 or V")
  }
  
  if (X.incl.X0 == TRUE & any(X[, 1] != 1)) {
    stop("1st column of X does not contain all 1's. Did you intend X.incl.X0 = FALSE?")
  }
  
  if (X.incl.X0 == FALSE & all(X[, 1] == 1)) {
    warning("1st column of X contains all 1's. Did you intend X.incl.X0 = TRUE?")
  }
  
  if (X.incl.X0 == FALSE) {
    X <- cbind(1, X)
  }
  
  for (b in 1:length(B)) {
    if (ncol(X) != length(B[[b]])) {
      stop("X and B should have the same number of predictors.")
    }
  }
  
  if (is.matrix(X) == FALSE) {
    warning("X should be matrix. Attempting conversion.")
    X <- as.matrix(X)
  }
  
  # store exp(XB) for each of v = 1, ..., V
  exp.v.old <- NULL
  for (v in 1:length(B)) {
    exp.vv <- exp(X %*% B[[v]])
    exp.v <- cbind(exp.v.old, exp.vv)
    # print(exp.v)
    exp.v.old <- exp.v
  }
  
  # reference group is last level
  if (length(B) == V - 1) {
    exp.v <- cbind(exp.v, 1)
  }
  
  #print(exp.v)

  
  # store probabilities
  p <- matrix(NA, nrow = nrow(exp.v), ncol = V)
  for (n in 1:nrow(exp.v)) {
    for (v in 1:V) {
      p[n, v] <- exp.v[n, v] / sum(exp.v[n, ])
    }
  }
  
  return(p)
}












