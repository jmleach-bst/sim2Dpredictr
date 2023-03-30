#' Build and Take the Cholesky Decomposition of a Covariance Matrix
#'
#' The function first builds a correlation matrix with correlation.builder, 
#' converts that matrix to a covariance matrix if necessary, and then takes 
#' the Cholesky decomposition of the matrix using either base R or the R 
#' package \code{spam}. Note that \code{spam} is particularly effective when
#' the matrix is sparse.
#'
#' @param matrix.type Determines whether to build a covariance matrix, 
#' \code{"cov"}, or a precision matrix, \code{"prec"}. See 
#' \code{\link{correlation_builder}{sim2Dpredictr}} and
#' \code{\link{precision_builder}{sim2Dpredictr}} for more details.
#' @inheritParams correlation_builder
#' @inheritParams precision_builder
#' @param use.spam If \code{use.spam = TRUE} then use tools from the R package
#' \code{spam}; otherwise, base R functions are employed. For large dimension
#' MVN with sparse correlation structure, \code{spam} is recommended; 
#' otherwise, base R may be faster. Defaults to \code{FALSE}.
#' @param triangle Determine whether to output an upper 
#' (\code{triangle = "upper"}) or lower (\code{triangle = "lower"}) triangular
#' matrix.
#' @param sigma Specify the desired standard deviations; the default is 1, in
#' which case the Cholesky decomposition is of a correlation matrix. If 
#' \code{sigma} != 1, then the Cholesky decomposition is of a covariance 
#' Matrix.
#'   \itemize{
#'     \item If sigma is a vector then length(sigma) must be equal to the
#'     total number of locations, i.e. \eqn{(n.row * n.col) by (n.row * n.col)}.
#'     \item sigma can take any scalar value when specifying common standard
#'     deviation.
#'   }
#' @param return.cov,return.prec Logical. When \code{TRUE}, also return the 
#' covariance or precision matrix, respectively. This is recommended when 
#' using \code{spam} to generate draws from the MVN.
#' @param print.R,print.S,print.Q Logical. When \code{TRUE}, then print the 
#' correlation, covariance, or precision matrix before taking the Cholesky 
#' decomposition. If \code{sigma} = 1, then S = R.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Banerjee:2015}{sim2Dpredictr}
#'
#' \insertRef{Ripley:1987}{sim2Dpredictr}
#'
#' \insertRef{Rue:2001}{sim2Dpredictr}
#'
#' \insertRef{spam}{sim2Dpredictr}
#' @return Matrix of dimension (n.row x n.col) x (n.row x n.col). If either 
#' \code{return.cov} or \code{return.prec} is \code{TRUE}, then returns a 
#' list where the first element is the covariance or precision matrix, and the
#' second element is the Cholesky factor.
#' @examples
#'
#' ## Use R package spam for Cholesky decomposition
#' R_spam <- chol_s2Dp(im.res = c(3, 3), matrix.type = "prec",
#'                     use.spam = TRUE, neighborhood = "ar1",
#'                     triangle = "upper")
#'
#' ## Use base R for Cholesky decomposition
#' R_base <- chol_s2Dp(corr.structure = "ar1", 
#'                     im.res = c(3, 3), rho = 0.15,
#'                     neighborhood = "round", 
#'                     r = 3, use.spam = FALSE)
#'
#' ## Specify standard deviations instead of default of sigma = 1.
#' R_sd <- chol_s2Dp(corr.structure = "ar1",
#'                   im.res = c(3, 3), rho = 0.15,
#'                   neighborhood = "round", r = 3, 
#'                   sigma = runif(9, 1.1, 4))
#' \dontrun{
#' ## Print options ON
#' R_pr_on <- chol_s2Dp(corr.structure = "ar1", 
#'                      im.res = c(3, 3), rho = 0.15,
#'                      sigma = 1:9, neighborhood = "round", 
#'                      r = 3, print.R = TRUE, print.S = TRUE)
#'  }
#'
#' @export
chol_s2Dp <- function(matrix.type = "cov", im.res, use.spam = FALSE,
                     corr.structure = "ar1", rho = NULL, phi = NULL,
                     tau = 1, alpha = 0.75, corr.min = NULL,
                     neighborhood = "none", w = NULL, h = NULL, r = NULL,
                     print.R = FALSE, print.S = FALSE, print.Q = FALSE,
                     sigma = 1, triangle = "upper",
                     print.all = FALSE, round.d = FALSE,
                     return.cov = TRUE, return.prec = TRUE){
  p <- prod(im.res)

  # use correlation matrix
  if (matrix.type == "cov") {
    # For unit variance, the covariance matrix is the correlation matrix
    if (corr.structure %in% c("ar1", "CS")) {
      if (is.null(rho) == TRUE) {
        stop("You must supply a value for rho (between -1 and 1).")
      }
      if (rho == 0) {
        R <- diag(p)
      } else {
        R <- correlation_builder(corr.structure = corr.structure,
                                 rho = rho, phi = phi,
                                 corr.min = corr.min,
                                 im.res = im.res, round.d = round.d,
                                 w = w, h = h, r = r, 
                                 neighborhood = neighborhood,
                                 print.all = print.all)
      }
    } else if (corr.structure %in% c("exponential", "gaussian")) {
      if (is.null(phi) == TRUE) {
        stop("You must supply a value for phi (phi > 0).")
      }
      R <- correlation_builder(corr.structure = corr.structure,
                               corr.min = corr.min,
                               rho = rho, phi = phi,
                               im.res = im.res, round.d = round.d,
                               w = w, h = h, r = r, 
                               neighborhood = neighborhood,
                               print.all = print.all)
    }
    if (print.R == TRUE) {
      cat("Below lies the correlation matrix. \n")
      print(R)
    }
    if ( identical(sigma, 1) | identical(sigma, rep(1, p)) ){S <- R}
    # if not using unit variance then need separate covariance matrix
    else {
      D <- diag(sigma, nrow = p, ncol = p)
      S <- D %*% R %*% D
    }
    if (print.S == TRUE) {
      cat("Below lies the covariance matrix. \n")
      print(S)
    }
    if(identical(S, diag(prod(im.res)))) {
      Rc <- S
      L <- S
    }
    else {
      if(use.spam == FALSE) {
        Rc <- chol(S)
        L <- t(Rc)
      }
      if(use.spam == TRUE) {
        S <- spam::as.spam(S)
        Rc <- spam::chol.spam(S)
        L <- spam::t(Rc)
      }
    }
    if (triangle == "upper") {
      if (return.cov == TRUE) {
        return(list(S = S, R = Rc))
      } else {
        return(Rc)
      }
    }
    if (triangle == "lower") {
      if (return.cov == TRUE) {
        return(list(S = S, L = L))
      } else {
        return(L)
      }
    }
  }

  # use precision matrix
  if (matrix.type == "prec") {
    if (neighborhood == "none") {
      neighborhood <- "ar1"
    }
    Q <- precision_builder(im.res = im.res, neighborhood = neighborhood,
                           tau = tau, alpha = alpha, w = w, h = h, r = r)
    if (print.Q == TRUE) {
      cat("Below lies the precision matrix. \n")
      print(Q)
    }
    if (print.S == TRUE) {
      cat("Below lies the covariance matrix. \n")
      print(solve(Q))
    }
    if (use.spam == FALSE) {
      Rc <- chol(Q)
      L <- t(Rc)
    }
    if (use.spam == TRUE) {
      Q <- spam::as.spam(Q)
      Rc <- spam::chol.spam(Q)
      L <- spam::t(Rc)
    }
    if (triangle == "upper") {
      if (return.prec == TRUE) {
        return(list(Q = Q, R = Rc))
      } else {
        return(Rc)
      }
    }
    if (triangle == "lower") {
      if (return.prec == TRUE) {
        return(list(Q = Q, L = L))
      } else {
        return(L)
      }
    }
  }
}

