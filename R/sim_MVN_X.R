#' Simulate Spatially Correlated MVN Data
#'
#' Takes N draws from a Multivariate Normal (MVN) distribution using either
#' base R or the R package \code{spam}. This function requires the Cholesky
#' decomposition of the desired covariance matrix.
#'
#' @param N The number of draws to take from MVN; i.e., the number of subjects.
#' @param mu One of the following:
#'  \itemize{
#'     \item A single scalar value for common mean.
#'     \item A vector of length \code{nrow(R)} (equivalently \code{nrow(R)})
#'      containing means for the MVN.
#'  }
#' @param L,R \code{L} and \code{R} are lower and upper triangular matrices, respectively,
#'  and are the Choleskly factor(s) of the desired covariance matrix for the MVN.
#'  Obtain \code{L} or \code{R} via \code{chol_s2Dp()} with settings
#'  \code{triangle = "lower"} or \code{triangle = "upper"}, respectively.
#'  Specify either \code{L} or \code{R}, but NOT both.
#' @param use.spam Logical. If \code{use.spam = TRUE} then use tools from the R package \code{spam};
#' otherwise, base R functions are employed. For large dimension MVN with sparse correlation
#' structure, \code{spam} is recommended; otherwise, base R may be faster. Defaults to \code{FALSE}.
#' Requires either the covariance matrix \code{S} or precision matrix, \code{Q}, that corresponds to the
#' Cholesky factor.
#' @param use.MASS Logical. When \code{TRUE} draws X from MVN using \code{mvrnorm} from \code{MASS}.
#' Note that this requires specification of the covariance matrix, \code{S}. Specifying the precision
#' matrix instead may slow down the process for large dimensions. Recommended to use \code{spam} to
#' generate draws when specifying a precision matrix, \code{Q}.
#' @param S,Q A covariance or precision matrix respectively. These are for use with \code{spam},
#' and can be extracted from output of \code{\link[sim2Dpredictr]{chol_s2Dp}} after choosing
#' \code{return.cov = TRUE} or \code{return.prec = TRUE}, respectively.
#' @param X.categorical Default is \code{X.categorical = FALSE}. If
#' \code{X.categorical = TRUE} then thresholds are applied to categorize
#' each predictor/image value.
#' @param X.num.categories A scalar value denoting the number of categories
#' in which to divide the data.
#' @param X.category.type Tells R how to categorize the data. Options are
#'  \code{X.category.type = c("percentile", "manual")}.
#'  If \code{X.category.type = "percentile"} then the data are divided into
#'  percentiles based on \code{X.num.categories}; e.g. if \code{X.num.categories = 4}
#'  then the values are divided into quartiles, and values in Q1 equal 0, betwen Q1
#'  and Q2 equal 1, between Q2 and Q3 equal 2, and greater than Q3 equal 3.
#'  If \code{X.category.type = "manual"} then specify the cutoff points with
#'  \code{X.manual.thresh}.
#' @param X.percentiles A vector of percentiles to be used in thresholding when
#' \code{X.categorical = TRUE} and \code{X.category.type = "percentile"}. The length
#' of this vector should equal the number of categories minus one, and all values
#' should be between zero and one.
#' @param X.manual.thresh A vector containing the thresholds for categorizing the
#'  values; e.g. if \code{X.num.categories = 4} and \code{X.manual.thresh = c(-3, 1, 17)},
#'  then values less than -3 are set to 0, equal or greater than -3 and less than 1
#'  are set to 1, equal or greater than 1 but less than 17 are set to 2, and equal or
#'  greater than 17 are set to 3. Note that \code{length(X.manual.thresh)} must always
#'  equal \code{X.num.categories - 1}.
#' @param X.cat.names A vector of category names. If \code{X.cat.names = NULL} then
#'  then the initial integers assigned are left as the values; the names in
#'  \code{X.cat.names} are assinged in ascending order.
#' @examples
#' ## verify MVN with base R
#' set.seed(732)
#' Lex <- chol_s2Dp(corr.structure = "ar1", im.res = c(3, 3), rho = 0.25,
#'                  sigma = 1, use.spam = FALSE, corr.min = 0.02,
#'                  triangle = "lower", return.cov = TRUE)
#' XbR = sim_MVN_X(N = 1000, mu = 0, L = Lex$L)
#'
#' apply(XbR, 2, mean)
#' cov(XbR)
#' Lex$S
#'
#' ## verify MVN with \code{spam}
#' set.seed(472)
#' Rex <- chol_s2Dp(im.res = c(3, 3), matrix.type = "prec",
#'                 use.spam = TRUE, neighborhood = "ar1",
#'                 triangle = "upper", return.prec = TRUE)
#'
#' Xspam = sim_MVN_X(N = 1000, mu = 0, R = Rex$R, Q = Rex$Q)
#'
#' apply(Xspam, 2, mean)
#' solve(cov(Xspam))
#' as.matrix(Rex$Q)
#'
#' ## Categories
#' set.seed(832)
#' Xtest <- sim_MVN_X(N = 30, mu = 0, L = Lex$L,
#'                    X.categorical = TRUE,
#'                    X.num.categories = 3,
#'                    X.category.type = "percentile",
#'                    X.cat.names = c("A", "B", "C"))
#' Xtest
#'
#' @return Matrix of dimension \code{N} x \code{(nrow(L))} (or equivalently \code{N} x \code{(nrow(R))})
#' where each row is draw from MVN, and each column represents a different "variable";
#' e.g. location in an image.
#' @note This function requires the Cholesky decomposition of the desired covariance
#'  matrix for the MVN; this allows for using this function in simulating multiple
#'  datasets of \code{N} MVN draws while only taking the Cholesky decomposition of
#'  the covariance matrix once.
#' @importFrom stats rnorm quantile
#' @importFrom Rdpack reprompt
#' @importFrom MASS mvrnorm
#' @importFrom spam rmvnorm.spam rmvnorm.prec
#' @references
#' \insertRef{spam}{sim2Dpredictr}
#'
#' \insertRef{Ripley:1987}{sim2Dpredictr}
#'
#' \insertRef{Rue:2001}{sim2Dpredictr}
#' @export
sim_MVN_X <- function(N, mu = 0, L = NULL, R = NULL,
                      S = NULL, Q = NULL,
                      use.spam = FALSE, use.MASS = FALSE,
                      X.categorical = FALSE, X.num.categories = 2,
                      X.category.type = "percentile",
                      X.percentiles = NULL,
                      X.manual.thresh = NULL,
                      X.cat.names = NULL){

  # Predictable errors and warnings
  if ( (is.null(R) == TRUE) & (is.null(L) == TRUE) & (use.MASS == FALSE & is.null(S) == TRUE) ){
    stop("Function requires either specification of L or R unless using MASS.")
  }
  if ( (is.null(R) == FALSE) & (is.null(L) == FALSE) ){
    warning("Specify either L or R, not both. Function proceeding with R.")
    L <- NULL
  }

  if (class(R) == "spam.chol.NgPeyton" | class(L) == "spam.chol.NgPeyton") {
    use.spam <- TRUE
    if (is.null(Q) == FALSE) {
      Q <- spam::as.spam(Q)
    }
    if (is.null(S) == FALSE) {
      S <- spam::as.spam(S)
    }
  }

  if (use.spam == TRUE & (is.null(Q) == TRUE & is.null(S) == TRUE)) {
    stop("Require either covariance or precision matrix. \n")
  }

  if (use.spam == TRUE & (is.null(Q) == FALSE & is.null(S) == FALSE)) {
    warning("S and Q both specified. Dropping Q and keeping S. \n")
    Q <- NULL
  }

  # When using {spam} convert matrices to spam objects if not already.

  if (is.null(R) == FALSE & use.spam == TRUE) {
    if (class(R) != "spam.chol.NgPeyton") {
      R <- spam::as.spam(R)
    }
  }
  if (is.null(L) == FALSE & use.spam == TRUE) {
    if (class(L) != "spam.chol.NgPeyton") {
      L <- spam::as.spam(L)
    }
  }

  if (is.null(Q) == FALSE & use.spam == TRUE) {
    if (class(Q) != "spam") {
      Q <- spam::as.spam(Q)
    }
  }

  if (is.null(S) == FALSE & use.spam == TRUE) {
    if (class(S) != "spam") {
      S <- spam::as.spam(S)
    }
  }

  # errors/warnings regarding MASS
  if (use.MASS == TRUE & (is.null(S) == TRUE & is.null(Q) == TRUE)) {
    stop("MASS function mvrnorm requires covariance matrix, S. \n")
  }

  if (use.MASS == TRUE & is.null(S) == TRUE) {
    warning("Better to specify S when using MASS. \n")
    S <- solve(Q)
  }

  # Obtain number of parameters.
  LRSQ <- list(L = L, R = R, S = S, Q = Q)
  LRSQw <- LRSQ[which(lapply(LRSQ, is.null) == FALSE)]
  LRSQw1 <- LRSQw[1]
  p <- dim(LRSQw[[1]])[1]

  if (length(mu) == 1) {
    mu <- rep(mu, p)
  } else {
    if (length(mu) != p) {cat("Invalid mean vector length.", "\n")}
  }
  if (use.spam == FALSE) {
    if (is.null(L) == TRUE & use.MASS == FALSE) {
      L <- t(R)
    }
    if (use.MASS == TRUE) {
      X <- MASS::mvrnorm(n = N, mu = mu, Sigma = S)
    } else {
      Z <- matrix(rnorm(p * N), ncol = N, nrow = p)
      mu.mat <- matrix(rep(mu, N), nrow = length(mu), ncol = N)
      X <- mu.mat + L %*% Z
      X <- t(X)
    }
  }
  if (use.spam == TRUE) {
    if (is.null(R) == TRUE) {
      R <- spam::t(L)
    }
    if (is.null(S) == FALSE) {
      X <- spam::rmvnorm.spam(n = N, mu = mu, Rstruc = R, Sigma = S)
    } else {
      X <- spam::rmvnorm.prec(n = N, mu = mu, Rstruc = R, Q = Q)
    }
  }

  # Threshold to produce categories (if desired)
  if (X.categorical == TRUE) {

    # Determine/ensure that num.categories is a whole number.
    is.wholenumber <- function(x, tol = .Machine$double.eps^2)  abs(x - round(x)) < tol
    if (is.wholenumber(X.num.categories) == FALSE) {
      warning("Rounding the number of categories to the nearest whole number.")
      X.num.categories = round(X.num.categories)
    }

    if (X.category.type == "percentile") {
      if (is.null(X.percentiles) == FALSE) {
        if ( length(X.percentiles) != (X.num.categories - 1) ) {
          stop("The number of percentiles should equal the number categories minus 1.")
        } else {
          if (any(X.percentiles >= 1) | any(X.percentiles <= 0)) {
            stop("Percentiles must all be between zero and one.")
          } else {
            Xps <- quantile(X, probs = X.percentiles)
          }
        }
      } else {
        Xps <- quantile(X, probs = seq( 1 / X.num.categories, 1,
                                        1 / X.num.categories))
        }
      } else if (X.category.type == "manual") {
        if (length(X.manual.thresh) != (X.num.categories - 1)) {
          stop("The number manual threshold values should equal the number of categories minus 1.")
        } else {
          Xps <- c(X.manual.thresh)
        }
      }
    if (X.category.type == "percentile") {
      cat("Category boundaries: ", Xps, "\n")
    }
    Xps <- c(Xps, max(X))
    Xcat <- rep(NA, length(X))
    Xv <- as.vector(t(X))
    for (i in 1:X.num.categories) {
      Xcat[(Xv <= Xps[i]) & (is.na(Xcat) == TRUE)] <- (i - 1)
    }
    X <- matrix(Xcat, byrow = TRUE, nrow = N)
    if (is.null(X.cat.names) == FALSE) {
      if (length(X.cat.names) != X.num.categories) {
        stop("Must have the same number of category names as categories.")
      } else {
        for (i in 1:X.num.categories) {
          if (is.null(X.cat.names) == TRUE) {
            X.cat.names <- 0:(X.num.categories - 1)
            print(X.cat.names)
          }
          X[X == (i - 1)] <- X.cat.names[i]
        }
      }
    }
  }

  return(X)
}
