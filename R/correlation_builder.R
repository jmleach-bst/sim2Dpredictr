#' Build a Correlation Matrix for 2D Spatial Data
#'
#' This function "builds" a correlation matrix based on user specifications.
#'
#' @param corr.structure One of \code{"ar1"}, \code{exponential}, 
#' \code{gaussian}, or \code{"CS"}. Correlations between locations i and j 
#' are \code{rho}\eqn{^{d}} for \code{corr.structure = "ar1"}, 
#' \eqn{exp(-phi * d)} for \code{corr.structure = "exponential"},
#' \eqn{exp(-phi * d ^ 2)} for \code{corr.structure = "gaussian"}, and 
#' \code{rho} when \code{corr.structure = "CS"}. Note that \code{d} is the
#'  Euclidean distance between locations i and j.
#' @param im.res A vector defining the dimension of spatial data. The first 
#' entry is the number of rows and the second  entry is the number of columns.
#' @param rho This is the maximum possible correlation between locations i 
#' and j. For all i,j \code{rho} MUST be between -1 and 1.
#' @param phi A scalar value greater than 0 that determines the decay rate of
#' correlation. This argument is only utilized when \code{corr.structure 
#' \%in\% c("exponential", "gaussian")}.
#' @param round.d If \code{round.d = TRUE}, then d is rounded to the nearest 
#' whole number.
#' @param corr.min Scalar value to specify the minimum non-zero correlation. 
#' Any correlations below \code{corr.min} are set to 0. Especially for high 
#' image resolution using this option can result in a sparser covariance 
#' matrix, which may significantly speed up draws when using \code{spam}. 
#' This option is preferred to using \code{neighborhood} and associated 
#' arguments when the primary concern is to avoid very small correlations 
#' and improve computation efficiency. Default is \code{NULL}, which places
#' no restrictions on the correlations.
#' @param neighborhood Defines the neighborhood within which marginal 
#' correlations are non-zero. The default is \code{"none"}, which allows 
#' marginal correlations to extend indefinitely. \code{neighborhood = "round"}
#' defines a circular neighborhood about locations and 
#' \code{neighborhood = "rectangle"} defines a rectangular neighborhood about
#' locations. Note that this argument differs from that in 
#' \code{\link[sim2Dpredictr]{precision_builder}}, in which \code{neighborhood} defines conditional non-zero
#' correlations.
#' @param r If \code{neighborhood = "round"}, then if locations i,j are 
#' separated by distance \eqn{d \ge r}, the correlation between them is zero.
#' @param w,h If \code{neighborhood = "rectangle"} then w and h are the number
#' of locations to the left/right and above/below a location i that define 
#' its neighborhood. Any locations outside this neighborhood have have zero 
#' correlation with location i.
#' @param print.all If \code{print.all = TRUE}, then prints each correlation 
#' and allows you to check whether the correlations are as you intended. This
#' option is NOT recommended for large point lattices/images.
#' @importFrom matrixcalc is.positive.definite
#' @note Caution is recommended when using \code{corr.min} or 
#' \code{neighborhood} to set many correlations to 0, as not all 
#' specifications will result in a positive definite matrix. In particular, 
#' sharp drop-offs tend to result in non-positive definite matrices.
#' @return Returns \eqn{(nr*nc) by (nr*nc)} correlation matrix.
#' @examples
#' ## examples
#' correlation_builder(corr.structure = "ar1", im.res = c(3, 3), rho = 0.5,
#'                     neighborhood = "round", r = 6, print.all = TRUE)
#'
#' correlation_builder(corr.structure = "exponential", im.res = c(3, 3), 
#'                     phi = 0.5,
#'                     neighborhood = "round", r = 3, print.all = TRUE)
#'
#' correlation_builder(corr.structure = "CS", im.res = c(3, 3),
#'                     rho = 0.5, print.all = TRUE)
#'
#' ## no "true" zeros, but gets close
#' c.nr <- correlation_builder(corr.structure = "ar1", neighborhood = "none",
#'                     corr.min = NULL, im.res = c(15, 15), rho = 0.5)
#' length(c.nr[c.nr > 0])
#' min(c.nr)
#'
#' ## set corr.min gives many zero entries; sparser structure
#' c.r <- correlation_builder(corr.structure = "ar1", neighborhood = "none",
#'                     corr.min = 0.01, im.res = c(15, 15), rho = 0.5)
#' ## raw number > 0
#' length(c.r[c.r > 0])
#' ## proportion > 0
#' length(c.r[c.r > 0]) / length(c.nr)
#' @export
correlation_builder <- function(corr.structure = "ar1", im.res,
                                corr.min = NULL, neighborhood = "none",
                                rho = NULL, phi = NULL,
                                w = NULL, h = NULL, r = NULL,
                                print.all = FALSE, round.d = FALSE
){
  if (corr.structure %in% c("exponential", "gaussian")) {
    if (is.null(phi) == FALSE) {
      if (phi <= 0) {
        stop(paste0("phi must be > 0, but current value is ", phi))
      }
    } else {
      if (is.null(phi) == TRUE) {
        stop("You need to specify a value for phi. (Note: phi > 0)")
      }
    }

  }

  if (corr.structure %in% c("ar1", "CS")) {
    if (is.null(rho) == FALSE) {
      if ((rho < -1) | (rho > 1)) {
        stop(paste0("Correlation must be between -1 and 1, but rho = ", 
                    rho, "."))
      }
    } else {
      if (is.null(rho) == TRUE) {
        stop("You need to specify a value for rho. (Note: -1 < rho < 1)")
      }
    }
  }

  # Used to determine/ensure that w, h are whole numbers if not NULL
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^2)  abs(x - round(x)) < tol

  if ( !(corr.structure %in% c("ar1", "CS", "exponential", "gaussian")) ) {
    stop("Invalid correlation structure. \n
         Options are ar1, CS, exponential, or gaussian.")
  }
  if ( !(neighborhood %in% c("round", "rectangle", "none")) ) {
    stop("Invalid neighborhood: Choose none, round, or rectangle.")
  }

  if (neighborhood == "rectangle") {
    if ( is.null(w) == FALSE ) {
      if ( w < 0 ){
        stop("Invalid neighborhood dimensions! Choose w > 0.")
      } else if (is.wholenumber(w) == FALSE) {
        warning("rounding w to the nearest whole number.")
        w <- round(w)
      }
    } else {
      if (is.null(w) == TRUE) {
        stop("Must specify a value for w. (Note: w is a positive whole number)")
      }
    }

    if ( is.null(h) == FALSE ) {
      if ( h < 0 ){
        stop("Invalid neighborhood dimensions! Choose h > 0.")
      } else if (is.wholenumber(h) == FALSE) {
        warning("rounding h to the nearest whole number.")
        h <- round(h)
      }
    } else {
      if (is.null(h) == TRUE) {
        stop("Must specify a value for h. (Note: h is a positive whole number)")
      }
    }
  }

  if (neighborhood == "round") {
    if ( is.null(r) == FALSE ) {
      if ( r < 0 ){
        stop("Invalid neighborhood dimensions! Choose r > 0.")
      }
    } else {
      if (is.null(r) == TRUE) {
        stop("Must specify a value for r. (Note: r > 0)")
      }
    }
  }

  p <- prod(im.res)
  R.up <- diag(p)
  # need an odd workaround to impose i>=j
  K <- 1:im.res[1]
  V <- 1:im.res[2]

  for (i in 1:im.res[1]) {
    for (j in 1:im.res[2]) {
      for (k in K[K >= i]) {
        if (k == i) {
          for (v in V[V >= j]) {
            if (!(i == k & j == v)) {
              rho.ij.kv <- corr_fun(corr.structure = corr.structure, 
                                    im.res = im.res,
                                    rho = rho, phi = phi, 
                                    corr.min = corr.min,
                                    round.d = round.d, 
                                    neighborhood = neighborhood,
                                    w = w, h = h, r = r, 
                                    i = i, j = j, k = k, v = v)
              R.up[((i - 1) * im.res[2] + j),((k - 1) * im.res[2] + v)] <- rho.ij.kv
              if (print.all == TRUE){
                cat("Correlation b/t location", i, j, 
                    "and location", k, v, "is", rho.ij.kv, "\n")
              }
            }
          }
        } else {
          for (v in V) {
            if (!(i == k & j == v)) {
              rho.ij.kv <- corr_fun(corr.structure = corr.structure, 
                                    im.res = im.res,
                                    rho = rho, phi = phi, 
                                    corr.min = corr.min,
                                    round.d = round.d, 
                                    neighborhood = neighborhood,
                                    w = w, h = h, r = r, 
                                    i = i, j = j, k = k, v = v)
              R.up[((i - 1) * im.res[2] + j),((k - 1) * im.res[2] + v)] <- rho.ij.kv
              if (print.all == TRUE){
                cat("Correlation b/t location", i, j, 
                    "and location", k, v, "is", rho.ij.kv, "\n")
              }
            }
          }
        }
      }
    }
  }
  R <- t(R.up) + R.up - diag(rep(1, p))
  if (matrixcalc::is.positive.definite(R) == FALSE) {
    stop("Corr matrix is not positive definite. Rethink your parameter settings.")
  }
  return(R)
}
