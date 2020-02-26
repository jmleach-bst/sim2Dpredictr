#' Construct a Precision Matrix
#'
#' This function constructs the precision matrix for a Conditional Autoregression (CAR).
#'
#' @param im.res A vector defining the dimension of spatial data. The first entry is the
#' number of rows and the second  entry is the number of columns.
#' @param tau A vector containing precision parameters. If of length 1, then all precisions are assumed equal.
#' Otherwise the length of \code{tau} should equal the number of variables.
#' @param alpha A scalar value between 0 and 1 that defines the strength of correlations. Note that when
#' \code{alpha = 0} the data are independent and when \code{alpha = 1}, the joint distribution is
#' the improper Intrinsic Autoregression (IAR), which cannot be used to generate data. Note also that  while
#' \code{alpha} does control dependence it is not interpretable as a correlation.
#' @param neighborhood Defines the neighborhood within which conditional correlations are non-zero.
#' This differs from use in \code{\link[sim2Dpredictr]{correlation_builder}}, where the neighborhood defines non-zero marginal
#' correlations. The default is \code{"ar1"}, which creates a neighborhood where the spatial locations
#' directly above, below, left, and right of a location are included in the neighborhood. More complicated
#' neighborhoods can be specified by \code{neighborhood = "round"}, which defines a circular neighborhood
#' about each location, and \code{neighborhood = "rectangle"}, which defines a rectangular neighborhood
#' about each location.
#' @param weight Determines how weights are assigned. \code{"distance"} assigns weights as the inverse of
#' Euclidean distance times a constant, \code{phi}. \code{"binary"} assigns weights to 1 for neighbors and 0 otherwise.
#' @param phi When \code{weight = "distance"} a constant by which to multiply the inverse of Euclidean distance.
#' Defaults to 1, and must exceed 0.
#' @param r If \code{neighborhood = "round"}, then locations i,j are separated by
#' distance \eqn{d \ge r} are conditionally independent.
#' @param w,h If \code{neighborhood = "rectangle"} then \code{w} and \code{h} are the number of locations
#' to the left/right and above/below a location i that define its neighborhood. Any locations
#' outside this neighborhood are conditionally independent of the specified location.
#' @param digits.Q Determines the number of digits to round entries in the precision matrix. Default is 10.
#' @importFrom matrixcalc is.positive.definite
#' @details This formulation of the CAR model is based on a formulation found in \insertCite{Banerjee:2015}{sim2Dpredictr}
#' where the joint distribution of the of the conditionally specified random variables are assumed to be
#' \eqn{N(0, [diag(tau^2)(D - alpha W)] ^ {-1})} and all neighbors are weighted 1. When weights other than 1 are desired,
#' the joint distribution is \eqn{N(0, [diag(tau^2) D (I - alpha D^{-1}W)] ^ {-1})}, e.g. as in
#' \insertCite{Jin+Carlin+Banerjee:2005}{sim2Dpredictr}.
#' @return A (precision) matrix.
#' @examples
#'
#' precision_builder(im.res = c(5, 5), tau = 1, alpha = 0.75,
#'                   neighborhood = "ar1")
#'
#' ## binary weights
#' precision_builder(im.res = c(5, 5), tau = 1, alpha = 0.75,
#'                   neighborhood = "round", r = 3)
#'
#' ## weights based on distance
#' precision_builder(im.res = c(5, 5), tau = 1, alpha = 0.75,
#'                   weight = "distance", phi = 1,
#'                   neighborhood = "round", r = 3)
#'
#' precision_builder(im.res = c(5, 5), tau = 1, alpha = 0.75,
#'                   neighborhood = "rectangle", w = 2, h = 2)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Banerjee:2015}{sim2Dpredictr}
#'
#' \insertRef{Jin+Carlin+Banerjee:2005}{sim2Dpredictr}
#'
#' @export
precision_builder <- function(im.res, tau = 1, alpha = 0.75,
                              neighborhood = "ar1", weight = "binary",
                              phi = 1, r = NULL, w = NULL, h = NULL,
                              digits.Q = 10) {
  if (length(tau) != 1) {
    if (length(tau) != prod(im.res)) {
      stop("tau must be length 1 or length prod(im.res)!")
    }
  }

  if (phi < 0) {
    stop("phi must exceed 0.")
  }

  if (any(tau <= 0)) {
    stop("Precision, tau, must exceed 0.")
  }

  if (alpha < 0 | alpha >= 1) {
    stop("alpha must be at least 0 and strictly less than 1.")
  }
  # generate adjacency matrix
  W <- proximity_builder(im.res = im.res, neighborhood = neighborhood,
                         type = "full", weight = weight, phi = phi,
                         r = r, w = w, h = h)

  # Diagonal matrix containing the number of neighbors for each location
  D <- diag(apply(W, 1, function(x) length(x[x != 0]) ))

  if (weight == "binary") {
    if (length(tau) == 1) {
      Q <- tau ^ 2 * (D - alpha * W)
    } else {
      Q <- diag(tau ^ 2) %*% (D - alpha * W)
    }
  } else {
    if (length(tau) == 1){
      Dt <- tau ^ 2 * D
    } else {
      Dt <- diag(tau ^ 2) %*% D
    }
    B <- solve(D) %*% W
    I <- diag(x = 1, nrow = prod(im.res), ncol = prod(im.res))
    Q <- Dt %*% (I - alpha * B)
    Q <- round(Q, digits = digits.Q)
  }
  if (matrixcalc::is.positive.definite(Q) == FALSE) {
    stop("Q is not positive definite. Rethink your parameter settings.")
  }
  return(Q)
}


