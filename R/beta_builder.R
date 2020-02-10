#' Create a Parameter Vector from Lattice Locations
#'
#' Specify the locations in the lattice/image that have non-zero parameters
#'  as well as the values for those parameters, and the function creates
#'  the parameter vector that matches the correct locations in the design matrix.
#' @param row.index,col.index Vectors of row/columns indices for non-zero parameters.
#' If \code{index.type = "manual"}, each vector should contain specific coordinates.
#' If \code{index.type = "rectangle"}, each vector should specify rectangle length;
#' e.g. row.index = 1:3 means the rectangle's 'length' is from rows 1 to 3.
#' If \code{index.type = "ellipse"}, the arguments should be scalar values specifying
#' the row/column coordinates for the center of the ellipse.
#' If \code{index.type = "decay"}, the arguments should specify the row/column coordinates,
#' respectively, of the peak parameter value.
#' @param im.res A vector specifying the dimension/resolution of the image. The first entry is
#'  the number of 'rows' in the lattice/image, and the second entry is the number of
#'  'columns' in the lattic/image.
#' @param B0 is the "true" intercept value; default is 0.
#' @param B.values is a vector "true" parameter values for non-zero parameters.
#' The order of assignment is by row. If B.values argument is a single value,
#' then all non-zero parameters are assigned to that value, unless \code{decay.fn}
#' has been specified, in which case \code{B.values} is the "peak", and non-zero
#' parameters decay smoothly by distance.
#' @param index.type is one of index.type = c("manual", "rectangle", "ellipse", "decay")
#' \itemize{
#'     \item \code{index.type = "manual"} uses row.index and col.index arguemnts to specify
#'                   manually selected non-zero locations. This setting is good
#'                   for irregular shaped regions.
#'     \item \code{index.type = "rectangle"} uses row.index and col.index arguments to
#'                   specify a rectangular region of non-zero parameters.
#'     \item \code{index.type = "ellipse"} uses \code{w} and \code{h} argmuents to specify
#'                   elliptical region of non-zero parameters
#'     \item \code{index.type = "decay"} allows the user to specify a peak location with
#'           \code{row.index} and \code{col.index}, as with \code{index.type = "ellipse"}.
#'           However, the non-zero parameter values decay as a function of distance from
#'           the peak.
#' }
#' @param decay.fn An argument to specify the decay function of non-zero parameters decay
#' from the peak when \code{index.type = "decay"}. Options are "exponential" or "gaussian".
#' The rate of decay is given by \eqn{B.values * exp(-phi * d)} and
#' \eqn{B.values * exp(-phi * d ^ 2)} for "exponential" and "gaussian", respectively.
#' The default is \code{decay.fn = "gaussian"}. Note that \eqn{d} is the Euclidean distance
#' between the peak and a specified location, while \eqn{phi} is the rate of decay and is set
#' by the user with \code{phi}.
#' @param phi A scalar value greater than 0 that determines the decay rate of non-zero
#' parameters when \code{index.type = "decay"}. The default is \code{phi = 0.5}.
#' @param max.d When \code{index.type = "decay"}, \code{max.d} determines the maximum
#' Euclidean distance from the peak that is allowed to be non-zero; parameters for locations
#' further than \code{max.d} from the peak are set to zero. If this argument is not set by
#' the user then all parameter values are determined by the decay function.
#' @param w,h If index.type = "ellipse" then the width and height of the ellipse, respectively.
#' @param bayesian If \code{TRUE}, then parameters are drawn from distributions based on initial
#' \code{B.values} vector. Default is \code{FALSE}.
#' @param bayesian.dist When \code{bayesian = TRUE}, specifies the distribution of the
#' parameters. Options are \code{"gaussian"} and \code{uniform}.
#' @param bayesian.scale A list. When \code{bayesian = TRUE} and \code{bayesian.dist = "gaussian"},
#' specifies the sd for the distributions of parameters. When \code{bayesian = TRUE} and
#' \code{bayesian.dist = "uniform"}, specifies the width for the uniform distributions for the
#' parameters. The first entry should be one of \code{"unique", "binary"}. If \code{"unique"},
#' then the second entry in the list should be a vector with length equal to \code{B.values + 1}
#' with unique values for the sd's/widths, including B0. B0 can be set to a constant value by setting
#' the first position of \code{bayesian.scale[[2]]} to 0. If \code{"binary"}, then the second entry
#' in the list should be a 3-element vector whose first entry is the sd/width of B0, second entry the
#' sd/width of "non-zero" or "important" parameters, and the third entry is the sd/width of the "zero"
#' or "irrelevant"  parameters.
#' @param output.indices If \code{output.indices = TRUE}, then the first element of the
#' returned list contains the indices for the non-zero parameter locations (Default).
#' If \code{output.indices = FALSE}, then only the parameter vector is returned.
#' @return A list containing (1) indices for the locations of "true" non-zero parameters,
#' and (2) a parameter vector.
#' @note The order of the parameters is by row. That is, if the lattice/image is 4x4,
#' then parameters 1-4 make up the first row, 5-8 then second, and so forth.
#' @examples
#' ## example
#' Bex1 <- beta_builder(row.index = c(5, 5, 6, 6),
#'                     col.index = c(5, 6, 5, 6),
#'                     im.res = c(10, 10),
#'                     B0 = 0, B.values = rep(1, 4))
#'
#' ## True non-zero parameters are locations 45, 46, 55, 56 in B
#' ## i.e. locations (5, 5), (5, 6), (6, 5), (6, 6)
#'
#' ## Suppose that we index rows by i = 1, ... , I
#' ##                       cols by j = 1, ... , J
#'
#' ## The index for a parameter is given by  J * (i - 1) + j
#' ## In this example, I = 10, J = 10; Thus:
#'
#' ## (5, 5) -> 10 * (5 - 1) + 5 = 45
#' ## (5, 6) -> 10 * (5 - 1) + 6 = 46
#' ## (6, 5) -> 10 * (6 - 1) + 5 = 55
#' ## (6, 6) -> 10 * (6 - 1) + 6 = 45
#' Bex1
#' ## length 101 (includes B0 w/ 100 variable parameter values)
#' length(Bex1$B)
#'
#' ## example: index.type = "rectangle"
#' Bex2 <- beta_builder(row.index = 12:15, col.index = 6:19,
#'                     im.res = c(20, 20), B0 = 16,
#'                     B.values = 1:(length(12:15) * length(6:19)),
#'                     index.type = "rectangle")
#'
#' Bex2
#' matrix(Bex2$B[-1], nrow = 20, byrow = TRUE)
#'
#' ## example: index.type = "ellipse"
#' Bex3 <- beta_builder(row.index = 4, col.index = 5,
#'                     im.res = c(10, 10),
#'                     B0 = 16, B.values = 3,
#'                     index.type = "ellipse",
#'                     h = 5, w = 4)
#' Bex3
#' matrix(Bex3$B[-1], nrow = 10, byrow = TRUE)
#'
#' ## decaying parameter values
#' Bex4 <- beta_builder(row.index = 10, col.index = 20,
#'                      im.res = c(30, 30), B0 = 0, B.values = 10,
#'                      index.type = "decay", max.d = 7,
#'                      output.indices = FALSE)
#' inf_2D_image(B = Bex4, im.res = c(30, 30), binarize.B = FALSE)
#'
#' Bex5 <- beta_builder(row.index = 4, col.index = 5,
#'                      im.res = c(10, 10),
#'                      B0 = 16, B.values = 5,
#'                      index.type = "ellipse",
#'                      h = 5, w = 4,
#'                      bayesian = TRUE,
#'                      bayesian.dist = "gaussian",
#'                      bayesian.scale = list("binary", c(0, 1, 0.25)))
#'
#' inf_2D_image(B = Bex5$B, im.res = c(10, 10), binarize.B = FALSE)
#'
#' @export
beta_builder = function(row.index, col.index, im.res,
                        B0 = 0, B.values, index.type = "manual",
                        decay.fn = "gaussian", phi = 0.5, max.d = Inf,
                        h, w,
                        bayesian = FALSE, bayesian.dist = NULL,
                        bayesian.scale = NULL,
                        output.indices = TRUE) {

  # check for positive values of row.index, col.index
  if ( any(c(row.index, col.index) < 0) ) {
    stop("Index values must be > 0.")
  }

  # ensure values in row.index, col.index specify possible locations in image.
  if ( any(row.index > im.res[1]) | any(col.index > im.res[2]) ) {
    stop("Index values must not exceed the bounds of the image.")
  }

  # calculate exact coordinates for retangular region
  if (index.type == "rectangle") {
    is.sequential.incr <- function(x){
      all(diff(x) == 1)
    }

    if( (is.sequential.incr(row.index) == FALSE) |
        (is.sequential.incr(col.index) == FALSE) ) {
      stop("row.index and col.index must be sequentially increasing.")
    }

    ci <- c()
    ri <- rep(row.index[1], length(col.index))
    for( i in 2:length(row.index)) {
      ri.next <- rep(row.index[i], length(col.index))
      ri.new <- c(ri, ri.next)
      ri <- ri.new
    }

    ci <- rep(col.index, length(row.index))

    row.index <- ri
    col.index <- ci
  }

  # calculate exact coordinates for elliptical region
  if (index.type == "ellipse") {
    center <- c(row.index, col.index)
    vert <- h/2
    horz <- w/2
    ri <- c()
    ci <- c()
    for (i in 1:im.res[1]) {
      for (j in 1:im.res[2]) {
        if ( ((i - center[1])^2 / vert^2)
             + ((j - center[2])^2 / horz^2) <= 1 ) {
          ri <- c(ri, i)
          ci <- c(ci, j)
        }
      }
    }
    row.index <- ri
    col.index <- ci
  }

  if (index.type == "decay") {
    if (is.null(phi) == FALSE) {
      if (phi <= 0) {
        stop(paste0("phi must be > 0, but current value is ", phi))
      }
    }
    if (length(B.values) > 1) {
      stop("B.values should be of length 1 when using decay.")
    }

    Bmat <- matrix(NA, nrow = im.res[1], ncol = im.res[2])
    peak.coords <- c(row.index, col.index)
    row.index <- c()
    col.index <- c()

    for (i in 1:im.res[1]) {
      for (j in 1:im.res[2]) {
        dij <- sqrt((i - peak.coords[1]) ^ 2 + (j - peak.coords[2]) ^ 2)
        if (dij > max.d) {
          Bmat[i, j] <- 0
        } else {
          if (decay.fn == "exponential") {
            Bmat[i, j] <- B.values * exp(-phi * dij)
          }
          if (decay.fn == "gaussian") {
            Bmat[i, j] <- B.values * exp(-phi * dij ^ 2)
          }
          row.index <- c(row.index, i)
          col.index <- c(col.index, j)
        }
      }
    }

    B.values <- t(Bmat)[t(Bmat) > 0]
  }

  # ensure non-zero parameters go to the correct location of the
  # parameter vector.
  B.not0.m=matrix(c(row.index, col.index), ncol = 2)
  I = im.res[1]
  J = im.res[2]
  B.not0.v = c()
  for (i in 1:nrow(B.not0.m)) {
    B.not0.v[i] = J * (B.not0.m[i, 1] - 1) + B.not0.m[i, 2]
  }

  # expand B.values if all parameters have same value
  if ( (index.type %in% c("manual", "ellipse")) &
       (length(B.values) == 1) &
       (length(B.not0.v) > 1) ) {
    B.values <- rep(B.values, length(B.not0.v))
  }

  # Error if the number of parameter values is not equal the number
  # non-zero parameters
  if ( length(B.values) != length(B.not0.v) ) {
    stop("The number of parameter values must equal the number of non-zero locations.")
  } else {
    B = rep(0, prod(im.res))
    for (j in 1:length(B.not0.v)) {
      B[B.not0.v[j]] = B.values[j]
    }
    B = c(B0, B)
  }

  # bayesian parameter draws
  if (bayesian == TRUE) {
    # ensure valid distribution for parameters
    if (!(bayesian.dist %in% c("gaussian", "uniform"))) {
      stop("Distribution options are gaussian or uniform.")
    }
    # checks
    if (is.list(bayesian.scale) == FALSE) {
      stop("bayesian.scale should be a list. See documentation.")
    }
    if (length(bayesian.scale) != 2) {
      stop("bayesian.scale should be length 2. See documentation.")
    }
    if (!(bayesian.scale[[1]] %in% c("unique", "binary"))) {
      stop("bayesian.scale[[1]] options are unique and binary.")
    }

    if (bayesian.dist == "gaussian"){

      # unique sd's
      if (bayesian.scale[[1]] == "unique") {
        if (length(bayesian.scale[[2]]) != length(B)) {
          stop("When bayesian.scale[[1]] is unique, length(bayesian.scale[[2]] must be the same as length(B).")
        }
        sd.b <- bayesian.scale[[2]]
      }

      # 2 sd's - 1 for relevant and 1 for irrelevant parameters
      if (bayesian.scale[[1]] == "binary") {
        if (length(bayesian.scale[[2]]) != 3) {
          stop("When bayesian.scale[[1]] is binary, length(bayesian.scale[[2]] must be 3.")
        }
        sd.b <- rep(bayesian.scale[[2]][3], length(B))
        sd.b[B > 0] = bayesian.scale[[2]][2]
        sd.b[1] <- bayesian.scale[[2]][1]
      }

      # new vector for parameters
      B.bayes <- c()
      for (b in 1:length(B)) {
        B.bayes[b] <- rnorm(1, B[b], sd.b[b])
      }
    }

    if (bayesian.dist == "uniform"){

      # unique sd's
      if (bayesian.scale[[1]] == "unique") {
        if (length(bayesian.scale[[2]]) != length(B)) {
          stop("When bayesian.scale[[1]] is unique, length(bayesian.scale[[2]] must be the same as length(B).")
        }
        w.b <- bayesian.scale[[2]]
      }

      # 2 sd's - 1 for relevant and 1 for irrelevant parameters
      if (bayesian.scale[[1]] == "binary") {
        if (length(bayesian.scale[[2]]) != 3) {
          stop("When bayesian.scale[[1]] is binary, length(bayesian.scale[[2]] must be 3.")
        }
        width.b <- rep(bayesian.scale[[2]][3], length(B))
        width.b[B > 0] = bayesian.scale[[2]][2]
        width.b[1] <- bayesian.scale[[2]][1]
      }

      # new vector for parameters
      B.bayes <- c()
      for (b in 1:length(B)) {
        B.bayes[b] <- runif(1, B[b] - 0.5 * width.b[b], B[b] + 0.5 * width.b[b])
      }
    }

    B <- B.bayes
  }


  if (output.indices == FALSE){
    return(B)
  } else {
    out = list(Truth.Indices = B.not0.v, B = B)
    }
}


