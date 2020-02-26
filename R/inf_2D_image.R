#' Display Inference Results for 2D Predictors
#'
#' Provide graphics for spatial extent of predictor parameters, rejections
#' and/or the truth/falsity of the rejections.
#'
#' @return An image depicting the spatial extent of some image characteristic.
#'
#' @param im.res A vector defining the dimension of spatial data. The first entry is the
#' number of rows and the second  entry is the number of columns.
#' @inheritParams sample_FP_Power
#' @param grid.color Specify the color for the grid lines.
#' @param binarize.B Either \code{TRUE} (default) or \code{FALSE}. When
#' \code{binarize.B = TRUE} the parameter vector is converted to a binary vector
#' where 1 indicates non-zero parameter and 0 indicates zero-valued parameter.
#' @param n.colors Determines the number of colors in the printed image. Default is
#' \code{length(unique(B))}, but it is recommended to use trial and error to determine
#' the ideal setting for specific situations.
#' @param plot.title When \code{plot.title = TRUE} a title accompanies the output graph, and
#' \code{plot.title = FALSE} suppresses the title.
#' @importFrom grDevices terrain.colors
#' @importFrom graphics box grid image legend
#' @note If both \code{rejections} and \code{B} are specified then the function
#' provides an image with separate color each for:
#' \itemize{
#'     \item No rejection and \code{B[i] = 0} (i.e. True Negative).
#'     \item No rejection and \code{B[i] != 0} (i.e. False Negative).
#'     \item Rejection and \code{B[i] = 0} (i.e. False Positive).
#'     \item Rejection and \code{B[i] != 0} (i.e. True Positive).
#' }
#' @examples
#' ## parameter vector
#' Bex <- beta_builder(row.index = c(rep(5, 3), rep(6, 3), rep(7, 3)),
#'                     col.index = rep(c(5, 6, 7), 3),
#'                     im.res = c(10, 10), index.type = "manual",
#'                     B0 = 0, B.values = 1:9,
#'                     output.indices = FALSE)
#'
#' ## co-opt beta builder to get rejections
#' rejex <- beta_builder(row.index = c(rep(4, 3), rep(5, 3), rep(6, 3)),
#'                       col.index = rep(c(4, 5, 6), 3),
#'                       im.res = c(10, 10), index.type = "manual",
#'                       B0 = 0, B.values = rep(1, 9),
#'                       output.indices = FALSE)[-1]
#'
#' rejex.sm2 <- beta_builder(row.index = 5:6, col.index = 5:6,
#'                           im.res = c(10, 10),
#'                           B0 = 0, B.values = 1,
#'                           output.indices = FALSE)[-1]
#'
#' ## just B
#' inf_2D_image(B = Bex, im.res = c(10, 10))
#' ## just rejections
#' inf_2D_image(rejections = rejex, im.res = c(10, 10))
#'
#' ## both B and rejections
#' inf_2D_image(rejections = rejex, B = Bex, im.res = c(10, 10))
#' inf_2D_image(rejections = rejex.sm2, B = Bex, im.res = c(10, 10))
#'
#' ## larger dimension example
#' Bex2 <- beta_builder(row.index = 5:15, col.index = 16:20,
#'                      im.res = c(50, 50), B0 = 0,
#'                      B.values = 1:(length(5:15) * length(16:20)),
#'                      index.type = "rectangle",
#'                      output.indices = FALSE)
#' rejex2 <- beta_builder(row.index = 13:21, col.index = 30:41,
#'                        im.res = c(50, 50), B0 = 0,
#'                        B.values = rep(1, (length(13:21) * length(30:41))),
#'                        index.type = "rectangle",
#'                        output.indices = FALSE)[-1]
#' rejex3 <- beta_builder(row.index = 5:20, col.index = 16:30,
#'                        im.res = c(50, 50), B0 = 0,
#'                        B.values = rep(1, (length(5:20) * length(16:30))),
#'                        index.type = "rectangle",
#'                        output.indices = FALSE)[-1]
#' rejex4 <- beta_builder(row.index = 5:10, col.index = 16:17,
#'                        im.res = c(50, 50), B0 = 0,
#'                        B.values = rep(1, (length(5:10) * length(16:17))),
#'                        index.type = "rectangle",
#'                        output.indices = FALSE)[-1]
#' ## images
#' inf_2D_image(B = Bex2, im.res = c(50, 50))
#' inf_2D_image(B = Bex2, im.res = c(50, 50), binarize.B = FALSE)
#' inf_2D_image(rejections = rejex2, im.res = c(50, 50))
#'
#' ## No TP
#' inf_2D_image(rejections = rejex2, B = Bex2, im.res = c(50, 50))
#' ## ALL TP
#' inf_2D_image(rejections = Bex2[-1], B = Bex2, im.res = c(50, 50))
#' ## No FN
#' inf_2D_image(rejections = rejex3, B = Bex2, im.res = c(50, 50))
#' ## No FP, but FN
#' inf_2D_image(rejections = rejex4, im.res = c(50, 50))
#' inf_2D_image(B = Bex2, im.res = c(50, 50))
#' inf_2D_image(rejections = rejex4, B = Bex2, im.res = c(50, 50))
#' @export
inf_2D_image <- function(rejections = NULL, B = NULL, im.res,
                         test.statistic = NULL, reject.threshold = NULL,
                         binarize.B = TRUE, grid.color = "grey",
                         n.colors = length(unique(B)),
                         B.incl.B0 = TRUE, plot.title = TRUE) {

  # otherwise the images are screwy-looking
  rotate = function(x){
    t(apply(x, 2, rev))
  }

  if ( (B.incl.B0 == TRUE) & (is.null(B) == FALSE) ) {
    B = B[-1]
  }

  if (is.null(rejections) == TRUE &
      (is.null(test.statistic) == FALSE & is.null(reject.threshold) == FALSE)) {
    rejections <- make_rejection(B = B, reject.threshold = reject.threshold,
                                 test.statistic = test.statistic)
  }

  if( is.null(rejections) == FALSE ) {
    rej.mat <- matrix(rejections,
                      nrow = im.res[1],
                      ncol = im.res[2],
                      byrow = TRUE)
    if ( is.null(B) == TRUE ) {
      if (plot.title == FALSE) {
        image(rotate(rej.mat), col = c("white", "#CC79A7"), axes = FALSE)
      } else {
        image(rotate(rej.mat), col = c("white", "#CC79A7"), axes = FALSE,
              main = "Pink Indicates Rejected Locations")
      }
      box()
      grid(nx = im.res[1], ny = im.res[2],
           col = grid.color, lty = 1)
    }
  }


  if( is.null(B) == FALSE ) {
    if (binarize.B == TRUE) {
      B[B > 0] <- 1
    }
    B.mat <- matrix(B, nrow = im.res[1],
                    ncol = im.res[2],
                    byrow = TRUE)
    if ( is.null(rejections) == TRUE) {
      if (length(unique(B)) < 3) {
        if (plot.title == FALSE) {
          image(rotate(B.mat), col = c("white", "#56B4E9"), axes = FALSE)
        } else {
          image(rotate(B.mat), col = c("white", "#56B4E9"), axes = FALSE,
                main = "Blue Indicates Non-Zero Parameter Locations")
        }
        box()
        grid(nx = im.res[1], ny = im.res[2],
             col = grid.color, lty = 1)
      }
      else{
        image(rotate(B.mat),
              col = c("white", rev(terrain.colors(n = n.colors))),
                      axes = FALSE)
        box()
        grid(nx = im.res[1], ny = im.res[2],
             col = grid.color, lty = 1)
      }
    }
  }

  if( (is.null(rejections) == FALSE) & (is.null(B) == FALSE) ) {
    im.mat <- matrix(0, nrow = im.res[1],
                     ncol = im.res[2],
                     byrow = TRUE)
    # FN
    im.mat[(rej.mat == 0) & (B.mat != 0)] <- 1
    # FP
    im.mat[(rej.mat != 0) & (B.mat == 0)] <- 2
    # TP
    im.mat[(rej.mat != 0) & (B.mat != 0)] <- 3

    if (length(im.mat[im.mat == 3]) == 0 ) {
      im.col <- c("white", "#E69F00", "#D55E00")
    }
    else {
      im.col <- c("white", "#E69F00", "#D55E00", "#009E73")
    }
    if (plot.title == FALSE) {
      image(rotate(im.mat), col = im.col, axes = FALSE)
    } else {
      image(rotate(im.mat), col = im.col,
            axes = FALSE, main = "Hypothesis Testing Results")
    }
    box()
    grid(nx = im.res[1], ny = im.res[2], col = "black", lty = 1)
    legend("right", legend = c("TN", "FN", "FP", "TP"),
           fill = c("white", "#E69F00", "#D55E00", "#009E73"))
  }


}

