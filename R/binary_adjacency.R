#' Generate an Adjacency Matrix
#'
#' Generates an adjacency matrix where neighbors are weighted 1 and all others are weighted 0.
#'
#' @param im.res A vector defining the dimension of spatial data. The first entry is the
#' number of rows and the second  entry is the number of columns.
#' @param neighborhood Determines how to assign neighbor status to locations; i.e. 1 for neighbors, 0 otherwise.
#'  \code{type = "round"} assigns neighbor status to locations within radius \code{r}. \code{type = "ar1"}
#'  assigns 1 to locations directly above or beside. \code{type = "rectangle"} assigns neighbor status to
#'  locations within \code{w} units to the left or right and  \code{h} units up or down.
#' @param type Specifies either sparse (\code{type = "sparse"}) or full (\code{type = "full"}) adjacency matrix.
#' @param r,h,w When \code{neighborhood = "round"}, \code{r} specifies the radius within which other locations
#'  are neighbors. When \code{neighborhood = "rectangle"}, \code{w} and \code{h} specify the number of units to
#'  the left/right and above/below the location are to be counted as neighbors.
#' @param include.coords If \code{type = "sparse"} and \code{include.coords = TRUE}, then the coordinates of
#'  neighbors are returned along with their indices.
#' @param print.im Allows user to print the 2D "image" matrix with index labels to visually verify that the adjacency
#'   matrix is as expected.
#' @examples
#' ## adjacency matrix with sparse structure (i.e., 2 columns) and ar1 neighborhood
#' sp.ar1 <- binary_adjacency(im.res = c(4, 4),
#'                            neighborhood = "ar1",
#'                            type = "sparse")
#' ## adjacency matrix with full structure (i.e., prod(im.dim) rows & columns) and ar1 neighborhood
#' full.ar1 <- binary_adjacency(im.res = c(4, 4),
#'                              neighborhood = "ar1",
#'                              type = "full")
#'
#' ## adjacency matrix with sparse structure (i.e., 2 columns) and rectangle neighborhood
#' sp.rec <- binary_adjacency(im.res = c(4, 4),
#'                            neighborhood = "rectangle",
#'                             w = 1, h = 1,
#'                             type = "sparse")
#'
#' ## adjacency matrix with sparse structure (i.e., 2 columns) and round neighborhood
#' sp.round <- binary_adjacency(im.res = c(5, 5),
#'                              neighborhood = "round",
#'                              r = 2, type = "sparse")
#' @export
binary_adjacency <- function(im.res, neighborhood = "ar1",
                             type = c("sparse", "full"),
                             r = NULL, h = NULL, w = NULL,
                             include.coords = FALSE,
                             print.im = FALSE){
  .Deprecated("proximity_builder")
  # bad hack... but otherwise no "global binding"...
  nb.y <- c()

  # real stuff we need
  J = prod(im.res)
  row.id <-rep(1, im.res[2])
  for (i in 2:im.res[1]) {
    row.id <- c(row.id, rep(i, im.res[2]))
  }
  coords <- data.frame(index = 1:J,
                       row.id = row.id,
                       col.id = rep(c(1:im.res[2]), im.res[1]) )
  # 4 corners have 2 neighbors,
  # top/bottom row, left-/right-most columns have 3 neighbors,
  # remaining locations have 4 neighbors.
  # Currently un-used -> may include if current code is too slow.
  sparse.adj <- matrix(0, nrow = 4 * 2 +
                         3 * (2 * im.res[1] + 2 * im.res[2]) +
                         4 * (J - 4 + 2 * im.res[1] + 2 * im.res[2]),
                       ncol = 2)

  for (x in 1:im.res[1]) {
    for (y in 1:im.res[2]) {

      # neighbors are just 1 unit up/down, left/right of location [x, y]
      if (neighborhood == "ar1") {
        nb.xy <- data.frame(nb.x = c(x - 1, x + 1, x , x),
                            nb.y = c(y, y, y - 1, y + 1))
      }

      # neighbors are locations with `w` units left/right and `h` units up/down
      if (neighborhood == "rectangle") {
        # x-coord for neighbors
        xn <- seq(x - w, x + w)
        # y-coord for neighbors
        yn <- seq(y - h, y + h)

        # generate (x, y) neighbor coordinates
        nb.x <- rep(xn[1], length(yn))
        for (i in 2:length(xn)) {
          nb.x <- c(nb.x, rep(xn[i], length(yn)))
        }

        nb.xy <- data.frame(cbind(nb.x, nb.y = rep(yn, length(xn))))
        nb.xy <- dplyr::filter(nb.xy, !(nb.x == x & nb.y == y))
      }

      # neighbors are locations within radius r of location [x, y]
      if (neighborhood == "round") {
        nb.xy <- neighbors_by_dist(x = x, y = y, r = r,
                                   coords = coords,
                                   im.res = im.res)
        nb.xy <- nb.xy[, -3]
      }

      nb.xy <- dplyr::filter(nb.xy, dplyr::between(nb.x, 1, im.res[1]) & dplyr::between(nb.y, 1, im.res[2]))

      nb.index <- coords$index[coords$row.id == nb.xy$nb.x[1] & coords$col.id == nb.xy$nb.y[1]]
      for (i in 2:nrow(nb.xy)) {
        nb.index <- c(nb.index, coords$index[coords$row.id == nb.xy$nb.x[i] & coords$col.id == nb.xy$nb.y[i]])
      }
      nb.xy <- cbind(location.index = rep(coords$index[coords$row.id == x & coords$col.id == y], nrow(nb.xy)),
                     x = x, y = y,
                     nb.xy, nb.index)
      if (x == 1 & y == 1) {
        nb.xy.old <- nb.xy
      } else {
        nb.xy.old <- rbind(nb.xy.old, nb.xy)
      }
    }
  }

  # determine whether to return sparse or full matrix
  if (type == "sparse") {
    if (include.coords == FALSE) {
      nb.xy.old <- nb.xy.old[ , -c(2:5)]
    }
    return(nb.xy.old)
  } else {
    adj.mat <- matrix(0, nrow = J, ncol = J)
    for (i in 1:nrow(nb.xy.old)) {
      adj.mat[nb.xy.old$location.index[i], nb.xy.old$nb.index[i]] <- 1
    }
    return(adj.mat)
  }

  if (print.im == TRUE) {
    im <- matrix(1:J, byrow = TRUE,
                 nrow = im.res[1], ncol = im.res[2])
    print(im)
  }

}

