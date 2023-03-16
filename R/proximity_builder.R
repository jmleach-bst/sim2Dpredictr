#' Generate a Proximity Matrix
#'
#' Generates a proximity matrix where non-zero entries are the weights 
#' associated with neighbors, and zero entries are not neighbors.
#'
#' @param im.res A vector defining the dimension of spatial data. The first 
#' entry is the number of rows and the second  entry is the number of columns.
#' @param weight Determines how weights are assigned. \code{"distance"} 
#' assigns weights as the inverse of Euclidean distance times a constant, 
#' \code{phi}. \code{"binary"} assigns weights to 1 for neighbors and 0 
#' otherwise. 
#' @param phi When \code{weight = "distance"} a constant by which
#' to multiply the inverse of Euclidean distance. Defaults to 1.
#' @param neighborhood Determines how to assign neighbor status to locations;
#' i.e. 1 for neighbors, 0 otherwise. \code{type = "round"} assigns neighbor 
#' status to locations within radius \code{r}. \code{type = "ar1"} assigns
#' 1 to locations directly above or beside. \code{type = "rectangle"} assigns 
#' neighbor status to locations within \code{w} units to the left or right 
#' and  \code{h} units up or down.
#' @param type Specifies either sparse (\code{type = "sparse"}) or full 
#' (\code{type = "full"}) proximity matrix.
#' @param r,h,w When \code{neighborhood = "round"}, \code{r} specifies the 
#' radius within which other locations are neighbors. When 
#' \code{neighborhood = "rectangle"}, \code{w} and \code{h} specify the number
#' of units to the left/right and above/below the location are to be counted 
#' as neighbors.
#' @param include.coords If \code{type = "sparse"} and 
#' \code{include.coords = TRUE}, then the coordinates of neighbors are
#' returned along with their indices.
#' @param print.im Allows user to print the 2D "image" matrix with index 
#' labels to visually verify that the proximity matrix is as expected.
#' @return A (proximity) matrix.
#' @importFrom dplyr mutate
#' @examples
#' ## adjacency matrix with sparse structure (i.e., 2 columns) 
#' ## and ar1 neighborhood
#' sp.ar1 <- proximity_builder(im.res = c(3, 3),
#'                             weight = "binary",
#'                             neighborhood = "ar1",
#'                             type = "sparse")
#' ## adjacency matrix with full structure 
#' ## (i.e., prod(im.dim) rows & columns) and ar1 neighborhood
#' full.ar1 <- proximity_builder(im.res = c(3, 3),
#'                               weight = "binary",
#'                               neighborhood = "ar1",
#'                               type = "full")
#'
#' ## proximity matrix weighted by distance (sparse)
#' sp.rnd <- proximity_builder(im.res = c(3, 3),
#'                             weight = "distance",
#'                             neighborhood = "round", r = 2,
#'                             type = "sparse",
#'                             include.coords = TRUE)
#'
#' ## proximity matrix weighted by distance (full)
#' full.rnd <- proximity_builder(im.res = c(3, 3),
#'                               weight = "distance",
#'                               neighborhood = "round", r = 2,
#'                               type = "full")
#' @export
proximity_builder <- function(im.res, neighborhood = "ar1",
                              type = c("sparse", "full"),
                              weight = "binary", phi = 1,
                              r = NULL, h = NULL, w = NULL,
                              include.coords = FALSE,
                              print.im = FALSE) {
  # bad hack... but otherwise no "global binding"...
  nb.y <- c()

  # real stuff we need
  J <- prod(im.res)
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

  # determine neighbors
  for (x in 1:im.res[1]) {
    for (y in 1:im.res[2]) {

      # neighbors are just 1 unit up/down, left/right of location [x, y]
      if (neighborhood == "ar1") {
        nb.xy <- data.frame(nb.x = c(x - 1, x + 1, x , x),
                            nb.y = c(y, y, y - 1, y + 1))
      }

      # neighbors are locations with `w` units left/right and `h` 
      # units up/down
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

      nb.xy <- dplyr::filter(nb.xy, 
                             dplyr::between(nb.x, 1, im.res[1]) & 
                               dplyr::between(nb.y, 1, im.res[2]))

      nb.index <- coords$index[coords$row.id == nb.xy$nb.x[1] & 
                                 coords$col.id == nb.xy$nb.y[1]]
      for (i in 2:nrow(nb.xy)) {
        nb.index <- c(nb.index, 
                      coords$index[coords$row.id == nb.xy$nb.x[i] & 
                                     coords$col.id == nb.xy$nb.y[i]])
      }
      nb.xy <- cbind(location.index = rep(
        coords$index[coords$row.id == x & coords$col.id == y], 
        nrow(nb.xy)),
        x = x, y = y,
        nb.xy, nb.index)
      if (x == 1 & y == 1) {
        nb.xy.old <- nb.xy
      } else {
        nb.xy.old <- rbind(nb.xy.old, nb.xy)
      }
    }
  }

  # determine weights
  if (weight == "distance") {
    nb.xy.old <- dplyr::mutate(
      nb.xy.old,
      weights = phi / sqrt((x - nb.x) ^ 2 + (y - nb.y) ^ 2)
    )
  } else {
    nb.xy.old$weights <- 1
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
      adj.mat[nb.xy.old$location.index[i], 
              nb.xy.old$nb.index[i]] <- nb.xy.old$weights[i]
    }
    return(adj.mat)
  }

  if (print.im == TRUE) {
    im <- matrix(1:J, byrow = TRUE,
                 nrow = im.res[1], ncol = im.res[2])
    print(im)
  }

}
