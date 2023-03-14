#' Determine and store neighbors by Euclidean Distance Constraints
#'
#' @param x,y are the row and column coordinates, respectively.
#' @param coords A dataframe containing indices and coordinates for the image.
#' @param im.res A vector containing the number of rows and columns, respectively.
#' @param r A scalar value determining the radius within which other locations
#' are neighbors to the current location (x, y).
#' @param print.ring When \code{print.ring = TRUE}, each iteration is shown,
#' with corresponding information regarding the number of neighbors present
#' in each ring. This argument primarily exists to allow the user to test
#' whether the neighborhood structure specified is as desired.
#' @return A tibble whose first column contains x indices, second column contains y indices, and
#' third column denotes the current ring about a location.
#' @note  This function avoids testing all points for being with a certain distance
#'  in order to determine neighbor status of a given point by progressively widening
#'  a box around the point. Each iteration widens the box by an extra ring, and we
#'  only test points in the new ring. If at the end of testing a ring there are no
#'  new neighbors then we stop expanding the box and return the neighbors' coordinates.
#'  For computational efficiency, this function assumes that all arguments except the
#'  current point's coordinates have been specified.
#' @examples
#' ## Necessary pre-specified arguments required for the function to work.
#'
#' ## image resoluation + number of spatial predictors
#' im.res <- c(5, 5)
#' J <- prod(im.res)
#'
#' ## create predictor indices w/ coordinates
#' row.id <-rep(1, im.res[2])
#' for (i in 2:im.res[1]) {
#'  row.id <- c(row.id, rep(i, im.res[2]))
#' }
#' coords <- data.frame(index = 1:J,
#'                      row.id = row.id,
#'                      col.id = rep(c(1:im.res[2]), im.res[1]) )
#'
#' neighbors_by_dist(x = 2, y = 2, im.res = im.res, coords = coords, r = 2)
#'
#' @export
neighbors_by_dist <- function(x, y, coords, im.res, r,
                              print.ring = FALSE) {

  # location index
  li <- coords$index[coords$row.id == x & coords$col.id == y]
  # 1st ring
  ring <- 1
  # x-coord for ring
  xn <- seq(max(x - ring, 1), min(x + ring, im.res[2]))
  # y-coord for ring
  yn <- seq(max(y - ring, 1), min(y + ring, im.res[1]))

  # initialize storage
  for (i in 1:length(xn)) {
    for (j in 1:length(yn)) {
      if (sqrt((xn[i] - x) ^ 2 + (yn[j] - y) ^ 2) <= r
          & !(xn[i] == x & yn[j] == y)  ) {
        nb.ij <- tibble::tibble(nb.x = xn[i],
                        nb.y = yn[j],
                        ring = ring)
        if (exists("nb.ij.old") == TRUE) {
          nb.ij.old <- rbind(nb.ij.old, nb.ij)
        } else {
          nb.ij.old <- nb.ij
        }
      }
    }
  }

  # only throw this error on ring 1.
  if (exists("nb.ij.old") == FALSE) {
    stop("r is too small; (x, y) are lonely because they have no neighbors.")
  }

  num.nb.ring <- nrow(nb.ij.old)
  nb.xy.old <- nb.ij.old

  if (print.ring == TRUE) {
    cat("The number of neighbors in ring ", ring, " is ", num.nb.ring, "\n")
    print(nb.ij.old)
  }

  # update ring
  ring <- 2

  # now iterate
  while (num.nb.ring > 0 & !(x - ring < 1 & x + ring > im.res[1]
                             & y - ring < 1 & y + ring > im.res[2])) {

    # top and bottom row y-coordinates
    tb.row.y <- seq(max(1, y - ring), min(im.res[2], y + ring ))

    # top and bottom row x-coordinates (if those rows exist)
    if (x - ring > 0) {
      top.row.x <- rep(x - ring, length(tb.row.y) )
      ring.coords <- tibble::tibble(nb.x = top.row.x,
                            nb.y = tb.row.y)
    }
    if (x + ring <= im.res[1]) {
      bottom.row.x <- rep(x + ring, length(tb.row.y) )
      if (exists("ring.coords") == TRUE) {
        ring.coords <- rbind(ring.coords,
                             tibble::tibble(nb.x = bottom.row.x,
                                    nb.y = tb.row.y))
      } else {
        ring.coords <- tibble::tibble(nb.x = bottom.row.x,
                              nb.y = tb.row.y)
      }
    }

    # generate left/right columns of ring (if they exist)

    # column x-coordinates are the same for left/right
    c.x <- seq(max(1, x - ring + 1), min(im.res[1], x + ring - 1))

    # left column y-coordinates
    if (y - ring > 0) {
      lc.y <- rep(y - ring, length(c.x))

      if (exists("ring.coords") == TRUE) {
        ring.coords <- rbind(ring.coords,
                             tibble::tibble(nb.x = c.x,
                                    nb.y = lc.y))
      } else {
        ring.coords <- tibble::tibble(nb.x = c.x,
                              nb.y = lc.y)
      }
    }

    # right column y-coordinates
    if (y + ring <= im.res[2]) {
      rc.y <- rep(y + ring, length(c.x))

      if (exists("ring.coords") == TRUE) {
        ring.coords <- rbind(ring.coords,
                             tibble::tibble(nb.x = c.x,
                                    nb.y = rc.y))
      } else {
        ring.coords <- tibble::tibble(nb.x = c.x,
                              nb.y = rc.y)
      }
    }

    # evaluate whether euclidean distance between points <= r
    nb.ring <- rep(0, nrow(ring.coords))
    for (d in 1:length(nb.ring)) {
      if(sqrt((ring.coords$nb.x[d] - x) ^ 2 + (ring.coords$nb.y[d] - y) ^ 2) <= r) {
        nb.ring[d] <- 1
      }
    }

    if (sum(nb.ring) == 0){
      num.nb.ring <- 0
    } else {
      nb.ij <- ring.coords[nb.ring == 1, ]
      nb.ij$ring = rep(ring, nrow(nb.ij))
      nb.ij.old <- rbind(nb.ij.old, nb.ij)
      num.nb.ring <- sum(nb.ring)
    }

    if (print.ring == TRUE) {
      cat("The number of neighbors in ring ", ring, " is ", num.nb.ring, "\n")
      print(nb.ij.old)
    }

    # update ring
    ring <- ring + 1
    remove(ring.coords)

  }
  return(nb.ij.old)
}
