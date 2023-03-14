#' Convert a 2D Space to Grid Coordinates
#'
#' Input the limits of a 2D space and the desired image resolution, then
#' the function outputs the appropriate grid/lattice coordinates.
#'
#' @param im.res A vector specifying the dimension/resolution of the image. The first entry is
#'  the number of 'rows' in the lattice/image, and the second entry is the number of
#'  'columns' in the lattice/image.
#' @param xlim,ylim These are the 2D image limits. Defaults for both are \code{c(0, 1)}.
#' It is not recommended to alter these arguments unless changing the limits has a
#' specific practical utility.
#' @return A data frame whose first column is x-coordinates and whose second column is y-coordinates.
#' @export
generate_grid <- function(im.res, xlim = c(0, 1), ylim = c(0, 1)) {
  x.length <- xlim[2] - xlim[1]
  y.length <- ylim[2] - ylim[1]
  x.cell <- seq(xlim[1], xlim[2], x.length / im.res[2])
  y.cell <- seq(ylim[2], ylim[1], -y.length / im.res[1])
  x.coord <- c()
  for (i in 1:(length(x.cell) - 1)){
    x.coord[i] <- mean(c(x.cell[i], x.cell[i + 1]))
  }
  y.coord <- c()
  for (i in 1:(length(y.cell) - 1)){
    y.coord[i] <- mean(c(y.cell[i], y.cell[i + 1]))
  }
  x <- rep(x.coord, length(y.coord))
  y <- c()
  old.y <- 1
  for (i in 1:length(y.coord)) {
    y[old.y:(old.y + length(x.coord) - 1)] <- y.coord[i]
    old.y <- old.y + length(x.coord)
  }
  xy <- data.frame(cbind(x, y))
  
  return(xy[order(xy$x, xy$y), ])
}