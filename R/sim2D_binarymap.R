#' Generate a Binary Map via the Boolean Method
#'
#' Use a Homogenous Poisson Process to generate random "events", a uniform distribution
#' to generate circles of random radii about the events, and take the union to obtain
#' a random set. This is mapped onto a lattice to obtain a binary map.
#'
#' @inheritParams sim2D_RandSet_HPPP
#' @param im.res A vector specifying the dimension/resolution of the image. The first entry is
#' the number of 'rows' in the lattice/image, and the second entry is the number of
#' 'columns' in the lattice/image.
#' @param store.type One of \code{c("list", "matrix", "all")}. When \code{store.type = "list"},
#' the output is a list where each element is a matrix defining a subject image. If
#' \code{store.type = "matrix"}, then the images are vectorized by row and each row
#' of the output matrix contains an image vector for a single subject. 
#' @param output.randset Logical. When \code{TRUE}, stores the data frame of original draws from
#' the HPPP and and random radii from \code{sim2D_RandSet_HPPP()}. This data frame is stored in the
#' first element of the output list named \code{randset}. The second element of the output list is a 
#' list/matrix of the final subject images depending on \code{store.type} and named \code{images}. 
#' @return A list; each element is a matrix of zeroes and ones.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Cressie+Wikle:2011}{sim2Dpredictr}
#' @examples
#' bin_ims <- sim2D_binarymap(N = 5, im.res = c(10, 10), store.type = "list",
#'                            lambda = 50, sub.area = TRUE,
#'                            min.sa = c(0.10, 0.10), max.sa = c(0.5, 0.5),
#'                            radius.bounds.min.sa = c(0.015, 0.04),
#'                            radius.bounds.max.sa = c(0.041, 0.06))
#'
#' rotate = function(x){
#'   t(apply(x, 2, rev))
#' }
#'
#' for (i in 1:length(bin_ims)) {
#'   image(rotate(bin_ims[[i]]),
#'         col = c("white", "darkgreen"),
#'         axes = FALSE)
#'   box()
#'   grid(nx = 10, ny = 10, col = "black",
#'        lty = 1)
#' }
#' @export
sim2D_binarymap <- function(N, xlim = c(0, 1), ylim = c(0, 1), im.res,
                            radius.bounds = c(0.02, 0.1), lambda = 50,
                            random.lambda = FALSE, lambda.sd = 10,
                            lambda.bound = NULL,
                            prior = "gamma", sub.area = FALSE,
                            min.sa = c(0.1, 0.1), max.sa = c(0.3, 0.3),
                            radius.bounds.min.sa = c(0.02, 0.05),
                            radius.bounds.max.sa = c(0.08, 0.15),
                            print.subj.sa = FALSE, print.lambda = FALSE,
                            print.iter = FALSE,
                            store.type = "list", output.randset = FALSE) {

  # generate the random set
  rs <- sim2D_RandSet_HPPP(xlim = xlim, ylim = ylim, N = N, radius.bounds = radius.bounds,
                           lambda = lambda, lambda.sd = lambda.sd, lambda.bound = lambda.bound,
                           random.lambda = random.lambda, prior = prior, sub.area = sub.area,
                           min.sa = min.sa, max.sa = max.sa,
                           radius.bounds.min.sa = radius.bounds.min.sa,
                           radius.bounds.max.sa = radius.bounds.max.sa,
                           print.subj.sa = print.subj.sa, print.lambda = print.lambda,
                           print.iter = print.iter)

  # generate grid/lattice
  g <- generate_grid(im.res = im.res, xlim = xlim, ylim = ylim)

  # determine which lattice points are within the random set (subject-wise)
  if (store.type == "list") {
    out <- list()
  } else if (store.type == "matrix") {
    out <- matrix(NA, nrow = N, ncol = prod(im.res))
  }

  for (i in 1:N) {
    # binary map information
    ws <- within_area(grid.centers = g, radii = rs$radii[rs$subj == i],
                      event.xcoord = rs$xcoord[rs$subj == i],
                      event.ycoord = rs$ycoord[rs$subj == i])
    if (store.type == "list") {
      out[[i]] <- matrix(ws$in.set, byrow = TRUE, nrow = im.res[1])
    } else if (store.type == "matrix") {
      out[i, ] <- ws$in.set
    }
  }
  
  if (output.randset == TRUE) {
    return(list(randset = rs, images = out))
  } else {
    return(out)
  }
}
