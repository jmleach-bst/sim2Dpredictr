#' Determine rejections
#'
#' @inheritParams sample_FP_Power
#'
make_rejection <- function(B, reject.threshold, test.statistic) {
  rejections <- rep(0, length(B))
  if (reject.threshold[[2]] == "greater") {
    rejections[test.statistic > reject.threshold[[1]]] <- 1
  }
  if (reject.threshold[[2]] == "less") {
    rejections[test.statistic < reject.threshold[[1]]] <- 1
  }
  if (reject.threshold[[2]] == "2-tailed") {
    rt <- c(reject.threshold[[1]], -reject.threshold[[1]])
    rejections[test.statistic < min(rt) | test.statistic > max(rt)] <- 1
  }
  return(rejections)
}
