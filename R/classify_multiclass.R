#' Classify subjects based on predicted probabilities for each class
#' 
#' TBD
#' 
#' @param predicted.probs A matrix where the number of rows is equal to the number of subjects and the 
#' number of columns equals the number of categories. \code{predicted.probs[i, j]} contains the probability 
#' that subject \code{i} belongs to category \code{j}. 
#' @param category.names A vector containing the names of each category. The order of names should match 
#' the order of columns in \code{predicted.probs}; correspondingly, the length of the vector should equal 
#' the number of columns in \code{predicted.probs}.
#' @param keep.probs Logical. When \code{TRUE}, the output is data frame consisting of the information in 
#' \code{predicted.probs} with an additional column \code{predicted.class} that contains the predicted class 
#' for each subject. When \code{FALSE}, a vector of the predicted classes is returned.  
#' @inheritParams generate_multinom_probs 
#' @return Depending on the option selected for \code{keep.probs}, returns a data frame or vector.
#' @details TBD
#' @examples 
#' ## number of categories
#' vt <- 3
#'
#' ## covariate values
#' xt <- matrix(rnorm(10 * 2), ncol = 2, nrow = 10)
#' 
#' ## list of parameter vectors (over-parameterized model)
#' bu <- list(b1 = c(0, 0.25, 0.25),
#'            b2 = c(0, -0.25, -0.25),
#'            b3 = c(0, 0.25, -0.25))
#'
#' ## subject specific probabilities for each category (over-parameterized model)
#' prp <- generate_multinom_probs(V = vt, X = xt, B = bu)
#' 
#' classify_multiclass(predicted.probs = prp, category.names = c("A", "B", "C"))
#' 
#' ## generate predicted probabilities within function
#' classify_multiclass(predicted.probs = NULL, category.names = c("A", "B", "C"),
#'                     X = xt, B = bu)
#'                     
#' @export
classify_multiclass <- function(predicted.probs = NULL, category.names, keep.probs = TRUE,
                                B = NULL, X = NULL, X.incl.X0 = FALSE) {
  
  if (is.null(predicted.probs) == FALSE) {
    if (is.null(B) == FALSE | is.null(X) == FALSE) {
      warning("predicted.probs and (B and/or X) both supplied. Ignoring B and/or X.")
    }
  }
  
  if (is.null(predicted.probs) == TRUE) {
    if (is.null(B) == TRUE | is.null(X) == TRUE) {
      stop("Use must supply values for B and X.")
    }
    V <- length(category.names)
    predicted.probs <- generate_multinom_probs(
      V = V, B = B, X = X,
      X.incl.X0 = X.incl.X0
    )
  }
  
  # check that the number of categories matches everywhere
  if (ncol(predicted.probs) != length(category.names)) {
    stop("ncol(predicted.probs) must equal length(category.names)")
  }
  
  predicted.class <- c()
  for (i in 1:nrow(predicted.probs)) {
    predicted.class[i] <- category.names[which(
      predicted.probs[i, ] == max(predicted.probs[i, ])
      )]
  }
  
  df <- data.frame(
    cbind(predicted.probs,
          predicted.class)
  )
  colnames(df) <- c(category.names, 
                    "predicted.class")
  
  if (keep.probs == TRUE) {
    return(df)
  } else {
    return(df$predicted.class)
  }
}







