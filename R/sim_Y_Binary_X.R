#' Simulate Scalar Outcomes from Simulated Spatially Dependent Binary Predictors
#'
#' N spatially dependent binary design vectors are simulated using \code{sim2D_binarymap}.
#' These design vectors are used to then simulate scalar outcomes that have
#' one of Gaussian, Binomial, or Poisson distributions.
#' @inheritParams sim2D_binarymap
#' @param B A vector parameter values; i.e. "betas". Note that \code{length(B)}
#' must equal \code{p + 1 = n.row * n.col + 1}; e.g. for normal outcomes
#' \eqn{Y = XB + e} with \code{Y} a scalar outcome and \code{e} the random error.
#' @param rand.err A scalar for the random error variance when \code{dist = "gaussian"}.
#' @param dist The distribution of the scalar outcome.
#' \itemize{
#'     \item \code{dist = "gaussian"} has \eqn{Y = XB + e}, where
#'     \eqn{e ~ N(0, rand.err)}.
#'     \item \code{dist = "binomial"} is drawn from eqn{Bin(XB, XB(1-XB))}
#'     using \code{rbinom()} when \code{binary.method = "Traditional"}. If
#'     \code{binary.method = "Gaussian"}, then simulation is based on a
#'     cutoff using \code{binary.cutoff}.
#'     \item \code{dist = "poisson"} is drawn from \eqn{Poisson(XB)} using
#'     \code{rpois()}.
#' }
#' @param binomial.method One of \code{c("traditional", "gaussian manual", "gaussian percentile")}.
#' Only specified when \code{dist = "binomial"}, and determines
#' whether draws are directly taken from a binomial distribution or if draws are
#' taken from a Multivariate Normal Distribution (in the manner of \code{dist = "gaussian"})
#' and thresholds imposed to binarize the outcomes. \code{binomial.method = "gaussian manual"}
#' allows the user to specify specific values for categorizing outcomes. \code{binomial.method = "gaussian percentile"}
#' allows the user to specify percentiles for binarizing the data. Both approaches use
#' \code{Y.thresh} to specify the cutoff value(s). If \code{binomial.method = "gaussian percentile"}
#' and \code{Y.thresh = NULL} then the median is used as the threshold. If
#' \code{binomial.method = "gaussian manual"} and \code{Y.thresh = NULL}, then 0 is used as
#' the threshold. Default is \code{binomial.method = "traditional"}.
#' @param count.method One of \code{c("traditional", "rounding")}. When \code{count.method = "traditional"},
#' the outcomes are drawn sequentially using \code{rpois()}. When \code{count.method = "traditional"},
#' the outcomes are drawn from an MVN, then values less than or equal to 0 are set to 0, and all other
#' values are rounded to the nearest whole number.
#' @param Y.thresh When \code{binomial.method = "traditional"}
#' @param incl.subjectID When \code{incl.subjectID = TRUE} a column of subject indices
#' is generated.
#' \code{Y.thresh = NULL} (default). If \code{binomial.method = "gaussian manual"},
#' then \code{Y.thresh} should be any scalar real number; values equal or above
#'  this cutoff are assigned 1 and values below are assigned 0.
#'  If \code{binomial.method = "gaussian percentile"}, then values equal or above
#'  this percentile are assigned 1, and other wise 0; in this case values should
#'  be between 0 and 1. For example, if \code{Y.thresh = 0.9}, then the
#'  cutoff is the 90th percentile.
#' @param print.out If \code{print.out = TRUE} then print the following for
#'  each subject, indexed y: \itemize{
#'      \item \code{X[y] \%*\% B}
#'      \item \code{p[y]}, \code{lambda[y]} for Binomial, Poisson, respectively.
#'  }
#' This is useful to see the effect of image parameter selection and beta
#' parameter selection on distributional parameters for the outcome of interest.
#' @note Careful parameter selection, i.e. \code{B}, is necessary to ensure that
#' simulated outcomes are reasonable; in particular, counts arising from the Poisson
#' distribution can be unnaturally large.
#' @examples
#'
#' ## Define non-zero beta values
#' Bex <- beta_builder(row.index = c(3, 3, 4, 4), col.index = c(3, 4, 3, 4),
#'                     im.res = c(5, 5),
#'                     B0 = 0, B.values = rep(1, 4),
#'                     output.indices = FALSE)
#' ## Simulate Datasets
#' ## parameter values
#' Nex = 100
#' set.seed(28743)
#'
#' Gauss.ex <- sim_Y_Binary_X(N = Nex, B = Bex, dist = "gaussian", im.res = c(5, 5))
#' hist(Gauss.ex$Y)
#'
#' ## direct draws from binomial
#' Bin.ex <- sim_Y_Binary_X(N = Nex, B = Bex, im.res = c(5, 5),
#'                          dist = "binomial", print.out = TRUE)
#' table(Bin.ex$Y)
#'
#' ## manual cutoff
#' Bin.ex2 <- sim_Y_Binary_X(N = Nex, B = Bex, im.res = c(5, 5),
#'                           dist = "binomial",
#'                           binomial.method = "gaussian manual",
#'                           Y.thresh = 1.25)
#' table(Bin.ex2$Y)
#'
#' ## percentile cutoff
#' Bin.ex3 <- sim_Y_Binary_X(N = Nex, B = Bex, im.res = c(5, 5),
#'                           dist = "binomial",
#'                           binomial.method = "gaussian percentile",
#'                           Y.thresh = 0.75)
#' table(Bin.ex3$Y)
#'
#' Pois.ex <- sim_Y_Binary_X(N = Nex, B = Bex, im.res = c(5, 5),
#'                           dist = "poisson", print.out = TRUE)
#' mean(Pois.ex$Y)
#' quantile(Pois.ex$Y, probs = c(0, 0.1, 0.25, 0.45, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
#' hist(Pois.ex$Y)
#' @return A data frame where each row consists of a single subject's data.
#' Col 1 is the outcome, Y, and each successive column contains the subject
#' predictor values.
#' @importFrom stats rnorm rbinom quantile rpois
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Cressie+Wikle:2011}{sim2Dpredictr}
#'
#' \insertRef{Ripley:1987}{sim2Dpredictr}
#' @export
sim_Y_Binary_X = function(N, B, rand.err = 1, dist, incl.subjectID = TRUE,
                       binomial.method = "traditional",
                       count.method = "traditional",
                       Y.thresh = NULL, print.out = FALSE,
                       xlim = c(0, 1), ylim = c(0, 1), im.res,
                       radius.bounds = c(0.02, 0.1), lambda = 50,
                       random.lambda = FALSE, lambda.sd = 10,
                       lambda.bound = NULL,
                       prior = "gamma", sub.area = FALSE,
                       min.sa = c(0.1, 0.1), max.sa = c(0.3, 0.3),
                       radius.bounds.min.sa = c(0.02, 0.05),
                       radius.bounds.max.sa = c(0.08, 0.15),
                       print.subj.sa = FALSE, print.lambda = FALSE,
                       print.iter = FALSE){
  n <- NULL

  bin.options <- c("traditional",
                   "gaussian manual",
                   "gaussian percentile")

  if ( dist == "binomial" & !(binomial.method %in% bin.options) ){
    stop(paste0("Invalid selection; binomial.method must be one of ", bin.options))
  }

  # Obtain number of parameters.
  p = prod(im.res)

  if (length(B) != p + 1){stop("B must have length 1 + nrow(L)")}
  out.names=c("Y", "X0")
  for (i in 3:(p + 2)){
    out.names[i] = paste0("X", i - 2)
  }

  # generate predictors and outcome
  X0 = rep(1, N)
  X = sim2D_binarymap(N = N, xlim = xlim, ylim = ylim, im.res = im.res,
                      radius.bounds = radius.bounds, lambda = lambda,
                      random.lambda = random.lambda, lambda.sd = lambda.sd,
                      lambda.bound = lambda.bound,
                      prior = prior, sub.area = sub.area,
                      min.sa = min.sa, max.sa = max.sa,
                      radius.bounds.min.sa = radius.bounds.min.sa,
                      radius.bounds.max.sa = radius.bounds.max.sa,
                      print.subj.sa = print.subj.sa, print.lambda = print.lambda,
                      print.iter = print.iter,
                      store.type = "matrix")
  Xn = cbind(X0, X)

  if ( (dist == "gaussian") |
       (dist == "binomial" & binomial.method != "traditional") |
       (dist == "poisson" & count.method != "traditional") ) {
    if (length(rand.err) == 1){
      En = rnorm(N, 0, rand.err)
    }
    else{
      if (length(rand.err) != p){stop("Length of rand.err must be 1 or p + 1 = n.col*n.row+1")}
      En = rnorm(N, 0, rand.err[n])
    }
    Y = Xn %*% B + En
    if (print.out == TRUE) {
      print(Xn %*% B)
    }
  }

  if (dist == "binomial"){

    if (binomial.method == "traditional" & is.null(Y.thresh) == FALSE ) {
      warning(paste0("A cutoff value is not applicable when binomial.method = ", binomial.method))
    }

    if (binomial.method != "traditional" & is.null(Y.thresh) == TRUE ) {
      if (binomial.method == "gaussian percentile") {
        warning(paste0("binomial method = ", binomial.method," requires user specified cutoff value. Defaulting to the median."))
      }
      if (binomial.method == "gaussian manual") {
        warning(paste0("binomial method = ", binomial.method," requires user specified cutoff value. Defaulting to 0."))
      }
    }

    if (binomial.method == "traditional") {
      Y = c()
      XB = Xn %*% B
      # for logistic regression GLM must use either:
      # p = exp(XB[y]) / (1 + exp(XB[y]))
      # p = 1 / (1 + 1 / exp(XB[y]))
      for (y in 1:N) {
        if (print.out == TRUE) {
          cat("Subject ", y, " has XB = ", XB[y],
              " and prob of 'success' =", 1 / (1 + 1 / exp(XB[y])), "\n")
        }
        Y[y] = stats::rbinom(n = 1, size = 1, prob = 1 / (1 + 1 / exp(XB[y])))
      }
    } else {

      # prepare for thresholding
      Y0 <- Y
      Y <- c()
      # Manual Cutoff
      if (binomial.method == "gaussian manual") {

        # if no user specified value, the threshold is 0.
        if (is.null(Y.thresh) == TRUE) {
          Y.thresh = 0
        }

        Y[Y0 > Y.thresh] <- 1
        Y[Y0 <= Y.thresh] <- 0
      }

      # Percentile Cutoff
      if (binomial.method == "gaussian percentile") {

        # If no user specified value, the threshold is the median.
        if (is.null(Y.thresh) == TRUE) {
          Y.thresh = 0.50
        }

        if ( ( Y.thresh < 0) | (Y.thresh > 1) ) {
          stop(paste0("Percentiles must be between 0 and 1, but Y.thresh = ", Y.thresh))
        } else {
          perc.ct <- quantile(x = Y0, probs = Y.thresh, type = 3)
          cat("The ", 100 * Y.thresh, "th Percentile (threshold) is ", perc.ct, ".", "\n")
          Y[Y0 > perc.ct] <- 1
          Y[Y0 <= perc.ct] <- 0
        }
      }
    }
  }

  if (dist == "poisson"){

    if (count.method == "traditional") {
      # For Poisson GLM lambda = exp(XB[y])
      Y = c()
      XB = Xn %*% B
      for (y in 1:N) {
        if (print.out == TRUE) {
          cat("Subject ", y, " has XB =", XB[y],
              " and lambda =", exp(XB[y]), "\n")
        }
        Y[y] = stats::rpois(n = 1, lambda = exp(XB[y]))
      }
    } else if (count.method == "rounding") {
      Y[Y <= 0] <- 0
      Y <- round(Y)
    }
  }

  # create final dataset
  data.1 = data.frame(Y, Xn)
  colnames(data.1) = out.names
  if (incl.subjectID == TRUE) {
    data.1$subjectID = 1:N
  }
  return(data.1[ , -2])
}


