#' Simulate Scalar Outcomes from Simulated Spatially Correlated Predictors
#'
#' N spatially correlated design vectors are simulated from an MVN. These
#' design vectors are used to then simulate scalar outcomes that have
#' one of Gaussian, Binomial, or Poisson distributions.
#' @inheritParams sim_MVN_X 
#' @param B A vector parameter values; i.e. "betas". Note that \code{length(B)}
#' must equal \code{p + 1 = n.row * n.col + 1}; e.g. for normal outcomes
#' \eqn{Y = XB + e} with \code{Y} a scalar outcome and \code{e} the random error.
#' Note that when \code{dist = "multinomial"} then \code{B} should be a list with
#' length equal to \code{V - 1}, i.e., should contain parameter values associated
#' with all categories except the reference category.
#' @param rand.err A vector for the random error standard deviation when \code{dist = "gaussian"},
#' or thresholding is used to obtain non-Normal draws. Must have length 1 or length N.
#' @param dist The distribution of the scalar outcome.
#' \itemize{
#'     \item \code{dist = "gaussian"} has \eqn{Y = XB + e}, where
#'     \eqn{e ~ N(0, rand.err)}.
#'     \item \code{dist = "binomial"} is drawn from eqn{Bin(XB, XB(1-XB))}
#'     using \code{rbinom()} when \code{binary.method = "traditional"}. If
#'     \code{binary.method = "gaussian"}, then simulation is based on a
#'     cutoff using \code{binary.cutoff}.
#'     \item \code{dist = "multinomial"} is drawn from \code{sample()} using probabilities
#'     generated based on Chapter 6.1.3 of Agresti (2007). Threshold-based approaches are not
#'     currently supported.
#'     \item \code{dist = "poisson"} is drawn from \eqn{Poisson(XB)} using
#'     \code{rpois()}.
#' }
#' @param V A numeric value stating the number of categories desired when \code{dist = "multinomial"}.
#' @param threshold.method One of \code{"none", "manual", "percentile", "round"}.
#' When \code{"none"} draws from Binomial or Poisson distributions are taken subject-wise
#' using base \code{R} functions. For the remaining options, draws are first taken from a
#' Normal distribution and then thresholded. \code{"manual"} uses \code{Y.thresh} to manually
#' select a cutoff, \code{"percentile"} uses \code{Y.thresh} to select percentiles used to bin
#' outcomes, and \code{"round"} sets values equal or less than 0 to 0, and rounds all positive
#' values to the nearest whole number.
#' @param Y.thresh A manual value used to threshold when \code{threshold.method = "manual"}; values
#' equal or greater than the cutoff are assigned 1 and all others 0. When \code{threshold.method = "percentile"},
#' a percentile to use to bin outcomes.
#' @param incl.subjectID When \code{incl.subjectID = TRUE} a column of subject indices
#' is generated.
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
#' ## generate precision matrix and take Cholesky decomposition
#' Rpre <- chol_s2Dp(im.res = c(5, 5), matrix.type = "prec",
#'                   use.spam = TRUE, neighborhood = "ar1",
#'                   triangle = "upper", return.prec = TRUE)
#' ## Generate correlation matrix & take Cholesky decomposition
#' Rcov <- chol_s2Dp(corr.structure = "ar1", im.res = c(5, 5), rho = 0.5,
#'                   triangle = "upper",
#'                   use.spam = FALSE, neighborhood = "none")
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
#' ## with precision matrix
#' Gauss.exp <- sim_Y_MVN_X(N = Nex, B = Bex,
#'                          R = Rpre$R, Q = Rpre$Q,
#'                          dist = "gaussian")
#' hist(Gauss.exp$Y)
#'
#' ## with covariance matrix
#' Gauss.exc <- sim_Y_MVN_X(N = Nex, B = Bex,
#'                          R = Rcov$R, S = Rcov$S,
#'                          dist = "gaussian")
#' hist(Gauss.exc$Y)
#'
#' ## direct draws from binomial
#' Bin.ex <- sim_Y_MVN_X(N = Nex, B = Bex, R = Rcov$R, S = Rcov$S,
#'                       dist = "binomial", print.out = TRUE)
#' table(Bin.ex$Y)
#'
#' ## manual cutoff
#' Bin.ex2 <- sim_Y_MVN_X(N = Nex, B = Bex,
#'                        R = Rcov$R, S = Rcov$S,
#'                        dist = "binomial",
#'                        threshold.method = "manual",
#'                        Y.thresh = 1.25)
#' table(Bin.ex2$Y)
#'
#' ## percentile cutoff
#' Bin.ex3 <- sim_Y_MVN_X(N = Nex, B = Bex,
#'                        R = Rcov$R, S = Rcov$S,
#'                        dist = "binomial",
#'                        threshold.method = "percentile",
#'                        Y.thresh = 0.75)
#' table(Bin.ex3$Y)
#'
#' ## Poisson Example - note the large counts
#' Pois.ex <- sim_Y_MVN_X(N = Nex, B = Bex,
#'                        R = Rcov$R, S = Rcov$S,
#'                        dist = "poisson", print.out = TRUE)
#' mean(Pois.ex$Y)
#' quantile(Pois.ex$Y, probs = c(0, 0.1, 0.25, 0.45, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
#' hist(Pois.ex$Y)
#' @return A data frame where each row consists of a single subject's data.
#' Col 1 is the outcome, Y, and each successive column contains the subject
#' predictor values.
#' @importFrom stats rnorm rbinom quantile rpois
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{spam}{sim2Dpredictr}
#'
#' \insertRef{Ripley:1987}{sim2Dpredictr}
#'
#' \insertRef{Rue:2001}{sim2Dpredictr}
#' 
#' \insertRef{Agresti:2007}{sim2Dpredictr}
#' @export
sim_Y_MVN_X = function(N, B, L = NULL, R = NULL,
                       S = NULL, Q = NULL, use.spam = TRUE,
                       mu = 0, rand.err = 1,
                       dist = "gaussian", V = NULL,
                       incl.subjectID = TRUE,
                       threshold.method = "none",
                       Y.thresh = NULL,
                       X.categorical = FALSE, X.num.categories = 2,
                       X.category.type = "percentile", X.manual.thresh = NULL,
                       X.cat.names = NULL, print.out = FALSE){
  n <- NULL

  if ( X.categorical == TRUE & X.num.categories > 2) {
    stop("Current version of sim_Y_MVN_X requires X.num.categories <= 2.")
  }

  if (threshold.method == "round" & dist != "poisson") {
    stop("Rounding is only appropriate for count outcomes.")
  }

  # Predictable errors and warnings...
  if ( (is.null(R) == TRUE) & (is.null(L) == TRUE) ){
    stop("Function requires either specification of L or R.")
  }
  if ( (is.null(R) == FALSE) & (is.null(L) == FALSE) ){
    warning("Specify either L or R, not both. Function proceeding with L.")
    R <- NULL
  }

  # Obtain number of parameters.
  LRSQ <- list(L = L, R = R, S = S, Q = Q)
  LRSQw <- LRSQ[which(lapply(LRSQ, is.null) == FALSE)]
  LRSQw1 <- LRSQw[1]
  p <- dim(LRSQw[[1]])[1]

  # build and/or verify mu is acceptable
  if (length(mu) == 1) {
    mu <- rep(mu, p)
  } else {
    if (length(mu) != p) {cat("Invalid mean vector length.", "\n")}
  }
  
  if (dist == "multinomial") {
    for (v in 1:length(B)) {
      if (length(B[[v]]) != p + 1){stop("Each B must have length 1 + nrow(L)")}
    }
  } else {
    if (length(B) != p + 1){stop("B must have length 1 + nrow(L)")}
  }
  
  out.names=c("Y", "X0")
  for (i in 3:(p + 2)){
    out.names[i] = paste0("X", i - 2)
  }

  # generate predictors and outcome
  X0 <- rep(1, N)
  X <- sim_MVN_X(N = N, mu = mu, L = L, R = R,
                 S = S, Q = Q, use.spam = use.spam,
                 X.categorical = X.categorical,
                 X.num.categories = X.num.categories,
                 X.category.type = X.category.type,
                 X.manual.thresh = X.manual.thresh,
                 X.cat.names = X.cat.names)
  Xn <- cbind(X0, X)

  if ( (dist == "gaussian") |
       (dist %in% c("binomial", "poisson") & threshold.method != "none")) {

    if (!(length(rand.err) %in% c(1, N))){
      stop("Length of rand.err must be 1 or N")
    } else {
      En <-stats:: rnorm(N, 0, rand.err)
      Y <- Xn %*% B + En
      if (print.out == TRUE) {
        print(Xn %*% B)
      }
    }
  }

  if (dist == "binomial"){

    if (threshold.method == "none" & is.null(Y.thresh) == FALSE ) {
      warning(paste0("A cutoff value is not applicable when threshold.method = ", threshold.method))
    }

    if (threshold.method != "none" & is.null(Y.thresh) == TRUE ) {
      if (threshold.method == "percentile") {
        warning(paste0("threshold.method = ",
                       threshold.method,
                       " requires user specified cutoff value. Defaulting to the median."))
      }
      if (threshold.method == "manual") {
        warning(paste0("threshold.method = ",
                       threshold.method," requires user specified cutoff value. Defaulting to 0."))
      }
    }

    if (threshold.method == "none") {
      Y <- c()
      XB <- Xn %*% B
      # for logistic regression GLM must use either:
      # p = exp(XB[y]) / (1 + exp(XB[y]))
      # p = 1 / (1 + 1 / exp(XB[y]))
      for (y in 1:N) {
        if (print.out == TRUE) {
          cat("Subject ", y, " has XB = ", XB[y],
              " and prob of 'success' =", 1 / (1 + 1 / exp(XB[y])), "\n")
        }
        Y[y] <- stats::rbinom(n = 1, size = 1, prob = 1 / (1 + 1 / exp(XB[y])))
      }
    } else {

      # prepare for thresholding
      Y0 <- Y
      Y <- c()
      # Manual Cutoff
      if (threshold.method == "manual") {

        # if no user specified value, the threshold is 0.
        if (is.null(Y.thresh) == TRUE) {
          Y.thresh = 0
        }

        Y[Y0 > Y.thresh] <- 1
        Y[Y0 <= Y.thresh] <- 0
      }

      # Percentile Cutoff
      if (threshold.method == "percentile") {

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
  
  if (dist == "multinomial") {
    Y <- c()
    
    # generate subject-specific probabilities
    p.mn <- generate_multinom_probs(V = V, B = B, X = Xn, X.incl.X0 = TRUE)
    
    Y <- c()
    for (i in 1:nrow(p.mn)) {
      Y[i] <-  sample(x = 1:V, size = 1, prob = p.mn[i, ])
    }
  }

  if (dist == "poisson"){

    if (threshold.method == "none") {
      # For Poisson GLM lambda = exp(XB[y])
      Y <- c()
      XB <- Xn %*% B
      for (y in 1:N) {
        if (print.out == TRUE) {
          cat("Subject ", y, " has XB =", XB[y],
              " and lambda =", exp(XB[y]), "\n")
        }
        Y[y] = stats::rpois(n = 1, lambda = exp(XB[y]))
      }
    } else if (threshold.method == "rounding") {
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


