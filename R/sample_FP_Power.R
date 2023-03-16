#' Obtain Sample False Positive Rates and Power
#'
#' This function calculates sample FDR, FWER, and Power for large numbers of
#' predictors, given a vector of "true" parameter values and a vector of 
#' associated rejections. In the case that more than 1 predictor has a "true"
#' non-zero parameter, then Power is defined as the proportion/percentage of
#' those "true" parameters identified.
#'
#' @param rejections A binary vector; \code{rejection[i] = 1} means the null
#' hypothesis is rejected for parameter \code{B[i]}, whereas 
#' \code{rejection[i] = 0} means that the null hypothesis was not rejected
#'  for parameter \code{B[i]}.
#' @param test.statistic A vector of test statistics; e.g., t-statistics or
#' p-values that are used to determine whether or not to reject the null
#' hypothesis.
#' @param reject.threshold A list whose first element is the rejection 
#' criteria, e.g., the minimum t-statistic or maximum p-value for which to
#' reject the null hypothesis. The second element is one of 
#' \code{c("greater", "less", "2-tailed")}, which tell the function to reject
#' when the values in \code{test.statistic} are greater than or less than the
#' threshold, the test is a 2-tailed, respectively. In the latter case the 
#' function internally calculates the upper or lower threshold needed for the
#' 2-tailed test.
#' @param FP,TP Binary vectors of false positive and true positive indicators,
#' respectively. \code{FP[i] = 1} means the null hypothesis was incorrectly 
#' rejected, and \code{TP[i] = 1} means the null hypothesis was correctly 
#' rejected. If either argument is \code{NULL}, then these vectors are 
#' computed; this is the default setting.
#' @param B A vector of "true" parameter values. For inference purposes, this
#' can be a vector of actual parameter values, or a binary vector indicating
#' non-zero status.
#' @param B.incl.B0 If \code{B.incl.B0 = TRUE} then the first entry should be
#' the intercept, B0. \code{B.incl.B0 = FALSE} indicates that the first entry
#' of B is not an intercept.
#' @param full.summary If \code{full.summary = TRUE} then the total numbers of
#' rejections, false positives, true positives, and non-zero parameters are
#' output along with FDR, FWER, and Power; otherwise, only FDR, FWER, and 
#' Power are output.
#' @return A data frame with columns for sample FDR, FWER, and Power.
#' @note The default operating approach is that the null hypothesis is 
#' \code{B[i] = 0} for each parameter. If other hypotheses are being tested 
#' then \code{B} should be converted to a binary vector indicating whether 
#' the null hypothesis \emph{should} have been rejected.
#' @examples
#' ## example 1
#'
#' ## rejection vector
#' rej.ex <- c(0, 1, 1, 0, 0, 1, 0)
#' ## false positive vector
#' fp.ex  <- c(0, 0, 1, 0, 0, 0, 0)
#' ## true positive vector
#' tp.ex  <- c(0, 1, 0, 0, 0, 1, 0)
#' ## parameter vector
#' par.ex <- c(0, 4, 0, 0, 3, 9, 0)
#'
#' sample_FP_Power(rej.ex, fp.ex, tp.ex, par.ex, B.incl.B0 = FALSE)
#'
#' ## Function can calculate TP and FP vectors
#' sample_FP_Power(rejections = rej.ex, 
#'                 FP = NULL, TP = NULL, 
#'                 B = par.ex, B.incl.B0 = FALSE)
#'
#' ## example 2: sum(FP, TP) must equal sum(rejections) or
#' ## function stops execution
#'
#' rej.ex2 <- c(0, 1, 0, 0, 0, 1, 0)
#' fp.ex2  <- c(0, 0, 1, 0, 0, 0, 0)
#' tp.ex2  <- c(0, 1, 0, 0, 0, 1, 0)
#' par.ex2 <- c(0, 4, 0, 0, 3, 9, 0)
#'
#' \dontrun{sample_FP_Power(rej.ex2, 
#'                          fp.ex2, tp.ex2, par.ex2, 
#'                          B.incl.B0 = FALSE)}
#'
#' ## example 3: calculate rejections from vector of test statistics
#' zstat <- c(-0.5, 1.98, 2.01, 1.45, -1.99)
#' # 2-tailed
#' sample_FP_Power(test.statistic = zstat,
#'                 reject.threshold = list(1.96, "2-tailed"),
#'                 B = c(0, 0, 4, 1, -2), B.incl.B0 = FALSE)
#' # 1-tailed (upper)
#' sample_FP_Power(test.statistic = zstat,
#'                 reject.threshold = list(1.96, "greater"),
#'                 B = c(0, 0, 4, 1, -2), B.incl.B0 = FALSE)
#' ## p-value
#' sample_FP_Power(test.statistic = c(0.44, 0.04, 0.01, 0.06, 0.02 ),
#'                 reject.threshold = list(0.05, "less"),
#'                 B = c(0, 0, 4, 1, -2), B.incl.B0 = FALSE)

#' @export
sample_FP_Power <- function(rejections = NULL, 
                            FP = NULL, TP = NULL,
                            test.statistic = NULL, 
                            reject.threshold = NULL,
                            B = NULL, B.incl.B0 = TRUE, 
                            full.summary = FALSE) {

  # remove intercept from B
  if (B.incl.B0 == TRUE & is.null(B) == FALSE){
    B <- B[-1]
  }

  if (is.null(rejections) == TRUE &
      (is.null(test.statistic) == FALSE & is.null(reject.threshold) == FALSE)) {
    rejections <- make_rejection(B = B, reject.threshold = reject.threshold,
                                 test.statistic = test.statistic)
  }

  if ( (is.null(FP) == TRUE) | (is.null(TP) == TRUE) ) {
    FP <- rep(0, length(rejections))
    TP <- FP
    if (sum(rejections) != 0) {
      FP[(rejections == 1) & (B == 0)] <- 1
      TP[(rejections == 1) & (B != 0)] <- 1
    }
  }

  if (sum(FP, TP) != sum(rejections)) {
    stop("Input Inconsistency: sum(TP, FP) must equal sum(rejections)")
  }

  # FDR is proportion of true positives out of total positives
  FDP <- c()

  # FWER is either 1 or 0 for a given sample; i.e. there is
  # either 0 false positives (no error) or >0 false positives
  # (1 or more errors).
  FWE <- c()

  # Power, in this context, is the proportion of true positives
  # out of total truly significant voxels in a sample.
  Power <- c()
  if (sum(FP) == 0) {
    FDP <- 0
    FWE <- 0
  }
  if (sum(FP) > 0) {
    FWE <- 1
  }
  if (sum(rejections) > 0 & sum(FP) > 0) {
    FDP <- sum(FP) / sum(rejections)
  }
  if (sum(B != 0) > 0) {
    Power <- sum(TP) / sum(B != 0)
  }
  else {Power <- NA}
  if (full.summary == TRUE) {
    return(data.frame(NumNonzeroB = sum(B != 0),
                      NumReject = sum(rejections),
                      NumFP = sum(FP), NumTP = sum(TP),
                      FDP, FWE, Power))
  } else {
    return(data.frame(FDP, FWE, Power))
  }
}

