
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sim2Dpredictr

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jmleach-bst/sim2Dpredictr.svg?branch=master)](https://travis-ci.org/jmleach-bst/sim2Dpredictr)
<!-- badges: end -->

The goal of sim2Dpredictr is to facilitate straightforward simulation of
spatially dependent predictors (continuous or binary), which may then be
used to simulate continuous, binary, or count outcomes within a
(generalized) linear model framework. Continuous predictors are
simulated using Multivariate Normal (MVN) distributions with a focus on
specific correlation structures; alternatively, one can specify
conditional dependence via a precision matrix, specifically for a
Conditional Autoregressive (CAR) model. Tools are included for easily
constructing and taking the Cholesky decomposition of a covariance or
precision matrix with either base  or the  package , which makes this
process faster when the matrix is sparse. The Boolean Model and
thresholding of MVN’s are used to simulate spatially dependent binary
maps. The package also includes a tool for easily specifying a parameter
vector with spatially clustered non-zero elements. These simulation
tools are designed for, but not limited to, testing the performance of
variable selection methods when predictors are spatially correlated.

## Installation

You can install the latest version of sim2Dpredictr from
[github](https://github.com) with:

``` r
devtools::install_github("jmleach-bst/sim2Dpredictr")
```

## Example

A simple demonstration is as follows; suppose each subject has a
\(5 \times 5\) standardized continuous-valued predictor image, and a
binary outcome. We can generate a spatial cluster of non-zero parameter
values with `beta_builder()`, simulate and take the Cholesky
decomposition of a correlation (or covariance) matrix with
`chol_s2Dp()`, and generate both the images and outcomes with
`sim_Y_MVN_X()`.

``` r
library(sim2Dpredictr)
# Construct spatially clusterd non-zero parameters.
Bex <- sim2Dpredictr::beta_builder(row.index = c(1, 1, 2), 
                                   col.index = c(1, 2, 1),
                                   im.res = c(3, 3),
                                   B0 = 0, B.values = rep(1, 3))

# Construct and take Cholesky decomposition of correlation matrix.
Rex <- sim2Dpredictr::chol_s2Dp(corr.structure = "ar1", 
                                im.res = c(3, 3), rho = 0.5,
                                use.spam = TRUE)

# Simulate a dataset with spatially dependent design matrix and binary outcomes.
sim.dat <- sim2Dpredictr::sim_Y_MVN_X(N = 3, B = Bex$B, 
                                      R = Rex$R, S = Rex$S, 
                                      dist = "binomial")

sim.dat
#>   Y         X1        X2        X3        X4         X5         X6         X7
#> 1 1 0.01648842 0.5932870 2.3363894 1.1158671 -0.7113958  0.5699940 -0.1636371
#> 2 1 0.68162682 0.6137754 0.8938886 0.5550313  1.4689111 -0.0611457  0.3572808
#> 3 1 0.02529938 1.4692546 2.0659131 0.3617255  0.7413636  0.3820786 -0.4709804
#>           X8         X9 subjectID
#> 1 -0.3385803 0.70755570         1
#> 2  0.1441915 0.39455996         2
#> 3  1.2407625 0.06535728         3
```

Once the dependence framework and non-zero parameter vector is set,
`sim_Y_MVN_X()` can be used to draw as many datasets as necessary, upon
each of which variable selection methods are applied; summaries from
each analyzed dataset can be obtained and then used to evaluate variable
selection performance. The documentation provides details about how to
use these functions (and others) to create desired simulations, and a
detailed vignette is being written to provide further guidance.
