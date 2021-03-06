
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sim2Dpredictr

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jmleach-bst/sim2Dpredictr.svg?branch=master)](https://travis-ci.org/jmleach-bst/sim2Dpredictr)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/sim2Dpredictr)](https://cran.r-project.org/package=sim2Dpredictr)
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

The cleaned up version is available on
[CRAN](https://cran.r-project.org/package=sim2Dpredictr):

``` r
install.packages("sim2Dpredictr")
```

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
#>   Y         X1         X2         X3        X4         X5         X6        X7
#> 1 0  1.1786151 -0.2661864 -1.5314784 0.2846413 -0.1856100  0.4776866 0.2602763
#> 2 0 -0.4664315 -0.8259383 -1.3538080 0.1751565 -1.4733417 -0.8754491 0.6765135
#> 3 1  1.0768827  2.3713662  0.7175469 1.4566237  0.7241528  0.7297517 1.8309959
#>           X8        X9 subjectID
#> 1  0.3592724 0.5775858         1
#> 2 -1.0393505 0.6504225         2
#> 3 -1.5098455 1.7212279         3
```

Once the dependence framework and non-zero parameter vector is set,
`sim_Y_MVN_X()` can be used to draw as many datasets as necessary, upon
each of which variable selection methods are applied; summaries from
each analyzed dataset can be obtained and then used to evaluate variable
selection performance. The documentation provides details about how to
use these functions (and others) to create desired simulations, and a
detailed vignette is being written to provide further guidance.
