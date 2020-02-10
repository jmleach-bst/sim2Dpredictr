
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sim2Dpredictr

<!-- badges: start -->

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

You can install the released version of sim2Dpredictr from
[github](https://github.com) with:

``` r
library(devtools)
library(remotes)
install_github("jmleach-bst/sim2Dpredictr", force=T, build_vignettes=T)
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
#>   Y         X1         X2         X3         X4         X5         X6
#> 1 1  0.7566782  0.8031276  0.8451634  1.3631704  0.6313988  0.9044230
#> 2 0 -0.1832055 -0.1108454  0.9215691 -0.7816555 -0.5723892  0.3313123
#> 3 0 -0.7854501 -1.3878088 -0.8661513  0.4638652  0.1423327 -0.1853822
#>           X7         X8        X9 subjectID
#> 1  0.3117991 -0.1180066 0.2538337         1
#> 2 -1.2781437 -0.8254655 0.3106721         2
#> 3  0.1108675  1.8265644 0.6759272         3
```

Once the dependence framework and non-zero parameter vector is set,
`sim_Y_MVN_X()` can be used to draw as many datasets as necessary, upon
each of which variable selection methods are applied; summaries from
each analzed dataset can be obtained and then used to evaluate variable
selection performance. The documentation provides details about how to
use these functions (and others) to create desired simulations, and a
detailed vignette is being written to provide further guidance.
