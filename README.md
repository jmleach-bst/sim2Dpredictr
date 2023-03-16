
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sim2Dpredictr

<!-- badges: start -->
<!-- [![Travis build status](https://travis-ci.com/jmleach-bst/sim2Dpredictr.svg?branch=master)](https://travis-ci.com/jmleach-bst/sim2Dpredictr) -->

[![](https://www.r-pkg.org/badges/version/sim2Dpredictr?color=green)](https://cran.r-project.org/package=sim2Dpredictr)
[![](https://img.shields.io/badge/devel%20version-0.1.1-green.svg)](https://github.com/jmleach-bst/sim2Dpredictr)
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- [![CRAN checks](https://badges.cranchecks.info/summary/sim2Dpredictr.svg)](https://cran.r-project.org/web/checks/check_results_sim2Dpredictr.html) -->
<!-- [![R build status](https://github.com/jmleach-bst/sim2Dpredictr/workflows/R-CMD-CHECK/badge.svg)](https://github.com/jmleach-bst/sim2Dpredictr/actions) -->
<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/sim2Dpredictr)](https://cran.r-project.org/package=sim2Dpredictr) -->
[![R-CMD-check](https://github.com/jmleach-bst/sim2Dpredictr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmleach-bst/sim2Dpredictr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `sim2Dpredictr` is to facilitate straightforward simulation
of spatially dependent predictors (continuous or binary), which may then
be used to simulate continuous, binary, count, or categorical ($> 2$
categories) outcomes within a (generalized) linear model framework. A
real-world example is when using medical images to model/predict
(scalar) clinical outcomes; such a scenario motivated the development of
`sim2Dpredictr`, which was used to simulate data to evaluate the
performance of methods for high-dimensional data analysis and prediction
(Leach, Aban, and Yi 2022; Leach et al. 2022).

In the first step, we simulate the predictors, i.e., the $\mathbf{X}_i$
part of a GLM, $g(E[Y_i]) = \mathbf{X}_i\mathbf{\beta}$, where
$g(\cdot)$ is an appropriate link function.

Continuous predictors are simulated using Multivariate Normal (MVN)
distributions with a focus on specific correlation structures;
alternatively, one can specify conditional dependence via a precision
matrix, specifically for a Conditional Autoregressive (CAR) model. Tools
are included for easily constructing and taking the Cholesky
decomposition of a covariance or precision matrix with either base `R`
or the `R` package `spam`, which makes this process faster when the
matrix is sparse. The Boolean Model and thresholding of MVN’s are used
to simulate spatially dependent binary maps. The package also includes a
tool for easily specifying a parameter vector with spatially clustered
non-zero elements. These simulation tools are designed for, but not
limited to, testing the performance of variable selection methods when
predictors are spatially correlated.

In the second step we use the predictor vectors $\mathbf{X}_i$ to
generate scalar outcomes, i.e., the $Y_i$ part of the GLM. The default
approach is to use the inverse link function to define subject specific
means, $\mu_i = E[Y_i] = g^{-1}(\mathbf{X}_i\mathbf{\beta})$. For
Normally distributed outcomes, $g(\cdot)$ is the identity link, so
$\mu_i = \mathbf{X}_i\mathbf{\beta}$ and a separate dispersion
parameter, $\sigma^2$, specifies the variance. In general, the variance
may be a function of the mean, as well as some other dispersion
parameter. We can draw directly from the desired distributions using
(some function of) $\mu_i$ (and if necessary, dispersion, $\sigma^2$) as
the parameter(s) for that distribution to obtain outcomes, $Y_i$.
Alternatively, outcomes can be initially generated as continuous and
then a threshold applied to obtain binary/categorical data if using the
inverse link function is too computationally expensive.

## Installation

`sim2Dpredictr` is available on
[CRAN](https://cran.r-project.org/package=sim2Dpredictr):

``` r
install.packages("sim2Dpredictr")
```

You can install the latest version of `sim2Dpredictr` from
[GitHub](https://github.com) with:

``` r
devtools::install_github("jmleach-bst/sim2Dpredictr")
```

## Example

A simple demonstration is as follows; suppose each subject has a
$5 \times 5$ standardized continuous-valued predictor image, and a
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
#>   Y        X1         X2        X3         X4          X5         X6         X7
#> 1 0 0.5198678 -1.1035062 -0.805563 -1.3767652 -1.70525968 -0.6388318 -1.5858393
#> 2 0 0.1798315  0.5728877  1.463024 -0.9348662 -0.05594379  0.8592733  0.2619289
#> 3 1 0.1996273  1.7450479  2.500030  1.4882150  1.61903239  1.2503542  0.3439712
#>            X8          X9 subjectID
#> 1 -1.72102736  0.06016316         1
#> 2  0.08300277 -0.28904666         2
#> 3  0.63695327  0.85655952         3
```

Once the dependence framework and non-zero parameter vector is set,
`sim_Y_MVN_X()` can be used to draw as many data sets as necessary, upon
each of which variable selection methods are applied; summaries from
each analyzed data set can be obtained and then used to evaluate
variable selection performance. The documentation provides details about
how to use these functions (and others) to create desired simulations.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Leach2022a" class="csl-entry">

Leach, Justin M, Inmaculada Aban, and Nengjun Yi. 2022. “Incorporating
Spatial Structure into Inclusion Probabilities for Bayesian Variable
Selection in Generalized Linear Models with the Spike-and-Slab Elastic
Net.” *Journal of Statistical Planning and Inference* 217: 141–52.
<https://doi.org/10.1016/j.jspi.2021.07.010>.

</div>

<div id="ref-Leach2022b" class="csl-entry">

Leach, Justin M, Lloyd J Edwards, Rajesh Kana, Kristina Visscher,
Nengjun Yi, and Inmaculada Aban. 2022. “The Spike-and-Slab Elastic Net
as a Classification Tool in Alzheimer’s Disease.” *PLoS ONE* 17:
e0262367. <https://doi.org/10.1371/journal.pone.0262367>.

</div>

</div>
