
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sim2Dpredictr

<!-- badges: start -->

<!-- badges: end -->

The goal of sim2Dpredictr is to facilitate straightforward simulation of
spatially dependent predictors (continuous or binary), which may then be
used to simulate continuous, binary, or count outcomes within a
(generalized) linear model framework. Continuous predictors are
simulated using Gauss Markov Random Fields with a focus on specific
correlation structures; tools are included for easily constructing and
taking the Cholesky decomposition of a covariance matrix. The Boolean
Model is used to simulate spatially dependent binary maps. The package
also includes a tool for easily specifying a parameter vector with
spatially clustered non-zero elements. These simulation tools are
designed for, but not limited to, testing the performance of variable
selection methods when predictors are spatially correlated.

## Installation

You can install the released version of sim2Dpredictr from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sim2Dpredictr")
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
Bex <- beta_builder(row.index = c(3, 3, 4, 4), 
                    col.index = c(3, 4, 3, 4),
                    im.res = c(5, 5),
                    B0 = 0, B.values = rep(1, 4),
                    output.indices = "No")

# Construct and take Cholesky decomposition of correlation matrix.
Rex <- chol_s2Dp(corr.structure = "ar1", im.res = c(5, 5), rho = 0.5,
                 use.spam = "Yes", neighborhood = "round", r = 6)

# Simulate a dataset with spatially dependent design matrix and binary outcomes.
sim.dat <- Gauss.ex <- sim_Y_MVN_X(N = 10, B = Bex, R = Rex, 
                                   dist = "binomial")

sim.dat
#>    Y          X1         X2          X3          X4          X5
#> 1  0 -0.02061515 -1.5894632 -1.35559167  0.16484763 -0.26523381
#> 2  0 -0.43710475  0.5515873  0.65557146 -0.29795070 -1.51309368
#> 3  0 -1.19094081  0.7277853  0.66378189  0.22665260 -0.17825732
#> 4  1  1.21684167 -0.8262323 -2.14888691 -2.06738788 -3.23741986
#> 5  0  0.20291944  0.3080348  0.15841043 -0.05737876 -0.14313532
#> 6  0 -0.98804308 -0.2612674 -1.02490514 -1.35795956 -0.99930280
#> 7  1 -0.02165241 -0.7256593  1.75109874  0.16724044  0.06504147
#> 8  0  0.14961732 -0.2450088 -0.07224691 -1.74498776 -2.41541321
#> 9  0  2.06117251  1.0487435 -1.20116935 -0.28928825  0.22747426
#> 10 1  1.18684635 -0.3070706 -0.10425472 -0.32944192 -0.35790135
#>              X6           X7         X8          X9        X10         X11
#> 1   0.063798218 -2.101505134 -1.1737204  0.67995876 -0.4029497 -0.79952658
#> 2   0.154721604 -0.003466665  0.6805385 -1.02192100 -0.4028219  0.38923798
#> 3  -2.632438623 -0.731610447 -0.8164019  0.45248778  0.1704944 -0.28016817
#> 4   0.386062983  0.793145013 -0.4755736 -1.41307920 -1.3075042  1.76024910
#> 5  -0.334767730 -2.548124612 -1.3426384 -1.28011311 -1.3935061 -0.71437785
#> 6  -0.005078639 -1.113566217 -2.5962909 -3.17801226 -0.9762987 -0.75996737
#> 7  -2.165911165 -0.110272533  1.6413552  1.37636292 -0.1910851 -0.59411463
#> 8   1.060219499 -0.154291188 -0.8445898 -0.32893185 -0.9967608  0.00538278
#> 9   1.098371583  1.183811346 -1.4305157 -1.66025279 -1.2402013  1.32711796
#> 10  1.828325958  0.398185802 -1.2804263 -0.03882051  0.2173992 -0.25086508
#>           X12         X13         X14         X15        X16        X17
#> 1  -0.5452645  0.20961178  0.98752024 -0.13112465 -1.2976024 -0.6748704
#> 2   0.5298656 -0.53461958  0.73809818  0.62949047  0.4713020  1.5959354
#> 3  -1.0919398  0.20916085 -0.57722989 -0.33306245  0.9901808  0.7340450
#> 4   0.1738520 -0.06722856 -0.81154873 -1.96833602  0.7659767 -0.2512449
#> 5  -0.1928543 -1.47643027 -1.93661614 -0.80808206 -0.3898620 -1.2655171
#> 6  -1.0673832 -1.77616341 -1.64554901 -0.67380547 -0.5900710 -2.5482818
#> 7  -0.5811536  1.51768475  0.53734093  0.03918055 -1.4689740 -0.7993639
#> 8   1.6365247  0.03595908 -0.45098266 -1.48679757  2.5600802  0.2886936
#> 9  -0.1578330 -1.19189538 -1.08398430  0.05319404 -0.1104455  0.3501400
#> 10 -0.4183365  1.63082668  0.03781277 -0.22231414 -0.6630328  0.8234409
#>           X18         X19         X20        X21         X22         X23
#> 1   0.3872862  0.17466081  0.03230334 -1.4543758 -2.56053077 -0.15381477
#> 2  -0.2927720 -0.37455994  0.26270434  0.2083138  2.07469551  0.87097106
#> 3   0.9519440 -0.60383548 -0.49426494 -0.5173436 -1.39331988 -0.43529298
#> 4  -0.3137399 -0.02724486 -1.93371327  1.4731123  0.30658680  0.67993503
#> 5  -2.4064875 -0.71584311 -0.63541224  1.1598316 -1.59336909 -0.86197895
#> 6  -0.7272319 -1.41437884 -1.25377746 -0.9548261 -1.54367440  0.34064265
#> 7   0.5202346  0.55340477  2.39373633 -1.6589867 -0.54905412  0.09348630
#> 8   0.4771259 -0.76618159 -2.29425853  1.0699468  1.23570133  0.04138893
#> 9  -0.7018730 -1.05591837  1.51754001 -0.5255136 -0.04246015  1.14734700
#> 10  1.4292531  0.86024863 -0.22698654  1.0764390  0.68413304  1.24301906
#>           X24        X25 subjectID
#> 1   1.4695680  2.2471619         1
#> 2   0.4788780 -0.4656137         2
#> 3  -0.2934772 -0.8508825         3
#> 4   0.2317841 -0.7638719         4
#> 5  -0.4804482  0.6277547         5
#> 6   0.4428507 -1.5663854         6
#> 7   1.1586797  0.6246729         7
#> 8  -0.8493131 -0.9126043         8
#> 9  -0.3217513  0.8282605         9
#> 10  1.6363927  0.7151107        10
```

Once the dependence framework and non-zero parameter vector is set,
`sim_Y_MVN_X()` can be used to draw as many datasets as necessary, upon
each of which variable selection methods are applied; summaries from
each analzed dataset can be obtained and then used to evaluate variable
selection performance. The documentation provides details about how to
use these functions (and others) to create desired simulations, and a
detailed vignette is being written to provide further guidance.
"# sim2Dpredictr" 
