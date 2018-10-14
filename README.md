---
output: github_document
---

<!-- rmarkdown v1 -->

<!-- README.md is generated from README.Rmd. Please edit that file -->



# multilevelMatching

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/multilevelMatching)](https://cran.r-project.org/package=multilevelMatching)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Travis-CI Build Status](https://travis-ci.org/shuyang1987/multilevelMatching.svg?branch=master)](https://travis-ci.org/shuyang1987/multilevelMatching)
[![AppveyorCI Build status](https://ci.appveyor.com/api/projects/status/eu7vlcbu2j854cdo?svg=true)](https://ci.appveyor.com/project/BarkleyBG/multilevelmatching-3hh85)
[![Coverage status](https://codecov.io/gh/shuyang1987/multilevelMatching/branch/master/graph/badge.svg)](https://codecov.io/github/shuyang1987/multilevelMatching?branch=master)

### Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments 

Installation:


```r
devtools::install_github("shuyang1987/multilevelMatching")
```


### Visit the [package website](https://shuyang1987.github.io/multilevelMatching/)


# Description


This package implements methods to estimate causal effects from observational studies when there are 2+ distinct levels of treatment (i.e., "multilevel treatment") using matching estimators, as introduced in [Yang et al. (2016) Biometrics](https://doi.org/10.1111/biom.12505). Matching on covariates, and matching or stratification on modeled propensity scores, are made available. These methods require matching on only a scalar function of generalized propensity scores. For more information, see the Extended Description below or the main paper:

  - Yang, S., Imbens G. W., Cui, Z., Faries, D. E., & Kadziola, Z. (2016) Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments. *Biometrics*, 72, 1055-1065. https://doi.org/10.1111/biom.12505
  
    
Visit the [package website](https://shuyang1987.github.io/multilevelMatching/)


# Estimators available

- Matching on raw covariates: via `multiMatch()` and `multilevelMatchX()`
- Matching on estimated propensity scores: via `multiMatch()` and `multilevelGPSMatch()`
  - using ordinal logistic regression
  - using multinomial logistic regression
     - This method also provides two types of variance estimates
  - using user-provided propensity score values
     - This method does not provide variance estimates
- Stratification on propensity scores: via `multilevelGPSStratification()`

# Tutorial

This is a brief tutorial; an extended tutorial is provided in the vignette for [version 1.0.0](https://github.com/shuyang1987/multilevelMatching/releases/).
We will use the dataset provided with this package


```r
library(multilevelMatching)
simulated_data <- multilevelMatching::simulated_data
knitr::kable(head(simulated_data), digits = 2)
```



| outcome| treatment| covar1| covar2| covar3| covar4| covar5| covar6|
|-------:|---------:|------:|------:|------:|------:|------:|------:|
|   -5.13|         1|  -0.87|   0.24|   0.23|  -2.89|   0.21|      0|
|   -3.03|         1|   0.27|  -0.35|  -0.40|  -2.21|   0.07|      1|
|    3.05|         1|   1.42|   1.32|  -1.20|   0.06|   1.26|      1|
|   -6.09|         1|  -1.39|  -1.07|   1.12|  -2.36|   0.00|      0|
|   -2.46|         1|  -1.15|   0.95|   0.51|  -2.68|   0.07|      0|
|   -0.32|         1|   0.27|   0.42|  -0.45|   2.43|   0.60|      0|

We restructure the dataframe slightly, and use identifying names for the observations:


```r
outcome <- simulated_data$outcome
treatment <- simulated_data$treatment
covar_matrix <- as.matrix(
  simulated_data[ ,names(simulated_data) %in% paste0("covar", 1:6)]
)
identifying_names <- paste0(
  rep(letters[1:25],each = 12), rep(letters[1:25], 12)
)
names(treatment) <- identifying_names
```

## Matching on covariates


```r
set.seed(123)
fit <- multiMatch(
  Y = outcome,
  W = treatment,
  X = covar_matrix,
  match_on = "covariates"
)

fit
#> -------------- Causal estimates ---------------
#>         Param Trt1 Trt2   Estimate  Variance
#> 1 EY(2)-EY(1)    1    2 0.07927361 0.1792186
#> 2 EY(3)-EY(1)    1    3 0.86264929 0.1634754
#> 3 EY(3)-EY(2)    2    3 0.78337567 0.3221616
#> --- Matching on 'covariates' with M=1, J=1 ---
```

## Matching on the Estimated Generalized Propensity Score (GPS)

Propensity scores can be estimated with either of the following options

  - `match_on="multinom"` for multinomial logistic regression from `nnet::multinom()`
  - `match_on="polr"` for ordinal logistic regression from `MASS::polr()`
  - Or, estimated propensity scores can be supplied via the `X` argument when `match_on="existing"`
  

```r
match_on <- "multinom"
# match_on <- "polr" 

set.seed(123)
fit2 <- multiMatch(
  Y = outcome,
  W = treatment,
  X = covar_matrix,
  match_on = match_on,
  trimming = FALSE
)

fit
#> -------------- Causal estimates ---------------
#>         Param Trt1 Trt2   Estimate  Variance
#> 1 EY(2)-EY(1)    1    2 0.07927361 0.1792186
#> 2 EY(3)-EY(1)    1    3 0.86264929 0.1634754
#> 3 EY(3)-EY(2)    2    3 0.78337567 0.3221616
#> --- Matching on 'covariates' with M=1, J=1 ---
```


Please see the vignette for an extended tutorial.

# Extended Description

## Matching with 3 or more levels of treatment

In setting with where 3 or more levels of treatment (i.e., multilevel treatment), our goal is to estimate pairwise average treatment effects from a common population using matching methods.

This goal can not be acheived by matching one treatment with another one at a time, since the pairwise matched samples may differ from the target population systematically, and thus they are not compatitable. One implication is that from this approach, it is possible that treatment A is better than treatment B, treatment B is better than treatment C, and treatment C is better than treatment A. 

We focus on estimating the average values of potential outcomes for each treatment level by matching methods, which facilitate estimation of pairwise average treatment effects for a common population.

The estimation methods include generalized propensity score (GPS) matching, GPS stratification, matching with the full set of covariates, matching with the full set of GPS vector. Note that GPS matching and GPS straticication only require matching on a scalar function when estimating the average value of the potential outcome at a particular treatment level, which reduces the matching dimension to one, regardless of the number of covariates and the number of treatment levels. 

In order to ensure sufficient overlap, [Crump et al. (2009)](https://doi.org/10.1093/biomet/asn055)'s trimming method can be extended to this setting as well. 


# News

See [the News site](https://shuyang1987.github.io/multilevelMatching/news/index.html) for the changelog.

#### A note on `multiMatch()`

The `multiMatch()` function may return slightly different estimates than the original 2 matching functions in certain circumstances. We attempt to ensure that the functions implement are identical methods up to perhaps random number generation. Please file an issue if you have any questions or concerns.
