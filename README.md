---
  output: github_document
---
<!-- rmarkdown v1 -->

<!-- README.md is generated from README.Rmd. Please edit that file -->



# multilevelMatching

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/multilevelMatching)](https://cran.r-project.org/package=multilevelMatching)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build Status](https://travis-ci.org/shuyang1987/multilevelMatching.svg?branch=master)](https://travis-ci.org/shuyang1987/multilevelMatching)
[![AppveyorCI Build status](https://ci.appveyor.com/api/projects/status/eu7vlcbu2j854cdo?svg=true)](https://ci.appveyor.com/project/BarkleyBG/multilevelmatching-3hh85)
[![Coverage status](https://codecov.io/gh/shuyang1987/multilevelMatching/branch/master/graph/badge.svg)](https://codecov.io/github/shuyang1987/multilevelMatching?branch=master)

Propensity Score Matching and Subclassification in Observational Studies with Multi-level Treatments

### Visit the [package website](https://shuyang1987.github.io/multilevelMatching/)

## Description

In setting with Multi-level treatments, our goal is to estimate pairwise average treatment effects from a common population using matching methods.

This goal can not be acheived by matching one treatment with another one at a time, since the pairwise matched samples may differ from the target population systematically, and thus they are not compatitable. One implication is that from this approach, it is possible that treatment A is better than treatment B, treatment B is better than treatment C, and treatment C is better than treatment A. 

We focus on estimating the average values of potential outcomes for each treatment level by matching methods, which facilitate estimation of pairwise average treatment effects for a common population.

The estimation methods include generalized propensity score (GPS) matching, GPS stratification,
matching with the full set of covariates, matching with the full set of GPS vector. Note that GPS matching and GPS straticication only require matching on a scalar function when estimating the average value of the potential outcome at a particular treatment level, which reduces the matching dimension to one, regardless of the number of covariates and the number of treatment levels. 

In order to ensure sufficient overlap, Crump et al. (2009)'s trimming method can be extended to this setting as well. 

## Installation with `devtools`:


```r
devtools::install_github("shuyang1987/multilevelMatching")
```


## Use

The functions implementing matching methods are shown here:

| Version | |Covariate Matching Function | Propensity Score Matching Function |
|----|----|------|-----|
| `v0.1`+| | `multilevelMatchX()` | `multilevelGPSMatch()` |
| `v0.1.0.9000`+ | | `multiMatch()` | `multiMatch()` |

Note that these functions have slightly different syntax. The `multiMatch()` function was added to standardize some of the inputs, has more verbose output with S3 class `"multiMatch"`.

- To estimate by stratifying on the estimating propensity score, use `multilevelGPSStratification()`.

### Toy example


```r
## Outcome
Y <- c(102,105,120,130,100,80,94,108,96)
## Treatment level
W <- c(1,1,1,3,2,3,2,1,2)
## Baseline Covariates
X <- c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
```

#### Matching on covariates


```r
library(multilevelMatching)
set.seed(123)

multilevelMatchX(Y,W,X)
#> $tauestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>  -10.666667    6.666667   17.333333 
#> 
#> $varestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>    9.111111  615.580247  613.925926

multiMatch(Y,W,X, match_on = "covariates")
#> no unit IDs supplied; unit_ids will be assigned generically
#> Warning in prepareData(Y = Y, W = W, X = X, match_on = match_on, trimming
#> = trimming, : It is recommended that X is a matrix. X will be coerced to a
#> matrix.
#> -------------- Causal estimates ---------------
#>         Param Trt1 Trt2   Estimate   Variance
#> 1 EY(2)-EY(1)    1    2 -10.666667   9.111111
#> 2 EY(3)-EY(1)    1    3   6.666667 615.580247
#> 3 EY(3)-EY(2)    2    3  17.333333 613.925926
#> --- Matching on 'covariates' with M=1, J=1 ---
```

#### Matching on Estimated Propensity Score

Using multinomial logistic regression:

```r
## v0.1
multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
#> $tauestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>  -10.444444    6.666667   17.111111 
#> 
#> $varestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>    8.545953  616.913580  611.122085 
#> 
#> $varestimateAI2012
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>    8.302024  411.456234  434.247037 
#> 
#> $analysisidx
#> [1] 1 2 3 4 5 6 7 8 9

## v0.1.0.9000+
multiMatch(Y,W,X,trimming = 0, match_on = "multinom")
#> no unit IDs supplied; unit_ids will be assigned generically
#> Warning in prepareData(Y = Y, W = W, X = X, match_on = match_on, trimming
#> = trimming, : It is recommended that X is a matrix. X will be coerced to a
#> matrix.
#> -------------- Causal estimates ---------------
#>         Param Trt1 Trt2   Estimate   Variance VarianceAI2016
#> 1 EY(2)-EY(1)    1    2 -10.444444   8.545953       8.302024
#> 2 EY(3)-EY(1)    1    3   6.666667 616.913580     411.456234
#> 3 EY(3)-EY(2)    2    3  17.111111 611.122085     434.247037
#> --- Matching on 'multinom' with M=1, J=1 ---
```

Using ordinal logistic regression:

```r
## v0.1
multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
#> $tauestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>      -9.375       5.875      15.250 
#> 
#> $varestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>    7.794922  582.654297  576.304688 
#> 
#> $varestimateAI2012
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>    5.072057  383.848575  430.978089 
#> 
#> $analysisidx
#> 1 2 4 5 6 7 8 9 
#> 1 2 4 5 6 7 8 9

## v0.1.0.9000+
multiMatch(Y,W,X,trimming = 0, match_on = "multinom")
#> no unit IDs supplied; unit_ids will be assigned generically
#> Warning in prepareData(Y = Y, W = W, X = X, match_on = match_on, trimming
#> = trimming, : It is recommended that X is a matrix. X will be coerced to a
#> matrix.
#> -------------- Causal estimates ---------------
#>         Param Trt1 Trt2   Estimate   Variance VarianceAI2016
#> 1 EY(2)-EY(1)    1    2 -10.444444   8.545953       8.302024
#> 2 EY(3)-EY(1)    1    3   6.666667 616.913580     411.456234
#> 3 EY(3)-EY(2)    2    3  17.111111 611.122085     434.247037
#> --- Matching on 'multinom' with M=1, J=1 ---
```

#### Stratifying on Estimated Propensity Score

There are some additional arguments for using the stratification method


```r
NS <- 5 ## The number of strata to divide into
linearp <- FALSE ## Use subclassification, not linear prediction
nboot <- 10 ## Number of bootstrap samples for variance estimation
```


We will use the dataset provided with this package


```r
data("simulated_data", package = "multilevelMatching")
head(simulated_data)
#>      outcome treatment     covar1     covar2     covar3      covar4
#> 1 -5.1263661         1 -0.8655326  0.2365470  0.2296276 -2.89128973
#> 2 -3.0291972         1  0.2705180 -0.3528553 -0.4021112 -2.21465345
#> 3  3.0503040         1  1.4215891  1.3168653 -1.1991933  0.05568702
#> 4 -6.0918947         1 -1.3918079 -1.0725496  1.1201373 -2.36120362
#> 5 -2.4552555         1 -1.1465895  0.9467882  0.5100671 -2.67835507
#> 6 -0.3245296         1  0.2675347  0.4166909 -0.4467756  2.42632589
#>        covar5 covar6
#> 1 0.211508226      0
#> 2 0.070739054      1
#> 3 1.258638639      1
#> 4 0.003610748      0
#> 5 0.069593990      0
#> 6 0.604489077      0
```



```r
plot(
  x = simulated_data$treatment, 
  y = simulated_data$outcome, 
  xlab = "treatment level",
  ylab = "outcome",
  main = "Outcomes by treatment level"
)
```

![plot of chunk boxplot](man/figures/boxplot-1.png)


```r
multilevelGPSStratification(
  Y = simulated_data$outcome,
  W = simulated_data$treatment,
  X = as.matrix(simulated_data[,-(1:2)]),
  GPSM = "multinomiallogisticReg",
  NS = NS,
  linearp = linearp,
  nboot = nboot
)
#> $tauestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>  -0.2593016   0.2684324   0.5277340 
#> 
#> $varestimate
#> EY(2)-EY(1) EY(3)-EY(1) EY(3)-EY(2) 
#>   0.0752562   0.1793161   0.2744478
```


# News

See [the News site](https://shuyang1987.github.io/multilevelMatching/news/index.html) for the changelog (and updates since version 0.1).

#### A note on `multiMatch()`

The `multiMatch()` function may return slightly different estimates than the original 2 matching functions in certain circumstances. We attempt to ensure that the functions implement are identical methods up to perhaps random number generation. Please file an issue if you have any questions or concerns.


```r
fit <- multiMatch(Y,W,X,trimming = 0, match_on = "multinom")
#> no unit IDs supplied; unit_ids will be assigned generically
#> Warning in prepareData(Y = Y, W = W, X = X, match_on = match_on, trimming
#> = trimming, : It is recommended that X is a matrix. X will be coerced to a
#> matrix.

names(fit)
#>  [1] "results"           "analysis_idx"      "mu"               
#>  [4] "impute_mat"        "estimate_args"     "model"            
#>  [7] "propensity_scores" "impute_match_data" "impute_mat_sorted"
#> [10] "AI2016_args"

summary(fit)
#> ------------- Method arguments --------------
#> $match_on
#> [1] "multinom"
#> 
#> $model_options
#> $model_options$reference_level
#> [1] 1
#> 
#> 
#> $M_matches
#> [1] 1
#> 
#> $J_var_matches
#> [1] 1
#> 
#> $trt_levels
#> [1] 1 2 3
#> 
#> $N_per_trt
#> W
#> 1 2 3 
#> 4 3 2 
#> 
#> ------------- Causal estimates --------------
#>         Param Trt1 Trt2   Estimate   Variance VarianceAI2016
#> 1 EY(2)-EY(1)    1    2 -10.444444   8.545953       8.302024
#> 2 EY(3)-EY(1)    1    3   6.666667 616.913580     411.456234
#> 3 EY(3)-EY(2)    2    3  17.111111 611.122085     434.247037
#> ---------------------------------------------
```
