Propensity Score Matching and Subclassification in Observational Studies with Multi-level Treatments
=================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/multilevelMatching)](https://cran.r-project.org/package=multilevelMatching)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build Status](https://travis-ci.org/shuyang1987/multilevelMatching.svg?branch=master)](https://travis-ci.org/shuyang1987/multilevelMatching)
[![AppveyorCI Build status](https://ci.appveyor.com/api/projects/status/eu7vlcbu2j854cdo?svg=true)](https://ci.appveyor.com/project/BarkleyBG/multilevelmatching-3hh85)
[![Coverage status](https://codecov.io/gh/shuyang1987/multilevelMatching/branch/master/graph/badge.svg)](https://codecov.io/github/shuyang1987/multilevelMatching?branch=master)


# multilevelMatching

In setting with Multi-level treatments, our goal is to estimate pairwise average treatment effects from a common population using matching methods.

This goal can not be acheived by matching one treatment with another one at a time, since the pairwise matched samples may differ from the target population systematically, and thus they are not compatitable. One implication is that from this approach, it is possible that treatment A is better than treatment B, treatment B is better than treatment C, and treatment C is better than treatment A. 

We focus on estimating the average values of potential outcomes for each treatment level by matching methods, which facilitate estimation of pairwise average treatment effects for a common population.

The estimation methods include generalized propensity score (GPS) matching, GPS stratification,
matching with the full set of covariates, matching with the full set of GPS vector. Note that GPS matching and GPS straticication only require matching on a scalar function when estimating the average value of the potential outcome at a particular treatment level, which reduces the matching dimension to one, regardless of the number of covariates and the number of treatment levels. 

In order to ensure sufficient overlap, Crump et al. (2009)'s trimming method can be extended to this setting as well. 

### installation with `devtools`:

```S
devtools::install_github("shuyang1987/multilevelMatching")
```

### Use

- In the development version (v0.1.0.9000+), the `multiMatch()` function was introduced to combine the `multilevelMatchX()` and `multilevelGPSMatch()` functions. For stratification on the propensity score, use `multilevelGPSStratification`.

```S
X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y<-c(102,105,120,130,100,80,94,108,96)
W<-c(1,1,1,3,2,3,2,1,2)

multilevelMatchX(Y,W,X)
multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")

multiMatch(Y,W,X, match_on = "covariates")$results
multiMatch(Y,W,X,trimming = 0, match_on = "multinom")$results
multiMatch(Y,W,X,trimming = 1, match_on = "multinom")$results
```
```S
set.seed(111)
n    <- 5000*6
# X1-X3 3 MVN var 2, 1, 1, covars 1, -1, -.5
vars   <- c(2,1,1)
covars <- c(1,-1,-.5)
mu     <- c(0,0,0)
tau    <- 1
Sigma <- diag(vars)
Sigma[2,1] <- Sigma[1,2] <- covars[1]
Sigma[3,1] <- Sigma[1,3] <- covars[2]
Sigma[3,2] <- Sigma[2,3] <- covars[3]
trt1 <- 100; trt1
trt2 <- 100; trt2
trt3 <- 100; trt3
# draw Xs
X13 <- mvrnorm(n,mu=mu,Sigma=Sigma, empirical = FALSE)
X1 <- X13[,1]
X2 <- X13[,2]
X3 <- X13[,3]
X4 <- runif(n,-3,3)
X5 <- rchisq(n, df=1)
X6 <- rbinom(n,size=1,prob=.5)

xb2 <- 0.1*(X1^2+X2+X3+X4+X5+X6)
xb3 <- 0.1*(X1+X2^2+X3^2+X4+X5+X6)
exb2<-exp(xb2)
exb3<-exp(xb3)
pi1<-1/(1+exp(xb2)+exp(xb3))
pi2<-exp(xb2)/(1+exp(xb2)+exp(xb3))
pi3<-exp(xb3)/(1+exp(xb2)+exp(xb3))
pi<-cbind(pi1,pi2,pi3)
apply(pi,2,mean)

W<-matrix(NA,n,4)
colnames(W)   <- c("W1","W2","W3","W")
for(kk in 1:n){
    W[kk,1:3]<-rmultinom(1, 1, prob = pi[kk,])
}

sim.dat <- data.frame(W,X1,X2,X3,X4,X5,X6)
trt1.keep <- sample(which(sim.dat$W1==1),trt1,replace=FALSE)
trt2.keep <- sample(which(sim.dat$W2==1),trt2,replace=FALSE)
trt3.keep <- sample(which(sim.dat$W3==1),trt3,replace=FALSE)
sim.dat <- sim.dat[c(trt1.keep,trt2.keep,trt3.keep),]
sim.dat[,"W"]<-sim.dat[,"W1"]+2*sim.dat[,"W2"]+3*sim.dat[,"W3"]
sim.dat[,"W"]<-as.factor(sim.dat[,"W"])
W <- sim.dat[,"W"]
X <- as.matrix(sim.dat[,names(sim.dat)[-c(1:4)]])
X1 <- X[,"X1"]; X2 <- X[,"X2"]; X3 <- X[,"X3"]; X4 <- X[,"X4"]; X5 <- X[,"X5"];X6 <- X[,"X6"]

# outcome: treatment effect is zero
u  <- rnorm(nrow(X))
# ouctome (linear)
Y <- 	(W==1)*(  X1 +   X2 +   X3 +   X4 +    X5-1 +     X6-0.5)+
(W==2)*(2*X1 + 3*X2 +   X3 + 2*X4 + 2*(X5-1) + 2*(X6-0.5))+
(W==3)*(3*X1 +   X2 + 2*X3 -   X4 -   (X5-1) -   (X6-0.5))+u

match1  <- multilevelMatchX(Y,W,X)
match1b <- multiMatch(Y,W,X, match_on="covariates")
match2  <- multilevelGPSMatch(Y,W,X,Trimming=FALSE,GPSM="multinomiallogisticReg") 
match2b  <- multiMatch(Y,W,X,Trimming=FALSE,match_on="multinom") 
match3  <- multilevelGPSMatch(Y,W,X,Trimming=TRUE,GPSM="multinomiallogisticReg") 
match3b  <- multiMatch(Y,W,X,Trimming=TRUE,match_on="multinom") 

match1$results
match1b$results
match2$results
match2b$results
match3$results
match3b$results

strat1  <- multilevelGPSStratification(Y,W,X,NS=10,GPSM="multinomiallogisticReg",linearp=0,nboot=50)
strat1$results
```

# News

See [NEWS.md](NEWS.md) for changes since version 0.1.

#### A note on `multiMatch()`

The `multiMatch()` function may return slightly different estimates than the original 2 matching functions in certain circumstances. We attempt to ensure that the functions implement are identical methods up to perhaps random number generation. Please file an issue if you have any questions or concerns.




