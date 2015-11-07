Propensity Score Matching and Subclassification in Observational Studies with Multi-level Treatments
=================
In setting with Multi-level treatments, our goal is to estimate pairwise average treatment effects from a common population using matching methods.

This goal can not be acheived by matching one treatment with another one at a time, since the pairwise matched samples may differ from the target population systematically, and thus they are not compatitable. One implication is that from this approach, it is possible that treatment A is better than treatment B, treatment B is better than treatment C, and treatment C is better than treatment A. 

We focus on estimating the average values of potential outcomes for each treatment level by matching methods, which facilitate estimation of pairwise average treatment effects for a common population.

The estimation methods include generalized propensity score (GPS) matching, GPS stratification,
matching with the full set of covariates, matching with the full set of GPS vector. Note that GPS matching and GPS straticication only require matching on a scalar function when estimating the average value of the potential outcome at a particular treatment level, which reduces the matching dimension to one, regardless of the number of covariates and the number of treatment levels. 

In order to ensure sufficient overlap, Crump et al. (2009)'s trimming method can be extended to this setting as well. 

### install
with `devtools`:

```S
devtools::install_github("shuyang1987/multilevelMatching")
```

### use
There are only three functions in this package. 

```S
X<-c(5.5,10.6,3.1,8.7,5.1,10.2,9.8,4.4,4.9)
Y<-c(102,105,120,130,100,80,94,108,96)
W<-c(1,1,1,3,2,3,2,1,2)
multilevelMatchX(Y,W,X)
multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")
```

`multilevelMatchX()`, `multilevelGPSMatch()`, `multilevelGPSStratification()` make super awesome illustrations. 


