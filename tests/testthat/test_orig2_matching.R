
context("test_orig2 match + stratificaiton output from v0.1.0")

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
trt1 <- 100#; trt1
trt2 <- 100#; trt2
trt3 <- 100#; trt3
# draw Xs
X13 <- MASS::mvrnorm(n,mu=mu,Sigma=Sigma, empirical = FALSE)
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
# apply(pi,2,mean)

W<-matrix(NA,n,4)
colnames(W)   <- c("W1","W2","W3","W")
for(kk in 1:n){
  W[kk,1:3]<-rmultinom(1, 1, prob = pi[kk,])
}

sim_dat <- data.frame(W,X1,X2,X3,X4,X5,X6)
trt1_keep <- sample(which(sim_dat$W1==1),trt1,replace=FALSE)
trt2_keep <- sample(which(sim_dat$W2==1),trt2,replace=FALSE)
trt3_keep <- sample(which(sim_dat$W3==1),trt3,replace=FALSE)
sim_dat <- sim_dat[c(trt1_keep,trt2_keep,trt3_keep),]
sim_dat[,"W"]<-sim_dat[,"W1"]+2*sim_dat[,"W2"]+3*sim_dat[,"W3"]
sim_dat[,"W"]<-as.factor(sim_dat[,"W"])
W <- sim_dat[,"W"]
X <- as.matrix(sim_dat[,names(sim_dat)[-c(1:4)]])
X1 <- X[,"X1"]; X2 <- X[,"X2"]; X3 <- X[,"X3"]; X4 <- X[,"X4"]; X5 <- X[,"X5"];X6 <- X[,"X6"]

# outcome: treatment effect is zero
u  <- rnorm(nrow(X))
# ouctome (linear)
Y <- 	(W==1)*(  X1 +   X2 +   X3 +   X4 +    X5-1 +     X6-0.5)+
  (W==2)*(2*X1 + 3*X2 +   X3 + 2*X4 + 2*(X5-1) + 2*(X6-0.5))+
  (W==3)*(3*X1 +   X2 + 2*X3 -   X4 -   (X5-1) -   (X6-0.5))+u

id_vals <- paste0("unitID", 1:length(W))
names(W) <- id_vals
# match1<-multilevelMatchX(Y,W,X)
# match2<-multilevelGPSMatch(Y,W,X,Trimming=FALSE,GPSM="multinomiallogisticReg")
# match3<-multilevelGPSMatch(Y,W,X,Trimming=TRUE,GPSM="multinomiallogisticReg")
# match4<-multilevelGPSStratification(Y,W,X,NS=10,GPSM="multinomiallogisticReg",linearp=0,nboot=50)

test_that(
  "matchX returns same output as original",  {
    set.seed(22)
    fit <- multilevelMatchX(Y,W,X)


    fit_orig <-
      structure(
        list(
          tauestimate = structure(
            c(0.0792736112950251,
              0.862649285149429, 0.783375673854404),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.179218593368605,
              0.163475380955354, 0.322161565872657),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          )
        ),
        .Names = c("tauestimate", "varestimate")
      )

    if (names(fit)[[1]] == "tauestimate"){

    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )


    } else {
      expect_equal(
        fit$results$Estimate,
        fit_orig$tauestimate,
        tol=1e-5,
        check.names = FALSE
      )

      expect_equal(
        fit$results$Variance,
        fit_orig$varestimate,
        tol=1e-5,
        check.names = FALSE
      )
    }

  }
)

test_that(
  "GPSMatch on MLR returns same output as original",  {
    set.seed(22)
    fit <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="multinomiallogisticReg")
    set.seed(22)
    fit2 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="multinomiallogisticReg")


    fit_orig <- structure(list(

      tauestimate = structure(
        c(-0.745849911303934,
          0.560519194422989, 1.30636910572692),
        .Names = c("EY(2)-EY(1)",

                   "EY(3)-EY(1)", "EY(3)-EY(2)")
      ), varestimate = structure(
        c(0.69749075129289,
          0.485657947932583, 0.708962013743806),
        .Names = c("EY(2)-EY(1)",
                   "EY(3)-EY(1)", "EY(3)-EY(2)")
      ), varestimateAI2012 = structure(
        c(0.434057084711074,
          0.228778686381817, 0.422646595835199),
        .Names = c("EY(2)-EY(1)",
                   "EY(3)-EY(1)", "EY(3)-EY(2)")
      ), analysisidx = 1:300), .Names = c("tauestimate",
                                          "varestimate",
                                          "varestimateAI2012",
                                          "analysisidx"))

    fit2_orig <-
      structure(
        list(
          tauestimate = structure(
            c(-1.05644469547957, 1.0011894314108,
              2.05763412689036),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.792431016304647, 0.445723871851535,
              0.831026403846409),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)",
                       "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(0.487639607313493,
              0.21046030940101, 0.462406885897202),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          analysisidx = structure(
            as.integer(1:300)[-c(192, 207)],
            .Names = id_vals[-c(192, 207)]
          )
        ),
        .Names = c(
          "tauestimate",
          "varestimate",
          "varestimateAI2012",
          "analysisidx"
        )
      )
    if (names(fit)[[1]] == "tauestimate") {



    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit$varestimateAI2012,
      fit_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit$analysisidx,
      fit_orig$analysisidx,
      tol=1e-5
    )


    expect_equal(
      fit2$tauestimate,
      fit2_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$varestimate,
      fit2_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit2$varestimateAI2012,
      fit2_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit2$analysisidx,
      fit2_orig$analysisidx,
      tol=1e-5
    )



    } else {
      expect_equal(
        fit$results$Estimate,
        check.attributes = FALSE,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$results$Variance,
        check.names = FALSE,
        fit_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit$results$VarianceAI2012,
        check.names = FALSE,
        fit_orig$varestimateAI2012,
        tol=1e-5
      )

      ## new code outputs NULL for analysis_idx when no trimming
      expect_failure(
        expect_equal(
          fit$analysis_idx$indices_kept,
          # check.names = FALSE,
          fit_orig$analysisidx,
          tol=1e-5
        )
      )


      expect_equal(
        fit2$results$Estimate,
        check.names = FALSE,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$Variance,
        check.names = FALSE,
        fit2_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit2$results$VarianceAI2012,
        check.names = FALSE,
        fit2_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit2$analysis_idx$indices_kept,
        fit2_orig$analysisidx,
        tol=1e-5
      )
    }

  }
)

test_that(
  "GPSMatch on POLR returns same output as original",  {
    set.seed(22)
    fit <- multilevelGPSMatch(Y,W,X,Trimming=0,GPSM="ordinallogisticReg")
    set.seed(22)
    fit2 <- multilevelGPSMatch(Y,W,X,Trimming=1,GPSM="ordinallogisticReg")


    fit_orig <-
      structure(
        list(
          tauestimate = structure(
            c(-0.730380965999445,
              0.371149118665736, 1.10153008466518),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.726144991733452,

              0.596600530693367, 0.974334076537271),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(NA,
              NA, NA),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          analysisidx = 1:300
        ),
        .Names = c("tauestimate",
                   "varestimate",
                   "varestimateAI2012",
                   "analysisidx"))
    fit2_orig <-


      structure(
        list(
          tauestimate = structure(
            c(-0.769329621053784,
              0.444571526855366, 1.21390114790915),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.718662789168262,
              0.468335581580708, 0.897831146664908),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimateAI2012 = structure(
            c(NA,
              NA, NA),
            .Names = c("EY(2)-EY(1)", "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          analysisidx = structure(
            as.integer(1:300)[-c(192, 207)],
            .Names = id_vals[-c(192, 207)]
          )
        ),
        .Names = c(
          "tauestimate",
          "varestimate",
          "varestimateAI2012",
          "analysisidx"
        )
      )

    if (names(fit)[[1]] == "tauestimate"){


    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit$varestimateAI2012,
      fit_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit$analysisidx,
      fit_orig$analysisidx,
      tol=1e-5
    )


    expect_equal(
      fit2$tauestimate,
      fit2_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$varestimate,
      fit2_orig$varestimate,
      tol=1e-5
    )
    expect_equal(
      fit2$varestimateAI2012,
      fit2_orig$varestimateAI2012,
      tol=1e-5
    )

    expect_equal(
      fit2$analysisidx,
      fit2_orig$analysisidx,
      tol=1e-5
    )



    } else {
      expect_equal(
        fit$results$Estimate,
        fit_orig$tauestimate,
        check.attributes = FALSE,
        tol=1e-5
      )

      expect_equal(
        fit$results$Variance,
        check.names = FALSE,
        fit_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit$results$VarianceAI2012,
        check.names = FALSE,
        fit_orig$varestimateAI2012,
        tol=1e-5
      )

      ## new code outputs NULL for analysis_idx when no trimming
      expect_failure(
        expect_equal(
          fit$analysis_idx$indices_kept,
          # check.names = FALSE,
          fit_orig$analysisidx,
          tol=1e-5
        )
      )


      expect_equal(
        fit2$results$Estimate,
        check.names = FALSE,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$Variance,
        check.names = FALSE,
        fit2_orig$varestimate,
        tol=1e-5
      )
      expect_equal(
        fit2$results$VarianceAI2012,
        check.names = FALSE,
        fit2_orig$varestimateAI2012,
        tol=1e-5
      )

      expect_equal(
        fit2$analysis_idx$indices_kept,
        fit2_orig$analysisidx,
        tol=1e-5
      )
    }

  }
)


test_that(
  "GPSStratification on POLR returns same output as original",  {
    set.seed(22)
    fit <- multilevelGPSStratification(
      Y,W,X, GPSM = "multinomiallogisticReg",
      NS=4, linearp = 0, nboot=3
    )
    set.seed(22)
    fit2 <- multilevelGPSStratification(
      Y,W,X, GPSM = "multinomiallogisticReg",
      NS=4, linearp = 1, nboot=3
    )
    set.seed(22)
    fit3 <- multilevelGPSStratification(
      Y,W,X, GPSM = "ordinallogisticReg",
      NS=4, linearp = 0, nboot=3
    )

    fit_orig <-

      structure(
        list(
          tauestimate = structure(
            c(-0.117854433238482,
              0.264508979403283, 0.382363412641764),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.135509387708196,
              0.0055455916082906, 0.0865332894765468),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          )
        ),
        .Names = c("tauestimate", "varestimate")
      )



    fit2_orig <-
      structure(
        list(
          tauestimate = structure(
            c(-0.408025014163206,
              0.376717569336057, 0.784742583499263),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          ),
          varestimate = structure(
            c(0.012145964464371,
              0.0700560871820945, 0.0894561274872827),
            .Names = c("EY(2)-EY(1)",
                       "EY(3)-EY(1)", "EY(3)-EY(2)")
          )
        ),
        .Names = c("tauestimate", "varestimate")
      )


    fit3_orig <-
    structure(
      list(
        tauestimate = structure(
          c(-0.0697955536457575,
            0.56457647939882, 0.634372033044577),
          .Names = c("EY(2)-EY(1)",
                     "EY(3)-EY(1)", "EY(3)-EY(2)")
        ),
        varestimate = structure(
          c(0.369160126983748,
            0.0762443762621385, 0.143100046612838),
          .Names = c("EY(2)-EY(1)",
                     "EY(3)-EY(1)", "EY(3)-EY(2)")
        )
      ),
      .Names = c("tauestimate", "varestimate")
    )


    if (names(fit)[[1]] == "tauestimate"){

    expect_equal(
      fit$tauestimate,
      fit_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit$varestimate,
      fit_orig$varestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$tauestimate,
      fit2_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit2$varestimate,
      fit2_orig$varestimate,
      tol=1e-5
    )


    expect_equal(
      fit3$tauestimate,
      fit3_orig$tauestimate,
      tol=1e-5
    )

    expect_equal(
      fit3$varestimate,
      fit3_orig$varestimate,
      tol=1e-5
    )



    } else {
      expect_equal(
        fit$results$Estimate,
        check.names = FALSE,
        fit_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit$results$Variance,
        check.names = FALSE,
        fit_orig$varestimate,
        tol=1e-5
      )

      expect_equal(
        fit$results$VarianceAI2012,
        rep(NA, 3)
      )
      expect_null(
        fit$analysis_idx$indices_kept
      )
      ## new code outputs NULL for analysis_idx when no trimming
      # expect_failure(
        # expect_equal(
        #   fit$analysis_idx$indices_kept,
        #   # check.names = FALSE,
        #   fit_orig$analysisidx,
        #   tol=1e-5
        # )
      # )


      expect_equal(
        fit2$results$Estimate,
        check.names = FALSE,
        fit2_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$Variance,
        check.names = FALSE,
        fit2_orig$varestimate,
        tol=1e-5
      )

      expect_equal(
        fit2$results$VarianceAI2012,
        rep(NA, 3)
      )
      expect_null(
        fit2$analysis_idx$indices_kept
      )
      # expect_equal(
      #   fit2$results$VarianceAI2012,
      #   check.names = FALSE,
      #   fit2_orig$varestimateAI2012,
      #   tol=1e-5
      # )

      # expect_equal(
      #   fit2$analysis_idx$indices_kept,
      #   fit2_orig$analysisidx,
      #   tol=1e-5
      # )

      expect_equal(
        fit3$results$Estimate,
        check.names = FALSE,
        fit3_orig$tauestimate,
        tol=1e-5
      )

      expect_equal(
        fit3$results$Variance,
        check.names = FALSE,
        fit3_orig$varestimate,
        tol=1e-5
      )
    }

  }
)
