
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


simulated_data <- sim.dat[, -(1:2)]

names(simulated_data) <- c(
  "outcome",
  "treatment",
  paste0("covar", 1:6)
)
simulated_data$outcome <- Y
row.names(simulated_data) <- NULL

save(
  simulated_data,
  file = rprojroot::find_package_root_file("data", "simulated_data.RData")
)
