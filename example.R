######################################################
# This is an example of KMD test 
#      for the joint NIE of multiple mediators
######################################################

source("KMDtest.R")

# load packages
library(mvtnorm)
library(GMMAT)

# Generat Example Data
set.seed(100)
nsubs <- 200
p <- 20 # number of mediators
q <- 5 # number of confounders

## expit function for logit
getp <- function(x)
{
  exp(x)/(1+exp(x))
}

beta0 <- 1
beta1 <- 2
#gamma.c <- c(0.6, 1, -1, 1, -0.6, -0.6)
sd.y <- 1     ## error term y
sd.m <- 1     ## error term M
sd.c <- 1     ## error term C
betas_c <- c(1,3,-1,2,1)
betas_m <- numeric(p)
betas_m[1:10] <- c(1,0,1,0,-1,0,1,0,0,1)

## set up coefficients for M models
alpha0 <- rep(c(1,0,-1,0), 5) ##1*200
alphas <- c(1, rep(0,19))
alphac <- matrix(0, nrow = q, ncol = p)
alphac[1,1:10] <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
alphac[2,1:10] <- c(1, 0, 0, 0, 1, 1, 1, 0, 0, 0)
alphac[3,1:10] <- c(0, 1, 0, 0, 1, 0, 0, 1, 1, 0)
alphac[4,1:10] <- c(0, 0, 1, 0, 0, 1, 0, 1, 0, 1)
alphac[5,1:10] <- c(0, 0, 0, 1, 0, 0, 1, 0, 1, 1)

# set up covariance matrix
rho <- c(0.1, 0.3) # correlation for blocks
Jm1m1 <- matrix(data=rep(1, p*p), nrow=p, ncol=p)
Im1 <- diag(rep(1, p))
sigma.m <- (1-rho[1])*Im1 + rho[1]*Jm1m1

Jm2m2 <- matrix(data=rep(1, q*q), nrow=q, ncol=q)
Im2 <- diag(rep(1, q))
sigma.c <- (1-rho[2])*Im2 + rho[2]*Jm2m2

hfn <- function(z)
{
  #  4*(z[1]) - (2*z[2]^2) + exp(-z[3])*z[4] - 8*sin(z[2])*cos(z[3])
  #  + log(2*exp(z[4])+5*exp(z[3])) + z[5]
  1*z[6] + 1*z[8] - 1*z[10] + 1*z[12] + 1*z[15]
}
int.all <- rep(1, nsubs)
C <- rmvnorm(nsubs, sigma=sigma.c*sd.c)
p.all <- getp(rep(0.6,nsubs))
A <- sapply(p.all, rbinom, n=1, size=1)

rand.em <- rmvnorm(nsubs, sigma=sigma.m*sd.m)
M <- rep(1,nsubs) %*% t(alpha0) + s.all %*% t(alphas) +
  c.all %*% alphac + rand.em

rand.ey <- rnorm(nsubs, 0, sd.y)
y <- int.all*beta0 + A*beta1 + M %*% betas_m + C %*% betas_c + rand.ey # 1000x1

dat.all <- data.frame(cbind(1:nsubs, y, A, C, M))
names.str <- c("id", "Y", "A", paste0("C", 1:q), paste0("M", 1:p))
names(dat.all) <- names.str

## perform KMD test
KMDtest(y=y, A=A, M=M, X=C, knls="linear")

##-- End of the Code