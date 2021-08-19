##################################################################
### Kernel Machine Difference (KMD) Test for Mediation Effect 
### Detection from Multiple Mediators
##################################################################
#'
#' Fit the kernel machine regression models for the KMD test.
#'
#' @export
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param A a vector of exposure data of length \code{n}.
#' @param M an \code{n}-by-\code{p} matrix of mediator data. Each row represents an observation and each column represents a mediator.
#' @param X an \code{n}-by-\code{q} matrix of covariate data where each row represents an observation and each column represents a covariate.
#' @param knls an option of the kernel choice. Currently implemented for \code{linear} and \code{gaussian}. Data adaptive kernel selection is not available in this version.
#' @return a list of summary results with point estimates, test statistic and p-value for both NIE and NDE 
#' @references Shen, J, Schwartz J, Baccarelli AA. and Lin X (2021+) Causal Mediation Analysis of Multiple Mediators Using Kernel Machine Regression Bayesian Kernel Machine Regression. unpublished.
#' @import utils

### loading packages and functions
source("lnKN.R")
source("gsKN.R")

#library(data.table)
#library(GMMAT)

KMDtest <- function (y, A, M, X = NULL, knls = "linear", tol = 1e-7, taumin = 1e-7, taumax = 1e3) 
{
  ## Input Check
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  stopifnot (is.numeric(A), length(A) == length(y), anyNA(A) == FALSE)
  if (!inherits(M, "matrix"))  M <- as.matrix(M)
  stopifnot (is.numeric(M), nrow(M) == length(y), anyNA(M) == FALSE)
  if (!inherits(X, "matrix"))  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE) 

  dat.all <- data.frame(cbind(y,A))
  names(da.all) <- c("y", "A")
  dat.mc <- as.matrix(cbind(M,C))
  mat.mc <- scale(dat.mc, center = TRUE, scale = TRUE)
  mat.c <- scale(C, center = TRUE, scale = TRUE)
  x.mat <- as.matrix(cbind(rep(1,length(y)), A))

  if (!knl %in% c("linear", "gaussian")) {
    stop("family", family, "not yet implemented; must specify either 'linear' or 'gaussian'")
  }
  if (family == "linear") {
    message("Fitting model with linear kernels")
    ful.cov <- ln_kn(mat.mc)
    red.cov <- ln_kn(mat.c)
  }
  if (family == "gaussian") {
    message("Fitting model with gaussian kernels")
    ful.cov <- gs_kn(mat.mc, 1)
    red.cov <- gs_kn(mat.c, 1)
  }

  ## full outcome model
  modf <- glmmkin(y ~ A, data = dat.all, kins = ful.cov,
                   family = gaussian(link = "identity"), tol = tol,
                   taumin = taumin, taumax = taumax)

  ful.sigma2 <- modf$theta[1]
  if (ful.sigma2 < 1e-7) {
    message("zero sigma in the full outcome model!")
  } else {
    if (modf$theta[2] < 1e-7) {
      ful.tau <- 1e-7
    } else {
      ful.tau <- modf$theta[2]
    }
    ful.lami <- ful.tau/ful.sigma2 ## inverse lambda
    v.mat <- solve(diag(nsubs) + (ful.cov * ful.lami))
    modf$alpha <- solve(t(x.mat) %*% v.mat %*% x.mat)%*%(t(x.mat) %*% v.mat %*% y)
    ful.betas <- modf$alpha[2]
    ful.alpha <- ful.lami * (v.mat %*% (y - (x.mat)%*%(modf$alpha)))
    ful.h <- t(t(ful.alpha) %*% ful.cov)
    res.f <- y - (x.mat)%*%(modf$alpha) - ful.h
    Amat.f <- v.mat %*% ((ful.cov*ful.lami) + (x.mat)%*%solve(t(x.mat) %*% v.mat %*% x.mat)%*% t(x.mat) %*% v.mat)
    tr.f <- matrix.trace(Amat.f)
    kmAIC <- nsubs * log(t(res.f) %*% res.f) + 2 * tr.f
  }

  ## reduced outcome model
  rx.mat <- cbind(rep(1,length(y)), A)
  modr <- glmmkin(y ~ A, data = dat.all, kins = red.cov,
                   family = gaussian(link = "identity"), tol = tol,
                   taumin = taumin, taumax = taumax)
  red.sigma2 <- modr$theta[1]
  if (red.sigma2 < 1e-7) {
    message("zero sigma in the reduced outcome model!")
  } else {
    if (modr$theta[2] < 1e-7) {
      red.tau <- 1e-7
    } else {
      red.tau <- modr$theta[2]
    }
    red.lami <- red.tau/red.sigma2
    rv.mat <- solve(diag(nsubs) + (red.cov*red.lami))
    modr$alpha <- solve(t(rx.mat) %*% rv.mat %*% rx.mat)%*%(t(rx.mat) %*% rv.mat %*% y)
    red.betas <- modr$alpha[2]
    red.alpha <- red.lami * (rv.mat %*% (y - (rx.mat)%*%(modr$alpha)))
    red.h <- t(t(red.alpha) %*% red.cov)
    res.r <- y - (rx.mat)%*%(modr$alpha) - red.h
    Amat.r <- rv.mat %*% ((red.cov*red.lami) + (rx.mat)%*%solve(t(rx.mat) %*% rv.mat %*% rx.mat)%*% t(rx.mat) %*% rv.mat)
    tr.r <- matrix.trace(Amat.r)
    kmAIC.r <- nsubs * log(t(res.r) %*% res.r) + 2 * tr.r
  }

  res.f2 <- y - (x.mat)%*%(modf$alpha)
  res.r2 <- y - (rx.mat)%*%(modr$alpha)

  Af <- (t(x.mat) %*% v.mat %*% x.mat)
  Ar <- t(rx.mat) %*% rv.mat %*% rx.mat
  P.f <- solve(Af) %*% t(x.mat) %*% v.mat
  P.r <- solve(Ar) %*% t(rx.mat) %*% rv.mat
  P.diff <- P.r - P.f
  kmfr.var.i <- (P.r %*% t(P.r))[2,2] * max(r.sig2-f.sig2,0) + (P.diff %*% t(P.diff))[2,2] * f.sig2
  kmfr.var.d <- (P.f %*% t(P.f))[2,2] * f.sig2

  # Wald test
  km.nde <- f.beta[2]
  tmpt.fr.d <- km.nde/sqrt(kmfr.var.d)
  kmfr.p.d <- pnorm(abs(tmpt.fr.d), lower.tail = F, log.p = FALSE)*2
  kmfr.rej.d <- as.numeric(kmfr.p.d < 0.05)

  km.nie <- r.beta[2] - f.beta[2]
  tmpt.fr <- km.nie/sqrt(kmfr.var.i)
  kmfr.p.i <- pnorm(abs(tmpt.fr), lower.tail = F, log.p = FALSE)*2
  kmfr.rej.i <- as.numeric(kmfr.p.i < 0.05)

  ## summary of results
  res <- list(NIE_est=km.nie, NIE_se=sqrt(kmfr.var.i), T_kmd=tmpt.fr, P_kmd=kmfr.p.i,
              NDE_est=km.nde, NDE_se=sqrt(kmfr.var.d), T_nde=tmpt.fr.d, P_kmd=kmfr.p.d)
  return(res)
}