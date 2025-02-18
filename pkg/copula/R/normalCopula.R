## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


normalCopula <- function(param = NA_real_, dim = 2L, dispstr = "ex") {
    if(!is.numeric(param)) storage.mode(param) <- "double" # for NA, keeping attributes!
    stopifnot((pdim <- length(param)) >= 1)
    if(pdim == 1 && is.na(param)) ## extend it (rho)
	pdim <- length(param <- rep(param, length.out = npar.ellip(dim, dispstr)))
  new("normalCopula", # validRho() in ./Classes.R checks 'dispstr' and more:
      dispstr = dispstr,
      dimension = as.integer(dim),
      parameters = param,
      param.names = paste("rho", 1:pdim, sep="."),
      param.lowbnd = lowbnd.rho.ellip(dim, dispstr, pdim),
      param.upbnd = rep(1, pdim),
      fullname = "<deprecated slot>", # "Normal copula family"
      getRho = function(obj) obj@parameters)
}


rnormalCopula <- function(n, copula, ...)
    pnorm(rmvnorm(n, sigma = getSigma(copula)))

pmvnormAlgo <- function(dim, x, checkCorr = FALSE, ...) {
    if(dim <= 3 && !anyNA(x) && (!any(xI <- x == Inf) || all(xI)))
        TVPACK(...)
    else if(dim <= 5)
        Miwa(checkCorr=checkCorr, ...)
    else
        GenzBretz(...)
}

pnormalCopula <- function(u, copula, algorithm = NULL, keepAttr=FALSE, ...)
{
  dim <- copula@dimension
  ## stopifnot(is.matrix(u), ncol(u) == dim) # not needed, as called from pCopula()
  i.lower <- rep.int(-Inf, dim)
  sigma <- getSigma(copula)
  apply(qnorm(u), 1, function(x)
      if(anyNA(x)) NA_real_
      else { # algorithm depending on 'x' ..
          if(is.null(algorithm))
              algorithm <- pmvnormAlgo(dim, x=x, ...)
          pmvnorm(lower = i.lower, upper = x,
                  sigma=sigma, algorithm=algorithm, keepAttr=keepAttr, ...)
      })
}

dnormalCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  r <- rep(NA_real_, nrow(u))
  ok <- !apply(u, 1, anyNA)
  x <- qnorm(u[ok, , drop=FALSE])
  ## work in log-scale [less over-/under-flow, then (maybe) transform:
  r[ok] <- dmvnorm(x, sigma = sigma, log=TRUE) - rowSums(dnorm(x, log=TRUE))
  ## now happens in dCopula(): -- dnormalCopula() not called directly by user
  ## if(any(out <- !is.na(u) & (u <= 0 | u >= 1)))
  ##   val[apply(out, 1, any)] <- -Inf
  if(log) r else exp(r)
}

if(FALSE)# Note this was mvtnorm::dmvnorm() up to .. 2014-03-..
    ## it had the big advantage of given 'NaN' in case of a non-pos.def. sigma
## MM: We could use the new version *when* that allows to pass
## chol(sigma) as an alternative to sigm
dmvnorm <- function (x, mean, sigma, log=FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
### MM: this is *really* poor : mahalanobis() computes solve(sigma) !
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE,
                                   only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    exp(logretval)
}



## as long we think we need print.copula(), we also need this:
print.normalCopula <- printEllipCopula
setMethod("show", signature("normalCopula"), function(object) printEllipCopula(object))


lambdaNormalCopula <- function(copula) {
  i01 <- as.numeric(copula@parameters == 1) ## == rho
  c(lower = i01, upper = i01)
}

setMethod("rCopula", signature("numeric", "normalCopula"), rnormalCopula)

setMethod("pCopula", signature("matrix", "normalCopula"), pnormalCopula)
setMethod("dCopula", signature("matrix", "normalCopula"), dnormalCopula)
## pCopula() and dCopula() *generic* already deal with non-matrix case!
## setMethod("pCopula", signature("numeric", "normalCopula"),pnormalCopula)
## setMethod("dCopula", signature("numeric", "normalCopula"),dnormalCopula)

## rho := copula@parameters
setMethod("tau", "normalCopula", function(copula, ...) 2 * asin(copula@parameters) /pi)
setMethod("rho", "normalCopula", function(copula, ...) asin(copula@parameters / 2) * 6 / pi)
setMethod("lambda", signature("normalCopula"), lambdaNormalCopula)

setMethod("iRho", "normalCopula", function(copula, rho, ...) sinpi(rho / 6) * 2)

