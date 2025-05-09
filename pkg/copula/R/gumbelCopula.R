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


AGumbel <- function(copula, w) .AGumbel(w, copula@parameters[1])
.AGumbel <- function(w, alpha) {
  r <- w
  w <- w[i <- which(!(bnd <- w == 0 | w == 1))]
  r[i] <- (w^alpha + (1 - w)^alpha)^(1/alpha)
  r[bnd] <- 1
  r
}

dAduGumbel <- function(copula, w) {
  alpha <- copula@parameters[1]

  ## deriv(expression((w^alpha + (1 - w)^alpha)^(1/alpha)), "w", hessian=TRUE):
  value <- eval(expression({
    .expr2 <- 1 - w
    .expr4 <- w^alpha + .expr2^alpha
    .expr5 <- 1/alpha
    .expr7 <- .expr5 - 1
    .expr8 <- .expr4^.expr7
    .expr9 <- alpha - 1
    .expr14 <- w^.expr9 * alpha - .expr2^.expr9 * alpha
    .expr15 <- .expr5 * .expr14
    .expr22 <- .expr9 - 1
    .value <- .expr4^.expr5
    .grad <- array(0, c(length(.value), 1L), list(NULL, "w"))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL, "w", "w"))
    .grad[, "w"] <- .expr8 * .expr15
    .hessian[, "w", "w"] <- .expr4^(.expr7 - 1) * (.expr7 * .expr14) *
        .expr15 + .expr8 * (.expr5 * (w^.expr22 * .expr9 * alpha +
        .expr2^.expr22 * .expr9 * alpha))
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha=alpha, w=w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1=der1, der2=der2)
}



gumbelCopula <- function(param = NA_real_, dim = 2L,
			 use.indepC = c("message", "TRUE", "FALSE"))
{
  stopifnot(length(param) == 1)
  if(!is.na(param) && param == 1) {
      use.indepC <- match.arg(use.indepC)
      if(!identical(use.indepC, "FALSE")) {
	  if(identical(use.indepC, "message"))
	      message("parameter at boundary ==> returning indepCopula()")
	  return( indepCopula(dim=dim) )
      }
  }

  ## get expressions of cdf and pdf
  cdfExpr <- function(d) {
    expr <- "( - log(u1))^alpha"
    for (i in 2:d) {
      cur <- paste0( "(-log(u", i, "))^alpha")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("exp(- (", expr, ")^ (1/alpha))")
    parse(text = expr)
  }

  pdfExpr <- function(cdf, d) {
    val <- cdf
    for (i in 1:d) {
      val <- D(val, paste0("u", i))
    }
    val
  }

  cdf <- cdfExpr((dim <- as.integer(dim)))
  pdf <- if (dim <= 6) pdfExpr(cdf, dim) # else NULL
  new("gumbelCopula",
      dimension = dim,
      parameters = param,
      exprdist = c(cdf = cdf, pdf = pdf),
      param.names = "alpha",
      param.lowbnd = 1,
      param.upbnd = Inf,
      fullname = "<deprecated slot>")# "Gumbel copula family; Archimedean copula; Extreme value copula"
}


rgumbelCopula <- function(n, copula) {
  ## frailty is stable(1,0,0) with 1/alpha
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if(is.na(alpha <- copula@parameters[1])) stop("parameter is NA")
  ## reduce to indepCopula
  if (alpha - 1 < .Machine$double.eps ^(1/3) )
      return(rCopula(n, indepCopula(dim=dim)))
  b <- 1/alpha
  ## stable (alpha = b, beta = 1, gamma = **), 0 < b < 1
  fr <- if(identical(getOption("copula:rstable1"), "rPosStable")) ## back compatible
      rPosStableS(n, b) else rstable1(n, alpha = b, beta = 1,
				      gamma = cospi2(b)^alpha, pm=1)
  ## now gumbel copula
  val <- matrix(runif(dim * n), nrow = n)
  psi(copula, - log(val) / fr)
  ## = copGumbel@psi(- log(val) / fr, alpha)
}


pgumbelCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1] # needed in eval():
  pmax(eval(copula@exprdist$cdf), 0)
}

dgumbelCopula <- function(u, copula, log=FALSE, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  dim <- copula@dimension
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  if(log) stop("'log=TRUE' not yet implemented")
  alpha <- copula@parameters[1] # needed in eval():
  eval(copula@exprdist$pdf)
}

if(FALSE) # not used anymore -- rather copGumbel@dacopula
dgumbelCopula.pdf <- function(u, copula, log=FALSE) {
  dim <- copula@dimension
  if (dim > 10) stop("Gumbel copula PDF not implemented for dimension > 10.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1]
  ## FIXME: improve log-case
  if(log)
    log(c(eval(gumbelCopula.pdf.algr[dim])))
  else  c(eval(gumbelCopula.pdf.algr[dim]))
}


iTauGumbelCopula <- function(copula, tau) {
    if (any(neg <- tau < 0)) {
        warning("For the Gumbel copula, tau must be >= 0. Replacing negative values by 0.")
        tau[neg] <- 0
    }
    1/(1 - tau)
}


gumbelRhoFun <- function(alpha) {
  theta <- .gumbelRho$trFuns$forwardTransf(alpha, .gumbelRho$ss)
  c(.gumbelRho$assoMeasFun$valFun(theta))
}

gumbeldRho <- function(alpha) {
  ss <- .gumbelRho$ss
  forwardTransf <- .gumbelRho$trFuns$forwardTransf
  forwardDer <- .gumbelRho$trFuns$forwardDer
  valFun <- .gumbelRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(valFun(theta, 1)) * forwardDer(alpha, ss)
}


iRhoGumbelCopula <- function(copula, rho) {
    if (any(neg <- rho < 0)) {
        warning("For the Gumbel copula, rho must be >= 0. Replacing negative values by 0.")
        rho[neg] <- 0
    }
  gumbelRhoInv <- approxfun(x = .gumbelRho$assoMeasFun$fm$ysmth,
                            y = .gumbelRho$assoMeasFun$fm$x, rule = 2)
  ss <- .gumbelRho$ss
  theta <- gumbelRhoInv(rho)
  .gumbelRho$trFuns$backwardTransf(theta, ss)
}


setMethod("rCopula", signature("numeric", "gumbelCopula"), rgumbelCopula)

setMethod("pCopula", signature("matrix", "gumbelCopula"),
	  ## was  pgumbelCopula
	  function(u, copula, ...) .pacopula(u, copGumbel, theta=copula@parameters))
setMethod("dCopula", signature("matrix", "gumbelCopula"),
	  ## was  dgumbelCopula.pdf
	  function (u, copula, log = FALSE, checkPar=TRUE, ...)
	  copGumbel@dacopula(u, theta=copula@parameters, log=log, checkPar=checkPar, ...))


setMethod("A", signature("gumbelCopula"), AGumbel)
setMethod("dAdu", signature("gumbelCopula"), dAduGumbel)

setMethod("iPsi", signature("gumbelCopula"),
	  function(copula, u) copGumbel@iPsi(u, theta=copula@parameters))
setMethod("psi", signature("gumbelCopula"),
	  function(copula, s) copGumbel@psi(s, theta=copula@parameters))

setMethod("diPsi", signature("gumbelCopula"),
	  function(copula, u, degree=1, log=FALSE, ...) {
    s <- if(log || degree %% 2 == 0) 1. else -1.
    s* copGumbel@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
})


setMethod("tau", "gumbelCopula", function(copula) 1 - 1/copula@parameters)

setMethod("rho", "gumbelCopula", function(copula) gumbelRhoFun(copula@parameters))

setMethod("lambda", signature("gumbelCopula"), function(copula, ...)
    c(lower = 0, upper = 2 - 2^(1/copula@parameters)))

setMethod("iTau", signature("gumbelCopula"), function(copula, tau, ...) iTauGumbelCopula(copula, tau))
setMethod("iRho", signature("gumbelCopula"), function(copula, rho, ...) iRhoGumbelCopula(copula, rho))

setMethod("dRho", "gumbelCopula", function(copula) gumbeldRho(copula@parameters))

setMethod("dTau", "gumbelCopula", function(copula) 1 / copula@parameters^2)
