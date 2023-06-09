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


indepCopula <- function(dim = 2L) {
    ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    uis <- paste0("u", 1:n)
    expr <- paste(uis, collapse="*")
    parse(text = expr)
  }

  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste0("u", i))
    }
    val
  }
  cdf <- cdfExpr((dim <- as.integer(dim)))
  pdf <- pdfExpr(cdf, dim)

  new("indepCopula",
             dimension = dim,
             exprdist = c(cdf=cdf, pdf=pdf),
             parameters = double(0),
             param.names = character(0),
             param.lowbnd = double(0),
             param.upbnd = double(0),
             fullname = "Independence copula")
}

setMethod("A", signature("indepCopula"), function(copula, w) rep.int(1, length(w)))

setMethod("rCopula", signature("numeric", "indepCopula"),
          function(n, copula) matrix(runif(n * copula@dimension), nrow = n))

pindepCopula <- function(u, copula, log.p=FALSE) {
  stopifnot(ncol(u) == copula@dimension)
  if(log.p) rowSums(log(u)) else apply(u, 1, prod)
}

setMethod("pCopula", signature("numeric", "indepCopula"),pindepCopula)
setMethod("pCopula", signature("matrix", "indepCopula"), pindepCopula)

dindepCopula <- function(u, copula, log=FALSE, ...) {
  stopifnot(ncol(u) == copula@dimension)
  rep.int(if(log) 0 else 1, nrow(u))
}

setMethod("dCopula", signature("numeric", "indepCopula"),dindepCopula)
setMethod("dCopula", signature("matrix", "indepCopula"), dindepCopula)


setMethod("tau", "indepCopula", function(copula, ...) 0)
setMethod("rho", "indepCopula", function(copula, ...) 0)
setMethod("lambda", "indepCopula",
          function(copula, ...)  c(lower=0, upper = 0))

## GETR
setMethod("describeCop", c("indepCopula", "character"),
          function(x, kind = c("short", "very short", "long")) {
    kind <- match.arg(kind)
    if(kind == "very short") # e.g. for show() which has more parts
        return("independence copula")
    ## else
    d <- dim(x)
    ch <- paste("independence copula, dim. d =", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = ch,
           stop("invalid 'kind': ", kind))
})
