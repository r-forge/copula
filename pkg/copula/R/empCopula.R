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


## DEPRECATED
Cn <- function(x, w) {
    .Deprecated("C.n")
    C.n(w, x)
}


### Empirical distribution function ############################################

##' @title Empirical distribution function with smoothing
##' @param x (m, d) matrix of evaluation points
##' @param X (n, d) matrix of data based on which the empirical distribution
##'        function is computed
##' @param smoothing smoothing method
##' @param offset scaling factor of the empirical distribution function (= n/(n + offset))
##' @param method method string ("C" for C code; "R" for R code)
##' @return empirical distribution function of X at x
##' @author Ivan Kojadinovic and Marius Hofert
.Fn <- function(x, X, smoothing = c("none", "beta", "checkerboard"),
                offset = 0, method = c("C", "R"))
## {
##     if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
##     if(!is.matrix(X)) X <- rbind(X, deparse.level=0L)
##     structure(class = "mvFn", x=x, X=X, ## <- so the result can be plotted/printed ..
##               .Fn(x, X, offset=offset, method=method))
## }
## .Fn <- function(x, X, offset=0, method=c("C", "R"))
{
    stopifnot(is.numeric(d <- ncol(X)), is.matrix(x), d == ncol(x))
    n <- nrow(X)
    if(d == 1)
	vapply(x, function(x.) sum(X <= x.), NA_real_) / (n + offset)
    else { # d > 1
	method <- match.arg(method)
        smoothing <- match.arg(smoothing)
        type <- which(smoothing == c("none", "beta", "checkerboard"))
	switch(method,
	       "C" = {
		   m <- nrow(x)
		   .C(Cn_C, # see ../src/empcop.c
		      as.double(X),
		      as.integer(n),
		      as.integer(d),
		      as.double(x),
		      as.integer(m),
		      ec=double(m),
		      as.double(offset),
                      as.integer(type))$ec
	       },
	       "R" = {
                   ## == apply(x, 1, function(x.) sum(colSums(t(X)<=x.)==d)/(n+offset) )
                   ## but vapply is slightly faster (says MH)
                   tX <- t(X)
                   switch(smoothing,
                          "none"={
                             vapply(1:nrow(x), function(k) sum(colSums(tX <= x[k,]) == d),
                                    NA_real_) / (n + offset)
                          },
                          "beta"={
                             stop("no R implementation for 'smoothing': ", smoothing)
                          },
                          "checkerboard"={
                             stop("no R implementation for 'smoothing': ", smoothing)
                          },
                          stop("wrong 'smoothing': ", smoothing)
                          )
               },
               stop("wrong 'method': ", method))
    }
}

##' @title Empirical distribution function
##' @param x (m, d) matrix of evaluation points
##' @param X (n, d) matrix of data based on which the empirical distribution
##'        function is computed
##' @param smoothing smoothing method
##' @param offset scaling factor of the empirical distribution function (= n/(n + offset))
##' @param method method string ("C" for C code; "R" for R code)
##' @return empirical distribution function of X at x
##' @author Ivan Kojadinovic and Marius Hofert
##' @note We should actually remove this and call the above function F.n()
##'       (thus also making 'smoothing' available to F.n())
F.n <- function(x, X, offset = 0, method = c("C", "R"))
    .Fn(x = x, X = X, smoothing = "none", offset = offset, method = method)


### Estimated partial derivatives of a copula ##################################

##' @title Estimated Partial Derivatives of a Copula Given the Empirical Copula
##' @param u (m, d)-matrix of evaluation points
##' @param U (n, d)-matrix of pseudo-observations
##' @param j.ind dimensions for which the partial derivatives should be estimated
##' @param b bandwidth in (0, 1/2) for the approximation
##' @param ... additional arguments passed to .Fn()
##' @return (m, length(j.ind))-matrix containing the estimated partial
##'         derivatives with index j.ind of the empirical copula of U at u
##' @author Marius Hofert
dCn <- function(u, U, j.ind = 1:d, b = 1/sqrt(nrow(U)), ...)
{
    ## Check
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    if(!is.matrix(U)) U <- rbind(U, deparse.level=0L)
    stopifnot((d <- ncol(U)) == ncol(u),
              0 <= u, u <= 1, 0 <= U, U <= 1,
              1 <= j.ind, j.ind <= d, 0 < b, b < 0.5)

    ## Functions to change the entry in the jth column of u
    ## See Remillard, Scaillet (2009) "Testing for equality between two copulas"
    adj.u.up <- function(x){
        x. <- x + b
        x.[x > 1-b] <- 1
        x.
    }
    adj.u.low <- function(x){
        x. <- x - b
        x.[x < b] <- 0
        x.
    }

    ## (v)apply to each index...
    m <- nrow(u)
    res <- vapply(j.ind, function(j){
        ## Build the two evaluation matrices for Cn (matrix similar to u with
        ## column j replaced). Note, this is slightly inefficient for those u < b
        ## but at least we can "matricize" the problem.
        ## 1) compute .Fn(u.up)
        u.up <- u
        u.up[,j] <- adj.u.up(u.up[,j])
        Cn.up <- .Fn(u.up, U, ...)
        ## 2) compute .Fn(u.low)
        u.low <- u
        u.low[,j] <- adj.u.low(u.low[,j])
        Cn.low <- .Fn(u.low, U, ...)
        ## 3) compute difference quotient
        (Cn.up - Cn.low) / (2*b)
    }, numeric(m))
    res <- pmin(pmax(res, 0), 1) # adjust (only required for small sample sizes)
    if(length(j.ind)==1) as.vector(res) else res
}


### Empirical copula ###########################################################

##' @title Empirical CDF and hence copula of X at x
##' @param x (m, d) matrix of evaluation points
##' @param X (n, d) matrix of pseudo-data based on which the empirical copula
##'        is computed
##' @param smoothing usual, beta or checkerboard empirical copula
##' @param offset scaling factor of the empirical distribution function (= n/(n + offset))
##' @param method method string ("C" for C code; "R" for R code)
##' @param ties.method passed to pobs()
##' @return empirical copula of X at x
##' @author Ivan Kojadinovic, Marius Hofert and Martin Maechler (C.n -> F.n; re-organisation)
##' @note See ../man/empcop.Rd for a nice graphical check with the Kendall function
C.n <- function(u, X, smoothing = c("none", "beta", "checkerboard"),
                offset = 0, method = c("C", "R"),
                ties.method = c("max", "average", "first", "last", "random", "min"))
{
    if(any(u < 0, 1 < u))
        stop("'u' must be in [0,1].")
    ties.method <- match.arg(ties.method)
    .Fn(u, pobs(X, ties.method = ties.method),
        offset = offset, method = method, smoothing = smoothing)
}


### Empirical copula class #####################################################

## Constructor
empCopula <- function(X, smoothing = c("none", "beta", "checkerboard"), offset = 0,
                      ties.method = c("max", "average", "first", "last", "random", "min"))
{
    stopifnot(is.matrix(X) || is.data.frame(X), 0 <= X, X <= 1, 0 <= offset, offset <= 1)
    smoothing <- match.arg(smoothing)
    ties.method <- match.arg(ties.method)
    ## Construct object
    new("empCopula",
	dimension = as.integer(ncol(X)),
	X = X, # 'parameters' so to say
        smoothing = smoothing, # a 'parameter' so to say
        offset = offset,
        ties.method = ties.method,
        ## require to pass validity as a 'copula' object ('exprdist' not required)
        parameters = double(0),
        param.names = character(0),
        param.lowbnd = double(0),
        param.upbnd = double(0),
	fullname = "<deprecated slot>")
}


### Methods ####################################################################

## describe method
setMethod(describeCop, c("empCopula", "character"), function(x, kind, prefix="", ...) {
    kind <- match.arg(kind)
    if(kind == "very short") # e.g. for show() which has more parts
        return(paste0(prefix, "Empirical copula"))
    ## else
    d <- dim(x)
    ch <- paste0(prefix, "Empirical copula, dim. d = ", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = ch,
           stop("invalid 'kind': ", kind))
})

## dCopula method
setMethod("dCopula", signature("matrix", "empCopula"),
	  function(u, copula, log = FALSE, ...) {
    if(copula@smoothing == "beta") {
        R <- apply(copula@X, 2, rank, ties.method = copula@ties.method) # (n, d) matrix of ranks
        n <- nrow(copula@X)
        if(!is.matrix(u)) u <- rbind(u) # (m, d) matrix of evaluation points
        m <- nrow(u)
        d <- copula@dimension
        f <- vapply(1:m, FUN = function(k) dbeta(u[k,], shape1 = R, shape2 = n+1-R, log = log, ...),
                    FUN.VALUE = matrix(NA_real_, nrow = n, ncol = d)) # (n, d, m) array containing all f_{n,R_{ij}}(u_{kj}) (see copula book)
        offset <- copula@offset
        apply(f, 3, function(f.) {
            if(log) {
                lx <- rowSums(f.)
                lsum(lx) - log(n + offset)
            } else {
                sum(apply(f., 1, prod)) / (n + offset)
            }
        })
    } else stop("Empirical copula only has a density for smoothing = 'beta'")
})

## rCopula method
setMethod("rCopula", signature("numeric", "empCopula"),
	  function(n, copula) copula@X[sample(1:nrow(copula@X), size = n), ])

## pCopula method
setMethod("pCopula", signature("matrix", "empCopula"),
	  function(u, copula, log.p = FALSE, ...) {
    res <- C.n(u, X = copula@X, smoothing = copula@smoothing,
               offset = copula@offset, ties.method = copula@ties.method, ...)
    if(log.p) log(res) else res
})

## Measures of association (unclear what their values are for smoothing = "beta" or "checkerboard"
## setMethod("tau", signature("empCopula"), function(copula) ) # unclear
## setMethod("rho", signature("empCopula"), function(copula) )
## setMethod("lambda", signature("empCopula"), function(copula) )