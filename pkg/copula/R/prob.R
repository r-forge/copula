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


### Computing probabilities of falling in hyperrectangles

##' @title Computing volumes of a function
##' @param FUN a function to compute the volume of
##' @param lower vector of lower endpoints of the volume
##' @param upper vector of upper endpoints of the volume
##' @param ... additional arguments passed to FUN
##' @return volume
##' @author Marius Hofert, Martin Maechler
##' @note MWE:
##' library(copula)
##' d <- 3
##' l <- rep(0.2, d)
##' u <- rep(0.4, d)
##' cc <- claytonCopula(2, dim = d)
##' prob(cc, l = l, u = u)
##' f <- function(x) pCopula(x, copula = cc)
##' volume(f, lower = l, upper = u)
volume <- function(FUN, lower, upper, ...) {
    d <- length(upper)
    stopifnot(is.numeric(lower), is.numeric(upper),
              length(lower) == d) # could even be called with lower > upper
    D <- 2^d
    m <- 0:(D - 1)
    ## digitsBase() from package 'sfsmisc' {slightly simplified} :
    ## Purpose: Use binary representation of 0:N
    ## Author: Martin Maechler, Date: Wed Dec 4 14:10:27 1991
    II <- matrix(0, nrow = D, ncol = d)
    for(i in d:1L) {
        II[,i] <- m %% 2L + 1L
        if (i > 1) m <- m %/% 2L
    }
    Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2] # signs for checkerboard system; the ("upper",...,"upper") case has +1; = c(2,2,...,2)
    x <- array(cbind(lower, upper)[cbind(c(col(II)), c(II))], dim = dim(II)) # evaluation points
    ## Computing the volume
    sum(Sign * FUN(x, ...))
}

##' @title Compute the probability P[l < U <= u] where U ~ copula x
##' @param x copula object
##' @param l d-vector of lower "integration" limits
##' @param u d-vector of upper "integration" limits
##' @return the probability that a random vector following the given copula
##'         falls in the hypercube with lower and upper corner l and u, respectively.
##' @author Marius Hofert, Martin Maechler
setGeneric("prob", function(x, l, u) standardGeneric("prob"))

setMethod("prob", signature(x = "Copula"),
          function(x, l, u) {
              d <- dim(x)
              stopifnot(is.numeric(l), is.numeric(u),
                        length(u) == d, d == length(l),
                        0 <= l, l <= u, u <= 1)
              ## if(d > 30)
	      ##     stop("prob() for copula dimensions > 30 are not supported (yet)")
              volume(function(z) pCopula(z, copula = x), lower = l, upper = u)
          })
