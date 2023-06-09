## Copyright (C) 2016 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

### keep using numeric for parameters but add an attribute for fixed
### ===> Entirely back compatible with previous copula package 'copula' objects

### Fix some parameters can be called before calling fitCopula #################

##' @title Fix a Subset of a Parameter Vector --> ../man/fixedPar.Rd
fixParam <- function(param, fixed = TRUE) {
    stopifnot(length(fixed) %in% c(1L, length(param)), is.logical(fixed))
    if(identical(fixed, FALSE)) attr(param, "fixed") <- NULL
    stopifnot(isTRUE(fixed) || length(param) == length(fixed), is.logical(fixed))
    if(!any(fixed)) attr(param, "fixed") <- NULL
    attr(param, "fixed") <- fixed
    param
}

##' @title Whether Each Component of a Parameter is free
##' @param param Numeric parameter vector, possibly with "fixed" attribute
##' @return A vector of logicals, TRUE = free
##' @author Jun Yan
isFree <- function(param)
    if (is.null(fixed <- attr(param, "fixed"))) rep(TRUE, length(param)) else !fixed

##' @title Return the attribute "fixed" even if it does not explicitly exist
##' @param param Numeric parameter vector, possibly with "fixed" attribute
##' @return A vector of logicals, FALSE = fixed
##' @author Ivan Kojadinovic
isFixedP <- function(param)
    if (is.null(fixed <- attr(param, "fixed"))) rep(FALSE, length(param)) else fixed


##' @title Whether or not the copula has "fixed" attr in parameters
##' @param copula A 'copula' object
##' @return TRUE if has, otherwise FALSE
##' @author Jun Yan
## unused
## hasFixedPar <- function(copula) !is.null(attr(copula, "fixed"))
### This to be used in place when npar is needed ###############################

##' @title Number of Free Parameters of a Vector
##' @param x A numeric parameter vector with possible attribute "fixed"
##' @return Length of free parameter
##' @author Jun Yan
nFree <- function(param) {
    length(param[isFree(param)]) # GETR
}

### get and set free/fixed parameters, needed in fitCopula #####################

##' @title Get the Free or Fixed Parameters of a Copula
##' @param copula 'copula' object
##' @param freeOnly logical: TRUE = free; FALSE = fixed
##' @param attr logical indicating if lower and upper bound attributes should be returned.
##' @param named logical indicating if parameter \code{\link{names}} should be returned.
##' @return A numeric vector of parameters with attributes
##'         param.names, param.lowbnd, and param.upbnd.
##' @author Jun Yan and Martin Maechler
setMethod("getParam", signature("copula", "logicalOrMissing", "logicalOrMissing",
                                "logicalOrMissing"),
          function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
    par <- copula@parameters
    if (length(par) == 0) return(par) ## no parameters (e.g., indepCopula)
    fixed <- attr(par, "fixed")
    sel <- if (!is.null(fixed) && freeOnly) !fixed else TRUE
    ## no selected parameter
    if (!any(sel)) ## = all(!sel), but faster
	numeric(0)
    else {
	par <- if(named) setNames(par[sel], copula@param.names[sel])
	       else par[sel]
	if(!attr)
	    par
	else ## 'attr' i.e., with the param.* as attributes :
	    structure(par,
		      param.lowbnd = copula@param.lowbnd[sel],
		      param.upbnd  = copula@param.upbnd[sel])
    }
})

## "FIXME": have the exported *generic* function  setTheta() (for years)
##  -----   whereas this is currently hidden :
##' @title Set or Modify Values of the Free Parameters of a Copula
##' @param copula a copula object
##' @param value a numeric vector to be set for the parameters
##' @return A copula object with parameters set to be param
##' @author Jun Yan
setMethod("freeParam<-", signature("copula", "numeric"), function(copula, value) { # GETR
    oldpar <- copula@parameters
    fixed <- attr(oldpar, "fixed")
    if (is.null(fixed) || !any(fixed)) {
        ## stopifnot(length(oldpar) == length(value))
        if (length(oldpar) != length(value)) # IK
            stop("the length of 'value' is not equal to the number of free parameters")
        copula@parameters[] <- value
    }
    else {
        sel <- !fixed
        ## stopifnot(sum(sel) == length(value))
        if (sum(sel) != length(value)) # IK
            stop("the length of 'value' is not equal to the number of free parameters")
        copula@parameters[sel] <- value
    }
    ## special operation for copulas with df parameters
    ## if (has.par.df(copula))
    ##     copula@df <- copula@parameters[length(copula@parameters)]
    ## if (validObject(copula)) copula
    ## else stop("Invalid copula object.")
    copula
})

##' @title Set or Modify "Fixedness" of Copula Parameters
##' @param copula a copula object
##' @param value a logical vector to be set for the fixed attribute
##' @return A copula object with parameters set to be param
##' @author Jun Yan
setMethod("fixedParam<-", signature("copula", "logical"), function(copula, value) { # GETR
    stopifnot(length(value) %in% c(1L, length(copula@parameters)))
    if (anyNA(copula@parameters[value])) stop("Fixed parameters cannot be NA.")
    attr(copula@parameters, "fixed") <-
	if(identical(value, FALSE) || !any(value)) NULL else value
    copula
})


## BEGIN GETR
## logical indicating which parameters are free
setMethod("free", signature("copula"), function(copula) isFree(copula@parameters))

## logical indicating which parameters are fixed
setMethod("fixed", signature("copula"), function(copula) isFixedP(copula@parameters))

## number of free parameters
setMethod("nFreeParam", signature("copula"), function(copula) nFree(copula@parameters))

## number of parameters
setMethod("nParam", signature("copula"), function(copula) length(copula@parameters))
## END GETR
