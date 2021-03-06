## Copyright (C) 2010 Marius Hofert and Martin Maechler
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

### TODO Should we try to be compatible to  'Copula' package ??

## This has advantage that arithmetic with scalars works "for free" already:
setClass("interval", contains =  "numeric", # of length 2
         representation(open  = "logical"),# of length 2
	 validity = function(object) {
	     if(length(rng <- object@.Data) != 2) "interval must be of length 2"
	     else if(length(object@open) != 2) "'open' must be of length 2"
	     else if(rng[2] < rng[1]) "'range[2]' must not be smaller than range[1]"
	     else TRUE
	 })

setClassUnion("maybeInterval", c("interval", "NULL"))

### Mother class of all (simple, *NON*-nested) Archimedean Copula Types
### for *any* dimension d
setClass("acopula",
	 representation(name = "character",
                        psi = "function",    # of (t, theta) -- the generator
                        psiInv = "function", # of (p, theta) -- psi_inverse: \psi^{-1}(p) = t
                        theta = "numeric", # value of theta or  'NA'  (for unspecified)
                        paraConstr = "function", # of (theta) ; constr(theta) |--> TRUE: "fulfilled"
                        ## when theta is one-dimensional, specifying the interval is more convenient:
                        paraInterval = "maybeInterval", # [.,.]  (.,.], etc ..
                        V0 = "function",	# of (n,theta) -- RNGenerator
                        tau= "function",	# of (theta)
                        tauInv = "function",    # of (tau)
                        lambdaL = "function",    # of (theta) lower bound  \lambda_l
                        lambdaLInv = "function", # of (lambda) - Inverse of \lambda_l
                        lambdaU = "function",    # of (theta)  - upper bound  \lambda_u
                        lambdaUInv = "function", # of (lambda) - Inverse of \lambda_u

                        ## Nesting properties if the child copulas are of the same family :
                        nestConstr = "function", # of (th0, th1) ; TRUE <==> "fulfilled"
                        V01= "function"	# of (V0,theta0,theta1)
                        ),
         prototype = prototype(theta = NA_real_),
	 validity = function(object) {
	     if (length(nm <- object@name) != 1 || nchar(nm) < 1)
		 return("'name' must be a string (of at least 1 character)")
             th <- tryCatch(validTheta(object), error = function(e) e)
             if(is(th, "error")) return(th $ message)

	     checkFun <- function(sName, nArgs, chkVectorizing=TRUE) {
		 f <- slot(object, sName)
		 if (length(formals(f)) < nArgs)
		     paste("slot '",sName,"' must be a function of at least ",nArgs,
			   " arguments", sep="")
		 else if(chkVectorizing && nArgs <= 2) {
		     ## test that the function can be called with NULL
		     r0 <- tryCatch(if(nArgs == 2) f(NULL, theta = th) else f(NULL),
                                    error = function(e) e)
		     if(is(r0, "error"))
			 sprintf("calling %s(NULL%s fails: %s", sName,
				 if(nArgs == 2) sprintf(", %g)", th) else ")",
				 r0$message)
		     else if(!identical(r0,
                                        {n0 <- numeric(0)
                                         if(nArgs == 2) f(n0, theta = th) else f(n0) }))
			 sprintf("%s(NULL ..) is not the same as %s(num(0) ..)",
				 sName, sName)
		     else TRUE
		 } else TRUE
	     } ## {checkFun}

             if(!isTRUE(tt <- checkFun("psi", 2)))	return(tt)
             if(!isTRUE(tt <- checkFun("psiInv", 2)))	return(tt)
             if(!isTRUE(tt <- checkFun("tau", 1)))	return(tt)

             if(!isTRUE(tt <- checkFun("paraConstr", 1, chkVect=FALSE))) return(tt)
             if(!isTRUE(tt <- checkFun("nestConstr", 2, chkVect=FALSE))) return(tt)

             ## ...
             ## ... (TODO)

             ## Check more :
	     if (object@psi(0, theta= 1/2) != 1)
		 return("psi(0, theta=1/2) != 1 -- seems an invalid generator")

             ## ....
             ## ....
             ## ....
             ## ....

	     ## 'else'	ok :
	     TRUE
	 })

## Construct 'paraConstr' slot automatically from 'paraInterval' :
setMethod("initialize", signature(.Object = "acopula"),
	  function(.Object, paraInterval, ...) {
	      if(!missing(paraInterval) && is(paraInterval, "interval")) {
		  .Object@paraConstr <- mkParaConstr(paraInterval)
		  .Object@paraInterval <- paraInterval
              }
              callNextMethod()
	  })


## A utility, not yet exported _ TODO ? _

setGeneric("validTheta", function(x, u) standardGeneric("validTheta"))
setMethod("validTheta", signature(x = "acopula"),
	  function(x) {
	      if(is((int <- x@paraInterval), "interval")) {
		  if(is.finite(d <- diff(as.numeric(int)))) # take mid-value in interval
		      int[1] + d/2
		  else if (all(iinf <- is.infinite(is.finite(int)))) ##	 (-Inf, Inf)
		      1/2
		  else { ## at least one end is finite
		      if(int[2] == Inf) int[1] + 1
		      else if(int[1] == -Inf) int[2] - 1
		      else stop("invalid paraInterval slot??")
		  }
	      } else  { ## need to use parameter-contraint function and "random search":
		  f.c <- x@paraConstr
		  for (th in c(1/2, (0:16)/16, round(exp(rnorm(20)),2)))
		      if(f.th <- f.c(th)) break
		  if(f.th) th else stop("could not find validTheta(<c>) for this copula")
	      }
	  })


### Nested Archimedean Copulas with *specified* dimension(s)
setClass("nacopula",
	 representation(copula = "acopula",
                        comp = "integer", # from 1:d -- of length in [0,d]
                        childCops = "list" #of nacopulas, possibly empty
                        ## TODO? nesting properties (V01,..) for specific mother-child relations
                        ),
         validity = function(object) {
             if(length(d <- dim(object)) != 1 || !is.numeric(d) || d <= 0)
                 return("invalid dim(.)")
             ic <- object@comp
             if((lc <- length(ic)) > d)
                 return("more components than the dimension")
	     if(!all("nacopula" == sapply(object@childCops, class)))
                 return("All 'childCops' elements must be 'nacopula' objects")
## FIXME: we need to recursively apply a "comp" function

             allC <- allComp(object)
             if(length(allC) != d)
                 return("must have d coordinates (from 'comps' and child copulas)")
             ##
             TRUE
         })


## FIXME: Maybe define a "strict_dim_nac" class which has the extra checks
## -----  for the "nacopula"s that appear as *SUB*-copulas, we do NOT want to
## require that their { components } are == { 1:d }
### Nested Archimedean Copulas with *specified* dimension(s)
setClass("outer_nacopula", contains = "nacopula",
         validity = function(object) {
             ## *Extra* checks in addition to those of "nacopula" :
             d <- dim(object)
             ic <- object@comp
             allC <- allComp(object)
             if(length(allC) != d)
                 return("must have d coordinates (from 'comps' and child copulas)")
             if(!all(sort(allC) == 1:d))
                 return(paste("The implicit coordinates are not identical to 1:d; instead\n  ",
                              paste(allC, collapse=", ")))
             ##
             TRUE
         })


## The dim() method is nicely defined  *recursive*ly :
setMethod("dim", signature(x = "nacopula"),
	  function(x) length(x@comp) + sum(unlist(lapply(x@childCops, dim))))

## also needed in validity above -- defined recursively, too :
allComp <- function(x) {
    stopifnot(is(x, "nacopula"))
    c(x@comp, unlist(lapply(x@childCops, allComp)))
}
