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

### pdf = dCopula()  ---  but then also pCopula()
### ===   ~~~~~~~~~                     ~~~~~~~~~

require(copula)
## library(fCopulae)

source(system.file("Rsource", "utils.R", package="copula", mustWork=TRUE))
## All non-virtual copula classes:
source(system.file("Rsource", "cops.R", package="copula", mustWork=TRUE))
## --> copcl, copObs, copBnds,  excl.2 , copO.2, copBnd.2

### preparation for a grid
n1 <- 17
n2 <- 21
eps <- .001 ## <- not going to very extremes
u1 <- seq(0, 1, length=n1); u1[1] <- eps; u1[n1] <- 1-eps
u2 <- seq(0, 1, length=n2); u2[1] <- eps; u2[n2] <- 1-eps

### d=2 ########################################################################

Exp.grid <- function(...)
    as.matrix(expand.grid(..., KEEP.OUT.ATTRS = FALSE))

umat <- Exp.grid(u1=u1, u2=u2)

## all copulas give the same tau except galambos, amh
tau <- 0.5

pdf("densCop_2d.pdf")

fCols <- colorRampPalette(c("red", "white", "blue"), space = "Lab")


## frankCopula
theta.fr <- iTau(frankCopula(), tau)
dcop <- matrix(dCopula(umat, frankCopula(param=theta.fr, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( frankCopula(%.4g) )", theta.fr))
round(dcop, 3)


## claytonCopula
(theta.cl <- iTau(claytonCopula(), tau))
stopifnot(all.equal(theta.cl, copClayton@iTau(tau), tolerance = 1e-13))
dcop <- matrix(dCopula(umat, claytonCopula(param=theta.cl, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf(" Density( claytonCopula(%.4g) )", theta.cl))
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( claytonCopula(%.4g) )", theta.cl))
round(dcop, 3)



## gumbelCopula
theta.gu <- iTau(gumbelCopula(), tau)
stopifnot(all.equal(theta.gu, copGumbel@iTau(tau), tolerance = 1e-13))
dcop <- matrix(dCopula(umat, gumbelCopula(param=theta.gu, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( gumbelCopula(%.4g) )", theta.gu))
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( gumbelCopula(%.4g) )", theta.gu))
round(dcop, 3)


## normalCopula
uB <- cbind(c(0:1, .5), (1:9)/10) # values at boundaries
fC <- dCopula(uB, normalCopula(0.55))
stopifnot(is.finite(fC), length(fC)==nrow(uB), fC[-3*(1:3)] == 0)
theta.n <- iTau(normalCopula(), tau)
dcop <- matrix(dCopula(umat, normalCopula(param=theta.n, dim = 2)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( normalCopula(%.4g) )", theta.n))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( normalCopula(%.4g) )", theta.n))
round(dcop, 3)


## tCopula
fC <- dCopula(uB, tCopula(0.55))
stopifnot(is.finite(fC), length(fC)==nrow(uB), fC[-3*(1:3)] == 0)
(theta.t. <- iTau(tCopula(df=10), tau))
dcop <- matrix(dCopula(umat, tCopula(param=theta.t., dim = 2, df=10)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( tCopula(%.4g, df = 10) )", theta.t.))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( tCopula(%.4g, df = 10) )", theta.t.))
round(dcop, 3)
## tCopula -- df=4
(theta <- iTau(tCopula(df=4), tau))
dcop <- matrix(dCopula(umat, tCopula(param=theta, dim = 2, df=4)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( tCopula(%.4g, df = 4) )", theta))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( tCopula(%.4g, df = 4) )", theta))
round(dcop, 3)

## galambosCopula
#
(theta <- iTau(galambosCopula(), tau))
stopifnot(all.equal(tau, tau(galambosCopula(theta)), tolerance = 1e-5))
dcop <- matrix(dCopula(umat, galambosCopula(param=theta)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( galambosCopula(%.4g) )", theta))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( galambosCopula(%.4g) )", theta))
round(dcop, 3)


## amhCopula
tau <- 0.3 ## -- to be in range for AMH
(theta <- iTau(amhCopula(), tau))
stopifnot(all.equal(tau, tau(amhCopula(theta)), tolerance = 1e-5))
dcop <- matrix(dCopula(umat, amhCopula(param=theta)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( amhCopula(%.4g) )", theta))
round(dcop, 3)

showProc.time()
dev.off()

### d > 2 ######################################################################

### for package fCopulae
## dcop.f <- darchmCopula(umat, alpha=2, type=1, alternative=T)
## round(matrix(dcop.f, n1, n2), 3)
## dcop.fCopulae <- darchmCopula(umat, alpha=theta.fr, type = 5, alternative = TRUE)
## round (matrix(dcop.fCopulae, n1, n2), 3)

## dim = 3, frankCopula
um3 <- Exp.grid(u1=u1, u2=u2, u3=u1)
dcE <- tryCatch(dcopula(um3), error=identity)
dcop <- dCopula(um3, frankCopula(param=theta.fr, dim = 3))
stopifnot(dcop >= 0, # no NA here
          inherits(dcE, "error"),
          grepl("defunct", conditionMessage(dcE), ignore.case=TRUE))
round(array(dcop, c(n1,n2,n1)), 3)

## dim = 4 --- fine as long we are "out of corners"
um4 <- Exp.grid(u1=u1, u2=u2, u3=u1, u4=rev(u2))
dcop <- dCopula(um4, frankCopula(param=theta.fr, dim = 4))
stopifnot(dcop >= 0,# no NA here
	  all.equal(dcop, copFrank@dacopula(um4, theta = theta.fr)))
## round(array(dcop, c(n1,n2,n1,n2)), 3)

showProc.time()


### --- now look at dCopula() and pCopula() for *all* copulas:

##' u "mostly" in [0,1] .. with exceptions that are even NA, NaN
mku <- function(n, fr = 1/20, sd = 1/10) {
    r <- runif(n) + rnorm(n, sd=sd)
    r[sample(n, ceiling(fr*n))] <- c(NA,NaN)
    r
}

## This is from ./moments.R --- keep in sync! ---
tau.s <- c(       -.1, 0, 0.05805, (1:2)/9, 0.3)
names(tau.s) <- paste0("tau=", sub("0[.]", ".", formatC(tau.s)))
tTau <- suppressWarnings(
    sapply(tau.s, function(tau)
               sapply(copObs, iTau, tau = tau)) )
tTau["joeCopula", "tau=-.1"] <- 1 # ugly hack

set.seed(12)
u <- matrix(mku(1000), ncol= 2)# d = 2 required in the copula objects


u.outside.01 <- function(u) apply(u, 1, function(x) any(x <= 0, 1 <= x))
u.out <- u.outside.01(u) # has NAs
u.ina <- apply(u, 1, function(x) any(is.na(x)))
u.OUT <- u.out & !u.ina  # no  NAs

for(cNam in names(copObs)) {
    cat(cNam, paste0(":\n", paste(rep("-",nchar(cNam)), collapse=""), "\n\n"))
    for(th in unique(tTau[cNam,])) {
        cat(sprintf(" -- theta: %8.5g\n", th))
        cop <- setPar(copObs[[cNam]], th)
try({### FIXME!
        fu <- dCopula(u, cop); fu.na <- fu[u.ina]
        Fu <- pCopula(u, cop); Fu.na <- Fu[u.ina]
        stopifnot(fu[u.OUT] == 0,
		  all.equal(fu, exp(dCopula(u, cop, log=TRUE)), tolerance = 1e-15),
                  0 <= Fu[!u.ina], Fu[!u.ina] <= 1,
                  is.na(fu.na) | (0 <= fu.na),
                  is.na(Fu.na) | (0 <= Fu.na & Fu.na <= 1))
})### FIXME
    }
}



showProc.time()
