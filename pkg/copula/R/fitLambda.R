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


##' @title Estimators of the Matrix of Tail Dependence Coefficients
##' @param u (n, d)-data matrix of pseudo-observations
##' @param method method string, available are:
##'        "Schmid.Schmidt": see Jaworksi et al. (2009, p. 231) and Schmid and
##'                          Schmidt (2007, Equations (5) and (17); bivariate
##'                          conditional Spearman's rho).
##'        "Schmidt.Stadtmueller": see Schmidt and Stadtmueller (2006, Equations
##'                                (3), (12), (16, left)); empirical lower tail copula
##'                                (outperformed EVT-based estimator in simulation
##'                                study)
##'        "t": Fitting a t copula to the whole data set and returning the implied
##'             tail dependence coefficient
##' @param p (small) probability (in (0,1)) used as 'cut-off' point
##' @param lower.tail logical indicating whether the lower or upper tail dependence
##'        coefficient is to be computed. If FALSE, the lower tail dependence
##'        coefficient of the 'flipped' data 1-u is computed.
##' @param verbose logical indicating whether a progress bar is displayed
##' @param ... additional arguments, for method = "t" passed to the
##'        underlying optimize()
##' @return estimate of the matrix of tail dependence coefficients
##' @author Marius Hofert
##' @note Formerly used heuristic for 'p': min(100/nrow(u), 0.1)
fitLambda <- function(u, method = c("Schmid.Schmidt", "Schmidt.Stadtmueller", "t"),
                      p = 1/sqrt(nrow(u)), # p = k/n with k -> Inf and k/n -> 0 as n -> Inf; popular choice (see also FRAPO::tdc): k = sqrt(n)
                      lower.tail = TRUE, verbose = FALSE, ...)
{
    ## Checks
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    d <- ncol(u)
    stopifnot(0 <= u, u <= 1, 0 <= p, p <= 1, length(p) == 1,
              d >= 2, is.logical(lower.tail), is.logical(verbose))

    ## Preliminaries for all methods
    method <- match.arg(method)
    u. <- if(lower.tail) u else 1-u # for upper tail dependence compute the lower tail dependence for the survival copula
    if(verbose) { # setup progress bar
        l <- 0 # counter
        pb <- txtProgressBar(max = choose(d, 2), style = if(isatty(stdout())) 3 else 1) # setup progress bar
        on.exit(close(pb)) # on exit, close progress bar
    }

    ## Compute Lambda
    switch(method,
    "Schmid.Schmidt" = {

        ## Compute lambda estimate for each pair
        Lam <- diag(1, nrow = d)
        pmu <- matrix(pmax(0, p-u.), ncol = d)
        for(i in 2:d) {
            for(j in 1:(i-1)) {
                Lam[i,j] <- mean(apply(pmu[, c(i, j)], 1, prod)) # \hat{Lambda}_{ij}
                if(verbose) {
                    l <- l + 1
                    setTxtProgressBar(pb, l) # update progress bar
                }
            }
        }
        int.over.Pi <- (p^2 / 2)^2
        int.over.M <- p^3 / 3
        Lam <- (Lam - int.over.Pi) / (int.over.M - int.over.Pi) # proper scaling

        ## Sanity adjustment
        Lam[Lam < 0] <- 0
        Lam[Lam > 1] <- 1

        ## Symmetrize and return
        ii <- upper.tri(Lam)
        Lam[ii] <- t(Lam)[ii]
        Lam

    },
    "Schmidt.Stadtmueller" = {

        ## Compute lambda estimate (n/k) * C_n(k/n, k/n) = sqrt(n) * C_n(1/sqrt(n), 1/sqrt(n)) = C_n(p,p)/p for each pair
        Lam <- diag(1, nrow = d)
        u.p <- u. <= p # p = k/n = sqrt(n)/n = 1/sqrt(n)
        for(i in 2:d) {
            for(j in 1:(i-1)) {
                Lam[i,j] <- mean(rowSums(u.p[, c(i, j)]) == 2) / p # \hat{Lambda}_{ij}
                if(verbose) {
                    l <- l + 1
                    setTxtProgressBar(pb, l) # update progress bar
                }
            }
        }

        ## Sanity adjustment
        Lam[Lam < 0] <- 0
        Lam[Lam > 1] <- 1

        ## Symmetrize and return
        ii <- upper.tri(Lam)
        Lam[ii] <- t(Lam)[ii]
        Lam

    },
    "t" = {

        ## -Log-likelihood for a t copula
        nLLt <- function(nu, P, u) {
            x <- qt(u, df = nu)
            -sum(dmvt(x, sigma = P, df = nu, log = TRUE) - rowSums(dt(x, df = nu, log = TRUE)))
        }

        ## Containers
        Lam <- diag(1, nrow = d) # matrix of pairwise estimated tail-dependence coefficients
        P   <- diag(1, nrow = d) # matrix of pairwise estimated correlation coefficients
        Nu  <- matrix(, nrow = d, ncol = d) # matrix of pairwise estimated degrees of freedom (NA on diagonal)

        ## Compute lambda for each pair (manually here, for performance reasons)
        ## Note: We use the approach of Mashal, Zeevi (2002) here
        for(i in 2:d) {
            for(j in 1:(i-1)) { # go over lower triangular matrix
                u.. <- u.[, c(i, j)]
                P. <- sin(cor.fk(u..)*pi/2) # 2 x 2 matrix
                rho <- P.[2,1]
                P[i,j] <- rho
                nu <- optimize(nLLt, interval = c(1e-4, 1e3), P = P., u = u.., ...)$minimum
                Nu[i,j] <- nu
                Lam[i,j] <- 2 * pt(-sqrt((nu + 1) * (1 - rho) / (1 + rho)), df = nu + 1)
                if(verbose) {
                    l <- l + 1
                    setTxtProgressBar(pb, l) # update progress bar
                }
            }
        }
        ## Symmetrize and return
        ii <- upper.tri(Lam)
        Lam[ii] <- t(Lam)[ii]
        P[ii]   <- t(P)[ii]
        Nu[ii]  <- t(Nu)[ii]
        list(Lambda = Lam, P = P, Nu = Nu)

    },
    stop("Wrong 'method'"))
}

