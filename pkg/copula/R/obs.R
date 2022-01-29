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


### Auxiliary transformations for copula observations ##########################

##' @title Pseudo-observations
##' @param x matrix of random variates to be converted to pseudo-observations
##' @param na.last passed to rank()
##' @param ties.method passed to rank()
##' @param lower.tail if FALSE, pseudo-observations when apply the empirical
##'        marginal survival functions are returned.
##' @return pseudo-observations (of the same dimensions as x)
##' @author Marius Hofert & Martin Maechler
pobs <- function(x, na.last = "keep",
		 ## formals(rank) works in pre-2015-10-15 and newer version of rank():
		 ties.method = eval(formals(rank)$ties.method),
		 lower.tail = TRUE)
{
    ties.method <- match.arg(ties.method)
    U <- if(!is.null(dim(x)))
	     apply(x, 2, rank, na.last=na.last, ties.method=ties.method) / (nrow(x)+1)
	 else
	     rank(x, na.last=na.last, ties.method=ties.method) / (length(x)+1)
    if(inherits(x, "zoo")) # incl "xts" (but no similar) -- FIXME? and use:
### if(is.object(x) && !isS4(x) && !is.data.frame(x)) ## "zoo", "xts" et al
	attributes(U) <- attributes(x)
    if(lower.tail) U else 1-U
}

##' @title Transform to (More) Uniform Margins
##' @param x d-column matrix containing the d-dimensional sample whose ranks are
##'        used to index standard uniforms according to different methods
##' @param n number of samples to be constructed.
##' @param method character string indicating the method to be used to
##'        sample. Available are:
##'        "beta": rank(x) is used to index a matrix of columnwise sorted
##'                standard uniforms of the same size as x. Then
##'                rows of this matrix are resampled with replacement.
##'                This is samples from the empirical beta copula of x.
##'        "rank.indexed.unif": Generate an (m,d)-matrix where m is at
##'                least nrow(x) and a multiple of nrow(x). One then
##'                generates m uniforms in each dimension, iterates over
##'                the m/nrow(x) many blocks of uniforms, sorts the block,
##'                of uniforms, indexes them by rank(x), rbind()s
##'                all such blocks of resulting data together, and
##'                returns the first n rows of data.
##' @param ties.method passed to the underlying rank().
##' @return (n, d)-matrix of marginally transformed random variates.
##' @author Marius Hofert
##' @note ruobs() is used to smooth empirical copula samples (e.g. for
##'       sampling from the empirical beta copula (but there ties will
##'       be produced the latest when n > nrow(x), due to resampling)
ruobs <- function(x, n = nrow(x), method = c("beta", "rank.indexed.unif"),
                  ties.method = "random")
{
    ## Checks
    if(!is.matrix(x)) x <- rbind(x)
    N <- nrow(x)
    d <- ncol(x)
    stopifnot(n >= 1)

    ## Compute componentwise ranks of x
    R <- apply(x, 2, rank, ties.method = match.arg(ties.method))

    ## Method switch
    method <- match.arg(method)
    switch(method,
           "beta" = {
               ## For empirical beta copula, as suggested by
               ## Segers, Sibuya, Tsukahara (2016, "The Empirical Beta Copula")

               ## Generate uniforms
               U <- matrix(runif(N * d), ncol = d) # only N * d

               ## Create dependence by componentwise sorting the uniforms
               ## and indexing them by the ranks of x
               U.sort <- apply(U, 2, sort) # note: no need to keep order as we resample anyways
               U.sort.ranked <- sapply(1:d, function(j) U.sort[R[,j],j]) # sorted U's indexed by the ranks of x

               ## Sampling (can generate ties of course)
               U.sort.ranked[sample.int(N, size = n, replace = TRUE),] # resampling from U.sort.ranked
           },
           "rank.indexed.unif" = {
               ## Careful. For n = N, Segers, Sibuya, Tsukahara (2016,
               ## "The Empirical Beta Copula", under (1.3)) says that this
               ## method is "not quite correct"; maybe that only applies
               ## to being a sample from a beta copula (?).

               ## Preliminaries
               nblk <- ceiling(n/N) # number of N-blocks so that m = nblk * N >= n and m is a multiple of N
               m <- nblk * N # need at least N and a multiple of N (so that indexing with R works)

               ## Generate uniforms
               ## Note: These could also be a Sobol sequence of length m
               ##       library(qrng)
               ##       U <- sobol(m, d = d, randomize = "digital.shift") # => generates clusters around each x[i,]
               ##       U <- matrix(replicate(d, sobol(m, randomize = "digital.shift")), ncol = d) # => same
               U <- matrix(runif(m * d), ncol = d) # number of rows now a multiple of N

               ## Create dependence by componentwise sorting the uniforms and
               ## indexing blocks of size N of them by the ranks of x
               ## Note: Using the following generates "small copula samples" concatenated
               ##       along the diagonal (=> comonotone scaled copula samples).
               ##       U.sort <- apply(U, 2, sort)
               ##       rex <- matrix(rep((0:(nblk-1))*N, each = N), nrow = m, ncol = d) # rank expansion
               ##       R.m <- rex + do.call(rbind, rep(list(R), nblk)) # rbind(0*N + R, 1*N + R, 2*N + R, ..., (nblk-1)*N + R)
               ##       U.sort.ranked <- sapply(1:d, function(j) U.sort[R.m[,j],j]) # sorted U's indexed by the *expanded* ranks of x
               blk.data <- lapply(1:nblk, function(k) {
                   ii <- N * (k-1) + 1:N # block indices
                   U.blk <- U[ii,]
                   ## Now the same as method = "beta" but for this block only
                   U.blk.sort <- apply(U.blk, 2, sort)
                   sapply(1:d, function(j) U.blk.sort[R[,j],j]) # U.blk.sort.ranked
               })
               U.blk.ranked <- do.call(rbind, blk.data)

               ## Sampling
               ## Shuffling not necessary
               ## U.shuffled <- U.blk.ranked [sample.int(m),] # shuffle the rows of the now dependent U's
               ## U.shuffled[1:n,] # grab out the first n
               U.blk.ranked[1:n,] # grab out the first n
           },
           stop("Wrong 'method'"))
}
