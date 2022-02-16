library(copula)

## Why using given ranks *blockwise* is *not* a good idea
rEmpCopula.via.block.indexed.runif <- function(n, copula)
{
    ## Preliminaries
    x <- copula@X
    dm <- dim(x)
    n.x <- dm[1] # size of original sample
    d <- dm[2]

    ## Ranks of x
    R <- apply(x, 2, rank, na.last = "keep", ties.method = "average") # (n.x, d)-matrix

    ## Block variables
    ## Example: n.x = 10, n = 17 => nblk = 2 many blocks, each of size n.x
    ##          (last one not fully used) => m = 20 rows can be filled with
    ##          nblk-many blocks of size n.x each.
    nblk <- ceiling(n/n.x) # smallest number of n.x-blocks so that m := nblk * n.x >= n
    m <- nblk * n.x # need m >= n.x and a multiple of n.x (so that indexing with ranks works)

    ## Generate m >= n uniforms
    ## Important: These *cannot* be sorted, but rather need to be sorted per block!
    U <- matrix(runif(m * d), ncol = d) # number of rows now a multiple of n.x

    ## Iterate over all blocks (no point in trying to do something more
    ## efficient, turns out to be less efficient *and* less readable)
    U.rank.lst <- lapply(1:nblk, function(k) {
        ii <- (k-1) * n.x + 1:n.x # indices of kth block
        U.k <- U[ii,] # (n.x, d)-matrix of U's from the kth block
        U.k.sort <- apply(U.k, 2, sort) # note again: sorting needs to be done here and not outside the blocks
        vapply(1:d, function(j) U.k.sort[R[,j],j], numeric(n.x)) # rank-index
    }) # nblk-list of (n.x * d)-vectors of rank-indexed sorted U's
    U.rank.mat <- do.call(rbind, U.rank.lst) # (m,d)-matrix (m = nblk * n.x) consisting of concatenated blocks
    ## Note:
    ## - matrix(unlist(U.rank.lst), ncol = d) would be *wrong*
    ## - Interestingly, the following idea did *not* work:
    ##   + shift ranks via
    ##     R.vec <- as.vector(R + rep((0:(d-1)) * n.x, each = n.x)) # (n.x * d)-vector of shifted ranks (we shift the jth column by (j-1) * n.x to be able to use R.vec to index flattened (n.x, d)-matrices)
    ##   + index sorted U's from block k by R.vec via
    ##     U.k.sort[R.vec] # (n.x * d)-vector of rank-indexed sorted U.k

    ## Grab out n <= m rows of U.rank.mat
    U.rank.mat[sample(1:m, size = n),] # no 'replace = TRUE' here
}

## Check the idea of blockwise indexing runifs for smoothing
n.x <- 2000
set.seed(271)
X. <- rCopula(n.x, copula = claytonCopula(2))

## n <= n.x => fine!
ec <- empCopula(X., smoothing = "block.indexed.runif")
U <- rCopula(n.x, copula = ec)
plot(U, pch = ".") # => clusters!

## n > n.x => clusters appear!
U. <- rCopula(100 * n.x, copula = ec)
plot(U., pch = ".") # bad!
