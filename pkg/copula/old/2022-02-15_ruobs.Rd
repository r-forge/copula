\name{ruobs}
\alias{ruobs}
\title{Transform To Uniform Margins}
\description{
  Compute certain (randomized, marginal) transforms to uniform margins.
}
\usage{
ruobs(x, n = nrow(x), method = c("beta", "rank.indexed.unif"),
      ties.method = "random")
}
\arguments{
  \item{x}{\eqn{d}-column \code{\link{matrix}()} containing the
    \eqn{d}-dimensional sample whose ranks are used to index standard
    uniform random variates.}
  \item{n}{number of samples to be constructed.}
  \item{method}{\code{\link{character}} string indicating the method
    to be used to index sample standard uniforms. Available are:
    \describe{
      \item{\code{"beta":}}{\code{rank(x, ...)} is used to index a
	matrix of columnwise sorted standard uniforms of the same size
	as \code{x}. Then rows of this matrix are resampled with
	replacement. This method is used for sampling from the empirical
	beta copula of \code{x}.}
      \item{\code{"rank.index.unif":}}{Generate an
       (\eqn{m},\eqn{d})-matrix where \eqn{m} is at
       least \code{nrow(x)} and a multiple of \code{nrow(x)}.
       Then \code{m} uniforms are generated in each dimension,
       one iterates over the \eqn{m}/\code{nrow(x)}-many blocks of
       uniforms, sorts the current block of uniforms, indexes it by
       \code{rank(x)}, \code{rbind()}s all such blocks of resulting
       data together, and returns the first \code{n} rows of.}
    }
  }
  \item{ties.method}{\code{\link{character}} string specifying how ranks
    should be computed if there are ties in any of the coordinate
    samples of \code{x}; passed to the underlying \code{\link{rank}}.}
}
\value{\eqn{n\times d}{(n x d)}\code{\link{matrix}} containing a
  sample constructed according to the given \code{method}.}
\details{
  This is a specialized function which should only be used by experts. The
  function is used, for example, to sample from empirical beta copulas.
}
%\author{Marius Hofert}
\seealso{
  \code{\link{pobs}()} for pseudo-observations.
}
\examples{
## Generate some data
set.seed(271)
ccop <- claytonCopula(iTau(0.5, copula = claytonCopula()))
N <- 100
x <- rCopula(N, copula = ccop)
plot(x)

## Smooth with sampel size n = N (the original sample size of x)
plot(ruobs(x))
plot(ruobs(x, method = "rank.indexed.unif"))

## For n > N
n <- 813
plot(ruobs(x, n = n)) # resampling
plot(ruobs(x, n = n, method = "rank.indexed.unif")) # indexing new uniforms
## The upper tail is not captured correctly (but then also neither
## in the original data x nor with method = "beta")
}