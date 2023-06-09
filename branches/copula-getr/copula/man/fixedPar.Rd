\name{fixParam}
\title{Fix a Subset of a Copula Parameter Vector}
\alias{fixParam}
\alias{fixedParam<-}
\description{
  It is sometimes useful to keep fixed some components of a copula
  parameter vector whereas the others are \dQuote{free} and will be
  estimated, e.g., by \code{\link{fitCopula}}.
}
\usage{
fixParam(param, fixed = TRUE)
fixedParam(copula) <- value
}
\arguments{
  \item{param}{numeric parameter vector}
  \item{fixed, value}{logical vector of the same length as \code{param}
    indicating for each component \code{param[j]} if it is (going to be)
    fixed or not.}
  \item{copula}{a \code{"\linkS4class{copula}"} object.}
}
\value{
  \code{fixParam(param)} returns a \code{\link{numeric}} vector with
  attribute \code{"fixed"}(a \code{\link{logical}}, either \code{TRUE}
  or vector of the same length as \code{param}) to indicate which
  components of \code{param} are to be held fixed or not.

  \code{fixedParam<-} returns a \code{"\linkS4class{copula}"} object
  with a partly fixed parameter (slot), i.e., corresponding to
  \code{fixParam()} above.
}
%\author{Jun Yan}
%% \details{
%% }
%% \references{
%% }
\seealso{
  \code{\link{fitCopula}} for fitting;
  t-copulas, \code{\link{tCopula}(*, df.fixed=TRUE)} now uses parameter
  fixing for \code{"df"}.
}
\examples{
nc1 <- normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE,FALSE)),
                    dispstr = "un")
nc1
nc13 <- nc12 <- nc1
fixedParam(nc12) <- c(TRUE, TRUE, FALSE) ; nc12
fixedParam(nc13) <- c(TRUE, FALSE, TRUE) ; nc13
set.seed(17); x <- rCopula(100, nc1)
summary(f.13 <- fitCopula(nc13, x, method = "ml"))
f.13@copula  # 'rho.2' is estimated; the others kept fixed

## Setting to 'FALSE' (removes the "fixed" parts completely):
nc0 <- nc13; fixedParam(nc0) <- FALSE
nc0
stopifnot(is.null(attributes(nc0@parameters)))
}
\keyword{manip}
\keyword{utilities}
