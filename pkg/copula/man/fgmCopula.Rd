\name{fgmCopula}
\alias{fgmCopula}
\title{Construction of a fgmCopula Class Object}
\description{
  Constructs a multivariate multiparameter Farlie-Gumbel-Morgenstern
  copula class object with its corresponding parameters and dimension.
}
\usage{
fgmCopula(param, dim = 2)
}
\arguments{
  \item{param}{a numeric vector specifying the parameter values.}
  \item{dim}{the dimension of the copula.}
}
\value{
  A Farlie-Gumbel-Morgenstern copula object of class
  \code{"\linkS4class{fgmCopula}"}.
}
\note{
  Note that a \eqn{d}-dimensional \code{fgmCopula} must have
  \code{npar}\eqn{= 2^d - d - 1} parameters.
  The verification of the validity of the parameter values is of high
  complexity and may not work for high dimensional copulas.

  The random number generation needs to be properly tested, especially
  for dimensions higher than 2.
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Nelsen, R. B. (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{Copula}}, \code{\link{copula-class}},
  \code{\link{fitCopula}}.
}
\examples{
## length(param) = #{parameters}  for  d-dimensional  FGM copula:
d <- 2:10; rbind(d, npar = 2^d - d - 1)
## d       2    3    4    5    6    7    8    9   10
## npar    1    4   11   26   57  120  247  502 1013

## a bivariate example
fgm.cop <- fgmCopula(1)
x <- rCopula(1000, fgm.cop)
cor(x, method = "kendall")
tau(fgm.cop)
cor(x, method = "spearman")
rho(fgm.cop)
persp  (fgm.cop, dCopula)
contour(fgm.cop, dCopula)

## a trivariate example with wrong parameter values
try(
 fgm2.cop <- fgmCopula(c(1,1,1,1), dim = 3)
) # Error: "Bad vector of parameters"

## a trivariate example with satisfactory parameter values
fgm2.cop <- fgmCopula(c(.2,-.2,-.4,.6), dim = 3)
fgm2.cop
}
\keyword{distribution}
\keyword{multivariate}
