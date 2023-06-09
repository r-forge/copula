\name{getAcop}
\alias{getAcop}
\title{Get "acopula" Family Object by Name}
\description{
  Get one of our "acopula" family objects (see \link{acopula-families} by name.
}
\usage{
getAcop(family, check=TRUE)
}
\arguments{
  \item{family}{either a \code{\link{character}} string, the short or longer form
    of the Archimedean family name (e.g., "Clayton" or simply "C"), see
    the \code{\link{acopula-families}} documentation,
    \cr
    or an \code{\linkS4class{acopula}} family object.}
  \item{check}{logical indicating if the class of the return value
    should be checked to be \code{"\linkS4class{acopula}"}.}
}
\value{an \code{"\linkS4class{acopula}"} family object, typically one of
  one of our predefined ones.
}
\author{Martin Maechler}
\seealso{
  Our predefined \code{\link{acopula-families}};
  the class definition \code{"\linkS4class{acopula}"}.
}
\examples{
getAcop("Gumbel")

stopifnot(identical(getAcop("J"), getAcop("Joe")),
          identical(getAcop("J"), copJoe))
}
\keyword{manip}
