\name{cloud2-methods}
\alias{cloud2}
\alias{cloud2-methods}
\alias{cloud2,matrix-method}
\alias{cloud2,data.frame-method}
\alias{cloud2,copula-method}
\alias{cloud2,mvdc-method}
\docType{methods}
\title{Methods `cloud2' in Package `copula'}
\description{
  Methods \code{cloud2()} to draw cloud plots of
  three-dimensional distributions from package \pkg{copula}.
}
\usage{
\S4method{cloud2}{matrix}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{cloud2}{data.frame}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{cloud2}{copula}(x, n,
      xlim = 0:1, ylim = 0:1, zlim = 0:1,
      xlab = expression(U[1]), ylab = expression(U[2]),
      zlab = expression(U[3]), \dots)
\S4method{cloud2}{mvdc}(x, n,
      xlim = NULL, ylim = NULL, zlim = NULL,
      xlab = expression(X[1]), ylab = expression(X[2]),
      zlab = expression(X[3]), \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{xlim, ylim, zlim}{The x-, y- and z-axis limits.}
  \item{xlab, ylab, zlab}{The x-, y- and z-axis labels.}
  \item{scales}{A \code{\link{list}} determining how the axes are drawn;
    see \code{\link{cloud}()}.}
  \item{par.settings}{See \code{\link{cloud}()}.}
  \item{n}{The sample size.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{cloud}()}.}
}
\section{Methods}{
  Cloud plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{cloud}()}.
}
%\author{Marius Hofert}
\examples{
## For 'matrix' objects
cop <- gumbelCopula(2, dim = 3)
n <- 1000
set.seed(271)
U <- rCopula(n, copula = cop)
cloud2(U, xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))

## For 'copula' objects
set.seed(271)
cloud2(cop, n = n) # same as above

## For 'mvdc' objects
mvNN <- mvdc(cop, c("norm", "norm", "exp"),
             list(list(mean = 0, sd = 1), list(mean = 1), list(rate = 2)))
cloud2(mvNN, n = n)
}
\keyword{methods}
\keyword{hplot}