\name{contourplot2-methods}
\alias{contourplot2}
\alias{contourplot2-methods}
\alias{contourplot2,matrix-method}
\alias{contourplot2,data.frame-method}
\alias{contourplot2,copula-method}
\alias{contourplot2,mvdc-method}
\docType{methods}
\title{Methods `contourplot2' in Package `copula'}
\description{
  Methods \code{contourplot2()} to draw contour plots of
  two dimensional distributions from package \pkg{copula}.
}
\usage{
\S4method{contourplot2}{matrix}(x, aspect = 1,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      xlab = NULL, ylab = NULL,
      region = TRUE, col.regions = gray(seq(0.4, 1, length.out = 100)),
      cuts = 16, labels = !region,
      scales = list(alternating = c(1,1), tck = c(1,0)), \dots)
\S4method{contourplot2}{data.frame}(x, aspect = 1,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      xlab = NULL, ylab = NULL,
      region = TRUE, col.regions = gray(seq(0.4, 1, length.out = 100)),
      cuts = 16, labels = !region,
      scales = list(alternating = c(1,1), tck = c(1,0)), \dots)
\S4method{contourplot2}{copula}(x, FUN, n.grid = 26, xlim = 0:1, ylim = 0:1,
      xlab = expression(u[1]), ylab = expression(u[2]), \dots)
\S4method{contourplot2}{mvdc}(x, FUN, n.grid = 26, xlim, ylim,
      xlab = expression(x[1]), ylab = expression(x[2]), \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{aspect}{The aspect ratio.}
  \item{xlim, ylim}{The x- and y-axis limits.}
  \item{xlab, ylab}{The x- and y-axis labels. If at least one is
    \code{NULL}, useful \code{xlab} and \code{ylab} are determined
    automatically; the behavior depends on the class of \code{x}.}
  \item{labels}{A \code{\link{logical}} indicating whether the contour
    lines are labeled; see \code{\link{contourplot}()}.}
  \item{region}{A \code{\link{logical}} indicating whether regions
    between contour lines should be filled as in a level plot; see
    \code{\link{contourplot}()}.}
  \item{col.regions}{The colors of the regions if \code{region = TRUE};
    see \code{\link{contourplot}()}.}
  \item{cuts}{The number of levels; see \code{\link{contourplot}()}.}
  \item{scales}{A \code{\link{list}} determining how the axes are drawn;
    see \code{\link{contourplot}()}.}
  \item{FUN}{The \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{The number of grid points used in each dimension. This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{contourplot}()}.}
}
\section{Methods}{
  Contourplot plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{contourplot}()}.
}
%\author{Marius Hofert}
\seealso{
  The \code{\link{contour-methods}} for drawing perspective plots via base graphics.
}
\examples{
## For 'matrix' objects
## The Frechet--Hoeffding bounds W and M
n.grid <- 26
u <- seq(0, 1, length.out = n.grid)
grid <- expand.grid("u[1]" = u, "u[2]" = u)
W <- function(u) pmax(0, rowSums(u)-1) # lower bound W
M <- function(u) apply(u, 1, min) # upper bound M
x.W <- cbind(grid, "W(u[1],u[2])" = W(grid)) # evaluate W on 'grid'
x.M <- cbind(grid, "M(u[1],u[2])" = M(grid)) # evaluate M on 'grid'
contourplot2(x.W) # contour plot of W
contourplot2(x.M) # contour plot of M

## For 'copula' objects
cop <- frankCopula(-4)
contourplot2(cop, pCopula) # the copula
contourplot2(cop, pCopula, xlab = "x", ylab = "y") # adjusting the labels
contourplot2(cop, dCopula) # the density

## For 'mvdc' objects
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
             list(list(mean = 0, sd = 1), list(mean = 1)))
xl <- c(-2, 2)
yl <- c(-1, 3)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, contour = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, region = FALSE, labels = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, region = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl,
             col.regions = colorRampPalette(c("royalblue3", "maroon3"), space="Lab"))
}
\keyword{methods}
\keyword{hplot}