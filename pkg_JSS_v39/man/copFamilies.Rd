\name{copFamilies}
\title{Specific Archimedean Copula Families ("acopula" Objects)}
\concept{Copula Family}
\alias{acopula-families}
\alias{copAMH}
\alias{copClayton}
\alias{copFrank}
\alias{copGumbel}
\alias{copJoe}
\docType{data}
\description{
  Specific Archimedean families (\code{"\linkS4class{acopula}"} objects)
  implemented in the package \pkg{nacopula}.
}
\usage{
copAMH
copClayton
copFrank
copGumbel
copJoe
}
\value{
  A \code{"\linkS4class{acopula}"} object.
}
\details{
  All these are objects of the formal class \code{"\linkS4class{acopula}"}.
  \describe{
    \item{\code{copAMH}:}{Archimedean family of Ali-Mikhail-Haq with
      parametric generator
      \deqn{\psi(t)=(1-\vartheta)/(\exp(t)-\vartheta),\	t\in[0,\infty],
      }{psi(t)=(1-theta)/(exp(t)-theta), t in [0,Inf],}
      with \eqn{\vartheta\in[0,1)}{theta in [0,1)}. The range of
      admissible Kendall's tau is [0,1/3). Note that the lower and upper
      tail-dependence coefficients are both zero, i.e., this copula
      family does not allow for tail dependence.}

     \item{\code{copClayton}:}{Archimedean family of Clayton with
       parametric generator
       \deqn{\psi(t)=(1+t)^{-1/\vartheta},\ t\in[0,\infty],}{
	      psi(t)=(1+t)^{-1/theta}, t in [0,Inf],}
       with \eqn{\vartheta\in(0,\infty)}{theta in (0,Inf)}. The range of
     admissible Kendall's tau, as well as that of the lower
     tail-dependence coefficient, is (0,1). Note that this copula does
     not allow for upper tail dependence.}

     \item{\code{copFrank}:}{Archimedean family of Frank with parametric
       generator
       \deqn{-\log(1-(1-e^{-\vartheta})\exp(-t))/\vartheta,\ t\in[0,\infty]}{
	     -log(1-(1-e^{-theta})exp(-t))/theta, t in [0,Inf],}
     with \eqn{\vartheta\in(0,\infty)}{theta in (0,Inf)}. The range of
     admissible Kendall's tau is (0,1). Note that this copula family
     does not allow for tail dependence.}

     \item{\code{copGumbel}:}{Archimedean family of Gumbel with
       parametric generator
       \deqn{\exp(-t^{1/\vartheta}),\ t\in[0,\infty]}{
	      exp(-t^{1/theta}), t in [0,Inf],}
	    with
     \eqn{\vartheta\in[1,\infty)}{theta in [1,Inf)}. The range of
     admissible Kendall's tau, as well as that of the upper
     tail-dependence coefficient, is [0,1). Note that this copula does
     not allow for lower tail dependence.}

     \item{\code{copJoe}:}{Archimedean family of Joe with parametric
       generator
       \deqn{1-(1-\exp(-t))^{1/\vartheta},\ t\in[0,\infty]}{
	     1-(1-exp(-t))^{1/theta}, t in [0,Inf],}
       with \eqn{\vartheta\in[1,\infty)}{theta in [1,Inf)}. The range of
       admissible Kendall's tau, as well as that of the upper
       tail-dependence coefficient, is [0,1). Note that this copula does
       not allow for lower tail dependence.}
  }
  Note that staying within one of these Archimedean families, all of
  them can be nested if two (generic) generator parameters
  \eqn{\vartheta_0}{theta.0}, \eqn{\vartheta_1}{theta.1} satisfy
  \eqn{\vartheta_0\le\vartheta_1}{theta.0 <= theta.1}.
}
\seealso{
  The class definition, \code{"\linkS4class{acopula}"}.
  \cr
  \code{\link{getAcop}} accesses these families
  \dQuote{programmatically}.
}
\author{Marius Hofert}
\references{
  Hofert, M. (2010a),
  Efficiently sampling nested Archimedean copulas,
  \emph{Computational Statistics & Data Analysis}, in press.

  Hofert, M. (2010b),
  \emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
  Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.

  Marius Hofert and Martin Maechler (2011)
  Nested Archimedean Copulas Meet R: The nacopula Package.
  \emph{Journal of Statistical Software} \bold{39}(9), 1--20.
  \url{http://www.jstatsoft.org/v39/i09/}.
}
\examples{
## Print a copAMH object and its structure
copAMH
str(copAMH)

## Show admissible parameters for a Clayton copula
copClayton@paraInterval

## Generate random variates from a Log(p) distribution via V0 of Frank
p <- 1/2
copFrank@V0(100, -log(1-p))

## Plot the upper tail-dependence coefficient as a function in the
## parameter for Gumbel's family
curve(copGumbel@lambdaU(x), xlim = c(1, 10), ylim = c(0,1), col = 4)

## Plot Kendall's tau as a function in the parameter for Joe's family
curve(copJoe@tau(x), xlim = c(1, 10), ylim = c(0,1), col = 4)

## ------- Plot psi() and tau() - and properties of all families ----

## The copula families currently provided:
(famNms <- ls("package:nacopula", patt="^cop[A-Z]"))

op <- par(mfrow=c(length(famNms),2), mar = .6+c(2,1.4,1,1), mgp = c(1.1, 0.4,0))
for(nm in famNms) { Cf <- get(nm)
   thet <- Cf@tauInv(0.3)
   curve(Cf@psi(x, theta = thet), 0, 5,
         xlab=expression(x), ylab="", ylim=0:1,
         main = substitute(list(NAM *"  "* psi(x, theta == TH), tau == 0.3),
                           list(NAM=Cf@name, TH=thet)))
   I <- Cf@paraInterval
   Iu <- pmin(10, I[2])
   ## if(I[2] == Inf) { Iu <- 1000; log <- "x" } else { Iu <- I[2]; log <- "" }
   curve(Cf@tau(x), I[1], Iu,
         xlab = substitute(theta \%in\% II, list(II = format(I))), ylab="",
         main=substitute(NAM * " " * tau(theta), list(NAM=Cf@name)))
}
par(op)

## Construct a bivariate Clayton copula with parameter theta
theta <- 2
C2 <- onacopula("Clayton", C(theta, 1:2))
C2@copula # is an "acopula" with specific parameter theta

curve(C2@copula@psi(x, C2@copula@theta),
      main = quote("Generator" ~~ psi ~~ " of Clayton A.copula"),
      xlab = quote(theta1), ylab = quote(psi(theta1)),
      xlim = c(0,5), ylim = c(0,1), col = 4)

## What is the corresponding Kendall's tau?
C2@copula@tau(theta) # 0.5

## What are the corresponding tail-dependence coefficients?
C2@copula@lambdaL(theta)
C2@copula@lambdaU(theta)

## Generate n pairs of random variates from this copula
U <- rnacopula(n = 1000, C2)
## and plot the generated pairs of random variates
plot(U, asp=1, main = "n = 1000 from  Clayton(theta = 2)")
}
\keyword{datasets}
