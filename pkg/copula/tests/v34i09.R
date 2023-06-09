#### The R code from
####  Kojadinovic, I. and Yan, J. (2010).
####  Modeling Multivariate Distributions with Continuous Margins
####  Using the copula R Package.
####  Journal of Statistical Software *34*(9), 1--20.
####
#### === http://www.jstatsoft.org/v34/i09/supp/2
##---------------------------------------------------------------------------
## The following changes are made here: instead of   'N = 1000'  (5x)
## we use   'N = N'  and set N here:
N <- 1000 # == original simulation for paper; we use smaller N below !
doExtras <- FALSE # (for ease of manual evaluation)
(doExtras <- copula:::doExtras())

N  <- if(doExtras) 100 else 12
N. <- if(doExtras) N else 5 # for really expensive parts
options(warn = 1)# [print as they happen]
## [plus  one other (doExtras) clause below]
##---------------------------------------------------------------------------


###################################################
### lossPrep
###################################################
library("copula")
data("loss")
myLoss <- subset(loss, censored==0, select=c("loss", "alae"))
str(myLoss) ## 1466 x 2 (integer)


###################################################
### breakTies
###################################################
set.seed(123)
pseudoLoss     <- sapply(myLoss, rank, ties.method="random") / (nrow(myLoss) + 1)
pseudoLoss.ave <- sapply(myLoss, rank                      ) / (nrow(myLoss) + 1)

###################################################
### lossPlot
###################################################
par(mfrow=c(1, 2), mgp=c(1.5, 0.5, 0), mar=c(3.5,2.5,0,0))
plot(pseudoLoss, sub="(a) random rank for ties")
plot(pseudoLoss.ave, sub="(b) average rank for ties")


###################################################
### lossIndep
###################################################
(nm <- sprintf("indTstSim_n=%d,p=%d,N=%d", nrow(pseudoLoss), p = 2, N = N))
sFil <- system.file("rData", paste0(nm, ".rds"), package = "copula")
## [not yet] - worth it?
empsamp <- if(file.exists(sFil)) {
               readRDS(sFil)
           } else {
               print(system.time(
                   iTS <- indepTestSim(nrow(pseudoLoss), p = 2, N = N, verbose = FALSE)))
               ## N=100: 16.5 sec [lynne 2015]
               iTS
           }

indepTest(pseudoLoss, empsamp)

###################################################
### lossGof
###################################################
gum1 <- gumbelCopula(1, use.indepC="FALSE")# not the indep.copula
system.time(gofGumb.pb <- gofCopula(gum1, pseudoLoss, estim.method="itau",
                                    simulation="pb", N = N, verbose = FALSE))
## N=100: 10.7 sec [lynne 2015]
gofGumb.pb
system.time(gofGumb.mult <- gofCopula(gum1, pseudoLoss, estim.method="itau",
                                      simulation="mult", N = N))
## N=100: 2.3 sec [lynne 2015]
gofGumb.mult


###################################################
### lossFit
###################################################
fitCopula(gumbelCopula(), pseudoLoss, method="itau")


###################################################
### repeat
###################################################
myAnalysis <- function(u, verbose=TRUE) {
  u.pseudo <- sapply(u, rank, ties.method="random") / (nrow(u) + 1)
  indTest <- indepTest(u.pseudo, empsamp)$global.statistic.pvalue

  pv.gof <- function(COP) { if(verbose) cat("gofC..(",class(COP),", ...) ")
      gofCopula(COP, u.pseudo, estim.method="itau", simulation="mult", N = N)$pvalue
                            if(verbose) cat("[done]\n") }

  gum1 <- gumbelCopula(1, use.indepC="FALSE")# not the indep.copula
  gof.g <- pv.gof(gum1)
  gof.c <- pv.gof(claytonCopula(1))
  gof.f <- pv.gof(frankCopula(1))
  gof.n <- pv.gof(normalCopula(0))
  gof.p <- pv.gof(plackettCopula(1))
  gof.t <- pv.gof(tCopula(0, df = 4, df.fixed = TRUE))

  fit.g <- fitCopula(gum1, u.pseudo, method="itau")
  c(indep = indTest, gof.g=gof.g, gof.c=gof.c, gof.f=gof.f, gof.n=gof.n, gof.t=gof.t, gof.p=gof.p,
    est = fit.g@estimate, se = sqrt(fit.g@var.est))
}

system.time(print(round(
if(doExtras) {
    ##was myReps <- t(replicate(100, myAnalysis(myLoss)))
    myReps <- t(replicate(12, myAnalysis(myLoss)))
    apply(myReps, 2, summary)
} else myAnalysis(myLoss), 3)))
###################################################
### lossGof3
###################################################
gofCopula(gum1, pseudoLoss.ave, estim.method="itau", simulation="mult", N = N)


###################################################
### srPrep
###################################################
data("rdj")
str(rdj)
apply(rdj[,2:4], 2, function(x) length(unique(x)))

###################################################
### pseudoobs
###################################################
pseudoSR <- pobs(rdj[,2:4])


###################################################
### srMultSerialIndepTest
###################################################
set.seed(123)
system.time(srMultSerialIndepTest <-
            multSerialIndepTest(rdj[,2:4]^2, lag.max=4, N = N, verbose = FALSE))
srMultSerialIndepTest
dependogram(srMultSerialIndepTest)


###################################################
### multIndepTest
###################################################
system.time(empsamp <- indepTestSim(nrow(pseudoSR), p = 3, N = N, verbose = FALSE))
srMultIndepTest <- indepTest(pseudoSR, empsamp)
srMultIndepTest
dependogram(srMultIndepTest)


###################################################
### srGof
###################################################
tC3 <- tCopula(c(0,0,0), dim=3, dispstr="un", df=5, df.fixed=TRUE)

if(FALSE) {
    ## Cannot work yet as pCopula() cannot accept t copulas with non-integer
    ## degrees of freedom - because mvtnorm::pmvt() cannot :
    system.time(srGof.t.pboo <- gofCopula(tC3, pseudoSR, N = N., estim.method="mpl"))
    srGof.t.pboo
}
system.time(srGof.t.mult <-
                gofCopula(tC3, pseudoSR, N = N, estim.method="mpl", simulation="mult"))
## 6.2 sec [lynne 2015]
srGof.t.mult


###################################################
### srFit
###################################################
(fm.tC3 <- fitCopula(tC3, pseudoSR, method="mpl"))
