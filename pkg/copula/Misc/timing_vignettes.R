## Idea, etc: By Marius Hofert

## Timing vignettes


setwd("../vignettes")
Sys.unsetenv("TEXINPUTS") # <<<! currently MM's TEXINPUTS give substantial slowdown (NAS file access)

## Timing command
timer <- function(x) {
    tm <- system.time(res <- tryCatch(x, error = function(e) e))
    if(inherits(res, "error")) NA else round(tm[["elapsed"]], 2)
}
do.rm <- !(Sys.info()[["user"]] %in% c("maechler", "<you>"))#

## Determine all vignettes
all <- system("ls", intern = TRUE)
(all.Rnw <- all[grepl("\\.Rnw$", all)])
(all.Rmd <- all[grepl("\\.Rmd$", all)])
(n <- length(all.vign <- c(all.Rnw, all.Rmd)))
if(n < 1)
    stop("Found no vignettes in the current working directory")

## Time measurement

Rcmd <- tools::Rcmd
basenameNoExt <- function(.) tools::file_path_sans_ext(basename(.))

tm <- numeric(n)
for(j in seq_along(all.Rnw)) { # iterate over all .Rnw
    cat("==> File ", fil <- all.Rnw[j],": ")
    tm[j] <- timer(Rcmd(paste("Sweave --pdf", fil)))
    if(do.rm) {
        all <- system(paste0("ls ",basenameNoExt(fil),"*"), intern = TRUE)
        file.remove(all[!grepl(fil, all, fixed=TRUE)])# keep *.Rnw file (and backups of that)
    }
    cat("\n \\--> finished ", fil, "\n-==-==-==-==-\n")
}
n1 <- length(all.Rnw)

## rmarkdown::render() :
## - Careful with environments: the code is executed in parent.frame(); see, ?render.
## - For example, if a global counter 'i' is used here, it is overwritten by wild_animals.Rmd.
## - If new.env() is used in render(), then dNAC.Rmd and AR_Clayton.Rmd fail with weird problems.
for(j in seq_along(all.Rmd)) { # iterate over all .Rmd
    cat("==> File ", fil <- all.Rmd[j],": ")
    tm[j+n1] <- timer(rmarkdown::render(fil))
    if(do.rm) file.remove(paste0(basenameNoExt(fil),".html"))
    cat("\n-==-==-==-==-\n finished ", fil, "\n")
}
if(do.rm) file.remove("Rplots.pdf")

## Results
timings <- data.frame("Vignette" = all.vign,
                      "Elapsed.sec" = tm)
sink("../Misc/timings.out", split = TRUE)
timings[order(timings[,"Elapsed.sec"], decreasing = TRUE), ]
sum(timings[,"Elapsed.sec"])
sink()

## 25.Jan. 2023: *with* slow TEXINPUTS,  no TEXINPUTS
##                  Vignette Elapsed.sec  Elapsed.sec
## 6                 GIG.Rmd       35.76        26.17
## 2        nacopula-pkg.Rnw       28.50        19.47
## 1         Frank-Rmpfr.Rnw       26.53        11.02
## 3        rhoAMH-dilog.Rnw       24.68         9.64
## 11 empiricial_copulas.Rmd       12.80         7.62
## 12 logL_visualization.Rmd        9.59         6.93
## 8                NALC.Rmd        9.50         5.59
## 13               qrng.Rmd        6.61         5.38
## 4        AC_Liouville.Rmd        6.26         4.58
## 14       wild_animals.Rmd        6.09         4.50
## 10               dNAC.Rmd        5.86         2.16
## 9        copula_GARCH.Rmd        3.22         1.82
## 7                HAXC.Rmd        2.00         1.47
## 5          AR_Clayton.Rmd        1.56         1.42
##                                ------       ------
sum(timings[,"Elapsed.sec"])  #   178.96       107.77
