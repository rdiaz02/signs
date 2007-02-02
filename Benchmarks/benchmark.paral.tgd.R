## Use the three data sets of bair & Tibshirani. "Predicting Cancer
## Patient Survival with Gene Expression Data"  in PLOS.
## The data were downloaded from the SignS example page
## http://signs.bioinfo.cnio.es/Examples/example.data.sets.SignS.tar.gz
## and read into R as follows:

## aml.covar <- read.table("aml.covar.txt", sep = "\t", header = TRUE,
##                         comment.char = "")[, -1]
## dlbcl.covar <- read.table("dlbcl.160.covar.txt", sep = "\t", header = TRUE,
##                         comment.char = "")[, -1]
## breast.covar <- read.table("breast.covar.txt", sep = "\t", header = TRUE,
##                         comment.char = "")[, -1]

## aml.surv <- scan("aml.surv.txt", what = double(0), sep = "\t")
## aml.surv <- aml.surv[-length(aml.surv)] ## trailing white turned into NA.

## dlbcl.surv <- scan("dlbcl.160.surv.txt", what = double(0), sep = "\t")
## dlbcl.surv <- dlbcl.surv[-length(dlbcl.surv)] ## trailing white turned into NA.

## breast.surv <- scan("breast.surv.txt", what = double(0), sep = "\t")
## breast.surv <- breast.surv[-length(breast.surv)] ## trailing white turned into NA.

## aml.event <- scan("aml.event.txt", what = double(0), sep = "\t")
## aml.event <- aml.event[-length(aml.event)] ## trailing white turned into NA.

## dlbcl.event <- scan("dlbcl.160.event.txt", what = double(0), sep = "\t")
## dlbcl.event <- dlbcl.event[-length(dlbcl.event)] ## trailing white turned into NA.

## breast.event <- scan("breast.event.txt", what = double(0), sep = "\t")
## breast.event <- breast.event[-length(breast.event)] ## trailing white turned into NA.

## save(file = "benchmark.data.sets.RData", list = ls())


load("benchmark.data.sets.RData")
library(SignS2)
library(Rmpi)
library(snow)
TheCluster <- makeCluster(60, "MPI")
clusterExport(TheCluster, c("lik1", "tgd1InternalSnow",
                            "tgdTrain", "tgdPieceInternalSnow"))



fParal <- function(dataset, epi = 5e-6,
                   maxstep = 5000,
                   checkEvery = 50000,
                   nfold = 10,
                   arrays = NULL,
                   genes = NULL) {
    print(gc())
    covar    <- t(get(paste(dataset, "covar", sep = "."))) ## they use subjects in rows
    survtime <- get(paste(dataset, "surv", sep = "."))
    status   <- get(paste(dataset, "event", sep = "."))

    if(!is.null(arrays)) {
        ## No, dlbcl shows no index trend
        covar <- covar[1:arrays, ]
        survtime <- survtime[1:arrays]
        status <- status[1:arrays]
    }
    if(!is.null(genes)) covar <- covar[, 1:genes]

    cat("\n\n Running one iteration \n\n")
    walltime <- unix.time({
        tmp <- tauBestP(covar, survtime, status, thres = c(0, 1),
                        epi = epi, maxiter = maxstep, nfold = nfold,
                        thresGrid = 6, checkEvery = checkEvery,
                        fitWithBest = FALSE)
        })[3]

    return(walltime)
}

nreps <- 5



s.20.40 <- replicate(nreps, fParal("dlbcl", arrays = 20, genes = 40))
s.40.40 <- replicate(nreps, fParal("dlbcl", arrays = 40, genes = 40))
s.80.40 <- replicate(nreps, fParal("dlbcl", arrays = 80, genes = 40))
s.100.40 <- replicate(nreps, fParal("dlbcl", arrays = 100, genes = 40))
s.120.40 <- replicate(nreps, fParal("dlbcl", arrays = 120, genes = 40))
    
    
s.40.20 <- replicate(nreps, fParal("dlbcl", arrays = 40, genes = 20))
s.40.80 <- replicate(nreps, fParal("dlbcl", arrays = 40, genes = 80))
s.40.160 <- replicate(nreps, fParal("dlbcl", arrays = 40, genes = 160))
s.40.320 <- replicate(nreps, fParal("dlbcl", arrays = 40, genes = 320))

save(list = ls(),
         file = "parallel.RData")



### Now, do the above with different number of slaves per node:
## 2, 6, 20, 60 (12 ?)











## ## what is the best checkEvery?
## paral.aml.100 <- replicate(nreps, fParal("aml", checkEvery = 100))
## paral.dlbcl.100 <- replicate(nreps, fParal("dlbcl", checkEvery = 100))
## paral.breast.100 <- replicate(nreps, fParal("breast", checkEvery = 100))
## paral.aml.500 <- replicate(nreps, fParal("aml", checkEvery = 500))
## paral.dlbcl.500 <- replicate(nreps, fParal("dlbcl", checkEvery = 500))
## paral.breast.500 <- replicate(nreps, fParal("breast", checkEvery = 500))
## paral.aml.1000 <- replicate(nreps, fParal("aml", checkEvery = 1000))
## paral.dlbcl.1000 <- replicate(nreps, fParal("dlbcl", checkEvery = 1000))
## paral.breast.1000 <- replicate(nreps, fParal("breast", checkEvery = 1000))
## paral.aml.ne <- replicate(nreps, fParal("aml", checkEvery = 50001))
## paral.dlbcl.ne <- replicate(nreps, fParal("dlbcl", checkEvery = 50001))
## paral.breast.ne <- replicate(nreps, fParal("breast", checkEvery = 50001))





## paral.dlbcl.subs <- sapply(rep(c(160, 80, 40, 20), 5),
##                          function(x) fParal("dlbcl", arrays = x, genes =3500))
## paral.dlbcl.gene <- sapply(rep(c(7000, 3500, 1750, 875), 5),
##                          function(x) fParal("dlbcl", arrays = 40, genes = x))


## paral.dlbcl.subs.ne <- sapply(rep(c(160, 80, 40, 20), 5),
##                          function(x) fParal("dlbcl", arrays = x, genes =3500, checkEvery = 50001))
## paral.dlbcl.gene.ne <- sapply(rep(c(7000, 3500, 1750, 875), 5),
##                          function(x) fParal("dlbcl", arrays = 40, genes = x, checkEvery = 50001))



## paral.aml <- replicate(nreps, fParal("aml", checkEvery = 1000))
