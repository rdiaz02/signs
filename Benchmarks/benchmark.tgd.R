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

source("gd1.R")

library(SignS)




fSeq <- function(dataset, epi = 5e-6,
                 maxstep = 10000,
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

    
    thresS <- seq(from = 0, to = 1, length.out = 6)

    cat("\n\n Running one iteration \n\n")
    walltime <- unix.time({
        for(ths in thresS) {
            tmp <- gdcvpl(covar, survtime, status, ths,
                          epi = epi, maxstep = maxstep, nfold = nfold)
        }
    })[3]

    return(walltime)
}



fParal <- function(dataset, epi = 5e-6,
                   maxstep = 10000,
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
                        thresGrid = 6, checkEvery = 50,
                        fitWithBest = FALSE)
        })[3]

    return(walltime)
}


nreps <- 5

serial.aml <- replicate(nreps, fSeq("aml"))
serial.dlbcl <- replicate(nreps, fSeq("dlbcl"))
serial.breast <- replicate(nreps, fSeq("breast"))

serial.dlbcl.subs <- sapply(rep(c(160, 80, 40, 20), 5),
                         function(x) fSeq("dlbcl", arrays = x, genes =3500))
serial.dlbcl.gene <- sapply(rep(c(7000, 3500, 1750, 875), 5),
                         function(x) fSeq("dlbcl", arrays = 40, genes = x))


paral.aml <- replicate(nreps, fParal("aml"))
paral.dlbcl <- replicate(nreps, fParal("dlbcl"))
paral.breast <- replicate(nreps, fParal("breast"))

paral.dlbcl.subs <- sapply(rep(c(160, 80, 40, 20), 5),
                         function(x) fParal("dlbcl", arrays = x, genes =3500))
paral.dlbcl.gene <- sapply(rep(c(7000, 3500, 1750, 875), 5),
                         function(x) fParal("dlbcl", arrays = 40, genes = x))
