rm(list = ls())

load("benchmark.data.sets.RData")
library(SignS2)
library(papply)

##library(snow)

fParal <- function(dataset, epi = 5e-6,
                   maxstep = 5000,
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
                        thresGrid = 6, 
                        fitWithBest = FALSE)
        })[3]

    return(walltime)
}

nreps <- 1


## the new code running sequentially

seq.20.20 <- replicate(nreps, fParal("dlbcl", arrays = 20, genes = 40))




### the parallel code
library(Rmpi)
library(papply)
mpiSpawnAll()

par.20.20 <- replicate(nreps, fParal("dlbcl", arrays = 20, genes = 40))




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

