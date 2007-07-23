load("benchmark.data.sets.RData")
library(SignS2)
library(papply)
library(survival)

##library(snow)

fParal <- function(dataset,
                   Minp = 0.1,
                   MaxIterationsCox = 200,
                   MaxSize = 100,
                   MinSize = 10,
                   MinCor = 0.8,
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
        all.res1 <- dStep1.parallel(covar, survtime, status,
                                    Minp, MaxIterationsCox)
        try(all.res3 <- fitDave.res1Given(covar, survtime, status,    
                                      res1 = all.res1,            
                                      Minp, MaxSize,              
                                      MinSize, MinCor,            
                                      MaxIterationsCox,
                                      plot = FALSE,
                                      interactive = FALSE))
        try(cvDaveRun <- cvDave.parallel3(x = covar, time = survtime,
                                      event = status,
                                      p = Minp, maxSize = MaxSize,
                                      minSize = MinSize,
                                      minCor = MinCor,
                                      MaxIterationsCox = MaxIterationsCox,
                                      nfold = nfold,
                                      mpiHosts = NULL))
    })[3]
    return(walltime)
}


library(Rmpi)
library(papply)
mpiSpawnAll()
