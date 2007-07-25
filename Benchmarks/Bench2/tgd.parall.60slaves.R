rm(list = ls())

system("lamhalt")
system("lamboot -v lamb-host-60slaves.def")


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


library(Rmpi)
library(papply)
mpiSpawnAll()



tgd.paral.60cpu.mini <- matrix(cbind(rep(NA, 10),
                                    narrays = c(rep(40, 5), 20, 40, 80, 100, 120),
                                    ngenes = c(20, 40, 80, 160, 320, rep(40, 5))),
                              ncol = 3)
                         
colnames(tgd.paral.60cpu.mini) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(tgd.paral.60cpu.mini)) {
    tgd.paral.60cpu.mini[i, 1] <- fParal("dlbcl", arrays = tgd.paral.60cpu.mini[i, 2],
                                   genes = tgd.paral.60cpu.mini[i, 3])
    save(file = "tgd.paral.60cpu.mini.RData", tgd.paral.60cpu.mini)
}





tgd.paral.60cpu <- matrix(cbind(rep(NA, 8),
                         narrays = c(rep(160, 4), 20, 40, 80, 100),
                         ngenes = c(1000, 2000, 4000, 6000, rep(7399, 4))),
                         ncol = 3)
                         
colnames(tgd.paral.60cpu) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(tgd.paral.60cpu)) {
    tgd.paral.60cpu[i, 1] <- fParal("dlbcl", arrays = tgd.paral.60cpu[i, 2],
                                   genes = tgd.paral.60cpu[i, 3])
    save(file = "tgd.paral.60cpu.RData", tgd.paral.60cpu)

}





