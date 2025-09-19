rm(list = ls())

ll <- length(commandArgs())
numero <- as.numeric(commandArgs()[ll])
cat("\n numero is ", numero, "\n")

system("lamhalt")
if(numero != 1)
    system(paste("lamboot -v lamb-host-", numero, "slaves.def", sep = ""))



library(SignS2)
library(papply)
library(survival)
load("dlbcl.benchmark.RData")


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


if(numero != 10) {
  tgd.paral.new <- matrix(cbind(rep(NA, 3),
                                narrays = rep(160, 3),
                                ngenes = c(12000, 24000, 48000)),
                          ncol = 3)
} else {
  tgd.paral.new <- matrix(cbind(rep(NA, 8),
                                narrays = c(rep(160, 7), 20, 40, 80, 100),
                                ngenes = c(1000, 2000, 4000, 6000, 12000, 24000, 48000,
                                  rep(7399, 4))),
                          ncol = 3)
}

colnames(tgd.paral.new) <- c("time", "narrays", "ngenes")

for(i in 1:nrow(tgd.paral.new)) {
    
    cat(" \n\n arrays = ", tgd.paral.new[i, 2],
        " genes = ", tgd.paral.new[i, 3], "\n")
    
    tgd.paral.new[i, 1] <- fParal("dlbcl", arrays = tgd.paral.new[i, 2],
                                genes = tgd.paral.new[i, 3])
    assign(paste("tgd.paral.new.", numero, "cpu", sep = ""),
           tgd.paral.new)
    save(file = paste("tgd.paral.new.", numero, "cpu.RData", sep = ""),
         list = paste("tgd.paral.new.", numero, "cpu", sep = ""))
}

assign(paste("tgd.paral.new.", numero, "cpu", sep = ""),
       tgd.paral.new)
save(file = paste("tgd.paral.new.", numero, "cpu.RData", sep = ""),
     list = paste("tgd.paral.new.", numero, "cpu", sep = ""))
