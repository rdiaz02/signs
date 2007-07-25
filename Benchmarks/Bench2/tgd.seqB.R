# run as: (with the param 0 or 1) nohup R --no-save --no-restore --slave --args 0 < benchmark.serial.tgd.R > serial.many.Rout &


ll <- length(commandArgs())
doone <- as.numeric(commandArgs()[ll])

load("benchmark.data.sets.RData")

source("gd1.R")


fSeq <- function(dataset, epi = 5e-6,
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


tgd.seqB <- matrix(cbind(rep(NA, 8),
                         narrays = c(rep(160, 4), 20, 40, 80, 100),
                         ngenes = c(1000, 2000, 4000, 6000, rep(7399, 4))),
                         ncol = 3)
                         
colnames(tgd.seqB) <- c("time", "narrays", "ngenes")


## with doone we controll where we start.

        tgd.seqB[4, 1] <- fSeq("dlbcl", arrays = tgd.seqB[4, 2],
                                genes = tgd.seqB[4, 3])
        save(file = "tgd.seqB.RData", tgd.seqB)
