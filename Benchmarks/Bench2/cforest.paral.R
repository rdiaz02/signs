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
library(MASS)
load("dlbcl.benchmark.RData")


##library(snow)

fParal <- function(dataset,
                   Minp = 0.15,
                   MaxIterationsCox = 200,
                   MaxSize = 100,
                   MinSize = 5,
                   MinCor = 0.65,
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
        cf.all <- my.cforest(covar, survtime, status, 200, NULL)
        cf.cv.output <- my.cforest.cv(covar, survtime, status, 200)
    })[3]
    return(walltime)
}


library(Rmpi)
library(papply)
mpiSpawnAll()




cforest.paral <- matrix(cbind(rep(NA, 8),
                                 narrays = c(rep(160, 7), 20, 40, 80, 100),
                                 ngenes = c(1000, 2000, 4000, 6000, 12000, 24000, 48000,
                                   rep(7399, 4))),
                           ncol = 3)
                         
colnames(cforest.paral) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(cforest.paral)) {
  cat(" \n\n arrays = ", cforest.paral[i, 2],
      " genes = ", cforest.paral[i, 3], "\n")
  
  cforest.paral[i, 1] <- fParal("dlbcl", arrays = cforest.paral[i, 2],
                             genes = cforest.paral[i, 3])
  
  assign(paste("cforest.paral.", numero, "cpu", sep = ""),
         cforest.paral)
  save(file = paste("cforest.paral.", numero, "cpu.RData", sep = ""),
       list = paste("cforest.paral.", numero, "cpu", sep = ""))
}

assign(paste("cforest.paral.", numero, "cpu", sep = ""),
       cforest.paral)
save(file = paste("cforest.paral.", numero, "cpu.RData", sep = ""),
     list = paste("cforest.paral.", numero, "cpu", sep = ""))




