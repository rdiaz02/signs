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
library(party)
library(mboost)
load("dlbcl.benchmark.RData")


##library(snow)

fParal <- function(Minp = 0.15,
                   MaxIterationsCox = 200,
                   MaxSize = 100,
                   MinSize = 5,
                   MinCor = 0.65,
                   nfold = 10,
                   arrays = NULL,
                   genes = NULL) {
    print(gc())

    if(is.null(arrays)) arrays <- nrow(dlbcl.covar.t)
    if(is.null(genes)) genes <- ncol(dlbcl.covar.t)

    covar <- dlbcl.covar.t[1:arrays, ]
    survtime <- dlbcl.surv[1:arrays]
    status <- dlbcl.event[1:arrays]
    covar <- covar[, 1:genes]

    cat("\n gc after creating objects \n")
    print(gc())
    
    cat("\n\n Running one iteration \n\n")
    walltime <- unix.time({
        gb.all <- my.glmboost(covar, survtime, status, NULL)
        rm(gb.all)
        cat("\n gc after gb.all \n")
        print(gc())
        gb.cv.output <- my.glmboost.cv(covar, survtime, status)
    })[3]
    rm(covar)
    rm(survtime)
    rm(status)
    rm(gb.cv.output)
    cat("\n gc after deleting \n")
    print(gc())
    return(walltime)
}

library(Rmpi)
library(papply)
mpiSpawnAll()

glmboost.paral <- matrix(cbind(rep(NA, 8),
                                 narrays = c(rep(160, 7), 20, 40, 80, 100),
                                 ngenes = c(1000, 2000, 4000, 6000, 12000, 24000, 48000,
                                   rep(7399, 4))),
                           ncol = 3)
                         
colnames(glmboost.paral) <- c("time", "narrays", "ngenes")

dlbcl.covar.t <- t(dlbcl.covar)
rm(dlbcl.covar)

for(i in (1:nrow(glmboost.paral))) {
  cat(" \n\n arrays = ", glmboost.paral[i, 2],
      " genes = ", glmboost.paral[i, 3], "\n")
  
  glmboost.paral[i, 1] <- fParal(arrays = glmboost.paral[i, 2],
                                 genes = glmboost.paral[i, 3])
  
  assign(paste("glmboost.paral.", numero, "cpu", sep = ""),
         glmboost.paral)
  save(file = paste("glmboost.paral.", numero, "cpu.RData", sep = ""),
       list = paste("glmboost.paral.", numero, "cpu", sep = ""))
}

assign(paste("glmboost.paral.", numero, "cpu", sep = ""),
       glmboost.paral)
save(file = paste("glmboost.paral.", numero, "cpu.RData", sep = ""),
     list = paste("glmboost.paral.", numero, "cpu", sep = ""))




