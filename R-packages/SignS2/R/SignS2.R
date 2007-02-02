###  Not here because I do not want Rmpi or CGIwithR executed in the slaves.
#### library(CGIwithR)
#### library(Rmpi)
#### library(survival)
#### library(combinat)
#### library(MASS)

require(GDD)
require(survival)
require(imagemap)

imClose <- function (im) {
    ## prevent all the "Closing PNG device ..."
    dev.off(im$Device)
}

tauBestP <- function(x, time, event, thres = c(0, 1),
                     epi = 5e-06, thresGrid = 6, 
                     maxiter = 5000, checkEvery = 50000, ## we do not want early stopping
                     nfold = 10, fitWithBest = TRUE) {
    ## allows early stopping, but therefore increases communicationg overhead

    thresGrid <- 6 ## number of values of thres tried;
    ## you can change this, but then be sure to change the "60" below
    ## in the calls to clusterApply

    ## everything is carefully chosen to use 60 slaves (or more).
    
    if((!is.vector(time)) | (!is.vector(event)))
        stop("Time and event should be vectors")
    n <- length(time)
    thresS <- seq(from = thres[1], to = thres[2], length.out = thresGrid)
    cvpl.mat <- array(NA, c(length(thresS), maxiter))
    cvindex <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    clusterParams <- expand.grid(1:nfold, thresS)
    totalProcs <- dim(clusterParams)[1]


    if(checkEvery >= maxiter) {
        print("WARNING: checkEvery >= maxiter. Early stopping disabled.")
        checkEvery <- maxiter
        earlyStop <- FALSE
    } else {
        earlyStop <- TRUE
    }
    
    xTheCluster <<- x
    tTheCluster <<- time
    eTheCluster <<- event
    epTheCluster <<- epi
    mitTheCluster <<- maxiter
    ceTheCluster <<- checkEvery
    cviTheCluster <<- cvindex
    clusterParams <<- clusterParams
    
    time.initial.export <-
        unix.time(clusterExport(TheCluster,
                                c("xTheCluster", "tTheCluster",
                                  "eTheCluster",
                                  "epTheCluster",
                                  "mitTheCluster","ceTheCluster",
                                  "cviTheCluster",
                                  "clusterParams")))
    cat("\n\n     Export to slaves took ", time.initial.export[3],
        " seconds\n")
    ## clean master too
    rm(list = c("xTheCluster", "tTheCluster", "eTheCluster",
       "epTheCluster","mitTheCluster", "ceTheCluster",
       "cviTheCluster", "clusterParams"), envir = .GlobalEnv)

    ts1 <- unix.time(clusterApply(TheCluster, 1:60, tgdSnowSetUp))
    cat("\n \n      Slaves setup took ", ts1[3], " seconds\n")
    
    clusterOutput <- list()
    t.r1 <-
        unix.time(
                  clusterOutput <-
                  clusterApply(TheCluster, 1:60,
                                 function(x) tgd1InternalSnow()))
    cat("\n \n      First run took ", t.r1[3], " seconds\n\n")

    tmp.cvpl <- matrix(unlist(clusterOutput),
                      ncol = checkEvery, byrow = TRUE)

    cvpllymtgd <- matrix(NA, nrow = thresGrid, ncol = maxiter)
    ## Set cvpl to NULL to make the stored
    ## and later transmitted object as small as possible
    iter <- 1
    for(th in 1:thresGrid) {
        cvpllymtgd[th, iter : (iter + checkEvery - 1)] <-
            apply(tmp.cvpl[clusterParams[, 2] == thresS[th], ], 2, sum)
    }
    
    ## I do not do anything with the pscore returned. They are not used
    ## anywhere else.

    if (!earlyStop) {
        anyStillRunning <- FALSE
        cvpl.mat <- cvpllymtgd/nfold
    } else {
        anyStillRunning <- TRUE
        runIts <- rep(1, 60)
        indexThresRunning <- 1:thresGrid
        iter <- iter + checkEvery
    }
        
    while(anyStillRunning) { ## we only enter here if earlyStop

        t.r3 <-
            unix.time(
                      clusterOutput <-
                      clusterApply(TheCluster,
                                   runIts,
                                   function(vv) tgdPieceInternalSnow(vv)))
        cat("\n \n .... additional runs took ", t.r3[3], " seconds\n")

        tmp.cvpl <- matrix(unlist(clusterOutput),
                           ncol = checkEvery, byrow = TRUE)
        
        for(th in indexThresRunning) {
            cvpllymtgd[th, iter : (iter + checkEvery - 1)] <-
                apply(tmp.cvpl[clusterParams[, 2] == thresS[th], ], 2, sum)
        }        

        cvpl.mat <- cvpllymtgd/nfold
        iter <- iter + checkEvery

        ## Which thresholds should be stopped?
        for(thh in indexThresRunning) {
            tmp.cvplscore <- cvpl.mat[thh, ]
            checkAt <- iter - round(checkEvery * (4:0)/4) - 1
            if(all(diff(tmp.cvplscore[checkAt]) > 0)) {
                cat("\n Exiting iterations:\n")
                cat("   CVPL getting monotonically worse\n")
                cat(paste("   at threshold = ", thresS[thh],
                          "; max steps = ", iter - 1 , "\n"))
                cat(paste("minimum = ", min(tmp.cvplscore[1:(iter-1)], na.rm = TRUE),
                          "at ", which.min(tmp.cvplscore[1:(iter -1)]), " steps.\n\n\n"))
                indexThresRunning <- setdiff(indexThresRunning, thh)
            } else {
                allMin <- min(tmp.cvplscore[1:(iter - 1)])
                sectionMin <- min(tmp.cvplscore[(iter - checkEvery):(iter - 1)])
                sectionCoef <-
                    lm(tmp.cvplscore[(iter - checkEvery):(iter - 1)] ~
                       c(1:checkEvery))$coeff[2]
                if ((allMin < sectionMin) & (sectionCoef > 0)) {
                    cat("\n Exiting iterations:")
                    cat("\n   CVPL above the minimum and getting worse (positive slope)\n")
                    cat(paste("   at threshold = ", thresS[thh],
                              "; max steps = ", iter - 1 , "\n"))
                    cat(paste("minimum = ", min(tmp.cvplscore[1:(iter - 1)], na.rm = TRUE),
                              "at ", which.min(tmp.cvplscore[1:(iter - 1)]), " steps.\n\n\n"))
                ## we are above the minimum, and going up
                    indexThresRunning <- setdiff(indexThresRunning, thh)
                }
            }
        }
        ## get back indices to run
        if(length(indexThresRunning)) {
            theseProcs <-
                as.vector(sapply(indexThresRunning,
                                 function(x) ((x - 1) * nfold +1): (x * nfold)))
            runIts <- rep(0, 60)
            runIts[theseProcs] <- 1
            ## paranoid checks
            if(length(unique(clusterParams[theseProcs, 2])) !=
               length(indexThresRunning)) stop("oops, should not be here [21]")
        } else anyStillRunning <- FALSE 
        if( (iter + checkEvery - 1) > maxiter) anyStillRunning <- FALSE
    } ## closes while(stillRunning)
    
    ## allow for possibl of several minima
    tmp <- which(cvpl.mat == min(cvpl.mat, na.rm = TRUE), arr.ind = TRUE)
    tmp <- tmp[nrow(tmp), ]
    thresBest <- thresS[tmp[1]]
    stepBest <- tmp[2]

    if(is.na(stepBest)) {
        caughtUserError("An error occured, probably related to numerical problems. Please let us know")
    }

    ## do only if asked for.
    if(fitWithBest) {
        tgd.alldata <- tgdTrain(x, time, event, thresBest, epi,
                                steps = stepBest)
    } else tgd.alldata <- NA

## To incorporate ROC, do it here:
    ## train best models for each, and then find ROC, with CV.

    betas <- tgd.alldata$beta
    return(list(betas = betas, nonZeroBetas = sum(betas != 0), 
                threshold = thresBest, epi = epi,
                step = stepBest, cvpl.mat = cvpl.mat,
                tgd.alldata = tgd.alldata,
                thres.loc = tmp[1]))
}


tauBestP.noearly <- function(x, time, event, thres = c(0, 1),
                     epi = 5e-06, thresGrid = 6, 
                     maxiter = 5000, checkEvery = 50,
                     nfold = 10) {
    ## checkEvery is not used as such
    checkEvery <- maxiter
    if((!is.vector(time)) | (!is.vector(event)))
        stop("Time and event should be vectors")
    n <- length(time)
    thresS <- seq(from = thres[1], to = thres[2], length.out = thresGrid)
    cvpl.mat <- array(NA, c(length(thresS), maxiter))
    cvindex <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    clusterParams <- expand.grid(1:nfold, thresS)
    totalProcs <- dim(clusterParams)[1]

    xTheCluster <<- x
    tTheCluster <<- time
    eTheCluster <<- event
    epTheCluster <<- epi
    mitTheCluster <<- maxiter
    ceTheCluster <<- checkEvery
    cviTheCluster <<- cvindex
    clusterParams <<- clusterParams
    
    time.initial.export <-
        unix.time(clusterExport(TheCluster,
                                c("xTheCluster", "tTheCluster",
                                  "eTheCluster",
                                  "epTheCluster",
                                  "mitTheCluster","ceTheCluster",
                                  "cviTheCluster",
                                  "clusterParams")))
    cat("\n\n     Export to slaves took ", time.initial.export[3],
        " seconds\n")
    ## clean master too
    rm(list = c("xTheCluster", "tTheCluster", "eTheCluster",
       "epTheCluster","mitTheCluster", "ceTheCluster",
       "cviTheCluster", "clusterParams"), envir = .GlobalEnv)

    ts1 <- unix.time(clusterApply(TheCluster, 1:60, tgdSnowSetUp))
    cat("\n \n      Slaves setup took ", ts1[3], " seconds\n")
    
    clusterOutput <- list()
    t.r1 <-
        unix.time(
                  clusterOutput <-
                  clusterApply(TheCluster, 1:60,
                                 function(x) tgd1InternalSnow()))
    cat("\n \n      First run took ", t.r1[3], " seconds\n\n")

    tmp.cvpl <- matrix(unlist(clusterOutput),
                      ncol = checkEvery, byrow = TRUE)

##    browser()
    cvpllymtgd <- matrix(NA, nrow = thresGrid, ncol = maxiter)
    ## Set cvpl to NULL to make the stored
    ## and later transmitted object as small as possible
    for(th in 1:thresGrid) {
        cvpllymtgd[th, ] <-
            apply(tmp.cvpl[clusterParams[, 2] == thresS[th], ], 2, sum)
    }        

    cvpl.mat <- cvpllymtgd/nfold

    ## have to allow for possibility of several minima; take largest
    tmp <-
        which(cvpl.mat == min(cvpl.mat, na.rm = TRUE), arr.ind = TRUE)
    tmp <- tmp[nrow(tmp), ]
    thresBest <- thresS[tmp[1]]
    stepBest <- tmp[2]

    if(is.na(stepBest)) {
        caughtUserError("An error occured, probably related to numerical problems. Please let us know")
    }

    tgd.alldata <- tgdTrain(x, time, event, thresBest, epi,
                          steps = stepBest)

## To incorporate ROC, do it here:
    ## train best models for each, and then find ROC, with CV.

    betas <- tgd.alldata$beta
    return(list(betas = betas, nonZeroBetas = sum(betas != 0), 
                threshold = thresBest, epi = epi,
                step = stepBest, cvpl.mat = cvpl.mat,
                tgd.alldata = tgd.alldata,
                thres.loc = tmp[1]))
}

      
tgdSnowSetUp <- function(index) {
    assign("cvindex", cviTheCluster, env = .GlobalEnv)
    assign("epi",  epTheCluster, env = .GlobalEnv)
    assign("steps",  ceTheCluster, env = .GlobalEnv)
    assign("maxiter",  mitTheCluster, env = .GlobalEnv)

    assign("thres",  clusterParams[index, 2], env = .GlobalEnv)
    assign("i",  clusterParams[index, 1], env = .GlobalEnv)
    assign("x.train",
           as.matrix(xTheCluster[cvindex != i, , drop = FALSE]),
           env = .GlobalEnv)
    assign("time.train", tTheCluster[cvindex != i], env = .GlobalEnv)
    assign("event.train",  eTheCluster[cvindex != i], env = .GlobalEnv)
    
    assign("x.test",
           as.matrix(xTheCluster[cvindex == i, , drop = FALSE]), env = .GlobalEnv)
    assign("time.test",  tTheCluster[cvindex == i], env = .GlobalEnv)
    assign("event.test",  eTheCluster[cvindex == i], env = .GlobalEnv)
    
    rm(list = c("xTheCluster", "tTheCluster", "eTheCluster",
       "epTheCluster","mitTheCluster", "ceTheCluster",
       "cviTheCluster", "clusterParams", "i"), envir = .GlobalEnv)
}


tgd1InternalSnow <- function(){

## x, time, event: covariate, time to event, and event indicator (1
## observed, 0 censored) for the training data.
## xtest, time.test, event.test: same as above, but for testing data.
## thres and epi: oh well, thres and epi.
## steps: the maximum number of steps to do in this run.
## maxiter: total number of possibe iterations (sum(steps) <= maxiter)
    
    ## like tgd, but sets size of vectors to maxiter
    ## and allows returning a list with the complete state.
    pt1 <- proc.time()[3]

    
    nm <- dim(x.train)
    n <<- nm[1]
    m <<- nm[2]
    nm1 <- dim(x.test)
    n1 <<- nm1[1]
    m1 <<- nm1[2]
    r <<- rank(time.train)
    gradient <- rep(NA, n)
    scores <- rep(NA, n1)
    beta <- beta1 <- rep(0, m)
    beta <- as.matrix(beta)
    cvpl <- rep(0, steps)
    
    warning(paste("I am using: threshold", thres))

    for(iteration in 1:steps) {
        beta <- beta1
        if(n1 == 1) scores <- sum(x.test * beta)
        else scores <- x.test %*% beta

        ita <- x.train%*% beta
     
        epita <- exp(ita)
        d <- rep(0, n)
        dono <- rep(0, n)
    
        for(i in 1:n) {
            d[i] <- sum(event.train[r == r[i]])
            dono[i] <- sum(epita[r >= r[i]])
        }
        risk <- d/dono
        culrisk <- rep(0, n)
        for(i in 1:n) {
            culrisk[i] <- sum(unique(risk[r <= r[i]]))
        }
        
        gradient <- event.train - epita * culrisk
        gra1 <- crossprod(x.train, gradient)
        gra1[abs(gra1) < thres * max(abs(gra1))] <- 0
        beta1 <- beta + epi * gra1
        
        ## cvpl is a somewhat separate thing here.
        lik1.train <- sum((ita - log(dono)) * event.train)
        cvpl[iteration] <- (lik1.train -
                            lik1(c(ita, scores),
                                 c(time.train,  time.test),
                                 c(event.train, event.test)))/length(time.test)

    }
    rm("beta1GlobalEnv", envir = .GlobalEnv)
    assign("beta1GlobalEnv", beta1, envir = .GlobalEnv)

    ## recall we are using clusterApply
    ## and we count on the persistence for future visits

    warning(paste("**** First run: spent ", proc.time()[3] - pt1, "seconds"))
    return(cvpl)
}


tgdPieceInternalSnow <- function(runIt) {
### like tgd, but continues the iterations where the other left.
### we are using clusterApply, so we count on the exact positions in the cluster.

    pt1 <- proc.time()[3]
    if(runIt) {
        cvpl <- rep(0, steps)
        scores <- rep(0, n1)
        gradient <- rep(0, n)
        
        beta1 <- beta1GlobalEnv
##         warning(paste("I am using: threshold", thres,
##                       " runIt ", runIt))

        for(iteration in 1:steps){
            beta <- beta1
            if(n1 == 1) scores <- sum(x.test * beta)
            else scores <- x.test %*% beta
            
            ita <- x.train %*% beta
            
            epita <- exp(ita)
            d <- rep(0, n)
            dono <- rep(0, n)
            for(i in 1:n) {
                d[i] <- sum(event.train[r == r[i]])
                dono[i] <- sum(epita[r >= r[i]])
            }
            risk <- d/dono
            culrisk <- rep(0, n)
            for(i in 1:n) {
                culrisk[i] <- sum(unique(risk[r <= r[i]]))
            }
            gradient <- event.train - epita * culrisk
            gra1 <- crossprod(x.train, gradient)
            gra1[abs(gra1) < thres * max(abs(gra1))] <- 0
            beta1 <- beta + epi * gra1
            
            
            lik1.train <- sum((ita - log(dono)) * event.train)
            cvpl[iteration] <- (lik1.train -
                                lik1(c(ita, scores),
                                     c(time.train,  time.test),
                                     c(event.train, event.test)))/length(time.test)
        }
        rm("beta1GlobalEnv", envir = .GlobalEnv)
        assign("beta1GlobalEnv", beta1, envir = .GlobalEnv)
        warning(paste("**** Other runs: spent ", proc.time()[3] - pt1, "seconds"))
        return(cvpl)
    } else {return(rep(NA, steps))}
}


lik1 <- function(score, time, event)
{
### return the partial likelihood for Cox model 	
### like the original, but with sapply; slightly faster than with for.    
    n <- length(score)
    r <- rank(time)
    epita <- exp(score)
    dono <- rep(0, n)
    dono <- sapply(r, function(x) sum(epita[r >= x]))
    lik <- sum((score - log(dono)) * event)
    return(lik)
}

cvTGDP <- function(x, time, event, thres = c(0, 1),
                 epi = 5e-06, thresGrid = 6, 
                 maxiter = 5000, checkEvery = 50,
                 nfold = 10) {

    ##     if(is.null(nfold)) {
    ##         nfold <- length(y)
    ##     } ## I've never tested this. Should work for leave-one-out but ....
    
    ## we asume at least as many mpi Rslaves as nfold.
    if(nfold > (mpi.comm.size() - 1))
        stop("nfold > number of mpi Rslaves")
    
    
    ## nfold is the number of folds for cross-validation, and this are
    ## also the number of folds for the internal cv. We could have two
    ## different parameters but ...its a pain.
    
    ## all other parameters are also both for the internal and external
    ## runs.
   
    n <- length(time)
    index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    OOB.scores <- rep(NA, n)
    clusterResults <- list()
    for(i in 1:nfold) {
        cat("\n\n **** Doing fold ", i, "    ******\n")
        clusterResults[[i]] <- tgdCVPred(x[index.select != i, , drop = FALSE],
                                         time[index.select != i],
                                         event[index.select != i],
                                         x[index.select == i, , drop = FALSE],
                                         thres = thres,
                                         epi = epi, thresGrid = thresGrid,
                                         checkEvery = checkEvery,
                                         maxiter = maxiter,
                                         nfold = nfold)
        OOB.scores[index.select == i] <-
            clusterResults[[i]]$scoresTest
        clusterResults[[i]]$scoresTest <- NULL
    }
    
    out <- list(clusterResults = clusterResults,
                OOB.scores = OOB.scores)
    class(out) <- "cvTGD"
    return(out)
}


tgdTrain <- function(x, time, event, thres, epi, steps){
### R.D.-U.: from tgd but I eliminate all testing data related
###    stuff and I return gradient, scores (for plotting later)
###    and betas.    

    
    
### x: n*m predictor matrix, with m predictors.
### time: survival time for n observations.
### event: censoring event for n observations.
### x, time, event together are training dataset
### xtest, time.test, event.test together are corresponding testing dataset
### thres: threshold value between 0 and 1
### epi: maximum factor (step size) scaling gradient for incrementing selected coefficients at each step
### steps: maximum number of threshold gradient descent iterations
### program outputs:
### gra: n*steps, is the gradient for the training score (x%*%beta) in each steps.
### scores: n1*steps, is the score(pred%*%beta) for testing dataset in each steps.
### lik: testing data's partial likelihood in each steps.
### cvpl: will be used to calculate the cross-validated partial likelihood in "tgdcvpl".

    x <- as.matrix(x)
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    r <- rank(time)
    gradient <- matrix(0, n, steps)
    beta <- beta1 <- rep(0, m)
    beta <- as.matrix(beta)
    f <- 1
    while(f <= steps) {
        beta <- beta1
        ita <- x %*% beta
        epita <- exp(ita)
        d <- rep(0, n)
        dono <- rep(0, n)
        for(i in 1:n) {
            d[i] <- sum(event[r == r[i]])
            dono[i] <- sum(epita[r >= r[i]])
        }
        risk <- d/dono
        culrisk <- rep(0, n)
        for(i in 1:n) {
            culrisk[i] <- sum(unique(risk[r <= r[i]]))
        }
        gradient[, f] <- event - epita * culrisk
        gra1 <- crossprod(x, gradient[, f])
        gra1[abs(gra1) < thres * max(abs(gra1))] <- 0
        beta1 <- beta + epi * gra1
        f <- f + 1
    }
    return(list(gradient = gradient, beta = beta1, scores = x %*% beta1))
}


KM.visualize <- function(scores, surv, event, ngroups = 2,
                         addmain = "(Overfitt)") {
### Plots like in Dave et al., 2004, p. 2164
### will make it more general, to > 2 groups, later; zz
    
    cutting <- quantile(scores,
                        probs = seq(from = 0, to = 1, length.out = ngroups + 1),
                        type = 8,
                        na.rm = TRUE)

    group <- rep(NA, length(scores))
    group[scores < cutting[2]] <- "High"
    group[scores >= cutting[2]] <- "Low"
    group <- factor(group)

    if(length(unique(na.omit(group))) < 2) {
        plot(x = c(0, 1), y = c(0, 1), 
         type = "n", axes = FALSE, xlab = "", ylab = "")
    	 box()     
    	 text(0.5, 0.7, "There is only one group (and/or")
    	 text(0.5, 0.5,
    	 "only a null model was fitted).")
         text(0.5, 0.3, "Nothing to plot, therefore.")
	 return(NULL)
    } else {
        test <- survdiff(Surv(surv, event) ~ group)
        pv <- (1 - pchisq(test$chisq, df = 1))
        
        plot(survfit(Surv(surv, event) ~ group),
             main = paste(addmain, "Survival curves: low vs. high model scores"),
             log = FALSE,
             col = c("Blue", "Red"),
             legend.text = c("High score", "Low score"),
         legend.pos = c(0.7 * max(surv), 0.9), legend.bty = "o",
             cex = 1.5, lwd = 1.5)
        if(pv > 1e-3) pv <- round(pv, 4)
        else if (pv > 1e-5) pv <- round(pv, 5)
        else if (pv > 1e-7) pv <- round(pv, 7)
        else pv <- "< 1e-7"
        
        title(sub = paste("Log-rank test p-value ", pv))
    }
}



KM.visualize4 <- function(scores, surv, event, ngroups = 4,
                         addmain = "(Overfitt)") {
### Plots like in Dave et al., 2004, p. 2164
### will make it more general, to > 2 groups, later; zz
### Man, is this a terrible hack!!!
    ngroups <- 4 ## options is irrelevant
    cutting <- quantile(scores,
                        probs = seq(from = 0, to = 1, length.out = ngroups + 1),
                        type = 8,
                        na.rm = TRUE)

    group <- rep(NA, length(scores))
  
    group[scores < cutting[2]] <- "Q1"
    group[(scores >= cutting[2]) & (scores < cutting[3])] <- "Q2"
    group[(scores >= cutting[3]) & (scores < cutting[4])] <- "Q3"
    group[scores >= cutting[4]] <- "Q4"
    group <- factor(group)

   if(length(unique(na.omit(group))) < 2) {
        plot(x = c(0, 1), y = c(0, 1), 
         type = "n", axes = FALSE, xlab = "", ylab = "")
    	 box()     
    	 text(0.5, 0.7, "There is only one group (and/or")
    	 text(0.5, 0.5,
    	 "only a null model was fitted).")
         text(0.5, 0.3, "Nothing to plot, therefore.")
	 return(NULL)
    } else { 
        
        test <- survdiff(Surv(surv, event) ~ group)
        pv <- (1 - pchisq(test$chisq, df = 3))
        
        plot(survfit(Surv(surv, event) ~ group),
         main = paste(addmain, "Survival curves: comparing four quartiles"),
             log = FALSE,
             lty = c(2, 1, 2, 1), 
         col = c("Blue", "Brown", "Red", "Cyan"),
             legend.text = c("Q1", "Q2", "Q3", "Q4"),
             legend.pos = c(0.7 * max(surv), 0.9), legend.bty = "o",
             cex = 1.5, lwd = 1.5)
        if(pv > 1e-3) pv <- round(pv, 4)
        else if (pv > 1e-5) pv <- round(pv, 5)
        else if (pv > 1e-7) pv <- round(pv, 7)
        else pv <- "< 1e-7"
        
        title(sub = paste("Log-rank test p-value ", pv))
    }
}





KM.visualize3 <- function(scores, surv, event, ngroups = 3,
                         addmain = "(Overfitt)") {
### Plots like in Dave et al., 2004, p. 2164
### will make it more general, to > 2 groups, later; zz
### Man, is this a terrible hack!!!
    ngroups <- 3 ## options is irrelevant
    cutting <- quantile(scores,
                        probs = seq(from = 0, to = 1, length.out = ngroups + 1),
                        type = 8,
                        na.rm = TRUE)

    group <- rep(NA, length(scores))
    
    group[scores < cutting[2]] <- "Q1"
    group[(scores >= cutting[2]) & (scores < cutting[3])] <- "Q2"
    group[scores >= cutting[3]] <- "Q3"
    group <- factor(group)

   if(length(unique(na.omit(group))) < 2) {
        plot(x = c(0, 1), y = c(0, 1), 
         type = "n", axes = FALSE, xlab = "", ylab = "")
    	 box()     
    	 text(0.5, 0.7, "There is only one group (and/or")
    	 text(0.5, 0.5,
    	 "only a null model was fitted).")
         text(0.5, 0.3, "Nothing to plot, therefore.")
	 return(NULL)
    } else { 
        test <- survdiff(Surv(surv, event) ~ group)
        pv <- (1 - pchisq(test$chisq, df = 2))
        
        plot(survfit(Surv(surv, event) ~ group),
             main = paste(addmain, "Survival curves: comparing three terciles"),
             log = FALSE,
             lty = c(2, 1, 2), 
             col = c("Blue", "Brown", "Red"),
             legend.text = c("Q1", "Q2", "Q3"),
             legend.pos = c(0.7 * max(surv), 0.9), legend.bty = "o",
             cex = 1.5, lwd = 1.5)
        if(pv > 1e-3) pv <- round(pv, 4)
        else if (pv > 1e-5) pv <- round(pv, 5)
        else if (pv > 1e-7) pv <- round(pv, 7)
        else pv <- "< 1e-7"
        
        title(sub = paste("Log-rank test p-value ", pv))
    }
}




plot.cvpl <- function(cvpl.mat, epi, thres = c(0, 1), thresGrid = 6) {

## Diagnostic plot for the output from tauBest

    ## Would be neat to add \tau in the legend, but I am not able to.

    
    thresS <- seq(from = thres[1], to = thres[2],
                  length.out = thresGrid)
    
    m.step <- max(apply(cvpl.mat, 1, function(x) sum(!is.na(x))))
    nus <- epi * (1:m.step)
    maxy <- max(cvpl.mat, na.rm = TRUE)
    miny <- min(cvpl.mat, na.rm = TRUE)

    rainbow.col <- rainbow(thresGrid)

    tmp1 <- round(col2rgb(1:(thresGrid + 1))/255)
    plotcolors <- rgb(tmp1[1,], tmp1[2,], tmp1[3,])
    matplot(x = nus,
            y = t(cvpl.mat[, 1:m.step]), xlab = expression(nu),
            ylab = "Cross validated partial likelihood",
            type = "l", col = plotcolors, lwd = 1.9,
            main = expression(paste("CV partial likelihood: effects of ",
                tau, " and ", nu)))
    legend(x = 0.1 * (8* max(nus) + 2* min(nus)),
           y = 0.1 * (9 * maxy + 1 * miny),
           legend = as.character(thresS), 
           col = plotcolors,
           lty = 1:thresGrid, title = expression(tau),
           cex = 1.1)
}

summaryTGDrun <- function(x, time, event, z, epi, thres = c(0, 1),
                          thresGrid = 6, plot = TRUE,
                          genesOut = TRUE, outfile = "genes.all.out")  {
    ## spit out the selected betas (i.e., for chose threshold)
    ## but also provide a table such as Table 1 in Gui & Li, Pac.Symp...
    ## with CVPL added.

    ## Note that, except for the absolute best, we need to call tgdTrain
    ## again.  zz: this could be parallelized
    
    thresS <- seq(from = thres[1], to = thres[2],
                  length.out = thresGrid)

    mins.at <- apply(z$cvpl.mat, 1, function(x) which.min(x))
    mins <- apply(z$cvpl.mat, 1, function(x) min(x, na.rm = TRUE ))
    cvpls <- z$cvpl.mat[matrix(c(1:thresGrid, mins.at), ncol = 2)]

    ## why is bestBetas.m a list and not an array?
    ## I think because initially I'd return only non-zero
    ## betas.
    bestBetas.m <- list()
    ##bestBetas.m[[z$thres.loc]] <- z$betas
    thres.do <- (1:thresGrid)[-z$thres.loc]

    tgdTrainSnow <- function(varArgs, x, time, event, epi)
        return(tgdTrain(x, time, event, varArgs[1], epi, varArgs[2])[[2]])
    
    variableArgs <- list()
    for(tt in thres.do) variableArgs[[tt]] <- c(thresS[tt], mins.at[tt])
    bestBetas.m.pre <- clusterApply(TheCluster, variableArgs,
                                    tgdTrainSnow, x, time, event, epi)

    bestBetas.m <- bestBetas.m.pre
    bestBetas.m[[z$thres.loc]] <- z$betas
    ##     browser()
##     if(z$thres.loc == 1) bestBetas.m <- c(list(z$betas), bestBetas.m.pre)
##     else if (z$thres.loc == thresGrid) bestBetas.m <-
##         c(bestBetas.m.pre, list(z$betas))
##     else bestBetas.m <- c(bestBetas.m.pre[1:(z$thres.loc - 1)],
##                           list(z$betas),
##                           bestBetas.m.pre[(z$thres.loc):thresGrid])
   
    n.betas.nozero <- lapply(bestBetas.m, function(x) sum(x != 0))
    outm <- data.frame(cbind(Threshold = thresS,
                             Minimum.CVPL.at = mins.at,
                             CVPL = cvpls,
                             Number.non.zero.coefficients = n.betas.nozero))

    bb <- z$betas[z$betas != 0]
    bb <- bb[order(abs(bb), decreasing = TRUE)]
    cat("\n\n Selected genes (ordered by decreasing value of their coefficient)\n")
    bb
    cat("\n\n\n Cross-validated partial likelihood and number of selected\n")
    cat(" genes for different thresholds (tau), with delta nu (epi) =", epi, ".\n")
    cat("\n ===============================================================\n\n")
    print(outm)

    if(plot) {
        webPNG(file = "fstdgrun.png", width = 1.5 * png.width,
               height = 2 * png.height,
               pointsize = png.pointsize,
               family = png.family)
        par(cex.axis = 0.75); par(cex.lab = 1.4); par(cex.main = 1.5)

        par(mfrow = c(3, 2)) ## z: change for more general setups
        for(ip in 1:thresGrid) {
            plot(bestBetas.m[[ip]], xlab = "Gene index",
                 ylab = "Coefficient",
                 type = "l",
                 main = paste("Threshold (tau) =", thresS[ip]))
        }
        dev.off()

        pdf(file = "fstdgrun.pdf", width = 1.5 * png.width,
               height = 2* png.height, onefile = FALSE)
        par(mfrow = c(3, 2)) ## z: change for more general setups
        for(ip in 1:thresGrid) {
            plot(bestBetas.m[[ip]], xlab = "Gene index",
                 ylab = "Coefficient",
                 type = "l",
                 main = paste("Threshold (tau) =", thresS[ip]))
        }
        dev.off()
    }
    if(genesOut) {
        tmp <- as.data.frame(bestBetas.m)
        colnames(tmp) <- thresS
        ## do I get rownames??
        write.table(file = outfile, tmp,
                    quote = FALSE, sep = "\t")
    }
    return(outm)
}
   
tgdCVPred <- function(x, time, event,
                     xtest, 
                     thres, epi, thresGrid,
                     checkEvery, maxiter,
                     nfold) {
    ## find best params by CV and predict on a new set.
    bestTrain <- tauBestP(x, time, event,
                         thres, epi, thresGrid,
                         maxiter, checkEvery,
                         nfold)
    
    return(list(scoresTest = xtest %*% bestTrain$betas,
                betas = bestTrain$betas,
                nonZeroBetas = bestTrain$nonZeroBetas, 
                threshold = bestTrain$threshold,
                step = bestTrain$step))
}

summary.cvTGD <- function(object, allDataObject, subjectNames) {

    cat("\n Out-of-bag scores\n\n")
    oobs <- matrix(object$OOB.scores, ncol = 1)
    rownames(oobs) <- subjectNames
    print(oobs)

    object <- object[[1]] ## don't need scores anymore. Simpler subsetting.

    ks <- length(object)
    
    ngenes <- unlist(lapply(object, function(x) x$nonZeroBetas))
    thresholds <- unlist(lapply(object, function(x) x$threshold))
    steps <- unlist(lapply(object, function(x) x$step))


    tmp.mat <- data.frame(Number.selected.genes = ngenes,
                          Optimal.Threshold = thresholds,
                          Optimal.Steps = steps)

    cv.names <- paste("CV.run.", 1:ks, sep = "")    
    rownames(tmp.mat) <- cv.names
    cat("\n\n Number of selected genes and parameters in cross-validation runs\n")
    cat("-------------------------------------------------------------------\n\n")
    print(tmp.mat)

    cat("\n\n Stability assessments \n")
    cat(    " ---------------------\n")
    cat("\n Genes selected in each of the cross-validation runs \n")

 
    for(i in 1:ks) {
        cat(paste("CV run  ", i, " (", ngenes[i], " genes selected):   ", sep = ""), "\n")
        print(rownames(object[[i]]$betas)[object[[i]]$betas != 0])
        cat("\n---\n")
    }

    tmp.genesSelected <- list()
    tmp.genesSelected[[1]] <- rownames(allDataObject$betas)[allDataObject$betas != 0]
    genesSelected.cv <- lapply(object, function(x)
                               rownames(x$betas)[x$betas != 0])
    tmp.genesSelected <- c(tmp.genesSelected, genesSelected.cv)

    shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 1))
    for(i in 1:(ks + 1)) {
        for(j in 1:(ks + 1)) { ## sure, need not be symmetric, but this is fast
            shared.genes[i, j] <-
                length(intersect(tmp.genesSelected[[i]],
                                 tmp.genesSelected[[j]]))
        }
    }

    ngenes <- c(allDataObject$nonZeroBetas, ngenes)
    prop.shared <- round(shared.genes/ngenes, 3)
    
    ngenesS <- paste("(", ngenes, ")", sep = "")
    colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", cv.names)
    rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", cv.names), ngenesS)
    
    options(width = 200)
    cat("\n\n Number of shared genes \n")
    print(as.table(shared.genes))
    
    cat("\n\n Proportion of shared genes (relative to row total) \n")
    print(as.table(prop.shared))
    
    options(width = 80)
    unlisted.genes.selected <- unlist(genesSelected.cv)
    
    in.all.data <-
        which(names(table(unlisted.genes.selected, dnn = NULL)) %in% tmp.genesSelected[[1]])
    cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
    print(sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE))
    cat("\n")
    
    
    cat("\n\n Gene frequencies in cross-validated runs \n\n")
    tmp.table <- sort(table(unlisted.genes.selected, dnn = NULL),
                      decreasing = TRUE)
    print(tmp.table)
    cat("\n")
}


    

############################################################
############################################################
############################################################

##################         Dave et al.  ####################

############################################################
############################################################
############################################################


coxph.fit.pomelo0 <- function (x, y, init = NULL,
                              control, method = "efron",  rownames = NULL) {
    warnStatus <- 0
        naindex <- which(is.na(x))
    if(length(naindex)) {
        x <- x[-naindex]
        y <- y[naindex, ]
    }
    x <- as.matrix(x) ## this ain't very efficient
    n <- nrow(y)
    if (is.matrix(x)) 
        nvar <- ncol(x)
    else if (length(x) == 0) 
        nvar <- 0
    else nvar <- 1
    time <- y[, 1]
    status <- y[, 2]
    sorted <- order(time)
    newstrat <- as.integer(rep(0, n))
    
    offset <- rep(0, n)
    weights <- rep(1, n)

    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])
    if (nvar == 0) {
        x <- as.matrix(rep(1, n))
        nullmodel <- TRUE
        nvar <- 1
        init <- 0
        maxiter <- 0
    }
    else {
        nullmodel <- FALSE
        maxiter <- control$iter.max
        if (!missing(init) && !is.null(init)) {
            if (length(init) != nvar) 
                stop("Wrong length for inital values")
        }
        else init <- rep(0, nvar)
    }
    coxfit <- .C("coxfit2", iter = as.integer(maxiter), as.integer(n), 
        as.integer(nvar), stime, sstat, x = x[sorted, ], as.double(offset[sorted] - 
            mean(offset)), as.double(weights), newstrat, means = double(nvar), 
        coef = as.double(init), u = double(nvar), imat = double(nvar * 
            nvar), loglik = double(2), flag = integer(1), double(2 * 
            n + 2 * nvar * nvar + 3 * nvar), as.double(control$eps), 
        as.double(control$toler.chol), sctest = as.double(method == 
            "efron"), PACKAGE = "survival")
    if (nullmodel) {
        score <- exp(offset[sorted])
        coxres <- .C("coxmart", as.integer(n), as.integer(method == 
            "efron"), stime, sstat, newstrat, as.double(score), 
            as.double(weights), resid = double(n), PACKAGE = "survival")
        resid <- double(n)
        resid[sorted] <- coxres$resid
        names(resid) <- rownames
        list(loglik = coxfit$loglik[1], linear.predictors = offset, 
            residuals = resid, method = c("coxph.null", "coxph"))
    }
    else {
        var <- matrix(coxfit$imat, nvar, nvar)
        coef <- coxfit$coef
        if (coxfit$flag < nvar) 
            which.sing <- diag(var) == 0
        else which.sing <- rep(FALSE, nvar)
        infs <- abs(coxfit$u %*% var)
        if (maxiter > 1) {
            if (coxfit$flag == 1000) { 
                warning("Ran out of iterations and did not converge")
                warnStatus <- 2
            }
            else {
                infs <- ((infs > control$eps) & infs > control$toler.inf * 
                  abs(coef))
                if (any(infs)) {
                  warning(paste("Loglik converged before variable ", 
                    paste((1:nvar)[infs], collapse = ","), "; beta may be infinite. "))
                  warnStatus <- 1
              }
            }
        }
        names(coef) <- dimnames(x)[[2]]
        lp <- c(x %*% coef) + offset - sum(coef * coxfit$means)
        score <- exp(lp[sorted])
        coef[which.sing] <- NA
        list(coefficients = coef, var = var, warnStatus = warnStatus,
             loglik = coxfit$loglik)
    }
}

dStep1.parallel <- function(x, time, event, p, MaxIterationsCox) { 
    res.mat <- matrix(NA, nrow = ncol(x), ncol = 6)
    sobject <- Surv(time,event)

    funpap3 <- function (x) {
        out1 <- coxph.fit.pomelo0(x, sobject, control = coxph.control(iter.max = 500))
        if(out1$warnStatus > 1) {
            return(c(0, NA,  out1$warnStatus))
        } else {
            sts <- out1$coef/sqrt(out1$var)
            return(c(out1$coef,
                     1- pchisq((sts^2), df = 1), 
                     out1$warnStatus))
        }
    }

    tmp <- matrix(unlist(papply(as.data.frame(x),
                             funpap3,
                             papply_commondata =list(sobject = sobject,
                             MaxIterationsCox = MaxIterationsCox))),
                             ncol = 3, byrow = TRUE)
    res.mat[, 1:2] <- tmp[, 1:2]
    res.mat[, 3] <- ifelse(res.mat[, 2] < p, 1, 0)
    res.mat[, 4] <- sign(res.mat[, 1]) * res.mat[, 3]
    res.mat[, 5] <- tmp[, 3]
    res.mat[, 6] <- p.adjust(tmp[, 2], method = "BH")
    res.mat[is.na(res.mat[, 2]), c(2, 6)] <- 999
    colnames(res.mat) <- c("coeff", "p.value", "keep", "pos.neg", "Warning", "FDR")
    return(res.mat)
}


# dStep1.parallel <- function(x, time, event, p, MaxIterationsCox) { 
#     ## Painfully slow!!!! use C++ code 
#     ## We use papply to parallelize; but what will happen when slaves run this???
#     res.mat <- matrix(NA, nrow = ncol(x), ncol = 4)
#     sobject <- Surv(time,event)

#     funpap3 <- function(x) {
#         rcph <- try(coxph(sobject ~ x,
#                           control = coxph.control(iter.max = MaxIterationsCox)))
#         if(class(rcph) == "try-error")
#             return(c(0, 99)) ## so, if error, an impossible value for p-value
#         else
#             return(c(rcph$coef, 1 - pchisq(rcph$wald, df = 1)))
#     }

#     res.mat[, 1:2] <-
#         matrix(unlist(papply(as.data.frame(x),
#                              funpap3,
#                              papply_commondata =list(sobject = sobject,
#                              MaxIterationsCox = MaxIterationsCox))),
#                              ncol = 2, byrow = TRUE)
    
#     res.mat[, 3] <- ifelse(res.mat[, 2] < p, 1, 0)
#     res.mat[, 4] <- sign(res.mat[, 1]) * res.mat[, 3]
    
#     colnames(res.mat) <- c("coeff", "p.value", "keep", "pos.neg")
#     return(res.mat)
# }


##  2. Cluster; independently for those with pos and neg beta.
dStep2 <- function(x, res.mat, maxSize, minSize,
                   minCor, plot,
                   interactive,
                   plotSizes = c(0.5, 1, 2)) {
##     basicSizeHeight <- 9
##     basicSizeWidth <-  19
    pdfheight <- 19
    pdfwidth <- 12.7
    height <- 1200
    width <- 800
    ps <- 12

    res.mat[is.na(res.mat[, 4]), 4] <- 0

    if(sum(res.mat[, 4] == 1) >= minSize) {
        pos.data <- x[, res.mat[, 4] == 1]
        pdok <- TRUE
    } else {
        pdok <- FALSE
        warning(paste("Not enough positive coeff. genes that",
                      "meet the p restrictions."))
    }
    if(sum(res.mat[, 4] == -1) >= minSize) {
        neg.data <- x[, res.mat[, 4] == -1]
        pnok <- TRUE
    } else {
        pnok <- FALSE
        warning(paste("Not enough negative coeff. genes that",
                      "meet the p restrictions."))
    }

    if((!pnok) & (!pdok)) {
	if(interactive) caughtUserError("No gene was above the minimal p threshold")
        else {
            return(NA)
        }
    }

    tn <- tp <- FALSE
    if(pdok) {
        pos.clus <- hclust(as.dist(1 -cor(pos.data)), method = "complete")
        pos.groups <- paste("P.", cutree(pos.clus, h = 1- minCor), sep = "")
        tpos <- table(pos.groups)
        pos.accept <- names(which((tpos >= minSize) & (tpos <= maxSize)))
        if(length(pos.accept)) {
            groupsPositive <- pos.groups[pos.groups %in% pos.accept]
            dataPositive <- pos.data[ , pos.groups %in% pos.accept]
            posGroups <- unique(groupsPositive)
            posMeanData <- matrix(NA, nrow = dim(dataPositive)[1],
                              ncol = length(posGroups))

            for(i in 1:length(posGroups)) {
                posMeanData[, i] <-
                    apply(dataPositive[, groupsPositive == posGroups[i]],
                          1, mean)
            }
            colnames(posMeanData) <- posGroups
            tp <- TRUE
        } else {
            tp <- FALSE
            warning(paste("No groups of positive coeff genes that",
                          "meet the p, minimum correlation and size restrictions."))
        }
    }

    if(pnok) {
        neg.clus <- hclust(as.dist(1 -cor(neg.data)), method = "complete")
        neg.groups <- paste("N.", cutree(neg.clus, h = 1- minCor), sep = "")
        npos <- table(neg.groups)
        neg.accept <- names(which((npos >= minSize) & (npos <= maxSize)))
        if(length(neg.accept)) {
            groupsNegative <- neg.groups[neg.groups %in% neg.accept]
            dataNegative <- neg.data[ , neg.groups %in% neg.accept]
            negGroups <- unique(groupsNegative)
            negMeanData <- matrix(NA, nrow = dim(dataNegative)[1],
                                  ncol = length(negGroups))

            for(i in 1:length(negGroups)) {
                negMeanData[, i] <-
                    apply(dataNegative[, groupsNegative == negGroups[i]],
                          1, mean)
            }
            colnames(negMeanData) <- negGroups
            tn <- TRUE
        } else {
            tn <- FALSE
            warning(paste("No groups of negative coeff. genes that",
                          "meet the p, minimum correlation and size restrictions."))
        }
    }

    if(!tn & !tp) {
        if(interactive) caughtUserError(paste("No groups that meet the p, minimum correlation",
                                              "and size restrictions."))
        else return(NA)
    }
    pdok <- tp & pdok
    pnok <- tn & pnok
    
    if(plot) {
       
        pdokf <- function(alllabels = FALSE) {
            pos.labels <- rep("                   ", ncol(pos.data))
            index.labels <- which(pos.groups %in% pos.accept)
            pos.labels[index.labels] <- colnames(pos.data)[index.labels]

            pos.dend <- as.dendrogram(pos.clus, hang = 0.001)
            par(mar = c(5, 2, 1, 8))
            plot(pos.dend, horiz = TRUE, xlab = "1 - correlation", leaflab = "none",
                 main = "Positive coefficients")
            abline(v = 1 - minCor, lty = 2, col = "blue")

            dfp <- data.frame(name = unlist(dendrapply(rev(pos.dend), nameLeave)),
                              height = unlist(dendrapply(rev(pos.dend), heightLeave)))
            dfp$pos.gr <- pos.groups[order.dendrogram(rev(pos.dend))]
            dfp$chosen.clus <- dfp$pos.gr %in% pos.accept
            dfp$y <- nrow(dfp) - as.numeric(rownames(dfp)) + 1
            rainbow.col <- rainbow(length(posGroups))
            for(i in 1:length(posGroups)) {
                dfpt <- dfp[dfp$pos.gr == posGroups[i], ]
                miny <- min(dfpt$y)
                maxy <- max(dfpt$y)
                axis(4, line = 5, at = c(miny, maxy), col = rainbow.col[i],
                     tick = TRUE, labels = FALSE, lw = 3)
                axis(4, line = 5, at = 0.5 * (miny + maxy),
                     col.axis = rainbow.col[i],
                     tick = FALSE, labels = posGroups[i], lw = 0,
                     cex.axis = 1.5)
                for(j in 1:nrow(dfpt))
                    text(dfpt$name[j], x = dfpt$height[j],
                         y = dfpt$y[j], col = rainbow.col[i],
                         font = 2, pos = 4)
            }
            if (alllabels) {
                dfpt <- dfp[dfp$chosen.clus == FALSE,]
		if(nrow(dfpt)) {
                for(j in 1:nrow(dfpt))
                    text(dfpt$name[j], x = dfpt$height[j],
                         y = dfpt$y[j], col = "black", cex = 0.8, pos = 4)
			 }
            }
            return(dfp)
        } ##</pdokf within plotting>

        dendmapp <- function(factor = 1, alllabels = FALSE) {
            psf <- ifelse(factor < 1, ps * factor, ps)
            nameIm <- paste("dend.P.factor", factor, ".alllabels", alllabels, sep = "")
            im1 <- imagemap3(nameIm, height = height * factor, width = width * factor,
                             ps = psf)
            dfp <- pdokf(alllabels = alllabels)
            
            for(np in 1:nrow(dfp)) {
                addRegion(im1) <- imRect(dfp[np, 2] + .030, dfp[np, 5] - 0.45,
                                         dfp[np, 2] - .1, dfp[np, 5] + 0.45,
                                         title = dfp[np, 1], alt = dfp[np, 1],
                                         href= linkGene2(dfp[np, 1]))
            }
            
            createIM(im1, file = paste("dend.P.factor", factor, ".alllabels",
                          alllabels, ".html", sep = ""))
            imClose(im1)
        }
        
            
        if (pdok) {
            for (plsz in plotSizes) {
                dendmapp(factor = plsz, alllabels = FALSE)
##                 pdf(file = paste("dend.P.factor", plsz, ".alllabels=FALSE",
##                     ".pdf", sep = ""),
##                     height = pdfheight * plsz,
##                     width = pdfwidth * plsz)
##                 pdokf(alllabels = FALSE)
##                dev.off()
            }
             for (plsz in plotSizes) {
                 dendmapp(factor = plsz, alllabels = TRUE)
##                  pdf(file = paste("dend.P.factor", plsz, ".alllabels = TRUE",
##                      ".pdf", sep = ""),
##                      height = pdfheight * plsz,
##                      width = pdfwidth * plsz)
##                  pdokf(alllabels = TRUE)
##                 dev.off()
             }
            
        } else {
            system("touch NoPositiveCluster")
##             for (plsz in plotSizes) {
##                 GDD(file = paste("dend.P.factor", plsz, ".alllabelsFALSE",
##                     ".png", sep = ""), w=width, h=height,
##                     type = "png", ps = ps)	  
##                 plot(x = c(0, 1), y = c(0, 1),
##                      type = "n", axes = FALSE, xlab = "", ylab = "")
##                 box()
##                 text(0.5, 0.7,
##                      "There are no genes with positive coefficients")
##                 text(0.5, 0.5, "that satisfy the p, minimum correlation")
##                 text(0.5, 0.3, "and size restrictions.")
##                 dev.off()
##             }
##             for (plsz in plotSizes) {
##                 GDD(file = paste("dend.P.factor", plsz, ".alllabelsTRUE",
##                     ".png", sep = ""), w=width, h=height,
##                     type = "png", ps = ps)	  
##                 plot(x = c(0, 1), y = c(0, 1),
##                      type = "n", axes = FALSE, xlab = "", ylab = "")
##                 box()
##                 text(0.5, 0.7,
##                      "There are no genes with positive coefficients")
##                 text(0.5, 0.5, "that satisfy the p, minimum correlation")
##                 text(0.5, 0.3, "and size restrictions.")
##                 dev.off()
##             }
        }


        pnokf <- function(alllabels = FALSE) {
            neg.labels <- rep("                   ", ncol(neg.data))
            index.labels <- which(neg.groups %in% neg.accept)
            neg.labels[index.labels] <- colnames(neg.data)[index.labels]

            neg.dend <- as.dendrogram(neg.clus, hang = 0.001)
            par(mar = c(5, 2, 1, 8))
            plot(neg.dend, horiz = TRUE, xlab = "1 - correlation", leaflab = "none",
                 main = "Negative coefficients")
            abline(v = 1 - minCor, lty = 2, col = "blue")

            dfp <- data.frame(name = unlist(dendrapply(rev(neg.dend), nameLeave)),
                              height = unlist(dendrapply(rev(neg.dend), heightLeave)))
            dfp$neg.gr <- neg.groups[order.dendrogram(rev(neg.dend))]
            dfp$chosen.clus <- dfp$neg.gr %in% neg.accept
            dfp$y <- nrow(dfp) - as.numeric(rownames(dfp)) + 1
            rainbow.col <- rainbow(length(negGroups))
            for(i in 1:length(negGroups)) {
                dfpt <- dfp[dfp$neg.gr == negGroups[i], ]
                miny <- min(dfpt$y)
                maxy <- max(dfpt$y)
                axis(4, line = 5, at = c(miny, maxy), col = rainbow.col[i],
                     tick = TRUE, labels = FALSE, lw = 3)
                axis(4, line = 5, at = 0.5 * (miny + maxy),
                     col.axis = rainbow.col[i],
                     tick = FALSE, labels = negGroups[i], lw = 0,
                     cex.axis = 1.5)
                for(j in 1:nrow(dfpt))
                    text(dfpt$name[j], x = dfpt$height[j],
                         y = dfpt$y[j], col = rainbow.col[i],
                         font = 2, pos = 4)
            }
            if (alllabels) {
                dfpt <- dfp[dfp$chosen.clus == FALSE,]
		if(nrow(dfpt)) {
                for(j in 1:nrow(dfpt))
                    text(dfpt$name[j], x = dfpt$height[j],
                         y = dfpt$y[j], col = "black", cex = 0.8, pos = 4)
			 }
            }
            return(dfp)
        } ##</pnokf within plotting>

        dendmapp <- function(factor = 1, alllabels = FALSE) {
            psf <- ifelse(factor < 1, ps * factor, ps)
            nameIm <- paste("dend.N.factor", factor, ".alllabels", alllabels, sep = "")
            im1 <- imagemap3(nameIm, height = height * factor, width = width * factor,
                             ps = psf)
            dfp <- pnokf(alllabels = alllabels)
            
            for(np in 1:nrow(dfp)) {
                addRegion(im1) <- imRect(dfp[np, 2] + .030, dfp[np, 5] - 0.45,
                                         dfp[np, 2] - .1, dfp[np, 5] + 0.45,
                                         title = dfp[np, 1], alt = dfp[np, 1],
                                         href= linkGene2(dfp[np, 1]))
            }
            
            createIM(im1, file = paste("dend.N.factor", factor, ".alllabels",
                          alllabels, ".html", sep = ""))
            imClose(im1)
        }

        if (pnok) {
            for (plsz in plotSizes) {
                dendmapp(factor = plsz, alllabels = FALSE)
##                 pdf(file = paste("dend.N.factor", plsz, ".alllabels=FALSE",
##                     ".pdf", sep = ""),
##                     height = pdfheight * plsz,
##                     width = pdfwidth * plsz)
##                 pnokf(alllabels = FALSE)
##                 dev.off()
            }
            for (plsz in plotSizes) {
                dendmapp(factor = plsz, alllabels = TRUE)
##                 pdf(file = paste("dend.N.factor", plsz, ".alllabels=TRUE",
##                     ".pdf", sep = ""),
##                     height = pdfheight * plsz,
##                     width = pdfwidth * plsz)
##                 pnokf(alllabels = TRUE)
##                 dev.off()
            }
        } else {
            system("touch NoNegativeCluster")
##             for (plsz in plotSizes) {
##                 GDD(file = paste("dend.N.factor", plsz, ".alllabelsFALSE",
##                     ".png", sep = ""), w=width, h=height,
##                     type = "png", ps = ps)	  
##                 plot(x = c(0, 1), y = c(0, 1),
##                      type = "n", axes = FALSE, xlab = "", ylab = "")
##                 box()
##                 text(0.5, 0.7,
##                      "There are no genes with negative coefficients")
##                 text(0.5, 0.5, "that satisfy the p, minimum correlation")
##                 text(0.5, 0.3, "and size restrictions.")
##                 dev.off()
##             }
##             for (plsz in plotSizes) {
##                 GDD(file = paste("dend.N.factor", plsz, ".alllabelsTRUE",
##                     ".png", sep = ""), w=width, h=height,
##                     type = "png", ps = ps)	  
##                 plot(x = c(0, 1), y = c(0, 1),
##                      type = "n", axes = FALSE, xlab = "", ylab = "")
##                 box()
##                 text(0.5, 0.7,
##                      "There are no genes with negative coefficients")
##                 text(0.5, 0.5, "that satisfy the p, minimum correlation")
##                 text(0.5, 0.3, "and size restrictions.")
##                 dev.off()
##             }            
        }
        
    } ##</ if plot>

    ## For predictions and results, which vars. correspond
    ## to which original genes

    ## The following two are the indices of the columns in the original (non
    ## divided data file)

    if(pdok) {
        posPositions <- which(res.mat[, 4] == 1)
        filteredPosPositions <- posPositions[pos.groups %in% pos.accept]
    }
    if(pnok) {
        negPositions <- which(res.mat[, 4] == -1)
        filteredNegPositions <- negPositions[neg.groups %in% neg.accept]
    }

    if(pnok & pdok)
        if(length(intersect(filteredPosPositions, filteredNegPositions)))
            stop("Non zero intersection between filtered pos and neg positions")
    if(!pnok) {
        negMeanData <- NULL
        groupsNegative <- NA
        filteredGroupsNegative <- NA
        filteredNegPositions <- NA
        negPositions <- NA
    }
    if(!pdok) {
        posMeanData <- NULL
        groupsPositive <- NA
        filteredGroupsPositive <- NA
        filteredPosPositions <- NA
        posPositions <- NA
    }
    
    return(list(md = cbind(posMeanData, negMeanData),
                filteredGroupsPositive = groupsPositive,
                filteredGroupsNegative = groupsNegative,
                filteredPosPositions = filteredPosPositions,
                filteredNegPositions = filteredNegPositions,
                posPositions = posPositions,
                negPositions = negPositions))
}


## 3. fit model

dStep3 <- function(res2, time, event, MaxIterationsCox) {
    if(all(is.na(res2))) return(NA)
    md <- res2$md
    sobject <- Surv(time, event)
    
    if(ncol(md) >= 2) {
        
        ## Select best two-genes model
        ncolmd <- ncol(md)
        modelsSizeTwo <- t(combn(1:ncolmd, 2))
        logliks <- rep(NA, nrow(modelsSizeTwo))
        for(i in 1:nrow(modelsSizeTwo)) {
            ## This is the logic: in some models, if we use a very large num.
            ## iterations, some scores end up being infinite. That sometimes
            ## can be "fixed" if we use a smaller number of iterations.
            ## If that does not work, then we set the value to NA
            trycox <- try(logliks[i] <-
                          coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                control = coxph.control(iter.max = MaxIterationsCox))$loglik[2])
            
            if(class(trycox) == "try-error")
                trycox <- try(logliks[i] <-
                              coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                    control = coxph.control(iter.max = 20))$loglik[2]) 
            
            if(class(trycox) == "try-error")
                trycox <- try(logliks[i] <-
                              coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                    control = coxph.control(iter.max = 10))$loglik[2]) 
            
            if(class(trycox) == "try-error")
                logliks[i] <- NA
        }
        ## zz: give warnings if any loglik is NA
        bestTwoGenes <- modelsSizeTwo[which.max(logliks), ]
        ## stepAIC
        mdf <- data.frame(md)
        attach(mdf) ## stepAIC not working otherwise
        trycox <- try(
                      bestTwoModel <-
                      coxph(eval(parse(text = paste("sobject ~",
                                       paste(colnames(mdf)[bestTwoGenes],
                                             sep = "", collapse = " + ")))),
                            control = coxph.control(iter.max = MaxIterationsCox))
                      )
        if(class(trycox) == "try-error")
            trycox <- try(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 20))
                          )
        if(class(trycox) == "try-error")
            trycox <- try(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 10))
                          )
        if(class(trycox) == "try-error")
            trycox <- try(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 5))
                          )
        if(class(trycox) == "try-error")
            trycox <- try(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 2))
                          )
        ## So we'll try stepwise; but if it doesn't work,
        ## we use forwards, with max size model  = num.variables.
        ## We could run into trouble, so we allow to decrease size
        ## of model, until we are back to bestTwoModel
        upperScope <- paste("~", paste(colnames(mdf), sep = "",
                                       collapse = " + "))
        tryaic <- try(
                      finalModel <- stepAIC(bestTwoModel, scope = upperScope,
                                            direction = "both")
                      )
        maxsize <- sum(event) - 2
        while(class(tryaic) == "try-error") {
            tryaic <-
                try(finalModel <- stepAIC(bestTwoModel,
                                          scope = upperScope,
                                          direction = "forward",
                                          steps = maxsize)
                    )
            maxsize <- maxsize - 1
            if(maxsize == 1) {
                finalModel <- bestTwoModel
                break
            }
        }
        
        
        ## Should we always be able to get something here? I don't
        ## think so.
        predictsFinalModel <- predict(finalModel, type = "lp")
        detach(mdf)
    } else {
        mdf <- data.frame(md)

        trycox <-
          try(
              finalModel <- coxph(sobject ~ ., data = mdf,
                                  control = coxph.control(iter.max = MaxIterationsCox)))
        if(class(trycox) == "try-error")
            trycox <-
            try(
                finalModel <- coxph(sobject ~ ., data = mdf,
                                    control = coxph.control(iter.max = 20)))
        if(class(trycox) == "try-error")
          trycox <-
            try(
                finalModel <- coxph(sobject ~ ., data = mdf,
                                    control = coxph.control(iter.max = 10)))
        if(class(trycox) == "try-error")
          trycox <-
            try(
                finalModel <- coxph(sobject ~ ., data = mdf,
                                    control = coxph.control(iter.max = 5)))
        if(class(trycox) == "try-error")
          trycox <-
            try(
                finalModel <- coxph(sobject ~ ., data = mdf,
                                    control = coxph.control(iter.max = 2)))
        
        if(class(trycox) == "try-error")
          predictsFinalModel <- rep(NA, length(time))
        else
          predictsFinalModel <- predict(finalModel, type = "lp")
    }
    
    out <- list(model = finalModel, scores = predictsFinalModel,
                clusterResults = res2)
    class(out) <- "fmDave"
    return(out)
}

dPredictNew <- function(res3, newdata) { ## returns predictions (scores)y
    ## given a model, get the newdata
    ## z is a list, like returned from dStep2
    if(all(is.na(res3))) return(NA)

    model <- res3$model
    
    if(! ("coefficients" %in% names(model))) 
	return(rep(NA, dim(newdata)[1]))
	## this is the Null
 	## model, only intercept
	
    z <- res3$clusterResults
    
    dataPositive <- newdata[, z$filteredPosPositions, drop = FALSE]
    if(ncol(dataPositive)) {
        posGroups <- unique(z$filteredGroupsPositive)
        posMeanData <- matrix(NA, nrow = dim(dataPositive)[1],
                              ncol = length(posGroups))
        for(i in 1:length(posGroups)) {
            posMeanData[, i] <-
                apply(dataPositive[, z$filteredGroupsPositive == posGroups[i], drop = FALSE],
                      1, mean)
        }
        colnames(posMeanData) <- posGroups
    } else {
        posMeanData <- NA
    }

    dataNegative <- newdata[, z$filteredNegPositions, drop = FALSE]
    if(ncol(dataNegative)) {
        negGroups <- unique(z$filteredGroupsNegative)
        negMeanData <- matrix(NA, nrow = dim(dataNegative)[1],
                              ncol = length(negGroups))
        for(i in 1:length(negGroups)) {
            negMeanData[, i] <-
                apply(dataNegative[, z$filteredGroupsNegative == negGroups[i], drop = FALSE],
                      1, mean)
        }
        colnames(negMeanData) <- negGroups
    } else {
        negMeanData <- NA
    }
    return(predict(model,
                   newdata = as.data.frame(cbind(posMeanData, negMeanData)),
                   type = "lp"))
}
    
fitDave.res1Given <- function(x, time, event, res1,
                              p, maxSize, 
                              minSize, minCor, MaxIterationsCox, plot,
                              interactive) {
    res2 <- dStep2(x, res1, maxSize, minSize, minCor, plot, interactive)
    res3 <- dStep3(res2, time, event, MaxIterationsCox)
    return(res3) ## i.e., an fmDave object returned
}    

DaveCVPred.res1Given <- function(x, time, event, res1,
                                 train.index, test.index,
                                 p, maxSize, 
                                 minSize, minCor, MaxIterationsCox,
                                 plot = FALSE,
                                 interactive = FALSE) {#OUT:scoresTest +
                                                       #fmDaveObject; to
    
    xtrain <- x[train.index, , drop = FALSE]
    xtest <- x[test.index, , drop = FALSE]
    timetrain <- time[train.index]
    eventtrain <- event[train.index]


    ## find best params by CV and predict on a new set.
    bestTrain <- fitDave.res1Given(xtrain, timetrain, eventtrain,
                                   res1 = res1,
                                   p = p, maxSize = maxSize,
                                   minSize = minSize, minCor = minCor,
                                   MaxIterationsCox = MaxIterationsCox,
                                   plot = FALSE,
                                   interactive =FALSE)

    testPred <- dPredictNew(res3 = bestTrain, newdata = xtest)
                
    return(list(scoresTest = testPred,
                fmDaveObject = bestTrain))
}

DaveCVPred.res1Given.InternalMPI <- function() {
## the following need to be passed with
##    mpi.bcast.Robj2slave
##    x, time, event, ,
##    p, maxSize, 
##    minSize, minCor, MaxIterationsCox,
##    index.select
##    for res1 supply the complete list

    res1 <- res1s[[foldNumber]]
    
    xtrain <- x[index.select != foldNumber, , drop = FALSE]
    xtest <- x[index.select == foldNumber, , drop = FALSE]
    timetrain <- time[index.select != foldNumber]
    eventtrain <- event[index.select != foldNumber]

    ## find best params by CV and predict on a new set.
    bestTrain <- fitDave.res1Given(xtrain, timetrain, eventtrain,
                                   res1 = res1,
                                   p = p, maxSize = maxSize,
                                   minSize = minSize, minCor = minCor,
                                   MaxIterationsCox = MaxIterationsCox,
                                   plot = FALSE,
                                   interactive =FALSE)

    testPred <- dPredictNew(res3 = bestTrain, newdata = xtest)
                
    return(list(scoresTest = testPred,
                fmDaveObject = bestTrain))
}

mpiSpawnAll <- function() {
    mpi.spawn.Rslaves(nslaves = mpi.universe.size())
    mpiMyCleanSetup()
}

mpiSpawnThis <- function(hosts) {
    mpi.spawn.Rslaves(nslaves = mpi.universe.size(), hosts = hosts)
    mpiMyCleanSetup()
}

mpiMyCleanSetup <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
    mpi.remote.exec(library(survival))
    mpi.remote.exec(library(combinat))
    mpi.remote.exec(library(SignS2))
    mpi.remote.exec(library(MASS))
}



DaveCVPred.res1Given.InternalMPI2 <- function(fnum) {
## to be used with papply
  ## the following need to be passed with
##    mpi.bcast.Robj2slave
##    x, time, event, ,
##    p, maxSize, 
##    minSize, minCor, MaxIterationsCox,
##    index.select
##    for res1 supply the complete list

    res1 <- res1s[[fnum]]
    
    xtrain <- x[index.select != fnum, , drop = FALSE]
    xtest <- x[index.select == fnum, , drop = FALSE]
    timetrain <- time[index.select != fnum]
    eventtrain <- event[index.select != fnum]

    ## find best params by CV and predict on a new set.
    bestTrain <- fitDave.res1Given(xtrain, timetrain, eventtrain,
                                   res1 = res1,
                                   p = p, maxSize = maxSize,
                                   minSize = minSize, minCor = minCor,
                                   MaxIterationsCox = MaxIterationsCox,
                                   plot = FALSE,
                                   interactive =FALSE)

    testPred <- dPredictNew(res3 = bestTrain, newdata = xtest)
                
return(list(scoresTest = testPred,
                fmDaveObject = bestTrain))
}



cvDave.parallel3 <- function(x, time, event,
                             p, maxSize,
                             minSize, minCor,
                             MaxIterationsCox,
                             nfold, mpiHosts) {

    if (mpi.comm.size(comm = 1) == 0) {
        mpiSpawnAll()
    } else { ## so mpi is running
        if ((mpi.comm.size(comm = 1) - 1) < mpi.universe.size()) {
            ## but few salves
            mpi.close.Rslaves()
            mpiSpawnAll()
        } else {
            mpiMyCleanSetup()
        }
    }
       
## we asume at least as many mpi Rslaves as nfold.
    ## mpiHosts used for the second spawning of slaves,
    ## only for the 10 fold runs.
    
    
    n <- length(time)
    index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    OOB.scores <- rep(NA, n)

    res1s <- list()

    ### Here we use papply, which is 
    argspapp <- list()
    cat("\n\n Computing gene-wise cox p-value\n")
    for(i in 1:nfold) {
        cat("\n  ....  fold ", i)
        xtr <- x[index.select != i, , drop = FALSE]
        ttr <- time[index.select != i]
        etr <- event[index.select != i]
        res1s[[i]] <- dStep1.parallel(xtr, ttr, etr, p = p,
                                      MaxIterationsCox = MaxIterationsCox)
    }

##    cat("\n\n Cleaning up MPI space, and setting up a new one\n\n")
    ## Clean previous stuff
##    mpi.close.Rslaves()
##    mpiSpawnAll()
    ## mpiSpawnThis(hosts = mpiHosts)
    cat("\n\n Cleaning up MPI space, and setting up again\n\n")	
    mpiMyCleanSetup()

    if(nfold > (mpi.comm.size() - 1))
        stop("nfold > number of mpi Rslaves")
    
    cat("\n\n Sending objects to MPI space\n\n")

    mpi.bcast.Robj2slave(x)
    mpi.bcast.Robj2slave(time)
    mpi.bcast.Robj2slave(event)
    mpi.bcast.Robj2slave(p)
    mpi.bcast.Robj2slave(maxSize)
    mpi.bcast.Robj2slave(minSize)
    mpi.bcast.Robj2slave(minCor)
    mpi.bcast.Robj2slave(MaxIterationsCox)
    mpi.bcast.Robj2slave(res1s)
    ## debug:
    global.res1s <<- res1s
    global.index.select <<- index.select
    mpi.bcast.Robj2slave(index.select)
##    mpi.bcast.cmd(foldNumber <- mpi.comm.rank())
##    mpi.bcast.Robj2slave(DaveCVPred.res1Given.InternalMPI)


    ## We use papply again
    cat("\n\n Computing the rest\n")
    tmp1 <- papply(as.list(1:nfold),
                 DaveCVPred.res1Given.InternalMPI2)
    cat("\n\n Cleaning up and closing MPI\n")
    try(mpi.close.Rslaves())
    
    for(i in 1:nfold) {
        OOB.scores[index.select == i] <-
            tmp1[[i]]$scoresTest
        tmp1[[i]]$scoresTest <- NULL  ## don't need this anymore
    }
    
    out <- list(cved.models = tmp1,
                OOB.scores = OOB.scores)
    class(out) <- "cvDave"
    return(out)
}




summary.cvDave <- function(object, allDataObject.selected, subjectNames, genenames, html = TRUE) {
    
    oobs <- matrix(object$OOB.scores, ncol = 1)
    rownames(oobs) <- subjectNames
    if (!html){
        cat("\n Out-of-bag scores\n\n")
        print(oobs)
    }
    if (html){
        cat("\n<h3>5.1. <a href=\"scores.oob.html\" target=\"scores_window\">View</a> out-of-bag scores.</h3>\n")
        cleanHTMLhead(file = "scores.oob.html",
                      title = "Linear predictor scores for out-of-bag data")
        write(paste("<TABLE frame=\"box\">\n",
                    "<tr><th>Subject/array</th> <th>Linear score</th></tr>\n"),
              file = "scores.oob.html", append = TRUE)
        wout <- ""
        for(i in 1:nrow(oobs)) {
            wout <- paste(wout, "\n <tr align=right>",
            "<td>", rownames(oobs)[i], "</td><td>", oobs[i], "</td></tr>\n")
        }
        wout <- paste(wout, "</TABLE>")
        write(wout, file = "scores.oob.html", append = TRUE)
        cleanHTMLtail(file = "scores.oob.html")


##         cat("<TABLE frame=\"box\">\n")
##         cat("<tr><th>Array/subject name</th> <th>Linear score</th></tr>\n")
##         for(i in 1:nrow(oobs)) {
##             cat("\n <tr align=right>")
##             cat("<td>", rownames(oobs)[i], "</td><td>", oobs[i], "</td></tr>\n")
##         }
##         cat("</TABLE>")
    }
        
    
    object <- object[[1]] ## don't need scores anymore. Simpler subsetting.

    allWarnings <- lapply(object, function(x) x$Warn)
    
    ks <- length(object)
    cv.names <- paste("CV.run.", 1:ks, sep = "")    


    if(html) {
        cat("\n\n <h3>5.2 Selected components, genes and coefficients in cross-validation runs</h3>\n")
        selectedInR <- list()
        for(cvr in 1:ks) {
            cat("\n\n <h4>CV run ", cvr, "</h4>\n")
            selectedInR[[cvr]] <- selectedSignatures(object[[cvr]]$fmDaveObject,
                                                     genenames = genenames, print = TRUE,
                                                     out = TRUE)
            cat("\n <hr align=\"left\" width=80>")
        }
    } else {
        cat("\n\n Selected components, genes and coefficients in cross-validation runs\n")
        cat("========================================================================\n\n")
        selectedInR <- list()
        for(cvr in 1:ks) {
            cat("\n\n CV run ", cvr, "\n")
            cat("------------------------\n")
            selectedInR[[cvr]] <- selectedSignatures(object[[cvr]]$fmDaveObject,
                                                     genenames = genenames, print = TRUE,
                                                     out = TRUE)
            cat("\n--------------------------------------------------------------------------------\n\n")
        }
    }


    
    tmp.genesSelected <- list()
    tmp.genesSelected[[1]] <- unlist(allDataObject.selected$Genes)
    genesSelected.cv <- lapply(selectedInR, function(x)
                               unlist(x$Genes))
    tmp.genesSelected <- c(tmp.genesSelected, genesSelected.cv)

    shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 1))
    for(i in 1:(ks + 1)) {
        for(j in 1:(ks + 1)) { ## sure, need not be symmetric, but this is fast
            shared.genes[i, j] <-
                length(intersect(tmp.genesSelected[[i]],
                                 tmp.genesSelected[[j]]))
        }
    }

    ngenes <- sapply(tmp.genesSelected, length)
    prop.shared <- round(shared.genes/ngenes, 3)
    
    ngenesS <- paste("(", ngenes, ")", sep = "")
    colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", cv.names)
    rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", cv.names), ngenesS)

    if(html) cat("<h3>5.3. Shared genes between runs</h3>\n")
    if(html) cat("<pre>")
    options(width = 200)
    cat("\n\n Number of shared genes \n")
    print(as.table(shared.genes))
    
    cat("\n\n Proportion of shared genes (relative to row total) \n")
    print(as.table(prop.shared))
    if(html) cat("</pre>")
    options(width = 80)
    unlisted.genes.selected <- unlist(genesSelected.cv)
    
    in.all.data <-
        which(names(table(unlisted.genes.selected, dnn = NULL)) %in% tmp.genesSelected[[1]])
    tmp1 <- sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE)
    if(!html) {
        cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
        print(tmp1)
        cat("\n")
    }
    if(html) {
        cat("\n<h3>5.4. Gene freqs. in cross-validated runs of genes selected in model with all data</h3> \n\n")
        cat("<TABLE frame=\"box\">\n")
        cat("<tr><th width=200>Gene name</th> <th width=200>Frequency</th></tr>\n")
        for(i in 1:length(tmp1)) {
            cat("\n <tr align=right>")
            cat("<td>", linkGene(names(tmp1)[i]), "</td><td>", tmp1[i], "</td></tr>\n")
        }
        cat("</TABLE>")    
    }
    tmp2 <- sort(table(unlisted.genes.selected, dnn = NULL),
                 decreasing = TRUE)
    if(!html) {
        cat("\n\n Gene frequencies in cross-validated runs \n\n")
        print(tmp2)
        cat("\n")
    }
    if(html){
        cat("\n\n<h3>5.5. Gene frequencies in cross-validated runs</h3> \n\n")
        cat("<TABLE frame=\"box\">\n")
        cat("<tr><th width=200>Gene name</th> <th width=200>Frequency</th></tr>\n")
        for(i in 1:length(tmp2)) {
            cat("\n <tr align=right>")
            cat("<td>", linkGene(names(tmp2)[i]), "</td><td>", tmp2[i], "</td></tr>\n")
        }
        cat("</TABLE>")
    }
}

selectedSignatures <- function(fmDave, genenames,
                               print = TRUE, out = TRUE,
                               html = TRUE) {

    if(is.na(fmDave)) {
        out <- NA
        if(print) print(NA)
        return(NA)
    }
        
    res2 <- fmDave$clusterResults
    selectedSign <- names(fmDave$model$coeff)
    signatures <- list()
    coeffs <- fmDave$model$coeff

    if(!length(selectedSign)) {
        signatures$component[1] <- 1
        signatures$Name[1] <- "Null model (only an intercept fitted, no components selected)"
        signatures$Coeffs[1] <- NA
        signatures$Genes[[1]] <- NA
        signatures$numGenes[1] <- 1
	}
	else {
    for(i in 1:length(selectedSign)) {
        signatures$component[i] <- i
        signatures$Name[i] <- selectedSign[i]
        signatures$Coeffs[i] <- coeffs[i]
        
        ind1 <- which(res2$filteredGroupsNegative ==
                      selectedSign[i])
        if(length(ind1)) {
            signatures$Genes[[i]] <- genenames[res2$filteredNegPositions[ind1]]
            signatures$numGenes[i] <- length(signatures$Genes[[i]])
        } else {
            ind2 <- which(res2$filteredGroupsPositive ==
                          selectedSign[i])
            if(length(ind2)) {
                signatures$Genes[[i]] <- genenames[res2$filteredPosPositions[ind2]]
                signatures$numGenes[i] <- length(signatures$Genes[[i]])
            } else {
                stop("No genes for this signature; weird!")
            }
        }
    }
   }
    signatures$totalComponents <- length(selectedSign)
    signatures$totalGenes <- sum(signatures$numGenes)

    if(print & (!html)) {
        total.length <- length(unlist(signatures$Genes)) +
            2* length(selectedSign) + 1
        ## OK, this is a horrible hack, but works for now
        c1 <- c2 <- c3 <- rep("", total.length)
        k <- 1
        c1[1] <- paste(rep("=", 14), sep = "", collapse = "")
        c3[1] <- paste(rep("=", 50), sep = "", collapse = "")
        c2[1] <- paste(rep("=", 11), sep = "", collapse = "")
        k <- 2
        for(i in 1:length(signatures$component)) {
            c1[k] <- signatures$Name[i]
            c2[k] <- round(signatures$Coeffs[i], 4)
            k <- k + 1
            for(j in 1:length(signatures$Genes[[i]])) {
                c3[k] <- signatures$Genes[[i]][j]
                k <- k + 1
            }
            c1[k] <- paste(rep("-", 14), sep = "", collapse = "")
            c3[k] <- paste(rep("-", 50), sep = "", collapse = "")
            c2[k] <- paste(rep("-", 11), sep = "", collapse = "")
            k <- k + 1
        }

        dffp <- data.frame(cbind("Component_Name" = c1,
                                 "Genes" = c3,
                                 "Coefficient" = c2))
        options(width = 120)
        cat(paste("\n\n Total of ", signatures$totalComponents,
                  "signature components selected and ", signatures$totalGenes,
                  "genes.\n"))
        cat("\n Selected signature components, genes and coefficients\n\n")
        print(dffp)
        options(width = 80)
    }
    if(print & html) {
        total.length <- length(unlist(signatures$Genes)) +
            2* length(selectedSign) + 1
        ## OK, this is a horrible hack, but works for now

        cat(paste("\n<p> Total of ", signatures$totalComponents,
                  "signature components selected and ", signatures$totalGenes,
                  "genes.</p>"))
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=150>Component name</th> <th width=200>Genes</th> <th width=80>Coefficient</th></tr>")
##         c1 <- c2 <- c3 <- rep("", total.length)
        k <- 1
##         c1[1] <- paste(rep("=", 14), sep = "", collapse = "")
##         c3[1] <- paste(rep("=", 50), sep = "", collapse = "")
##         c2[1] <- paste(rep("=", 11), sep = "", collapse = "")
        k <- 2
        for(i in 1:length(signatures$component)) {
            cat("\n<tbody>")
            cat("\n <tr align=center>\n")
            cat("<th> <b>", signatures$Name[i],
                "</b></th><th></th><th><b>",
                round(signatures$Coeffs[i], 4),
                "</b></th></tr>")
            ##             c1[k] <- signatures$Name[i]
            ##             c2[k] <- round(signatures$Coeffs[i], 4)
            k <- k + 1
            for(j in 1:length(signatures$Genes[[i]])) {
##                c3[k] <- signatures$Genes[[i]][j]
                k <- k + 1
                cat("\n<tr><td></td><td>",
                    linkGene(signatures$Genes[[i]][j]),
                    "</td><td></td></tr>")
                
            }
            cat("\n</tbody>\n")
##             c1[k] <- paste(rep("-", 14), sep = "", collapse = "")
##             c3[k] <- paste(rep("-", 50), sep = "", collapse = "")
##             c2[k] <- paste(rep("-", 11), sep = "", collapse = "")
            k <- k + 1
        }
        cat("\n</TABLE>\n")
        
##         dffp <- data.frame(cbind("Component_Name" = c1,
##                                  "Genes" = c3,
##                                  "Coefficient" = c2))
##         options(width = 120)
##         print(dffp)
##         options(width = 80)
    }

    if(out) return(signatures)
}





################################################################

#############       Utility functions

################################################################


linkGene <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else
        paste("<a href=\"http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&internal=0&org=",
              organism,"\" target=\"icl_window\" >",id,"</a>", sep = "")
## target=\"icl_window\"\
}
     


linkGene2 <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else 
        paste("http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&internal=0&org=",
              organism, sep = "")
}

imagemap3 <- function(filename,width=480,height=480,
                      title='Imagemap from R', ps = 12){

    GDD(file = paste(filename,".png",sep=''),w=width, h=height,
        type = "png", ps = ps)	  
	  
    im <- list()
    im$Device <- dev.cur()
    im$Filename=filename
    im$Height=height
    im$Width=width
    im$Objects <- list()
    im$HTML <- list()
    im$title <- title
    
    class(im) <- "imagemap"
    im
}

heightLeave <- function(x) {
    if (is.leaf(x))
        attr(x, "height")
    else
        NULL
}
nameLeave <- function(x) {
    if (is.leaf(x))
        attr(x, "label")
    else
        NULL
}











#### Why are we using papply? See debug.R in upper directory to find out.
#### Essentially, mpi by itslef will not deal well with errors even if we
#### use try.





### This is no longer used.
## old.cvDave.parallel3 <- function(x, time, event,
##                              p, maxSize,
##                              minSize, minCor,
##                              MaxIterationsCox,
##                              nfold, mpiHosts) {

##     if (mpi.comm.size(comm = 1) == 0) {
##         mpiSpawnAll()
##     } else { ## so mpi is running
##         if ((mpi.comm.size(comm = 1) - 1) < mpi.universe.size()) {
##             ## but few salves
##             mpi.close.Rslaves()
##             mpiSpawnAll()
##         } else {
##             mpiMyCleanSetup()
##         }
##     }
       
##     ## we asume at least as many mpi Rslaves as nfold.
##     ## mpiHosts used for the second spawning of slaves,
##     ## only for the 10 fold runs.
    
    
##     n <- length(time)
##     index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
##     OOB.scores <- rep(NA, n)

##     res1s <- list()

##     ### Here we use papply, which is 
##     argspapp <- list()
##     cat("\n\n Computing gene-wise cox p-value\n")
##     for(i in 1:nfold) {
##         cat("\n  ....  fold ", i)
##         xtr <- x[index.select != i, , drop = FALSE]
##         ttr <- time[index.select != i]
##         etr <- event[index.select != i]
##         res1s[[i]] <- dStep1.parallel(xtr, ttr, etr, p = p,
##                                       MaxIterationsCox = MaxIterationsCox)
##     }

##     cat("\n\n Cleaning up MPI space, and setting up a new one\n\n")
##     ## Clean previous stuff
##     mpi.close.Rslaves()
##     mpiSpawnThis(hosts = mpiHosts)

##     if(nfold > (mpi.comm.size() - 1))
##         stop("nfold > number of mpi Rslaves")
    
##     cat("\n\n Sending objects to MPI space\n\n")

##     mpi.bcast.Robj2slave(x)
##     mpi.bcast.Robj2slave(time)
##     mpi.bcast.Robj2slave(event)
##     mpi.bcast.Robj2slave(p)
##     mpi.bcast.Robj2slave(maxSize)
##     mpi.bcast.Robj2slave(minSize)
##     mpi.bcast.Robj2slave(minCor)
##     mpi.bcast.Robj2slave(MaxIterationsCox)
##     mpi.bcast.Robj2slave(res1s)
##     ## debug:
##     global.res1s <<- res1s
##     global.index.select <<- index.select
##     mpi.bcast.Robj2slave(index.select)
##     mpi.bcast.cmd(foldNumber <- mpi.comm.rank())
##     mpi.bcast.Robj2slave(DaveCVPred.res1Given.InternalMPI)

##     cat("\n\n Computing the rest\n")
##     tmp1 <- mpi.remote.exec(DaveCVPred.res1Given.InternalMPI())

##     ## debung
##     print("printing tmp1")
##     print(tmp1)
    
##     cat("\n\n Cleaning up and closing MPI\n")
##     ##try(mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir = .GlobalEnv)))
##  ##    cat("mpi.universe.size() ", mpi.universe.size())
## ##     cat("mpi.comm.size()", mpi.comm.size()) 
##     try(mpi.close.Rslaves())
    
##     for(i in 1:nfold) {
##         OOB.scores[index.select == i] <-
##             tmp1[[i]]$scoresTest
##         tmp1[[i]]$scoresTest <- NULL  ## don't need this anymore
##     }
    
##     out <- list(cved.models = tmp1,
##                 OOB.scores = OOB.scores)
##     class(out) <- "cvDave"
##     return(out)
## }









## summary.cvDave <- function(object, allDataObject.selected, subjectNames, genenames) {
##     cat("\n Out-of-bag scores\n\n")

##     cat("\n Out-of-bag scores\n\n")
##     oobs <- matrix(object$OOB.scores, ncol = 1)
##     rownames(oobs) <- subjectNames
##     print(oobs)


    
##     object <- object[[1]] ## don't need scores anymore. Simpler subsetting.

##     allWarnings <- lapply(object, function(x) x$Warn)
    
##     ks <- length(object)
##     cv.names <- paste("CV.run.", 1:ks, sep = "")    


    
##     cat("\n\n Selected components, genes and coefficients in cross-validation runs\n")
##     cat("========================================================================\n\n")
##     selectedInR <- list()
##     for(cvr in 1:ks) {
##         cat("\n\n CV run ", cvr, "\n")
##         cat("------------------------\n")
##         selectedInR[[cvr]] <- selectedSignatures(object[[cvr]]$fmDaveObject,
##                                                  genenames = genenames, print = TRUE,
##                                                  out = TRUE)
##         cat("\n--------------------------------------------------------------------------------\n\n")
##     }
    
##     tmp.genesSelected <- list()
##     tmp.genesSelected[[1]] <- unlist(allDataObject.selected$Genes)
##     genesSelected.cv <- lapply(selectedInR, function(x)
##                                unlist(x$Genes))
##     tmp.genesSelected <- c(tmp.genesSelected, genesSelected.cv)

##     shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 1))
##     for(i in 1:(ks + 1)) {
##         for(j in 1:(ks + 1)) { ## sure, need not be symmetric, but this is fast
##             shared.genes[i, j] <-
##                 length(intersect(tmp.genesSelected[[i]],
##                                  tmp.genesSelected[[j]]))
##         }
##     }

##     ngenes <- sapply(tmp.genesSelected, length)
##     prop.shared <- round(shared.genes/ngenes, 3)
    
##     ngenesS <- paste("(", ngenes, ")", sep = "")
##     colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", cv.names)
##     rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", cv.names), ngenesS)
    
##     options(width = 200)
##     cat("\n\n Number of shared genes \n")
##     print(as.table(shared.genes))
    
##     cat("\n\n Proportion of shared genes (relative to row total) \n")
##     print(as.table(prop.shared))
    
##     options(width = 80)
##     unlisted.genes.selected <- unlist(genesSelected.cv)
    
##     in.all.data <-
##         which(names(table(unlisted.genes.selected, dnn = NULL)) %in% tmp.genesSelected[[1]])
##     cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
##     print(sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE))
##     cat("\n")
    
    
##     cat("\n\n Gene frequencies in cross-validated runs \n\n")
##     tmp.table <- sort(table(unlisted.genes.selected, dnn = NULL),
##                       decreasing = TRUE)
##     print(tmp.table)
##     cat("\n")
## }

## ## why genenames? aren't these the columns of the data set??
## ## Probably allows to pass a vector with the names I want








## selectedSignatures <- function(fmDave, genenames,
##                                print = TRUE, out = TRUE) {

##     if(is.na(fmDave)) {
##         out <- NA
##         if(print) print(NA)
##         return(NA)
##     }
        
##     res2 <- fmDave$clusterResults
##     selectedSign <- names(fmDave$model$coeff)
##     signatures <- list()
##     coeffs <- fmDave$model$coeff

##     for(i in 1:length(selectedSign)) {
##         signatures$component[i] <- i
##         signatures$Name[i] <- selectedSign[i]
##         signatures$Coeffs[i] <- coeffs[i]
        
##         ind1 <- which(res2$filteredGroupsNegative ==
##                       selectedSign[i])
##         if(length(ind1)) {
##             signatures$Genes[[i]] <- genenames[res2$filteredNegPositions[ind1]]
##             signatures$numGenes[i] <- length(signatures$Genes[[i]])
##         } else {
##             ind2 <- which(res2$filteredGroupsPositive ==
##                           selectedSign[i])
##             if(length(ind2)) {
##                 signatures$Genes[[i]] <- genenames[res2$filteredPosPositions[ind2]]
##                 signatures$numGenes[i] <- length(signatures$Genes[[i]])
##             } else {
##                 stop("No genes for this signature; weird!")
##             }
##         }
##     }

##     signatures$totalComponents <- length(selectedSign)
##     signatures$totalGenes <- sum(signatures$numGenes)

##     if(print) {
##         total.length <- length(unlist(signatures$Genes)) +
##             2* length(selectedSign) + 1
##         ## OK, this is a horrible hack, but works for now
##         c1 <- c2 <- c3 <- rep("", total.length)
##         k <- 1
##         c1[1] <- paste(rep("=", 14), sep = "", collapse = "")
##         c3[1] <- paste(rep("=", 50), sep = "", collapse = "")
##         c2[1] <- paste(rep("=", 11), sep = "", collapse = "")
##         k <- 2
##         for(i in 1:length(signatures$component)) {
##             c1[k] <- signatures$Name[i]
##             c2[k] <- round(signatures$Coeffs[i], 4)
##             k <- k + 1
##             for(j in 1:length(signatures$Genes[[i]])) {
##                 c3[k] <- signatures$Genes[[i]][j]
##                 k <- k + 1
##             }
##             c1[k] <- paste(rep("-", 14), sep = "", collapse = "")
##             c3[k] <- paste(rep("-", 50), sep = "", collapse = "")
##             c2[k] <- paste(rep("-", 11), sep = "", collapse = "")
##             k <- k + 1
##         }

##         dffp <- data.frame(cbind("Component_Name" = c1,
##                                  "Genes" = c3,
##                                  "Coefficient" = c2))
##         options(width = 120)
##         cat(paste("\n\n Total of ", signatures$totalComponents,
##                   "signature components selected and ", signatures$totalGenes,
##                   "genes.\n"))
##         cat("\n Selected signature components, genes and coefficients\n\n")
##         print(dffp)
##         options(width = 80)
##     }
##     if(out) return(signatures)
## }


## Make dencrograms vertical!! Oh well, depends on dendrogram, whose API might change...
## Imagemaps and dendrograms: search in R
## Provide two "families": with all names and only with cluster names

## If all names in color, no need for boxes?
## Or just box, no fill-up?


##  2. Cluster; independently for those with pos and neg beta.

##  2. Cluster; independently for those with pos and neg beta.
## dStep2 <- function(x, res.mat, maxSize, minSize,
##                    minCor, plot,
##                    interactive) {
##     res.mat[is.na(res.mat[, 4]), 4] <- 0

##     if(sum(res.mat[, 4] == 1) >= minSize) {
##         pos.data <- x[, res.mat[, 4] == 1]
##         pdok <- TRUE
##     } else {
##         pdok <- FALSE
##         warning(paste("Not enough positive coeff. genes that",
##                       "meet the p restrictions."))
##     }
##     if(sum(res.mat[, 4] == -1) >= minSize) {
##         neg.data <- x[, res.mat[, 4] == -1]
##         pnok <- TRUE
##     } else {
##         pnok <- FALSE
##         warning(paste("Not enough negative coeff. genes that",
##                       "meet the p restrictions."))
##     }

##     if((!pnok) & (!pdok)) {
## 	if(interactive) caughtUserError("No gene was above the minimal p threshold")
##         else {
##             return(NA)
##         }
##     }

##     tn <- tp <- FALSE
##     if(pdok) {
##         pos.clus <- hclust(as.dist(1 -cor(pos.data)), method = "complete")
##         pos.groups <- paste("P.", cutree(pos.clus, h = 1- minCor), sep = "")
##         tpos <- table(pos.groups)
##         pos.accept <- names(which((tpos >= minSize) & (tpos <= maxSize)))
##         if(length(pos.accept)) {
##             groupsPositive <- pos.groups[pos.groups %in% pos.accept]
##             dataPositive <- pos.data[ , pos.groups %in% pos.accept]
##             posGroups <- unique(groupsPositive)
##             posMeanData <- matrix(NA, nrow = dim(dataPositive)[1],
##                               ncol = length(posGroups))

##             for(i in 1:length(posGroups)) {
##                 posMeanData[, i] <-
##                     apply(dataPositive[, groupsPositive == posGroups[i]],
##                           1, mean)
##             }
##             colnames(posMeanData) <- posGroups
##             tp <- TRUE
##         } else {
##             tp <- FALSE
##             warning(paste("No groups of positive coeff genes that",
##                           "meet the p, minimum correlation and size restrictions."))
##         }
##     }

##     if(pnok) {
##         neg.clus <- hclust(as.dist(1 -cor(neg.data)), method = "complete")
##         neg.groups <- paste("N.", cutree(neg.clus, h = 1- minCor), sep = "")
##         npos <- table(neg.groups)
##         neg.accept <- names(which((npos >= minSize) & (npos <= maxSize)))
##         if(length(neg.accept)) {
##             groupsNegative <- neg.groups[neg.groups %in% neg.accept]
##             dataNegative <- neg.data[ , neg.groups %in% neg.accept]
##             negGroups <- unique(groupsNegative)
##             negMeanData <- matrix(NA, nrow = dim(dataNegative)[1],
##                                   ncol = length(negGroups))

##             for(i in 1:length(negGroups)) {
##                 negMeanData[, i] <-
##                     apply(dataNegative[, groupsNegative == negGroups[i]],
##                           1, mean)
##             }
##             colnames(negMeanData) <- negGroups
##             tn <- TRUE
##         } else {
##             tn <- FALSE
##             warning(paste("No groups of negative coeff. genes that",
##                           "meet the p, minimum correlation and size restrictions."))
##         }
##     }

##     if(!tn & !tp) {
##         if(interactive) caughtUserError(paste("No groups that meet the p, minimum correlation",
##                                               "and size restrictions."))
##         else return(NA)
##     }
##     pdok <- tp & pdok
##     pnok <- tn & pnok
    
##     if(plot) {
## ####         if(dev.interactive()) {
## ####             op <- par(ask = TRUE, las = 1, mar = c(1, 5, 4, 1))
## ####         } else {
## ####             op <- par(las = 1, mar = c(1, 5, 4, 1))
## ####         }
## ####         on.exit(par(op))

## ###        par(las = 1, mar = c(3, 5, 4, 1))
        
##         ## works well with: pdf(file = "test1.pdf", height = 9, width = 19)

        
##         pdokf <- function() {
##             pl.p.groups <-
##                 as.numeric(sapply(pos.accept,
##                                   function(x) strsplit(x, "P.")[[1]][2]))
##             pos.labels <- rep("                   ", ncol(pos.data))
##             index.labels <- which(pos.groups %in% pos.accept)
##             pos.labels[index.labels] <- colnames(pos.data)[index.labels]
            
##             plot(pos.clus, hang = 0.001, sub = "", xlab = "",
##                  labels = pos.labels, cex = 0.7, main = "Positive coefficients",
##                  ylab = "1 - correlation")
##             abline(h = 1 - minCor, lty = 2, col = "blue")
##             rainbow.col <- rainbow(length(posGroups))
##             for(i in 1:length(posGroups)) {
##                 wthis.original <- which(pos.groups == posGroups[i])
##                 wthis.plot.order <- which(pos.clus$order %in% wthis.original)
##                 minx <- min(wthis.plot.order)
##                 maxx <- max(wthis.plot.order)
##                 polygon(x = c(minx, minx, maxx, maxx),
##                         y =c(0, 1 - minCor, 1 - minCor, 0),
##                         col = rainbow.col[i],
##                         density = 10,
##                         lty = 2)
##                 axis(1, line = 0, at = c(minx, maxx), col = rainbow.col[i],
##                      tick = TRUE, labels = FALSE, lw = 3)
##                 axis(1, line = 0, at = 0.5 * (maxx + minx),
##                      col.axis = rainbow.col[i],
##                      tick = FALSE, labels = posGroups[i], lw = 0,
##                      cex.axis = 1.5)
## ##                 segments(x0 = minx, y0= -0.15, x1 = maxx, y1 = -0.15,
## ##                          col = rainbow.col[i], lwd = 4, lty = 1)
## ##                 text(x = 0.5* (maxx + minx), y = -0.20,
## ##                      labels = posGroups[i])
##             }
##         } ##</pdok within plotting>

##         if (pdok) {
##             pdf(file = "ClusterPositiveCoeffs.pdf", height = 9, width = 19)
##             pdokf()
##             dev.off()
##             webPNG(file = "ClusterPositiveCoeffs.png", height = 9,
##                    width = 19, pointsize = png.pointsize,
##                    family = png.family)
##             par(cex.axis = 0.75); par(cex.lab = 1.4); par(cex.main = 1.5)
##             pdokf()
##             dev.off()
##         } else {
##             pdf(file = "ClusterPositiveCoeffs.pdf", height = 9, width = 9)
##             plot(x = c(0, 1), y = c(0, 1),
##                  type = "n", axes = FALSE, xlab = "", ylab = "")
##             box()
##             text(0.5, 0.7,
##             "There are no genes with positive coefficients")
##             text(0.5, 0.5, "that satisfy the p, minimum correlation")
##             text(0.5, 0.3, "and size restrictions.")
##             dev.off()

            
##             webPNG(file = "ClusterPositiveCoeffs.png", height = 9,
##                    width = 9, pointsize = png.pointsize,
##                    family = png.family)
##             plot(x = c(0, 1), y = c(0, 1),
##                  type = "n", axes = FALSE, xlab = "", ylab = "")
##             box()
##             text(0.5, 0.7,
##             "There are no genes with positive coefficients")
##             text(0.5, 0.5, "that satisfy the p, minimum correlation")
##             text(0.5, 0.3, "and size restrictions.")
##             dev.off()
##         }

##         pnokf <- function() {
##             pl.n.groups <-
##                 as.numeric(sapply(neg.accept,
##                                   function(x) strsplit(x, "N.")[[1]][2]))
##             neg.labels <- rep("                   ", ncol(neg.data))
##             index.labels <- which(neg.groups %in% neg.accept)
##             neg.labels[index.labels] <- colnames(neg.data)[index.labels]
            
##             ## need to remap from the original position to the position
##             ## a variable takes in the plot (its order in the clustering
##             ## object).
            
##             plot(neg.clus, hang = 0.001, sub = "", xlab = "",
##                  labels = neg.labels, cex = 0.7, main = "Negative coefficients",
##                  ylab = "1 - correlation")
##             abline(h = 1 - minCor, lty = 2, col = "blue")
            
##             rainbow.col <- rainbow(length(negGroups))
            
##             for(i in 1:length(negGroups)) {
##                 wthis.original <- which(neg.groups == negGroups[i])
##                 wthis.plot.order <- which(neg.clus$order %in% wthis.original)
##                 minx <- min(wthis.plot.order)
##                 maxx <- max(wthis.plot.order)
##                 polygon(x = c(minx, minx, maxx, maxx),
##                         y =c(0, 1 - minCor, 1 - minCor, 0),
##                         col = rainbow.col[i],
##                         density = 10,
##                         lty = 2)
##                 axis(1, line = 0, at = c(minx, maxx), col = rainbow.col[i],
##                      tick = TRUE, labels = FALSE, lw = 3)
##                 axis(1, line = 0, at = 0.5 * (maxx + minx),
##                      col.axis = rainbow.col[i],
##                      tick = FALSE, labels = negGroups[i], lw = 0,
##                      cex.axis = 1.5)
## ##                 segments(x0 = minx, y0= -0.15, x1 = maxx, y1 = -0.15,
## ##                          col = rainbow.col[i], lwd = 4, lty = 1)
## ##                 text(x = 0.5* (maxx + minx), y = -0.2 ,
## ##                      labels = negGroups[i])
##             }
##         } ## </ if(pnok) within plot

##         if (pnok) {
##             pdf(file = "ClusterNegativeCoeffs.pdf", height = 9, width = 19)
##             pnokf()
##             dev.off()
##             webPNG(file = "ClusterNegativeCoeffs.png", height = 9,
##                    width = 19, pointsize = png.pointsize,
##                    family = png.family)
##             par(cex.axis = 0.75); par(cex.lab = 1.4); par(cex.main = 1.5)
##             pnokf()
##             dev.off()
##         } else {
##             pdf(file = "ClusterNegativeCoeffs.pdf", height = 9, width = 9)
##             plot(x = c(0, 1), y = c(0, 1),
##                  type = "n", axes = FALSE, xlab = "", ylab = "")
##             box()
##             text(0.5, 0.7,
##             "There are no genes with negative coefficients")
##             text(0.5, 0.5, "that satisfy the p, minimum correlation")
##             text(0.5, 0.3, "and size restrictions.")
##             dev.off()

            
##             webPNG(file = "ClusterNegativeCoeffs.png", height = 9,
##                    width = 9, pointsize = png.pointsize,
##                    family = png.family)
##             plot(x = c(0, 1), y = c(0, 1),
##                  type = "n", axes = FALSE, xlab = "", ylab = "")
##             box()
##             text(0.5, 0.7,
##             "There are no genes with negative coefficients")
##             text(0.5, 0.5, "that satisfy the p, minimum correlation")
##             text(0.5, 0.3, "and size restrictions.")
##             dev.off()
##         }

        
##     } ##</ if plot>

##     ## For predictions and results, which vars. correspond
##     ## to which original genes

##     ## The following two are the indices of the columns in the original (non
##     ## divided data file)

##     if(pdok) {
##         posPositions <- which(res.mat[, 4] == 1)
##         filteredPosPositions <- posPositions[pos.groups %in% pos.accept]
##     }
##     if(pnok) {
##         negPositions <- which(res.mat[, 4] == -1)
##         filteredNegPositions <- negPositions[neg.groups %in% neg.accept]
##     }

##     if(pnok & pdok)
##         if(length(intersect(filteredPosPositions, filteredNegPositions)))
##             stop("Non zero intersection between filtered pos and neg positions")
##     if(!pnok) {
##         negMeanData <- NULL
##         groupsNegative <- NA
##         filteredGroupsNegative <- NA
##         filteredNegPositions <- NA
##         negPositions <- NA
##     }
##     if(!pdok) {
##         posMeanData <- NULL
##         groupsPositive <- NA
##         filteredGroupsPositive <- NA
##         filteredPosPositions <- NA
##         posPositions <- NA
##     }
    
##     return(list(md = cbind(posMeanData, negMeanData),
##                 filteredGroupsPositive = groupsPositive,
##                 filteredGroupsNegative = groupsNegative,
##                 filteredPosPositions = filteredPosPositions,
##                 filteredNegPositions = filteredNegPositions,
##                 posPositions = posPositions,
##                 negPositions = negPositions))
## }
