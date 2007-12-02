###  Not here because I do not want Rmpi or CGIwithR executed in the slaves.
#### library(CGIwithR)
#### library(Rmpi)
#### library(survival)
#### library(combinat)
#### library(MASS)

require(GDD)
#require(survival)
#require(imagemap) ### FIXME: include the needed code below
library(party)
library(mboost)



###############################################
##########                    #################
########## glmboost functions #################
##########                    #################
###############################################

my.glmboost <- function(x, time, event, newdata, mstop = 500) {

    if(!is.matrix(x)) x <- as.matrix(x)
    gb1 <- glmboost(x, Surv(time, event), family = CoxPH(),
                    control = boost_control(mstop = mstop,
                    center = TRUE))
    ## for 10-fold CV for risk; from help for cvrisk
    n <- nrow(x)
    k <- 10
    ntest <- floor(n / k)
    cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, n)), k - 1), 
                      rep(0, n * k - (k - 1) * (n + ntest))), nrow = n)
    gridf <- c(5, 10, 25, seq(from = 50, to = 500, length = 10))
    cvm <- cvrisk(gb1, folds = cv10f, grid = gridf)
    best.mstop <- mstop(cvm)
    cat("\n Best mstop is ", best.mstop, "\n")
    gb1 <- gb1[mstop(cvm)]
    selected.genes.rows <- which(abs(coef(gb1)) > 0)
    lsgenes <- length(selected.genes.rows)
    selected.genes.names <- names(selected.genes.rows)
    selected.genes.stats <- cbind(coef(gb1)[selected.genes.rows],
                                  rep(NA, lsgenes),
                                  rep(NA, lsgenes),
                                  rep(NA, lsgenes),
                                  rep(0, lsgenes),
                                  rep(NA, lsgenes))
    overfit_predicted_surv_time <- predict(gb1, newdata = x, type = "lp")
    if(!(is.null(newdata))) {
        pred.stime <- predict(gb1, newdata = newdata, type = "lp")
    } else {
        pred.stime <- NULL
    }
    return(list(selected.genes.stats = selected.genes.stats,
                selected.genes.rows = selected.genes.rows,
                selected.genes.names = selected.genes.names,
                selected.genes.number = lsgenes,
                glmboost_ob = gb1,
                predicted_surv_time = pred.stime,
                overfit_predicted_surv_time = overfit_predicted_surv_time))
}

my.glmboost.cv <- function(x, time, event, mstop = 500, nfold = 10) {
    n <- length(time)
    index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    OOB.scores <- rep(NA, n)

    my.glmboost.internal.MPI <- function(fnum) {
        ## to be used with papply
        cat("\n Doing fnum ", fnum, "\n") ## FIXME: debug
        xtrain <- x[index.select != fnum, , drop = FALSE]
        xtest <- x[index.select == fnum, , drop = FALSE]
        timetrain <- time[index.select != fnum]
        eventtrain <- event[index.select != fnum]
        retval <- my.glmboost(xtrain, timetrain, eventtrain, xtest)
        return(retval)
    }

    tmp1 <- papply(as.list(1:nfold),
                   my.glmboost.internal.MPI,
                   papply_commondata = list(x = x,
                   time = time, event = event, 
                   index.select = index.select))
    
    for(i in 1:nfold) {
        OOB.scores[index.select == i] <-
            tmp1[[i]]$predicted_surv_time
        tmp1[[i]]$predicted_surv_time <- NULL  ## don't need this anymore
    }

    out <- list(cved.models = tmp1,
                OOB.scores = OOB.scores)
    return(out)
}

    

###############################################
##########                    #################
########## cforest functions  #################
##########                    #################
###############################################


geneSelect <- function(x, sobject, numgenes) {
    ### Select the "best" (according to p-value from Cox) numgenes
    ## a modification of dStep1.serial
    
    MaxIterationsCox <- 200
    res.mat <- matrix(NA, nrow = numgenes, ncol = 6)
    funpap3 <- function (x) {
        out1 <-
            coxph.fit.pomelo0(x, sobject,
                              control = coxph.control(iter.max =  MaxIterationsCox))
        if(out1$warnStatus > 1) {
            return(c(0, NA, out1$warnStatus))
        } else {
            sts <- out1$coef/sqrt(out1$var)
            return(c(out1$coef, 1- pchisq((sts^2), df = 1), out1$warnStatus))
        }
    }
    tmp0 <- t(apply(x, 2, funpap3))
    to.keep <- order(tmp0[, 2])[1:numgenes]
    tmp <- tmp0[to.keep, ]
    res.mat[, 1:2] <- tmp[, 1:2]
    res.mat[, 3] <- NA # ifelse(res.mat[, 2] < p, 1, 0)
    res.mat[, 4] <- NA # sign(res.mat[, 1]) * res.mat[, 3]
    res.mat[, 5] <- tmp[, 3]
    res.mat[, 6] <- NA # p.adjust(tmp[, 2], method = "BH")
    res.mat[is.na(res.mat[, 2]), c(2, 6)] <- NA
    colnames(res.mat) <-  c("coeff", "p.value", "keep", "pos.neg", "Warning", "FDR")
    return(list(res.mat = res.mat, rows.to.keep = to.keep))
}

cf.median.pred.survtime <- function(object, newdata) {
    ## Obtain the median predicted survival time for a cforest survival object
    ## with data newdata
    tmp <- treeresponse(object, newdata = newdata)
    return(sapply(tmp, function(x) {x$time[which(x$surv < 0.50)[1]]}))
}


cf.mean.survtime <- function(object, newdata) {
    tmp <- treeresponse(object, newdata = newdata)
    f1 <- function(times, survs) {
        ## note that cforest returns estimated KM, but we do not know the max surv time.
        ## so we need to truncate
        diftime <- c(times[1], diff(times))
        probs <- c(1, survs[-length(survs)])
        return(sum(diftime * probs))
    }
    return(sapply(tmp, function(x) f1(x$time, x$surv)))
}
         

my.cforest <- function(x, time, event, ngenes, newdata) {
    ## Does "everything":
    ##   - select genes
    ##   - fit cforest model
    ##   - obtain predictions
    ##   (in the future maybe variable importances, but time consuming and not used now)
    sobject <- Surv(time, event)
    selected.genes <- geneSelect(x, sobject, ngenes)
    x1 <- data.frame(x[, selected.genes$rows.to.keep])
    x1$time <- time
    x1$event <- event
    cf1 <- cforest(Surv(time, event) ~ ., data = x1,
                   control = cforest_classical(ntree = 500))
    if(!(is.null(newdata))) {
        if(!(is.data.frame(newdata))) newdata <- data.frame(newdata)
        pred.stime <- cf.mean.survtime(cf1, newdata)
    } else {
        pred.stime <- NULL
    }
    if(!(is.data.frame(x))) x <- data.frame(x)
    overfit.pred.stime <- cf.mean.survtime(cf1, data.frame(x))
    return(list(selected.genes.stats = selected.genes$res.mat,
                selected.genes.rows = selected.genes$rows.to.keep,
                selected.genes.names = colnames(x)[selected.genes$rows.to.keep],
                selected.genes.number = length(selected.genes$rows.to.keep),
                cforest_ob = cf1,
                predicted_surv_time = pred.stime,
                overfit_predicted_surv_time = overfit.pred.stime))
}



my.cforest.cv <- function(x, time, event, ngenes, nfold = 10,
                              universeSize = 10) {
### Take care of MPI stuff
###     if (mpi.comm.size(comm = 1) == 0) {
###         mpiSpawnAll(universeSize)
###         cat("\n      MPI:  cond 1 \n")
###     } else { ## so mpi is running
###         if ((mpi.comm.size(comm = 1) - 1) < universeSize) {
###             ## but few salves
###             mpi.close.Rslaves()
###             mpiSpawnAll(universeSize)
###             cat("\n     MPI:  cond 2 \n")
###         } else {
###             mpiDelete()
###             cat("\n     MPI:  cond 3 \n")
###         }
###     }
    n <- length(time)
    index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    OOB.scores <- rep(NA, n)

    my.cforest.internal.MPI <- function(fnum) {
        ## to be used with papply
        ## the following need to be passed with
        ##    mpi.bcast.Robj2slave
        ##    x, time, event, ,
        ##    MaxIterationsCox,
        ##    index.select
        cat("\n Doing fnum ", fnum, "\n") ## FIXME: debug
        xtrain <- x[index.select != fnum, , drop = FALSE]
        xtest <- x[index.select == fnum, , drop = FALSE]
        timetrain <- time[index.select != fnum]
        eventtrain <- event[index.select != fnum]
        
        retval <- my.cforest(xtrain, timetrain, eventtrain, ngenes, xtest)
        return(retval)
    }

    tmp1 <- papply(as.list(1:nfold),
                   my.cforest.internal.MPI,
                   papply_commondata = list(x = x,
                   time = time, event = event, 
                   ngenes = ngenes, index.select = index.select))
    
    for(i in 1:nfold) {
        OOB.scores[index.select == i] <-
            tmp1[[i]]$predicted_surv_time
        tmp1[[i]]$predicted_surv_time <- NULL  ## don't need this anymore
    }

    out <- list(cved.models = tmp1,
                OOB.scores = OOB.scores)
    return(out)
}


###############################################
##########                    #################
########## util functions  #################
##########                    #################
###############################################


print.selected.genes <- function(object,
                                 idtype, organism) {
    ## writes to a file, and python generates a bunch of HTML tables;
    ##    these are then linked from the main HTML; code is in runAndCheck.py
    
    p.values.original <- data.frame(Names = object$selected.genes.names,
                                    p.value = object$selected.genes.stats[, 2],
                                    coeff = object$selected.genes.stats[, 1], 
                                    abs.coeff = abs(object$selected.genes.stats[, 1]),
                                    fdr = NA,
                                    Warning = object$selected.genes.stats[, 5])

###     if (any(is.na(p.values.original))) {
###         p.values.original[is.na(p.values.original)] <- 999.999
###         cat("Oh, oh, some nas in p.values.originalll",
###             file = "nas.in.p.values.WARN")
###     }

    write.table(file = "p.values.coeffs.txt",
                p.values.original, row.names = FALSE, 
                col.names = TRUE,
                quote = FALSE,
                sep = "\t")
    system(paste("/http/signs/cgi/order.html.py", idtype, organism)) ## call python to generate pre-sorted HTML tables
}


print.cv.results <- function(cvobject, allDataObject, subjectNames, html.level = 3, html = TRUE,
                             outfile = "results.txt",
                             pals_main = "Selected.genes.txt",
                             pals_all  = "Selected.and.CV.selected.txt") {
    sink(file = outfile)
    if(html) {
        oobs <- matrix(cvobject$OOB.scores, ncol = 1)
        rownames(oobs) <- subjectNames
        cat("\n <hr><h2>", html.level, ". Cross-validation runs</h2>\n", sep = "")
        cat("\n<h3>", html.level, ".1. <a href=\"scores.oob.html\" target=\"scores_window\">View</a> out-of-bag scores.</h3>\n", sep ="")
        cleanHTMLhead(file = "scores.oob.html",
                      title = "Predictor scores for out-of-bag data")
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
        
    } else {
        
        cat("\n Out-of-bag scores\n\n")
        oobs <- matrix(cvobject$OOB.scores, ncol = 1)
        rownames(oobs) <- subjectNames
        print(oobs)
    }
    
    cvobject <- cvobject[[1]] ## don't need scores anymore. Simpler subsetting.
    
    ks <- length(cvobject)
    
    ngenes <- unlist(lapply(cvobject, function(x) x$selected.genes.number))
    cv.names <- paste("CV.run.", 1:ks, sep = "")    
    ngenes <- matrix(ngenes, ncol = 1)
    rownames(ngenes) <- cv.names
    
    if(html) {
        cat("\n\n <h3>", html.level,
            ".2 Number of selected genes in cross-validation runs</h3>\n", sep = "")
        print(ngenes)
        
        cat("\n\n <h3>", html.level,
            ".3 Genes selected in each of the cross-validation runs</h3>\n", sep = "")
        for(i in 1:ks) {
            cat("\n\n <h4>CV run ", i, "</h4>\n")
            cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
            cat("<tr align=left><th width=200>Gene</th> </tr>")
            thesegenes <- cvobject[[i]]$selected.genes.names
            for(thegene in thesegenes) {
                cat("\n<tr><td>", linkGene(thegene), "</td></tr>")
            }
            cat("\n </TABLE>")
        }
    } else {
        cat("\n\n Number of selected genes and parameters in cross-validation runs\n")
        cat("-------------------------------------------------------------------\n\n")
        print(ngenes)
        
        cat("\n\n Stability assessments \n")
        cat(    " ---------------------\n")
        cat("\n Genes selected in each of the cross-validation runs \n")
        
        for(i in 1:ks) {
            cat(paste("CV run  ", i, " (", ngenes[i], " genes selected):   ", sep = ""), "\n")
            print(cvobject[[i]]$selected.genes.names)
            cat("\n---\n")
        }
    }
    
    tmp.genesSelected <- list()
    tmp.genesSelected[[1]] <- allDataObject$selected.genes.names
    genesSelected.cv <- lapply(cvobject, function(x) x$selected.genes.names)
    tmp.genesSelected <- c(tmp.genesSelected, genesSelected.cv)

    shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 1))
    for(i in 1:(ks + 1)) {
        for(j in 1:(ks + 1)) { ## sure, need not be symmetric, but this is fast
            shared.genes[i, j] <-
                length(intersect(tmp.genesSelected[[i]],
                                 tmp.genesSelected[[j]]))
        }
    }

    ngenes <- c(allDataObject$selected.genes.number, ngenes)
    prop.shared <- round(shared.genes/ngenes, 3)
    
    ngenesS <- paste("(", ngenes, ")", sep = "")
    colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", cv.names)
    rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", cv.names), ngenesS)
    
    options(width = 200)
    
    html.level <- html.level + 1
    if(html) {
        cat("\n\n <h2>", html.level, ". Stability assessments</h2>\n", sep = "")
        cat("\n\n <h3>", html.level, ".1 Number of shared genes</h3> \n", sep = "")
    } else {
        cat("\n\n Stability assessments \n")
        cat("\n\n Number of shared genes \n")
    }
    print(as.table(shared.genes))

    if(html) {
        cat("\n\n <h3>", html.level, ".2 Proportion of shared genes (relative to row total)</h3> \n", sep = "")
    } else {
        cat("\n\n Proportion of shared genes (relative to row total) \n")
    }
    print(as.table(prop.shared))
    
    options(width = 80)
    unlisted.genes.selected <- unlist(genesSelected.cv)
    
    in.all.data <-
        which(names(table(unlisted.genes.selected, dnn = NULL)) %in% tmp.genesSelected[[1]])

    if(html) {
        tmptmp <- sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE)
        rntmptmp <- rownames(tmptmp)
        cat("\n\n\n<h3>", html.level, ".3 Gene freqs. in cross-validated runs of genes selected in model with all data</h3> \n\n", sep = "")
        
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=200>Gene</th> <th width=50>Frequency</th></tr>")
            for(ii in 1:length(tmptmp)) {
                cat("\n<tr><td>", linkGene(rntmptmp[ii]), "</td>",
                    "<td>", tmptmp[ii], "</td></tr>")
            }
            cat("\n </TABLE>")
    } else {
        cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
        print(sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE))
    }

    if(html) {
        tmptmp <- sort(table(unlisted.genes.selected, dnn = NULL), decreasing = TRUE)
        rntmptmp <- rownames(tmptmp)
        cat("\n\n\n<h3>", html.level, ".4 Gene freqs. in cross-validated runs</h3> \n\n", sep = "")
        
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=200>Gene</th> <th width=50>Frequency</th></tr>")
            for(ii in 1:length(tmptmp)) {
                cat("\n<tr><td>", linkGene(rntmptmp[ii]), "</td>",
                    "<td>", tmptmp[ii], "</td></tr>")
            }
            cat("\n </TABLE>")
    } else {
        cat("\n\n\n Gene freqs. in cross-validated runs \n\n")
        print(sort(table(unlisted.genes.selected, dnn = NULL), decreasing = TRUE))
    }
    sink()

    #####    PaLS   #####
    writeForPaLS(tmp.genesSelected, pals_main, pals_all)
}

writeForPaLS <- function(genesSelected, pals_main, pals_all) {
    #####    PaLS   #####
    ## genesSelected is a list;
    ## first element has genes from run on all data;
    ## remaining components are each of the CV runs
    ## Works for lists without additional structure, such as
    ##    returned by TGD, cforest, etc. Not for Dave et al.
    ##    approach which has signatures.
    
    for(nr in (1:length(genesSelected))) {
        if(nr == 1) {
            cat("#Run_on_all_data\n", file = pals_main)
            cat("#Run_on_all_data\n", file = pals_all)
            for(gene in genesSelected[[nr]]) {
                cat(gene, "\n", file = pals_main, append = TRUE)
                cat(gene, "\n", file = pals_all, append = TRUE)
            }
        } else {
            cat("#CVRun_", nr - 1, sep = "", file = pals_all, append = TRUE)
            for(gene in genesSelected[[nr]]) {
                cat(gene, "\n", file = pals_all, append = TRUE)
            }
        }
    }
}

print.validation.results <- function(object, html.level = 5,
                                     outfile = "results.txt") {
    sink(file = outfile, append = TRUE)
    cat("\n <h2>", html.level, ". Validation data</h2>\n", sep = "")
    cat("\n\n <h3>", html.level,
        ".1. <a href=\"scores.validation.html\" target=\"vscores_window\">View</a>",
        "the scores for validation data.</h3>\n", sep = "")
    cleanHTMLhead(file = "scores.validation.html",
                  title = "Scores for validation data")
    write(paste("<TABLE frame=\"box\">\n",
                "<tr><th>Validation subject/array</th> <th>Linear score</th></tr>\n"),
          file = "scores.validation.html", append = TRUE)
    wout <- ""
    valpred <- object$predicted_surv_time
    subjnames <- names(valpred)
    for(i in 1:length(valpred)) {
        wout <- paste(wout, "\n <tr align=right>",
                      "<td>", subjnames[i], "</td><td>", valpred[i], "</td></tr>\n")
    }
    wout <- paste(wout, "</TABLE>")
    write(wout, file = "scores.validation.html", append = TRUE)
    cleanHTMLtail(file = "scores.validation.html")
    sink()
}

kmplots <- function(cv.scores, overfitt.scores, Time, Event) {
    gdd.width <- png.width <- 480
    gdd.height <- png.height <- 410
    pdf(file = "kmplot-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1
    dev.off()
    pdf(file = "kmplot-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2
    dev.off()
    GDD(file = "kmplot-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1
    dev.off()
    GDD(file = "kmplot-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2
    dev.off()
    pdf(file = "kmplot4-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize4(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1.4
    dev.off()
    pdf(file = "kmplot4-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize4(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.4
    dev.off()
    GDD(file = "kmplot4-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize4(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1.4
    dev.off()
    GDD(file = "kmplot4-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize4(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.4
    dev.off()


    pdf(file = "kmplot3-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize3(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    pdf(file = "kmplot3-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize3(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()
    GDD(file = "kmplot3-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(cv.scores,  Time,
                 Event,  addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    GDD(file = "kmplot3-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(overfitt.scores,  Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()

}

kmplots.validation <- function(scores, validationTime, validationEvent) {
    valpred <- scores
    gdd.width <- png.width <- 480
    gdd.height <- png.height <- 410
    pdf(file = "kmplot-validation.pdf", width = png.width,
        height = png.height)
    KM.visualize(valpred, validationTime,
                 validationEvent, addmain = NULL)
    dev.off()
    GDD(file = "kmplot-validation.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize(valpred, validationTime,                         
                 validationEvent, addmain = NULL)
    dev.off()
    
    pdf(file = "kmplot4-validation.pdf", width = png.width,
        height = png.height)
    KM.visualize4(valpred, validationTime,
                  validationEvent, addmain = NULL)
    dev.off()
    GDD(file = "kmplot4-validation.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize4(valpred, validationTime,                         
                  validationEvent, addmain = NULL)
    dev.off()
    
    
    pdf(file = "kmplot3-validation.pdf", width = png.width,
        height = png.height)
    KM.visualize3(valpred, validationTime,
                  validationEvent, addmain = NULL)
    dev.off()
    GDD(file = "kmplot3-validation.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(valpred, validationTime,                         
                  validationEvent,  addmain = NULL)
    dev.off()
}



imClose <- function (im) {
    ## prevent all the "Closing PNG device ..."
    dev.off(im$Device)
}


#### Several major changes:

##   - no early stopping anymore for tauBestP: communication way too expensive.
##     look at repos if you want it.

##   - No usage of snow. Now only papply.


tgdSingle <- function(x, time, event, unique.thres = 1, epi = 5e-06,
                      maxiter = 5000,
                      nfold = 10, fitWithBest = TRUE){
    tauBestP(x, time, event, thres = c(unique.thres, unique.thres),
             epi, thresGrid = 1, maxiter,  nfold, fitWithBest)
}


## to check tgdSingle
## arrays <- 40
## genes <- 30
## covar <- dlbcl.covar[1:arrays, ]
## survtime <- dlbcl.surv[1:arrays]
## status <- dlbcl.status[1:arrays]
## covar <- covar[, 1:genes]
## uu <- tgdSingle(covar, survtime, status, unique.thres = 0.9, maxiter = 10000)



tauBestP <- function(x, time, event, thres = c(0, 1),
                     epi = 5e-06, thresGrid = 6, 
                     maxiter = 5000, 
                     nfold = 10,
                     fitWithBest = TRUE) {
    ## checkEvery is not used as such
##     checkEvery <- maxiter
    if((!is.vector(time)) | (!is.vector(event)))
        stop("Time and event should be vectors")
    n <- length(time)
    nfold <- min(nfold, n)
    thresS <- seq(from = thres[1], to = thres[2], length.out = thresGrid)
    cvpl.mat <- array(NA, c(length(thresS), maxiter))
    cvindex <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    clusterParams <- expand.grid(1:nfold, thresS)
    totalProcs <- dim(clusterParams)[1]

##     xTheCluster <<- x
##     tTheCluster <<- time
##     eTheCluster <<- event
##     epTheCluster <<- epi
##     mitTheCluster <<- maxiter
##     ceTheCluster <<- checkEvery
##     cviTheCluster <<- cvindex
##     clusterParams <<- clusterParams
    
##     time.initial.export <-
##         unix.time(clusterExport(TheCluster,
##                                 c("xTheCluster",
##                                   "tTheCluster",
##                                   "eTheCluster",
##                                   "epTheCluster",
##                                   "mitTheCluster",
##                                   "ceTheCluster",
##                                   "cviTheCluster",
##                                   "clusterParams")))
##     cat("\n\n     Export to slaves took ", time.initial.export[3],
##         " seconds\n")
##     ts1 <- unix.time(clusterApply(TheCluster, 1:60, tgdSnowSetUp))
##     cat("\n \n      Slaves setup took ", ts1[3], " seconds\n")
    
    clusterOutput <- list()
    t.r1 <-
        unix.time(
                  clusterOutput <-
                  papply(as.list(1:totalProcs),
                         nodeRun,
                         papply_commondata = list(
                         clusterParams = clusterParams,
                         x = x,
                         time = time,
                         event = event,
                         epi = epi,
                         steps = maxiter,
                         cvindex = cvindex)))

    tmp.cvpl <- matrix(unlist(clusterOutput),
                      ncol = maxiter, byrow = TRUE)

    cvpllymtgd <- matrix(NA, nrow = thresGrid, ncol = maxiter)
    ## Set cvpl to NULL to make the stored
    ## and later transmitted object as small as possible
    for(th in 1:thresGrid) {
        cvpllymtgd[th, ] <-
            apply(tmp.cvpl[clusterParams[, 2] == thresS[th], ], 2, sum)
    }        

    cvpl.mat <- cvpllymtgd/nfold

    ## have to allow for possibility of several minima; take largest
##    browser()
    tmp <-
        which(cvpl.mat == min(cvpl.mat, na.rm = TRUE), arr.ind = TRUE)
    tmp <- tmp[nrow(tmp), ]
    thresBest <- thresS[tmp[1]]
    stepBest <- tmp[2]

    if(is.na(stepBest)) {
        caughtUserError("An error occured, probably related to numerical problems. Please let us know")
    }

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



      
## tgdSnowSetUp <- function(index) {
##     assign("cvindex", cviTheCluster, env = .GlobalEnv)
##     assign("epi",  epTheCluster, env = .GlobalEnv)
##     assign("steps",  ceTheCluster, env = .GlobalEnv)
##     assign("maxiter",  mitTheCluster, env = .GlobalEnv)

##     assign("thres",  clusterParams[index, 2], env = .GlobalEnv)
##     assign("i",  clusterParams[index, 1], env = .GlobalEnv)
##     assign("x.train",
##            as.matrix(xTheCluster[cvindex != i, , drop = FALSE]),
##            env = .GlobalEnv)
##     assign("time.train", tTheCluster[cvindex != i], env = .GlobalEnv)
##     assign("event.train",  eTheCluster[cvindex != i], env = .GlobalEnv)
    
##     assign("x.test",
##            as.matrix(xTheCluster[cvindex == i, , drop = FALSE]), env = .GlobalEnv)
##     assign("time.test",  tTheCluster[cvindex == i], env = .GlobalEnv)
##     assign("event.test",  eTheCluster[cvindex == i], env = .GlobalEnv)
    
##     rm(list = c("xTheCluster", "tTheCluster", "eTheCluster",
##        "epTheCluster","mitTheCluster", "ceTheCluster",
##        "cviTheCluster", "clusterParams", "i"), envir = .GlobalEnv)
## }

nodeRun <- function(index) {
    i <- clusterParams[index, 1]
    thres <- clusterParams[index, 2]
    x.train <- as.matrix(x[cvindex != i, , drop = FALSE])
    x.test <- as.matrix(x[cvindex == i, , drop = FALSE])
    time.train <- time[cvindex != i]
    time.test <- time[cvindex == i]
    event.train <- event[cvindex != i]
    event.test <- event[cvindex == i]
    tgd1InternalSnow(x.train, time.train, event.train,
                     x.test, time.test, event.test,
                     thres, epi, steps)
}




tgd1InternalSnow <- function(x.train, time.train, event.train,
                             x.test, time.test, event.test,
                             thres,
                             epi,
                             steps){
##    on.exit(browser())
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
    n <- nm[1]
    m <- nm[2]
    nm1 <- dim(x.test)
    n1 <- nm1[1]
    m1 <- nm1[2]
    r <- rank(time.train)
    gradient <- rep(NA, n)
    scores <- rep(NA, n1)
    beta <- beta1 <- rep(0, m)
    beta <- as.matrix(beta)
    cvpl <- rep(0, steps)
    
##    cat(paste("I am using: threshold", thres))

    for(iteration in 1:steps) {
        beta <- beta1
        if(n1 == 1) {
            scores <- sum(x.test * as.vector(beta))
        }
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
    ## what are these for?
    ##     rm("beta1GlobalEnv", envir = .GlobalEnv)
    ##     assign("beta1GlobalEnv", beta1, envir = .GlobalEnv)

    ## recall we are using clusterApply
    ## and we count on the persistence for future visits

##    cat(paste("**** First run: spent ", proc.time()[3] - pt1, "seconds"))
    return(cvpl)
}


## tgdPieceInternalSnow <- function(runIt) {
## ### like tgd, but continues the iterations where the other left.
## ### we are using clusterApply, so we count on the exact positions in the cluster.

##     pt1 <- proc.time()[3]
##     if(runIt) {
##         cvpl <- rep(0, steps)
##         scores <- rep(0, n1)
##         gradient <- rep(0, n)
        
##         beta1 <- beta1GlobalEnv
## ##         warning(paste("I am using: threshold", thres,
## ##                       " runIt ", runIt))

##         for(iteration in 1:steps){
##             beta <- beta1
##             if(n1 == 1) scores <- sum(x.test * beta)
##             else scores <- x.test %*% beta
            
##             ita <- x.train %*% beta
            
##             epita <- exp(ita)
##             d <- rep(0, n)
##             dono <- rep(0, n)
##             for(i in 1:n) {
##                 d[i] <- sum(event.train[r == r[i]])
##                 dono[i] <- sum(epita[r >= r[i]])
##             }
##             risk <- d/dono
##             culrisk <- rep(0, n)
##             for(i in 1:n) {
##                 culrisk[i] <- sum(unique(risk[r <= r[i]]))
##             }
##             gradient <- event.train - epita * culrisk
##             gra1 <- crossprod(x.train, gradient)
##             gra1[abs(gra1) < thres * max(abs(gra1))] <- 0
##             beta1 <- beta + epi * gra1
            
            
##             lik1.train <- sum((ita - log(dono)) * event.train)
##             cvpl[iteration] <- (lik1.train -
##                                 lik1(c(ita, scores),
##                                      c(time.train,  time.test),
##                                      c(event.train, event.test)))/length(time.test)
##         }
##         rm("beta1GlobalEnv", envir = .GlobalEnv)
##         assign("beta1GlobalEnv", beta1, envir = .GlobalEnv)
##         warning(paste("**** Other runs: spent ", proc.time()[3] - pt1, "seconds"))
##         return(cvpl)
##     } else {return(rep(NA, steps))}
## }


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
                 maxiter = 5000, 
                 nfold = 10) {

    ##     if(is.null(nfold)) {
    ##         nfold <- length(y)
    ##     } ## I've never tested this. Should work for leave-one-out but ....
    
    ## we asume at least as many mpi Rslaves as nfold.
    ## eh? why?
##    if(nfold > (mpi.comm.size() - 1))
##        stop("nfold > number of mpi Rslaves")
    
    
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



plot.cvpl <- function(cvpl.mat, epi, thres = c(1, 1), thresGrid = 1) {
##plot.cvpl <- function(cvpl.mat, epi, thres = c(0, 1), thresGrid = 6) {
   
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
            y = t(cvpl.mat[, 1:m.step, drop = FALSE]), xlab = expression(nu),
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


plot.cvpl.single  <- function(cvpl.mat, epi, thres) {
   ## for a single threshold
## Diagnostic plot for the output from tauBest
   
    m.step <- max(apply(cvpl.mat, 1, function(x) sum(!is.na(x))))
    nus <- epi * (1:m.step)
    maxy <- max(cvpl.mat, na.rm = TRUE)
    miny <- min(cvpl.mat, na.rm = TRUE)

    plot(x = nus,
         y = t(cvpl.mat[, 1:m.step, drop = FALSE]),
         xlab = paste(expression(nu), " ( = Delta nu * number iterations)"),
         ylab = "Cross validated partial likelihood",
         type = "l", col = "black", lwd = 1.9,
         main = expression(paste("CV partial likelihood: effects of nu")))
}



summaryTGDrun <- function(x, time, event, z, epi, thres = c(0, 1),
                          thresGrid = 6, plot = TRUE,
                          genesOut = TRUE, outfile = "genes.all.out",
                          html = TRUE)  {
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

    tgdTrainPap <- function(varArgs)
        return(tgdTrain(xdatasn, timedatasn, eventdatasn,
                        varArgs[1], epidatasn, varArgs[2])[[2]])
    
    variableArgs <- list()
    for(tt in thres.do) variableArgs[[tt]] <- c(thresS[tt], mins.at[tt])

    bestBetas.m.pre <- papply(variableArgs,
                              tgdTrainPap,
                              papply_commondata = list(
                              xdatasn = x,
                              timedatasn = time,
                              eventdatasn = event,
                              epidatasn = epi))
##     bestBetas.m.pre <- clusterApply(TheCluster, variableArgs,
##                                     tgdTrainSnow, x, time, event, epi)

    
    
    
    bestBetas.m <- bestBetas.m.pre
    bestBetas.m[[z$thres.loc]] <- z$betas

    n.betas.nozero <- lapply(bestBetas.m, function(x) sum(x != 0))
    outm <- data.frame(cbind(Threshold = thresS,
                             Minimum.CVPL.at = mins.at,
                             CVPL = cvpls,
                             Number.non.zero.coefficients = n.betas.nozero))

    bb <- z$betas[z$betas != 0]
    names(bb) <- rownames(z$betas)[z$betas != 0]
    bb <- bb[order(abs(bb), decreasing = TRUE)]
    rnnbb <- names(bb)

    if(html) {
##        browser()
        cat("\n\n<h3> 4.1 Selected genes (ordered by decreasing value of their coefficient)</h3>\n")
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=200>Gene</th> <th width=80>Coefficient</th></tr>")
        for(iii in 1:length(bb)) {
            cat("\n <tr align=center>\n")
            cat("<td>", linkGene(rnnbb[iii]), "</td><td>",
                bb[iii],
                "</td></tr>")
        }
        cat("\n </TABLE>\n")
     
##         cat("\n\n\n Cross-validated partial likelihood and number of selected\n")
##         cat(" genes for different thresholds (tau), with delta nu (epi) =", epi, ".\n")
##         cat("\n ===============================================================\n\n")
##         print(outm)

        
    } else {       
        cat("\n\n Selected genes (ordered by decreasing value of their coefficient)\n")
        print(bb)
##         cat("\n\n\n Cross-validated partial likelihood and number of selected\n")
##         cat(" genes for different thresholds (tau), with delta nu (epi) =", epi, ".\n")
##         cat("\n ===============================================================\n\n")
##         print(outm)
    }
    if(plot) {

        GDD(file = "fstdgrun.png", w = png.width,
            h = png.height,
            ps = png.pointsize,
            type = "png")
        par(cex.axis = 0.75); par(cex.lab = 1.4); par(cex.main = 1.5)

        for(ip in 1:thresGrid) {
            plot(bestBetas.m[[ip]], xlab = "Gene index",
                 ylab = "Coefficient",
                 type = "l",
                 main = paste("Threshold (tau) =", thresS[ip]))
        }
        dev.off()

        pdf(file = "fstdgrun.pdf", width =  png.width,
               height =  png.height, onefile = FALSE)
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
                      maxiter,
                      nfold) {
    ## find best params by CV and predict on a new set.
    bestTrain <- tauBestP(x, time, event,
                          thres, epi, thresGrid,
                          maxiter, nfold)
    
    return(list(scoresTest = xtest %*% bestTrain$betas,
                betas = bestTrain$betas,
                nonZeroBetas = bestTrain$nonZeroBetas, 
                threshold = bestTrain$threshold,
                step = bestTrain$step))
}

summary.cvTGD <- function(object, allDataObject, subjectNames, html = TRUE,
                          pals_main = "Selected.genes.txt",
                          pals_all  = "Selected.and.CV.selected.txt") {

    if(html) {
        oobs <- matrix(object$OOB.scores, ncol = 1)
        rownames(oobs) <- subjectNames

        cat("\n<h3>4.1. <a href=\"scores.oob.html\" target=\"scores_window\">View</a> out-of-bag scores.</h3>\n")
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

    } else {
        
        cat("\n Out-of-bag scores\n\n")
        oobs <- matrix(object$OOB.scores, ncol = 1)
        rownames(oobs) <- subjectNames
        print(oobs)
    }
    
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


    if(html) {
        cat("\n\n <h3>4.2 Number of selected genes and parameters in cross-validation runs</h3>\n")
        print(tmp.mat)
        
        cat("\n\n <h3>4.3 Genes selected in each of the cross-validation runs</h3>\n")
        for(i in 1:ks) {
            cat("\n\n <h4>CV run ", i, "</h4>\n")
            cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
            cat("<tr align=left><th width=200>Gene</th> </tr>")
            thesegenes <- rownames(object[[i]]$betas)[object[[i]]$betas != 0]
            for(thegene in thesegenes) {
                cat("\n<tr><td>", linkGene(thegene), "</td></tr>")
            }
            cat("\n </TABLE>")
        }
    } else {
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


    if(html) {
        cat("\n\n <h2>5. Stability assessments</h2>\n")
        cat("\n\n <h3>5.1 Number of shared genes</h3> \n")
    } else {
        cat("\n\n Stability assessments \n")
        cat("\n\n Number of shared genes \n")
    }
    print(as.table(shared.genes))

    if(html) {
        cat("\n\n <h3>5.2 Proportion of shared genes (relative to row total)</h3> \n")
    } else {
        cat("\n\n Proportion of shared genes (relative to row total) \n")
    }
    print(as.table(prop.shared))
    
    options(width = 80)
    unlisted.genes.selected <- unlist(genesSelected.cv)
    
    in.all.data <-
        which(names(table(unlisted.genes.selected, dnn = NULL)) %in% tmp.genesSelected[[1]])

    if(html) {
        tmptmp <- sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE)
        rntmptmp <- rownames(tmptmp)
        cat("\n\n\n<h3> 5.3 Gene freqs. in cross-validated runs of genes selected in model with all data</h3> \n\n")
        
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=200>Gene</th> <th width=50>Frequency</th></tr>")
            for(ii in 1:length(tmptmp)) {
                cat("\n<tr><td>", linkGene(rntmptmp[ii]), "</td>",
                    "<td>", tmptmp[ii], "</td></tr>")
            }
            cat("\n </TABLE>")
    } else {
        cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
        print(sort(table(unlisted.genes.selected, dnn = NULL)[in.all.data], decreasing = TRUE))
    }

    if(html) {
        tmptmp <- sort(table(unlisted.genes.selected, dnn = NULL), decreasing = TRUE)
        rntmptmp <- rownames(tmptmp)
        cat("\n\n\n<h3> 5.4 Gene freqs. in cross-validated runs</h3> \n\n")
        
        cat("\n <TABLE  frame=\"box\" rules=\"groups\">\n")
        cat("<tr align=left><th width=200>Gene</th> <th width=50>Frequency</th></tr>")
            for(ii in 1:length(tmptmp)) {
                cat("\n<tr><td>", linkGene(rntmptmp[ii]), "</td>",
                    "<td>", tmptmp[ii], "</td></tr>")
            }
            cat("\n </TABLE>")
    } else {
        cat("\n\n\n Gene freqs. in cross-validated runs \n\n")
        print(sort(table(unlisted.genes.selected, dnn = NULL), decreasing = TRUE))
    }

    writeForPaLS(tmp.genesSelected, pals_main, pals_all)
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



#### At least in our set up, doing dStep1 sequentially is faster than
#### parallel. These are some timings:
####  The "4x" means with 4 times the number of genes (by rbind)

## data set          1 node           sequential          62 nodes
## ---------------------------------------------------------------------
## aml                4.97             4.791              6.66
## aml  4x           19.182           18.909             24.604
## breast             3.351                               3.236
## breast 4x         13.106                              12.896
## dlbcl              6.388            6.347             10.233 
## dlbcl 4x          25.208           24.452             40.619
## dlbcl 8x          50.327           49.303             78.267 


dStep1.serial <- function(x, time, event, p, MaxIterationsCox) {
    res.mat <- matrix(NA, nrow = ncol(x), ncol = 6)
    sobject <- Surv(time, event)
    ## cat("\n Starting dStep1.serial at ", date(), " \n\n"); ptm <- proc.time()
    funpap3 <- function (x) {
        out1 <-
            coxph.fit.pomelo0(x, sobject,
                              control = coxph.control(iter.max =  MaxIterationsCox))
        if(out1$warnStatus > 1) {
            return(c(0, NA, out1$warnStatus))
        } else {
            sts <- out1$coef/sqrt(out1$var)
            return(c(out1$coef, 1- pchisq((sts^2), df = 1), out1$warnStatus))
        }
    }
    tmp <- t(apply(x, 2, funpap3))
    res.mat[, 1:2] <- tmp[, 1:2]
    res.mat[, 3] <- ifelse(res.mat[, 2] < p, 1, 0)
    res.mat[, 4] <- sign(res.mat[, 1]) * res.mat[, 3]
    res.mat[, 5] <- tmp[, 3]
    res.mat[, 6] <- p.adjust(tmp[, 2], method = "BH")
    res.mat[is.na(res.mat[, 2]), c(2, 6)] <- 999
    colnames(res.mat) <-  c("coeff", "p.value", "keep", "pos.neg", "Warning", "FDR")
    ##  cat("\n Finished dStep1.serial at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")
    return(res.mat)
}





dStep1.parallel <- function(x, time, event, p, MaxIterationsCox) { 
    res.mat <- matrix(NA, nrow = ncol(x), ncol = 6)
    sobject <- Surv(time, event)
    cat("\n Starting dStep1.parallel at ", date(), " \n\n"); ptm <- proc.time()
    funpap3 <- function (x) {
        out1 <-
            coxph.fit.pomelo0(x, sobject,
                              control = coxph.control(iter.max = MaxIterationsCox))
        if(out1$warnStatus > 1) {
            return(c(0, NA,  out1$warnStatus))
        } else {
            sts <- out1$coef/sqrt(out1$var)
            return(c(out1$coef,
                     1- pchisq((sts^2), df = 1), 
                     out1$warnStatus))
        }
    }
    funpap4 <- function(xmat) apply(xmat, 2, funpap3)
    
    ## split data into the right number of groups for parallelization
    nparalGroups <- (mpi.comm.size(comm = 1) - 1)
    if(nparalGroups > 1) {
        paralGroups <- as.numeric(factor(cut(1:ncol(x),
                                             nparalGroups, labels = FALSE)))
    } else {
        nparalGroups <- 1
        paralGroups <- rep(1, ncol(x))
    }

    datalist <- list()
    for(ng in 1:nparalGroups)
        datalist[[ng]] <- x[, ng == paralGroups]
    
    tmp <- matrix(unlist(papply(datalist,
                                funpap4,
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
    cat("\n Finished dStep1.parallel at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")
    return(res.mat)
}



dStep1.parallel.old <- function(x, time, event, p, MaxIterationsCox) { 
    res.mat <- matrix(NA, nrow = ncol(x), ncol = 6)
    sobject <- Surv(time,event)
    cat("\n Starting dStep1.parallel at ", date(), " \n\n"); ptm <- proc.time()
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
    cat("\n Finished dStep1.parallel at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")
    return(res.mat)
}


##  2. Cluster; independently for those with pos and neg beta.
dStep2.old <- function(x, res.mat, maxSize, minSize,
                   minCor, plot,
                   interactive,
                   plotSizes = c(0.5, 1, 2)) {
    cat("\n       Starting dStep2 at ", date(), " \n\n"); ptm <- proc.time()
    
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
       
#####         pdokf <- function(alllabels = FALSE) {
#####             pos.labels <- rep("                   ", ncol(pos.data))
#####             index.labels <- which(pos.groups %in% pos.accept)
#####             pos.labels[index.labels] <- colnames(pos.data)[index.labels]

#####             pos.dend <- as.dendrogram(pos.clus, hang = 0.001)
#####             par(mar = c(5, 2, 1, 8))
#####             plot(pos.dend, horiz = TRUE, xlab = "1 - correlation", leaflab = "none",
#####                  main = "Positive coefficients")
#####             abline(v = 1 - minCor, lty = 2, col = "blue")

#####             dfp <- data.frame(name = unlist(dendrapply(rev(pos.dend), nameLeave)),
#####                               height = unlist(dendrapply(rev(pos.dend), heightLeave)))
#####             dfp$pos.gr <- pos.groups[order.dendrogram(rev(pos.dend))]
#####             dfp$chosen.clus <- dfp$pos.gr %in% pos.accept
#####             dfp$y <- nrow(dfp) - as.numeric(rownames(dfp)) + 1
#####             rainbow.col <- rainbow(length(posGroups))
#####             for(i in 1:length(posGroups)) {
#####                 dfpt <- dfp[dfp$pos.gr == posGroups[i], ]
#####                 miny <- min(dfpt$y)
#####                 maxy <- max(dfpt$y)
#####                 axis(4, line = 5, at = c(miny, maxy), col = rainbow.col[i],
#####                      tick = TRUE, labels = FALSE, lw = 3)
#####                 axis(4, line = 5, at = 0.5 * (miny + maxy),
#####                      col.axis = rainbow.col[i],
#####                      tick = FALSE, labels = posGroups[i], lw = 0,
#####                      cex.axis = 1.5)
#####                 for(j in 1:nrow(dfpt))
#####                     text(dfpt$name[j], x = dfpt$height[j],
#####                          y = dfpt$y[j], col = rainbow.col[i],
#####                          font = 2, pos = 4)
#####             }
#####             if (alllabels) {
#####                 dfpt <- dfp[dfp$chosen.clus == FALSE,]
##### 		if(nrow(dfpt)) {
#####                 for(j in 1:nrow(dfpt))
#####                     text(dfpt$name[j], x = dfpt$height[j],
#####                          y = dfpt$y[j], col = "black", cex = 0.8, pos = 4)
##### 			 }
#####             }
#####             return(dfp)
#####         } ##</pdokf within plotting>





        
        
#####         dendmappP <- function(factor = 1, alllabels = FALSE) {
#####             psf <- ifelse(factor < 1, ps * factor, ps)
#####             nameIm <- paste("dend.P.factor", factor, ".alllabels", alllabels, sep = "")
#####             im1 <- imagemap3(nameIm, height = height * factor, width = width * factor,
#####                              ps = psf)
#####             dfp <- pdokf(alllabels = alllabels)
            
#####             for(np in 1:nrow(dfp)) {
#####                 addRegion(im1) <- imRect(dfp[np, 2] + .030, dfp[np, 5] - 0.45,
#####                                          dfp[np, 2] - .1, dfp[np, 5] + 0.45,
#####                                          title = dfp[np, 1], alt = dfp[np, 1],
#####                                          href= linkGene2(dfp[np, 1]))
#####             }
            
#####             createIM(im1, file = paste("dend.P.factor", factor, ".alllabels",
#####                           alllabels, ".html", sep = ""))
#####             imClose(im1)
#####         }
        
            
        if (pdok) {
            for (plsz in plotSizes) {
                dendmappP(factor = plsz, alllabels = FALSE)
            }
             for (plsz in plotSizes) {
                 dendmappP(factor = plsz, alllabels = TRUE)
             }
            
        } else {
            system("touch NoPositiveCluster")
        }




        

#####         pnokf <- function(alllabels = FALSE) {
#####             neg.labels <- rep("                   ", ncol(neg.data))
#####             index.labels <- which(neg.groups %in% neg.accept)
#####             neg.labels[index.labels] <- colnames(neg.data)[index.labels]

#####             neg.dend <- as.dendrogram(neg.clus, hang = 0.001)
#####             par(mar = c(5, 2, 1, 8))
#####             plot(neg.dend, horiz = TRUE, xlab = "1 - correlation", leaflab = "none",
#####                  main = "Negative coefficients")
#####             abline(v = 1 - minCor, lty = 2, col = "blue")

#####             dfp <- data.frame(name = unlist(dendrapply(rev(neg.dend), nameLeave)),
#####                               height = unlist(dendrapply(rev(neg.dend), heightLeave)))
#####             dfp$neg.gr <- neg.groups[order.dendrogram(rev(neg.dend))]
#####             dfp$chosen.clus <- dfp$neg.gr %in% neg.accept
#####             dfp$y <- nrow(dfp) - as.numeric(rownames(dfp)) + 1
#####             rainbow.col <- rainbow(length(negGroups))
#####             for(i in 1:length(negGroups)) {
#####                 dfpt <- dfp[dfp$neg.gr == negGroups[i], ]
#####                 miny <- min(dfpt$y)
#####                 maxy <- max(dfpt$y)
#####                 axis(4, line = 5, at = c(miny, maxy), col = rainbow.col[i],
#####                      tick = TRUE, labels = FALSE, lw = 3)
#####                 axis(4, line = 5, at = 0.5 * (miny + maxy),
#####                      col.axis = rainbow.col[i],
#####                      tick = FALSE, labels = negGroups[i], lw = 0,
#####                      cex.axis = 1.5)
#####                 for(j in 1:nrow(dfpt))
#####                     text(dfpt$name[j], x = dfpt$height[j],
#####                          y = dfpt$y[j], col = rainbow.col[i],
#####                          font = 2, pos = 4)
#####             }
#####             if (alllabels) {
#####                 dfpt <- dfp[dfp$chosen.clus == FALSE,]
##### 		if(nrow(dfpt)) {
#####                 for(j in 1:nrow(dfpt))
#####                     text(dfpt$name[j], x = dfpt$height[j],
#####                          y = dfpt$y[j], col = "black", cex = 0.8, pos = 4)
##### 			 }
#####             }
#####             return(dfp)
#####         } ##</pnokf within plotting>

#####         dendmappN <- function(factor = 1, alllabels = FALSE) {
#####             psf <- ifelse(factor < 1, ps * factor, ps)
#####             nameIm <- paste("dend.N.factor", factor, ".alllabels", alllabels, sep = "")
#####             im1 <- imagemap3(nameIm, height = height * factor, width = width * factor,
#####                              ps = psf)
#####             dfp <- pnokf(alllabels = alllabels)
            
#####             for(np in 1:nrow(dfp)) {
#####                 addRegion(im1) <- imRect(dfp[np, 2] + .030, dfp[np, 5] - 0.45,
#####                                          dfp[np, 2] - .1, dfp[np, 5] + 0.45,
#####                                          title = dfp[np, 1], alt = dfp[np, 1],
#####                                          href= linkGene2(dfp[np, 1]))
#####             }
            
#####             createIM(im1, file = paste("dend.N.factor", factor, ".alllabels",
#####                           alllabels, ".html", sep = ""))
#####             imClose(im1)
#####         }

        if (pnok) {
            for (plsz in plotSizes) {
                dendmappN(factor = plsz, alllabels = FALSE)
            }
            for (plsz in plotSizes) {
                dendmappN(factor = plsz, alllabels = TRUE)
            }
        } else {
            system("touch NoNegativeCluster")
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
    cat("\n Finished dStep2 at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")    
    return(list(md = cbind(posMeanData, negMeanData),
                filteredGroupsPositive = groupsPositive,
                filteredGroupsNegative = groupsNegative,
                filteredPosPositions = filteredPosPositions,
                filteredNegPositions = filteredNegPositions,
                posPositions = posPositions,
                negPositions = negPositions))
}



dStep2 <- function(x, res.mat, maxSize, minSize,
                   minCor, plot,
                   interactive,
                   plotSizes = c(0.5, 1, 2)) {
    cat("\n       Starting dStep2 at ", date(), " \n\n"); ptm <- proc.time()
    
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
        if(pdok & (! pnok)) {
            system("touch NoNegativeCluster")

            datalist <- list()
            plotSizes2 <- c(plotSizes, plotSizes)
            aLabs <- c(rep(TRUE, length(plotSizes)),
                       rep(FALSE, length(plotSizes)))
            for(jj in 1:length(plotSizes2)) {
                datalist[[jj]] <- list()
                datalist[[jj]]$dir <- getwd()
                datalist[[jj]]$factor <- plotSizes2[jj]
                datalist[[jj]]$alllabels <- aLabs[jj]
                datalist[[jj]]$pn.groups <- pos.groups
                datalist[[jj]]$pn.data <- pos.data
                datalist[[jj]]$pn.accept <- pos.accept
                datalist[[jj]]$pn.clus <- pos.clus
                datalist[[jj]]$main <- "Positive coefficients"
                datalist[[jj]]$pnGroups <- posGroups
                datalist[[jj]]$theName <- "dend.P.factor"
                datalist[[jj]]$minCor <- minCor
            }
            tmp <- papply(datalist, function(z) wrapDendmapp(z))

        } else if ((! pdok) & pnok) {
            system("touch NoPositiveCluster")

            datalist <- list()
            plotSizes2 <- c(plotSizes, plotSizes)
            aLabs <- c(rep(TRUE, length(plotSizes)),
                       rep(FALSE, length(plotSizes)))
            for(jj in 1:length(plotSizes2)) {
                datalist[[jj]] <- list()
                datalist[[jj]]$dir <- getwd()
                datalist[[jj]]$factor <- plotSizes2[jj]
                datalist[[jj]]$alllabels <- aLabs[jj]
                datalist[[jj]]$pn.groups <- neg.groups
                datalist[[jj]]$pn.data <- neg.data
                datalist[[jj]]$pn.accept <- neg.accept
                datalist[[jj]]$pn.clus <- neg.clus
                datalist[[jj]]$main <- "Negative coefficients"
                datalist[[jj]]$pnGroups <- negGroups
                datalist[[jj]]$theName <- "dend.N.factor"
                datalist[[jj]]$minCor <- minCor
            }
            tmp <- papply(datalist, function(z) wrapDendmapp(z))

            
        } else if (pdok & pnok) {
            datalist <- list()
            plotSizes2 <- c(plotSizes, plotSizes)
            aLabs <- c(rep(TRUE, length(plotSizes)),
                       rep(FALSE, length(plotSizes)))
            for(jj in 1:length(plotSizes2)) {
                datalist[[jj]] <- list()
                datalist[[jj]]$dir <- getwd()
                datalist[[jj]]$factor <- plotSizes2[jj]
                datalist[[jj]]$alllabels <- aLabs[jj]
                datalist[[jj]]$pn.groups <- neg.groups
                datalist[[jj]]$pn.data <- neg.data
                datalist[[jj]]$pn.accept <- neg.accept
                datalist[[jj]]$pn.clus <- neg.clus
                datalist[[jj]]$main <- "Negative coefficients"
                datalist[[jj]]$pnGroups <- negGroups
                datalist[[jj]]$theName <- "dend.N.factor"
                datalist[[jj]]$minCor <- minCor

            }
            plotSizes2 <- c(plotSizes, plotSizes)
            aLabs <- c(rep(TRUE, length(plotSizes)),
                       rep(FALSE, length(plotSizes)))
            for(jj in 1:length(plotSizes2)) {
                datalist[[jj + length(plotSizes2)]] <- list()
                datalist[[jj + length(plotSizes2)]]$dir <- getwd()
                datalist[[jj + length(plotSizes2)]]$factor <- plotSizes2[jj]
                datalist[[jj + length(plotSizes2)]]$alllabels <- aLabs[jj]
                datalist[[jj + length(plotSizes2)]]$pn.groups <- pos.groups
                datalist[[jj + length(plotSizes2)]]$pn.data <- pos.data
                datalist[[jj + length(plotSizes2)]]$pn.accept <- pos.accept
                datalist[[jj + length(plotSizes2)]]$pn.clus <- pos.clus
                datalist[[jj + length(plotSizes2)]]$main <- "Positive coefficients"
                datalist[[jj + length(plotSizes2)]]$pnGroups <- posGroups
                datalist[[jj + length(plotSizes2)]]$theName <- "dend.P.factor"
                datalist[[jj]]$minCor <- minCor
           }

            tmp <- papply(datalist, function(z) wrapDendmapp(z))
           
        } else {
            stop("We should never get here!!! Plot error ")
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
    cat("\n Finished dStep2 at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")    
    return(list(md = cbind(posMeanData, negMeanData),
                filteredGroupsPositive = groupsPositive,
                filteredGroupsNegative = groupsNegative,
                filteredPosPositions = filteredPosPositions,
                filteredNegPositions = filteredNegPositions,
                posPositions = posPositions,
                negPositions = negPositions))
}



pdnokf <- function(pn.groups,  pn.data, pn.accept, pn.clus,
                   main, pnGroups, minCor,
                   alllabels = FALSE) {
    pn.labels <- rep("                   ", ncol(pn.data))
    index.labels <- which(pn.groups %in% pn.accept)
    pn.labels[index.labels] <- colnames(pn.data)[index.labels]
    pn.dend <- as.dendrogram(pn.clus, hang = 0.001)
    par(mar = c(5, 2, 1, 8))
    plot(pn.dend, horiz = TRUE, xlab = "1 - correlation", leaflab = "none",
         main = main)
    abline(v = 1 - minCor, lty = 2, col = "blue")
    dfp <- data.frame(name = unlist(dendrapply(rev(pn.dend), nameLeave)),
                      height = unlist(dendrapply(rev(pn.dend), heightLeave)))
    dfp$pn.gr <- pn.groups[order.dendrogram(rev(pn.dend))]
    dfp$chosen.clus <- dfp$pn.gr %in% pn.accept
    dfp$y <- nrow(dfp) - as.numeric(rownames(dfp)) + 1
    rainbow.col <- rainbow(length(pnGroups))
    cat("\n ... main is ", main, "\n")
    for(i in 1:length(pnGroups)) {
        dfpt <- dfp[dfp$pn.gr == pnGroups[i], ]
        miny <- min(dfpt$y) ## by not setting na.rm = TRUE we would bomb if wrong set
        maxy <- max(dfpt$y)
        cat("\n ......... pnGroup ", i,"\n")
        axis(4, line = 5, at = c(miny, maxy), col = rainbow.col[i],
             tick = TRUE, labels = FALSE, lw = 3)
        axis(4, line = 5, at = 0.5 * (miny + maxy),
             col.axis = rainbow.col[i],
             tick = FALSE, labels = pnGroups[i], lw = 0,
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


dendmapp <- function(factor, alllabels, 
                     pn.groups,  pn.data, pn.accept, pn.clus,
                     main, pnGroups, minCor, theName = "dend.P.factor") {
    ps <- 12
    height <- 1200
    width <- 800
    psf <- ifelse(factor < 1, ps * factor, ps)
    nameIm <- paste(theName, factor, ".alllabels", alllabels, sep = "")
    im1 <- imagemap3(nameIm, height = height * factor, width = width * factor,
                     ps = psf)
    dfp <- pdnokf(pn.groups = pn.groups,
                  pn.data = pn.data,
                  pn.accept = pn.accept,
                  pn.clus = pn.clus,
                  main = main,
                  pnGroups = pnGroups,
                  minCor   = minCor,
                  alllabels = alllabels)
    
    for(np in 1:nrow(dfp)) {
        addRegion(im1) <- imRect(dfp[np, 2] + .030, dfp[np, 5] - 0.45,
                                 dfp[np, 2] - .1, dfp[np, 5] + 0.45,
                                 title = dfp[np, 1], alt = dfp[np, 1],
                                 href= linkGene2(dfp[np, 1]))
    }
    
    createIM(im1, file = paste(theName, factor, ".alllabels",
                  alllabels, ".html", sep = ""))
    imClose(im1)
}


wrapDendmapp <- function(li) {
    ## for papply usage
    setwd(li$dir)
    dendmapp(li$factor, li$alllabels, li$pn.groups, li$pn.data, li$pn.accept,
             li$pn.clus, li$main, li$pnGroups, li$minCor, li$theName)
}


## 3. fit model

tryCatch2 <- function(x) tryCatch(x, error = function(e) "Error")

dStep3 <- function(res2, time, event, MaxIterationsCox) {
    
    ## this is a terrible hack, but I am getting scoping problems in stepAIC
    assign("..___MaxIterationsCox", MaxIterationsCox, env = .GlobalEnv)
    
    cat("\n ..... Starting dStep3 at ", date(), " \n\n"); ptm <- proc.time()
    if(all(is.na(res2))) return(NA)
    md <- res2$md
    sobject <- Surv(time, event)
    
    if(ncol(md) >= 2) {
        
        ## Select best two-genes model
        ncolmd <- ncol(md)
        modelsSizeTwo <- t(combn(1:ncolmd, 2))
        logliks <- rep(NA, nrow(modelsSizeTwo))
        cat("\n                  Number of models of size two: ", nrow(modelsSizeTwo), "\n")
        for(i in 1:nrow(modelsSizeTwo)) { ## could parallelize here??
            ## This is the logic: in some models, if we use a very large num.
            ## iterations, some scores end up being infinite. That sometimes
            ## can be "fixed" if we use a smaller number of iterations.
            ## If that does not work, then we set the value to NA
            trycox <- tryCatch(logliks[i] <-
                          coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                control = coxph.control(iter.max = MaxIterationsCox))$loglik[2],
                          error = function(e) "Error")
            cat("\n              tried iter.max = ", MaxIterationsCox)

            if((length(trycox) == 1) & (trycox == "Error")) {
                cat("\n              failed with iter.max = ", MaxIterationsCox)

                trycox <- tryCatch2(logliks[i] <-
                              coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                    control = coxph.control(iter.max = 20))$loglik[2])
                cat("\n              tried iter.max = ", 20)

            }
            if((length(trycox) == 1) & (trycox == "Error")) {
                trycox <- tryCatch2(logliks[i] <-
                              coxph(sobject ~ md[, modelsSizeTwo[i, ]],
                                    control = coxph.control(iter.max = 10))$loglik[2])
                cat("\n              tried iter.max = ", 10)
            }
            if((length(trycox) == 1) & (trycox == "Error"))
                logliks[i] <- NA
        }
        ## zz: give warnings if any loglik is NA
        bestTwoGenes <- modelsSizeTwo[which.max(logliks), ]
        ## stepAIC
        mdf <- data.frame(md)
        attach(mdf) ## stepAIC not working otherwise
        thisenv <- environment()
        cat("\n                 Fitting bestTwoModel\n")
        trycox <- tryCatch2(                      
                      bestTwoModel <-
                      coxph(eval(parse(text = paste("sobject ~",
                                       paste(colnames(mdf)[bestTwoGenes],
                                             sep = "", collapse = " + ")))),
                            control = coxph.control(iter.max = ..___MaxIterationsCox))
                      )
        cat("\n              tried with iter.max = ", ..___MaxIterationsCox)

        if((length(trycox) == 1) & (trycox == "Error")) {
            cat("\n              failed with iter.max = ", ..___MaxIterationsCox)
            trycox <- tryCatch2(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 20))
                          )
            cat("\n              tried iter.max = ", 20)

        }
        if((length(trycox) == 1) & (trycox == "Error")) {
            trycox <- tryCatch2(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 10))
                          )
            cat("\n              tried iter.max = ", 10)
        }
        if((length(trycox) == 1) & (trycox == "Error")) {
            trycox <- tryCatch2(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 5))
                          )
            cat("\n              tried iter.max = ", 5)
        }
        if((length(trycox) == 1) & (trycox == "Error")) {
            trycox <- tryCatch2(
                          bestTwoModel <-
                          coxph(eval(parse(text = paste("sobject ~",
                                           paste(colnames(mdf)[bestTwoGenes],
                                                 sep = "", collapse = " + ")))),
                                control = coxph.control(iter.max = 2))
                          )
            cat("\n              tried iter.max = ", 2)
        }
        ## So we'll try stepwise; but if it doesn't work,
        ## we use forwards, with max size model  = num.variables.
        ## We could run into trouble, so we allow to decrease size
        ## of model, until we are back to bestTwoModel

        cat("\n                   Fitting finalModel\n")
        upperScope <- paste("~", paste(colnames(mdf), sep = "",
                                       collapse = " + "))
        tryaic <- tryCatch2(
                      finalModel <- stepAIC(bestTwoModel, scope = upperScope,
                                            direction = "both")
                      )
        maxsize <- sum(event) - 2
        while((length(tryaic) == 1) & (tryaic == "Error")) {
            cat("\n                stepwise failed; trying forward. Maxsize = ", maxsize, "\n")
            tryaic <-
                tryCatch2(finalModel <- stepAIC(bestTwoModel,
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
    } else { ## so we do not need to look over several models
        cat("\n                    Not looking over several models\n")
        mdf <- data.frame(md)

        for (mitercox in c(MaxIterationsCox, 20, 10, 5, 2)) {
            trycox <-
                tryCatch2(
                    finalModel <- coxph(sobject ~ ., data = mdf,
                                        control = coxph.control(iter.max = mitercox)))
            if((length(trycox) == 1) & (trycox == "Error")) break
            cat("\n                         failed with iter.max = ", mitercox)
        }
##         if(class(trycox) == "try-error")
##             trycox <-
##             try(
##                 finalModel <- coxph(sobject ~ ., data = mdf,
##                                     control = coxph.control(iter.max = 20)))
##         if(class(trycox) == "try-error")
##           trycox <-
##             try(
##                 finalModel <- coxph(sobject ~ ., data = mdf,
##                                     control = coxph.control(iter.max = 10)))
##         if(class(trycox) == "try-error")
##           trycox <-
##             try(
##                 finalModel <- coxph(sobject ~ ., data = mdf,
##                                     control = coxph.control(iter.max = 5)))
##         if(class(trycox) == "try-error")
##           trycox <-
##             try(
##                 finalModel <- coxph(sobject ~ ., data = mdf,
##                                     control = coxph.control(iter.max = 2)))
        
        if((length(trycox) == 1) & (trycox == "Error"))
          predictsFinalModel <- rep(NA, length(time))
        else
          predictsFinalModel <- predict(finalModel, type = "lp")
    }
    
    out <- list(model = finalModel, scores = predictsFinalModel,
                clusterResults = res2)
    class(out) <- "fmDave"
    cat("\n ..... Finishing dStep3 at ", date(), "; took ", (proc.time() - ptm)[3], " \n\n")
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
    cat("\n Starting fitDave.res1Given at ", date(), " \n\n"); ptm <- proc.time()
    res2 <- dStep2(x, res1, maxSize, minSize, minCor, plot, interactive)
    res3 <- dStep3(res2, time, event, MaxIterationsCox)
    cat("\n Ended fitDave.res1Given at ", date(), "; took ", (proc.time() - ptm)[3], " seconds \n\n")
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

mpiSpawnAll <- function(numslaves = NULL) {
    if (is.null(numslaves)) numslaves <- mpi.universe.size()
    mpi.spawn.Rslaves(nslaves = numslaves)
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
    mpi.remote.exec(library(GDD))
    mpi.remote.exec(library(R2HTML))
    mpi.remote.exec(library(party))
    mpi.remote.exec(library(mboost))
#    Mpi.remote.exec(library(imagemap))
}

mpiDelete <- function() {
    mpi.remote.exec(rm(list = ls(env = .GlobalEnv), envir =.GlobalEnv))
}




cvDave.parallel3 <- function(x, time, event,
                             p, maxSize,
                             minSize, minCor,
                             MaxIterationsCox,
                             nfold,
                             universeSize = 10) {
    cat("\n Starting cvDave.parallel3 at ", date(), " \n\n"); ptm <- proc.time()
    
    if (mpi.comm.size(comm = 1) == 0) {
        mpiSpawnAll(universeSize)
        cat("\n      cvDave.parallel3:  cond 1 \n")
    } else { ## so mpi is running
        if ((mpi.comm.size(comm = 1) - 1) < universeSize) {
            ## but few salves
            mpi.close.Rslaves()
            mpiSpawnAll(universeSize)
            cat("\n     cvDave.parallel3:  cond 2 \n")
        } else {
            mpiDelete()
            cat("\n     cvDave.parallel3:  cond 3 \n")
        }
    }
       
    n <- length(time)
    index.select <- sample(rep(1:nfold, length = n), n, replace = FALSE)
    OOB.scores <- rep(NA, n)

    
    cat("\n\n Computing gene-wise cox p-value\n")
    ##    res1s <- list()
    ##     for(i in 1:nfold) {
    ##         cat("\n  ....  fold ", i)
    ##         xtr <- x[index.select != i, , drop = FALSE]
    ##         ttr <- time[index.select != i]
    ##         etr <- event[index.select != i]
    ##         res1s[[i]] <- dStep1.parallel(xtr, ttr, etr, p = p,
    ##                                       MaxIterationsCox = MaxIterationsCox)
    ##     }
    f00 <- function(i) {
        xtr <- x[index.select != i, , drop = FALSE]
        ttr <- time[index.select != i]
        etr <- event[index.select != i]
        return(dStep1.serial(xtr, ttr, etr, p = p,
                               MaxIterationsCox = MaxIterationsCox))
    }
    res1s <- papply(as.list(1:nfold),
                    f00,
                    papply_commondata = list(
                    x = x,
                    time = time,
                    event = event,
                    p = p,
                    MaxIterationsCox = MaxIterationsCox,
                    index.select = index.select))
        
    cat("\n\n Cleaning up MPI slaves\n\n")	
    mpiDelete()

##     if(nfold > (mpi.comm.size() - 1))
##         stop(paste("nfold > number of mpi Rslaves; mpi.comm.size = ",
##                    mpi.comm.size() ))
    
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
##     ##    global.res1s <<- res1s
##     ##    global.index.select <<- index.select
##     mpi.bcast.Robj2slave(index.select)
##     ##    mpi.bcast.cmd(foldNumber <- mpi.comm.rank())
##     ##    mpi.bcast.Robj2slave(DaveCVPred.res1Given.InternalMPI)

##     ## We use papply again
    cat("\n\n Computing the rest\n")


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

    tmp1 <- papply(as.list(1:nfold),
                   DaveCVPred.res1Given.InternalMPI2,
                   papply_commondata = list(x = x,
                   time = time, event = event, p = p,
                   maxSize = maxSize, index.select = index.select,
                   minSize = minSize, minCor = minCor,
                   MaxIterationsCox = MaxIterationsCox,
                   res1s = res1s))
    ##    cat("\n\n Cleaning up and closing MPI\n")
    ##    try(mpi.close.Rslaves())
    
    for(i in 1:nfold) {
        OOB.scores[index.select == i] <-
            tmp1[[i]]$scoresTest
        tmp1[[i]]$scoresTest <- NULL  ## don't need this anymore
    }
    
    out <- list(cved.models = tmp1,
                OOB.scores = OOB.scores)
    class(out) <- "cvDave"
    cat("\n Ending cvDave.parallel3 at ", date(), "; it took ", (proc.time() - ptm)[3], " \n\n")
    ## For debugging, mainly of lam logs
    system("tar -czvf lam.logs.tar.gz *.log")
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
