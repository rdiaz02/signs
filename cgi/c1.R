library(survival)

coxph.fit.pomelo0 <- function (x, y, init = NULL,
                              control, method = "effron",  rownames = NULL) {
    warnStatus <- 0
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


#####  zzz: meter en signs2.R, en la libreríoa.
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
                             ncol = 2, byrow = TRUE)
    res.mat[, 1:2] <- tmp[, 1:2]
    res.mat[, 3] <- ifelse(res.mat[, 2] < p, 1, 0)
    res.mat[, 4] <- sign(res.mat[, 1]) * res.mat[, 3]
    res.mat[, 5] <- tmp[, 3]
    res.mat[, 6] <- p.adjust(tmp[, 2], method = "BH")
    res.mat[is.na(res.mat[, 2]), c(2, 6)] <- 999
    colnames(res.mat) <- c("coeff", "p.value", "keep", "pos.neg", "Warning", "FDR")
    return(res.mat)
}










## what follows is for pomelo
cox.parallel <- function(x, time, event, p, MaxIterationsCox = 200) { 
    ## Painfully slow!!!! use C++ code 
    ## We use papply to parallelize; but what will happen when slaves run this???
    res.mat <- matrix(NA, nrow = ncol(x), ncol = 6)
    sobject <- Surv(time,event)
    
    funpap3 <- function (x) {
        out1 <- coxph.fit.pomelo0(x, sobject, control = coxph.control(iter.max = MaxIterationsCox))
        if(out1$warnStatus > 1) {
            return(c(NA, NA, out1$warnStatus))
        } else {
            sts <- out1$coef/sqrt(out1$var)
            return(c(out1$coef,
                     1- pchisq((sts^2), df = 1), 
                     out1$warnStatus))
        }
    }
    
    tmp <- matrix(unlist(papply(as.data.frame(x),
                                funpap3,
                                papply_commondata =list(sobject = sobject))),
                  ncol = 2, byrow = TRUE)
    res.mat[, 1:3] <- tmp
    res.mat[, 4] <- p.adjust(tmp[, 2], method = "BH")
    colnames(res.mat) <- c("coeff", "p.value", "Warning", "FDR")
    return(res.mat)
}

rescox <- cox.paralel(x, time, event, p, MaxIterationsCox = 200)  
p.values.original <- rescox[, c(2, 4, 1)]
p.values.original <- data.frame(Row = 1:length(geneNames),
                                Names = geneNames, p.values.original, 
                                abs.coeff = abs(p.values.original[, 2]),
                                Warning = all.res1[, 3])

p.values.original <-
    p.values.original[order(p.values.original$FDR, -p.values.original$abs.coeff), ]
cat(rep("\n", 13), file = "multest_parallel.res")
    ## we have: Name, p.value, adj.p.value, coef, abs.coef, warnings
write.table(file = "multest_parallel.res",
            p.values.original, row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t",
            append = TRUE)

## how exactly shoudl we do this in Pomelo II? The C++ produces a table that is then read.
## Can we just read that? Yes: And then, modify generate.table
