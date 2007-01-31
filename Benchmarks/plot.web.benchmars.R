## First, run wrapWebBenchmarks.py
## That will generate two 

rr <- function(x) rep(x, x)

designMat <- expand.grid(sample = 1:5,
                         users = unlist(sapply(c(1, 2, 10, 20), rr) ))

designMat <- designMat[order(designMat[, 2], designMat[, 1]), ]

breast.results <- scan("breast.web.bnchmk.txt", sep ="\t", what = double(0))
dlbcl.results <- scan("dlbcl.web.bnchmk.txt", sep ="\t", what = double(0))


breast.means <-
    as.vector(tapply(breast.results, list(designMat[, 1], designMat[, 2]), mean))

dlbcl.means <-
    as.vector(tapply(dlbcl.results, list(designMat[, 1], designMat[, 2]), mean))

mdm <- expand.grid(nrep = 1:5, narray = c(1, 2, 10, 20))

par(mfrow = c(1, 2))
par(cex.lab = 1.5); par(cex.axis = 1.5); par(cex.main = 1.5)
boxplot(breast.means ~ mdm[, 2], xlab = "Number of simultaneous users",
        ylab = "User wall time", main = "Breast data set")
boxplot(dlbcl.means ~ mdm[, 2], xlab = "Number of simultaneous users",
        ylab = "User wall time", main = "DLBCL data set")
