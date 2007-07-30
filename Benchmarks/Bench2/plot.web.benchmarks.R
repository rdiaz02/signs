## First, run wrapWebBenchmarks1.py, which calls fl-run-test to do the timings.
## The script rb.sh does that. However, note that for every run we get one value
## for one user, and 20 for 20 users. So I did a few runs by hand, and placed them
## in file 1.1.1

nusers <- c(rep(10, 10), rep(15, 15), rep(1, 5), rep(5, 5))


### TGD

d1 <- scan("breast.web.tgd.10_1.txt", sep = "\t",
               what = double(0), 10)
d2 <- scan("breast.web.tgd.15_1.txt", sep = "\t",
               what = double(0), 15)
d3 <- scan("breastA.web.tgd.txt", sep = "\t",
               what = double(0), 10)
breast.tgd <- c(d1, d2, d3)

d1 <- scan("dlbcl.web.tgd.10_1.txt", sep = "\t",
               what = double(0), 10)
d2 <- scan("dlbcl.web.tgd.15_1.txt", sep = "\t",
               what = double(0), 15)
d3 <- scan("dlbclA.web.tgd.txt", sep = "\t",
               what = double(0), 10)
dlbcl.tgd <- c(d1, d2, d3)



### FCMS

d1 <- scan("../dlbcl.web.tests.reedition.10_1.txt", sep = "\t",
               what = double(0), 10)
d2 <- scan("../dlbcl.web.tests.reedition.15_1.txt", sep = "\t",
               what = double(0), 15)
d3 <- scan("../dlbclA.web.tests.reedition.txt", sep = "\t",
               what = double(0), 10)
dlbcl.fcms <-  c(d1, d2, d3)



d1 <- scan("../breast.web.tests.reedition.10_1.txt", sep = "\t",
               what = double(0), 10)
d2 <- scan("../breast.web.tests.reedition.15_1.txt", sep = "\t",
               what = double(0), 15)
d3 <- scan("../breastA.web.tests.reedition.txt", sep = "\t",
               what = double(0), 10)
breast.fcms <-  c(d1, d2, d3)






postscript(file = "fcms.bench-web.eps", height = 8, width = 12,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")
par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.2)

boxplot(breast.fcms ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "Breast data set (78 arrays x 4751 genes)")


boxplot(dlbcl.fcms ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "DLBCL data set (160 arrays x 7399 genes)")

dev.off()



postscript(file = "tgd.bench-web.eps", height = 8, width = 12,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")
par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.2)

boxplot(breast.tgd ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "Breast data set (78 arrays x 4751 genes)")


boxplot(dlbcl.tgd ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "DLBCL data set (160 arrays x 7399 genes)")

dev.off()
