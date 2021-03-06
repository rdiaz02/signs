## First, run wrapWebBenchmarks1.py, which calls fl-run-test to do the timings.
## The script rb.sh does that. However, note that for every run we get one value
## for one user, and 20 for 20 users. So I did a few runs by hand, and placed them
## in file 1.1.1

rr <- function(x) rep(x, x)

breast.results <- scan("breast.web.tests.txt", sep = "\t",
                       what = double(0), 40)
dlbcl.results <- scan("dlbcl.web.tests.txt", sep = "\t",
                       what = double(0), 40)
nusers <- c(rep(1, 5), rr(5), rr(10), rr(20))

breast.results.re <- scan("breast.web.tests.reedition.txt", sep = "\t",
                       what = double(0), 40)
dlbcl.results.re <- scan("dlbcl.web.tests.reedition.txt", sep = "\t",
                       what = double(0), 40)
nusers2 <- c(rep(1, 5), rr(5), rr(10), rr(20))

dl2 <- c(scan('dlbcl.web.tests.reedition.15_1.txt', sep = "\t", what = double(0), 15),
         scan('dlbcl.web.tests.reedition.50_1.txt', sep = "\t", what = double(0), 50))

br2 <- c(scan('breast.web.tests.reedition.15_1.txt', sep = "\t", what = double(0), 15),
         scan('breast.web.tests.reedition.50_1.txt', sep = "\t", what = double(0), 50))

nusers3 <- c(rr(15), rr(50))

nusers <- c(nusers, nusers2, nusers3)

dlbcl.results <- c(dlbcl.results, dlbcl.results.re, dl2)
breast.results <- c(breast.results, breast.results.re, br2)


postscript(file = "bench-web.eps", height = 8, width = 12,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")
par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.2)

boxplot(breast.results ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "Breast data set (78 arrays x 4751 genes)")


boxplot(dlbcl.results ~ nusers, xlab = "Number of simultaneous users",
        ylab = "User wall time (seconds)",
        main = "DLBCL data set (160 arrays x 7399 genes)")

dev.off()
