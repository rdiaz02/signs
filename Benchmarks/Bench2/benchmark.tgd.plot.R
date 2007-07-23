load("tgd.paral.1cpu.RData")
load("tgd.paral.2cpu.RData")
load("tgd.paral.20cpu.RData")
load("tgd.paral.60cpu.RData")

tgd.paral.out <- data.frame(rbind(tgd.paral.1cpu,
                       tgd.paral.2cpu,
                       tgd.paral.20cpu,
                       tgd.paral.60cpu))
tgd.paral.out$CPUS <- c(rep(1, 8), rep(2, 8), rep(20, 8), rep(60, 8))

load("tgd.seq.RData")
load("tgd.seq2.RData")
load("tgd.seqA.RData")
load("tgd.seqB.RData")
load("tgd.seqC.RData")

tgd.seq.out <- tgd.seq
tgd.seq.out[3, 1] <- tgd.seqA[3, 1]
tgd.seq.out[4, 1] <- tgd.seqB[4, 1]
tgd.seq.out[8, 1] <- tgd.seqC[8, 1]
tgd.seq.out[5:7, 1] <- tgd.seq2[5:7, 1]
tgd.seq.out <- data.frame(tgd.seq.out)

tgd.seq.out$CPUS <- NA
tgd.all <- rbind(tgd.seq.out, tgd.paral.out)
miny <- min(tgd.all[, 1])
maxy <- max(tgd.all[, 1])


tgd.all$farrays <- factor(tgd.all$narrays)


postscript(file = "bench.tgd.eps", height = 9.6, width = 14.4,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")

par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.5)
par(cex.lab = 1.1)
par(cex.main = 0.9)

plot(time ~ as.numeric(farrays),
     ylim = c(miny, maxy),
     xlab = "Number of arrays (samples)",
     ylab = "User wall time (seconds)",
     type = "b",
     log = "y",
     data = subset(tgd.all, is.na(CPUS) & ngenes == 7399),
     axes = FALSE,
     main = "Effect of number of arrays (number of genes = 7399)",
     lwd = 2)
axis(2)
box()
axis(1, at = c(1, 2, 3, 4), labels = c("20", "40", "80", "100"))

points(time ~ as.numeric(farrays),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 1)  & (ngenes == 7399)),
       lwd = 2,
       col = "blue")
points(time ~ as.numeric(farrays),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 2)  & (ngenes == 7399)),
       lwd = 2,
       col = "green")
points(time ~ as.numeric(farrays),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 20)  & (ngenes == 7399)),
       lwd = 2,
       col = "orange")
points(time ~ as.numeric(farrays),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 60)  & (ngenes == 7399)),
       lwd = 2,
       col = "red")

legend(1, 40000, c("Sequential", "Parall. 1 CPU",
                   "Parall. 2 CPUs", "Parall. 20 CPUs",
                   "Parall. 60 CPUs"),
       col = c("black", "blue", "green", "orange", "red"),
       pch = 21,
       lty = 1,
       lwd = 2,
       cex = 0.75)

plot(time ~ as.numeric(as.factor(ngenes)),
     ylim = c(miny, maxy),
     xlab = "Number of genes",
     ylab = "User wall time (seconds)",
     type = "b",
     log = "y",
     data = subset(tgd.all, is.na(CPUS) & narrays == 160),
     axes = FALSE,
     main = "Effect of number of genes (number of arrays = 160)",
     lwd = 2)
axis(2)
box()
axis(1, at = c(1, 2, 3, 4), labels = c(1000, 2000, 4000, 6000))

points(time ~ as.numeric(as.factor(ngenes)),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 1)  & (narrays == 160)),
       lwd = 2,
       col = "blue")
points(time ~ as.numeric(as.factor(ngenes)),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 2)  & (narrays == 160)),
       lwd = 2,
       col = "green")
points(time ~ as.numeric(as.factor(ngenes)),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 20)  & (narrays == 160)),
       lwd = 2,
       col = "orange")
points(time ~ as.numeric(as.factor(ngenes)),
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 60)  & (narrays == 160)),
       lwd = 2,
       col = "red")

dev.off()












arrays <- c(rep(c(20, 40, 80, 100, 120), rep(5, 5)), rep(40, 20))
genes <- c(rep(40, 25), c(rep(c(20, 80, 160, 320), rep(5, 4))))

dm <- cbind(arrays = arrays, genes = genes)
dm <- data.frame(rbind(dm, dm, dm, dm, dm, dm))

dm$type <- c(rep(c("P_2_slaves/node", "P_6_slaves/node", "P_12_slaves/node",
                   "P_20_slaves/node", "P_60_slaves/node", "Serial"),
                 rep(45, 6)))

dm$time <- NA


times <- NULL

load('parallel.2pernode.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40,
           s.40.20, s.40.80, s.40.160, s.40.320)
load('parallel.6pernode.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40,
           s.40.20, s.40.80, s.40.160, s.40.320)
load('parallel.12pernode.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40,
           s.40.20, s.40.80, s.40.160, s.40.320)
load('parallel.20.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40,
           s.40.20, s.40.80, s.40.160, s.40.320)
load('parallel.60pernode.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40,
           s.40.20, s.40.80, s.40.160, s.40.320)
load('serial.arrays.RData')
times <- c(times, s.20.40, s.40.40, s.80.40, s.100.40, s.120.40)
load('serial.genes.RData')
times <- c(times, s.40.20, s.40.80, s.40.160, s.40.320)

dm$time <- times
timings.tgd <- dm
save(timings.tgd, file = "timings.tgd.RData")


## some idea of measurement variability
tapply(dm$time, list(dm$type, dm$arrays, dm$genes), sd)

load("timings.tgd.RData")
tapply(log(timings.tgd$time),
       list(timings.tgd$type, timings.tgd$arrays, timings.tgd$genes), sd)

tapply(timings.tgd$time,
       list(timings.tgd$type, timings.tgd$arrays, timings.tgd$genes), mean)


### I probably should only plot means, o.w., too busy a plot.
##  I could use tapply, but a pain to reorganize later

geom.mean <- function(x) exp(mean(log(x)))

so <- seq(1, 270, by = 5)
timings.means <- data.frame(arrays = rep(NA, 54),
                            genes = rep(NA, 54),
                            type = rep(NA, 54),
                            time = rep(NA, 54))
for(i in 1:54) {
    j <- so[i]
    timings.means[i, c(1, 2, 3)] <- timings.tgd[j, c(1, 2, 3)]
    si <- j + (0:4)
    timings.means[i, 4] <- geom.mean(timings.tgd[si, 4])
}

## change in array number

posA <- which(timings.means$genes == 40)


pxyA <- function(type) {
    poses <- which((timings.means$genes == 40) & (timings.means$type == type))
    return(cbind(timings.means[poses, c(1, 4)]))
}

pxyB <- function(type) {
    poses <- which((timings.means$arrays == 40) & (timings.means$type == type))
    oo <- order(timings.means[poses, 2])
    return(cbind(timings.means[poses, c(2, 4)][oo, ]))
}

postscript(file = "bench.tgd.eps", height = 9.6, width = 14.4,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")
par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.5)
par(cex.lab = 1.1)
par(cex.main = 0.9)
plot(pxyA("Serial"), ylim = c(10, 5700),
     type = "b", lwd = 2, col = "black",
     xlab = "Number of arrays (samples)",
     ylab = "User wall time (seconds)", log = "y",
     xlim = c(20, 170), xaxt = "n",
     main = "Effect of number of arrays (number of genes = 40)"
     )
axis(1, at = c(20, 40, 80, 100, 120),
     labels = TRUE)

points(pxyA("P_2_slaves/node"),
       type = "b", col = "blue", log = "y", lwd = 2)
points(pxyA("P_6_slaves/node"),
       type = "b", col = "green", log = "y", lwd = 2)
points(pxyA("P_12_slaves/node"),
       type = "b", col = "orange", log = "y", lwd = 2)
points(pxyA("P_20_slaves/node"),
       type = "b", col = "brown", log = "y", lwd = 2)
points(pxyA("P_60_slaves/node"),
       type = "b", col = "red", log = "y", lwd = 2)

text(60, 5000, "Original sequential code")
text(60, 700, "Parallelized code") ##, vfont = c("sans serif", "italic"))

text(cbind(4, 0) + pxyA("P_2_slaves/node")[5, ], "2 slaves/node", col = "blue", adj = 0)
text(cbind(4, 0) + pxyA("P_6_slaves/node")[5, ], "6 slaves/node", col = "green", adj = 0)
text(cbind(4, 0) + pxyA("P_12_slaves/node")[5, ], "12 slaves/node", col = "orange", adj = 0)
text(cbind(4, 0) + pxyA("P_20_slaves/node")[5, ], "20 slaves/node", col = "brown", adj = 0)
text(cbind(4, 0) + pxyA("P_60_slaves/node")[5, ], "60 slaves/node", col = "red", adj = 0)




par(las = 1)
plot(pxyB("Serial"), ylim = c(15, 2000),
     type = "b", lwd = 2, col = "black",
     xlab = "Number of genes",
     ylab = "User wall time (seconds)", log = "y",
     xlim = c(20, 470), xaxt = "n",
     main = "Effect of number of genes (number of arrays = 40)"
     )
axis(1, at = c(20, 40, 80, 160, 320),
     labels = TRUE)

points(pxyB("P_2_slaves/node"),
       type = "b", col = "blue", log = "y", lwd = 2)
points(pxyB("P_6_slaves/node"),
       type = "b", col = "green", log = "y", lwd = 2)
points(pxyB("P_12_slaves/node"),
       type = "b", col = "orange", log = "y", lwd = 2)
points(pxyB("P_20_slaves/node"),
       type = "b", col = "brown", log = "y", lwd = 2)
points(pxyB("P_60_slaves/node"),
       type = "b", col = "red", log = "y", lwd = 2)

text(140, 2000, "Original sequential code")
text(150, 320, "Parallelized code") ##, vfont = c("sans serif", "italic"))

text(cbind(9, 0) + pxyB("P_2_slaves/node")[5, ], "2 slaves/node", col = "blue", adj = 0)
text(cbind(9, 0) + pxyB("P_6_slaves/node")[5, ], "6 slaves/node", col = "green", adj = 0)
text(cbind(9, 0) + pxyB("P_12_slaves/node")[5, ], "12 slaves/node", col = "orange", adj = 0)
text(cbind(9, 0) + pxyB("P_20_slaves/node")[5, ], "20 slaves/node", col = "brown", adj = 0)
text(cbind(9, 0) + pxyB("P_60_slaves/node")[5, ], "60 slaves/node", col = "red", adj = 0)

dev.off()
