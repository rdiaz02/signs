load("tgd.paral.1cpu.RData")
load("tgd.paral.2cpu.RData")
load("tgd.paral.20cpu.RData")
load("tgd.paral.60cpu.RData")

load("tgd.paral.new.1cpu.RData")
load("tgd.paral.new.2cpu.RData")
load("tgd.paral.new.20cpu.RData")
load("tgd.paral.new.60cpu.RData")
load("tgd.paral.new.10cpu.RData")


tgd.paral.out <- data.frame(rbind(tgd.paral.1cpu,
                                  tgd.paral.2cpu,
                                  tgd.paral.20cpu,
                                  tgd.paral.60cpu,
                                  tgd.paral.new.1cpu,
                                  tgd.paral.new.2cpu,
                                  tgd.paral.new.20cpu,
                                  tgd.paral.new.60cpu,
                                  tgd.paral.new.10cpu))
tgd.paral.out$CPUS <- c(rep(1, 8), rep(2, 8), rep(20, 8), rep(60, 8),
                        rep(1, 3), rep(2, 3), rep(20, 3), rep(60, 3),
                        rep(10, 11))

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

tgd.all <- tgd.all[order(tgd.all$CPUS, tgd.all$narrays, tgd.all$ngenes), ]
tgd.all$ratio <- tgd.all$time/tgd.all$time[1:11]
tgd.all$ratio <- 1/tgd.all$ratio

postscript(file = "bench.tgd.eps", height = 9.6, width = 14.4,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")

par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.3)
par(cex.lab = 1)
par(cex.main = 0.9)
par(cex.axis = 0.8)

plot(time ~ narrays,
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
axis(1, at = c(20, 40, 80, 100))

points(time ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 1)  & (ngenes == 7399)),
       lwd = 2,
       col = "blue")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 2)  & (ngenes == 7399)),
       lwd = 2,
       col = "green"
)
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 10)  & (ngenes == 7399)),
       lwd = 2,
       col = "violet")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 20)  & (ngenes == 7399)),
       lwd = 2,
       col = "orange")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 60)  & (ngenes == 7399)),
       lwd = 2,
       col = "red")

legend(20, 55000, c("Original Sequential",
                   "New sequential",
                   "Parall. 2 CPUs",
                   "Parall. 10 CPUs",
                   "Parall. 20 CPUs",
                   "Parall. 60 CPUs"),
       col = c("black", "blue", "green", "violet", "orange", "red"),
       pch = 21,
       lty = 1,
       lwd = 2,
       cex = 0.75)

plot(time ~ ngenes,
     ylim = c(miny, maxy),
     xlab = "Number of genes",
     ylab = "User wall time (seconds)",
     type = "b",
     log = "xy",
     xlim = c(1000, 48000),
     data = subset(tgd.all, is.na(CPUS) & narrays == 160),
     axes = FALSE,
     main = "Effect of number of genes (number of arrays = 160)",
     lwd = 2)
axis(2)
box()
axis(1, at = c(1000, 2000, 4000, 6000, 12000, 24000, 48000))

points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 1)  & (narrays == 160)),
       lwd = 2,
       col = "blue")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 2)  & (narrays == 160)),
       lwd = 2,
       col = "green")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 10)  & (narrays == 160)),
       lwd = 2,
       col = "violet")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 20)  & (narrays == 160)),
       lwd = 2,
       col = "orange")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 60)  & (narrays == 160)),
       lwd = 2,
       col = "red")

dev.off()












postscript(file = "bench.tgd.eps", height = 9.6, width = 14.4,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")

miny <- min(tgd.all$ratio)
maxy <- max(tgd.all$ratio)

par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.3)
par(cex.lab = 1)
par(cex.main = 0.9)
par(cex.axis = 0.8)

plot(ratio ~ narrays,
     ylim = c(miny, maxy),
     xlab = "Number of arrays (samples)",
     ylab = "Fold increase in speed (relative to 1 CPU)",
     type = "b",
     log = "y",
     data = subset(tgd.all, is.na(CPUS) & ngenes == 7399),
     axes = FALSE,
     main = "Effect of number of arrays (number of genes = 7399)",
     lwd = 2)
axis(2)
box()
axis(1, at = c(20, 40, 80, 100))

points(ratio ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 1)  & (ngenes == 7399)),
       lwd = 2,
       col = "blue")
points(ratio ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 2)  & (ngenes == 7399)),
       lwd = 2,
       col = "green"
)
points(ratio ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 10)  & (ngenes == 7399)),
       lwd = 2,
       col = "violet")
points(ratio ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 20)  & (ngenes == 7399)),
       lwd = 2,
       col = "orange")
points(ratio ~ narrays,
       type = "b", log = "y",
       data = subset(tgd.all, (CPUS == 60)  & (ngenes == 7399)),
       lwd = 2,
       col = "red")

tmp <- subset(tgd.all, (CPUS == 1)  & (ngenes == 7399))
text(y = tmp$ratio + 0.1, x = tmp$narrays,
     labels = round(tmp$time),
     log = "y")
     


plot(ratio ~ ngenes,
     ylim = c(miny, maxy),
     xlab = "Number of genes",
     ylab = "Fold increase in speed (relative to 1 CPU)",
     type = "b",
     log = "xy",
     xlim = c(1000, 48000),
     data = subset(tgd.all, is.na(CPUS) & narrays == 160),
     axes = FALSE,
     main = "Effect of number of genes (number of arrays = 160)",
     lwd = 2)
axis(2)
box()
axis(1, at = c(1000, 2000, 4000, 6000, 12000, 24000, 48000))

points(ratio ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 1)  & (narrays == 160)),
       lwd = 2,
       col = "blue")
points(ratio ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 2)  & (narrays == 160)),
       lwd = 2,
       col = "green")
points(ratio ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 10)  & (narrays == 160)),
       lwd = 2,
       col = "violet")
points(ratio ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 20)  & (narrays == 160)),
       lwd = 2,
       col = "orange")
points(ratio ~ ngenes,
       type = "b", log = "xy",
       data = subset(tgd.all, (CPUS == 60)  & (narrays == 160)),
       lwd = 2,
       col = "red")


tmp <- subset(tgd.all, (CPUS == 1)  & (narrays == 160))
text(y = tmp$ratio + 0.1, x = tmp$ngenes,
     labels = round(tmp$time),
     log = "y")


legend(6200, 0.99, c("Original Sequential",
                   "New sequential",
                   "Parall. 2 CPUs",
                   "Parall. 10 CPUs",
                   "Parall. 20 CPUs",
                   "Parall. 60 CPUs"),
       col = c("black", "blue", "green", "violet", "orange", "red"),
       pch = 21,
       lty = 1,
       lwd = 2,
       cex = 0.75)


dev.off()


