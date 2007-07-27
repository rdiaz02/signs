load("fcms.paral.1cpu.RData")
load("fcms.paral.2cpu.RData")
load("fcms.paral.20cpu.RData")
load("fcms.paral.60cpu.RData")
load("fcms.paral.10cpu.RData")


fcms.paral.out <- data.frame(rbind(fcms.paral.1cpu,
                                   fcms.paral.2cpu,
                                   fcms.paral.10cpu,
                                   fcms.paral.20cpu,
                                   fcms.paral.60cpu))
fcms.paral.out$CPUS <- c(rep(1, 11), rep(2, 11), rep(10, 11),
                         rep(20, 11), rep(60, 11))
fcms.all <- fcms.paral.out
miny <- min(fcms.all[, 1])
maxy <- max(fcms.all[, 1])

fcms.all$farrays <- factor(fcms.all$narrays)





postscript(file = "bench.fcms.eps", height = 9.6, width = 14.4,
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
     data = subset(fcms.all, (CPUS == 1) & (ngenes == 7399)),
     axes = FALSE,
     main = "Effect of number of arrays (number of genes = 7399)",
     lwd = 2,
     col = "blue")
axis(2)
box()
axis(1, at = c(20, 40, 80, 100))

## points(time ~ as.numeric(farrays),
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 1)  & (ngenes == 7399)),
##        lwd = 2,
##        col = "blue")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(fcms.all, (CPUS == 2)  & (ngenes == 7399)),
       lwd = 2,
       col = "green"
)
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(fcms.all, (CPUS == 10)  & (ngenes == 7399)),
       lwd = 2,
       col = "violet")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(fcms.all, (CPUS == 20)  & (ngenes == 7399)),
       lwd = 2,
       col = "orange")
points(time ~ narrays,
       type = "b", log = "y",
       data = subset(fcms.all, (CPUS == 60)  & (ngenes == 7399)),
       lwd = 2,
       col = "red")

legend(20, 5000, c("1 CPU",
                   "Parall. 2 CPUs",
                   "Parall. 10 CPUs",
                   "Parall. 20 CPUs",
                   "Parall. 60 CPUs"),
       col = c("blue", "green", "violet", "orange", "red"),
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
     data = subset(fcms.all, (CPUS == 1) & (narrays == 160)),
     axes = FALSE,
     main = "Effect of number of genes (number of arrays = 160)",
     lwd = 2,
     col = "blue")
axis(2)
box()
#par(las = 3)
axis(1, at =  c(1000, 2000, 4000, 6000, 12000, 24000, 48000))

## points(time ~ ngenes,
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 1)  & (narrays == 160)),
##        lwd = 2,
##        col = "blue")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(fcms.all, (CPUS == 2)  & (narrays == 160)),
       lwd = 2,
       col = "green")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(fcms.all, (CPUS == 10)  & (narrays == 160)),
       lwd = 2,
       col = "violet")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(fcms.all, (CPUS == 20)  & (narrays == 160)),
       lwd = 2,
       col = "orange")
points(time ~ ngenes,
       type = "b", log = "xy",
       data = subset(fcms.all, (CPUS == 60)  & (narrays == 160)),
       lwd = 2,
       col = "red")

dev.off()
