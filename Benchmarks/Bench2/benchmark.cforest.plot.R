load("cforest.paral.1cpu.RData")
load("cforest.paral.2cpu.RData")
load("cforest.paral.20cpu.RData")
load("cforest.paral.60cpu.RData")
load("cforest.paral.10cpu.RData")


cforest.paral.out <- data.frame(rbind(cforest.paral.1cpu,
                                   cforest.paral.2cpu,
                                   cforest.paral.10cpu,
                                   cforest.paral.20cpu,
                                   cforest.paral.60cpu))
cforest.paral.out$CPUS <- c(rep(1, 11), rep(2, 11), rep(10, 11),
                         rep(20, 11), rep(60, 11))
cforest.all <- cforest.paral.out
miny <- min(cforest.all[, 1])
maxy <- max(cforest.all[, 1])

cforest.all$farrays <- factor(cforest.all$narrays)

cforest.all$ratio <- cforest.all$time/cforest.all$time[1:11]
cforest.all$ratio <- 1/cforest.all$ratio




flab <- function(cpus, nge, narr, cex = 0.8) {
    if (is.null(narr)) {
        tmp <- subset(cforest.all, (CPUS == cpus)  & (ngenes == nge))
        text(y = 0.1 + tmp$ratio, x = tmp$narrays,
             labels = round(tmp$time), cex = cex)
    }
    else {
        tmp <- subset(cforest.all, (CPUS == cpus)  & (narrays == narr))
        text(y = 0.1 + tmp$ratio, x = tmp$ngenes,
             labels = round(tmp$time), cex = cex)
    }
}


miny <- min(cforest.all$ratio)
maxy <- max(cforest.all$ratio)

postscript(file = "bench.cforest.eps", height = 9.6, width = 14.4,
           horizontal = FALSE,
           onefile = FALSE, paper = "special")

par(mfrow = c(1, 2))
par(las = 1)
par(cex = 1.3)
par(cex.lab = 1)
par(cex.main = 0.9)
par(cex.axis = 0.8)
par(mgp = c(2.5, 1, 0))
par(mar = c(4, 4, 3, 2)) 

plot(ratio ~ narrays,
     ylim = c(miny, maxy),
     xlab = "Number of arrays (samples)",
     ylab = "Fold increase in speed (relative to 1 CPU)",
     type = "b",
#     log = "y",
     data = subset(cforest.all, (CPUS == 1) & (ngenes == 7399)),
     axes = FALSE,
     main = "Effect of number of arrays (number of genes = 7399)",
     lwd = 2,
     col = "blue")
axis(2)
box()
axis(1, at = c(20, 40, 80, 100))

points(ratio ~ narrays,
       type = "b", 
       data = subset(cforest.all, (CPUS == 2)  & (ngenes == 7399)),
       lwd = 2,
       col = "green"
)
points(ratio ~ narrays,
       type = "b", 
       data = subset(cforest.all, (CPUS == 10)  & (ngenes == 7399)),
       lwd = 2,
       col = "violet")
points(ratio ~ narrays,
       type = "b", 
       data = subset(cforest.all, (CPUS == 20)  & (ngenes == 7399)),
       lwd = 2,
       col = "orange")
points(ratio ~ narrays,
       type = "b", 
       data = subset(cforest.all, (CPUS == 60)  & (ngenes == 7399)),
       lwd = 2,
       col = "red")


flab(1, 7399, NULL)
flab(2, 7399, NULL)
flab(10, 7399, NULL)

legend(60, 3, c("1 CPU",
                   "Parall. 2 CPUs",
                   "Parall. 10 CPUs",
                   "Parall. 20 CPUs",
                   "Parall. 60 CPUs"),
       col = c("blue", "green", "violet", "orange", "red"),
       pch = 21,
       lty = 1,
       lwd = 2,
       cex = 0.75)

plot(ratio ~ ngenes,
     ylim = c(miny, maxy),
     xlab = "Number of genes",
     ylab = "Fold increase in speed (relative to 1 CPU)",
     type = "b",
     log = "x",
     data = subset(cforest.all, (CPUS == 1) & (narrays == 160)),
     axes = FALSE,
     main = "Effect of number of genes (number of arrays = 160)",
     lwd = 2,
     col = "blue")
axis(2)
box()
#par(las = 3)
axis(1, at =  c(1000, 2000, 4000, 6000, 12000, 24000, 48000))

points(ratio ~ ngenes,
       type = "b", 
       data = subset(cforest.all, (CPUS == 2)  & (narrays == 160)),
       lwd = 2,
       col = "green")
points(ratio ~ ngenes,
       type = "b", 
       data = subset(cforest.all, (CPUS == 10)  & (narrays == 160)),
       lwd = 2,
       col = "violet")
points(ratio ~ ngenes,
       type = "b", 
       data = subset(cforest.all, (CPUS == 20)  & (narrays == 160)),
       lwd = 2,
       col = "orange")
points(ratio ~ ngenes,
       type = "b", 
       data = subset(cforest.all, (CPUS == 60)  & (narrays == 160)),
       lwd = 2,
       col = "red")


flab(1, NULL, 160)
flab(2, NULL, 160)
flab(10, NULL, 160)

dev.off()











## postscript(file = "bench.fcms.eps", height = 9.6, width = 14.4,
##            horizontal = FALSE,
##            onefile = FALSE, paper = "special")

## par(mfrow = c(1, 2))
## par(las = 1)
## par(cex = 1.3)
## par(cex.lab = 1)
## par(cex.main = 0.9)
## par(cex.axis = 0.8)

## plot(time ~ narrays,
##      ylim = c(miny, maxy),
##      xlab = "Number of arrays (samples)",
##      ylab = "User wall time (seconds)",
##      type = "b",
##      log = "y",
##      data = subset(fcms.all, (CPUS == 1) & (ngenes == 7399)),
##      axes = FALSE,
##      main = "Effect of number of arrays (number of genes = 7399)",
##      lwd = 2,
##      col = "blue")
## axis(2)
## box()
## axis(1, at = c(20, 40, 80, 100))

## ## points(time ~ as.numeric(farrays),
## ##        type = "b", log = "y",
## ##        data = subset(fcms.all, (CPUS == 1)  & (ngenes == 7399)),
## ##        lwd = 2,
## ##        col = "blue")
## points(time ~ narrays,
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 2)  & (ngenes == 7399)),
##        lwd = 2,
##        col = "green"
## )
## points(time ~ narrays,
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 10)  & (ngenes == 7399)),
##        lwd = 2,
##        col = "violet")
## points(time ~ narrays,
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 20)  & (ngenes == 7399)),
##        lwd = 2,
##        col = "orange")
## points(time ~ narrays,
##        type = "b", log = "y",
##        data = subset(fcms.all, (CPUS == 60)  & (ngenes == 7399)),
##        lwd = 2,
##        col = "red")

## legend(20, 5000, c("1 CPU",
##                    "Parall. 2 CPUs",
##                    "Parall. 10 CPUs",
##                    "Parall. 20 CPUs",
##                    "Parall. 60 CPUs"),
##        col = c("blue", "green", "violet", "orange", "red"),
##        pch = 21,
##        lty = 1,
##        lwd = 2,
##        cex = 0.75)

## plot(time ~ ngenes,
##      ylim = c(miny, maxy),
##      xlab = "Number of genes",
##      ylab = "User wall time (seconds)",
##      type = "b",
##      log = "xy",
##      data = subset(fcms.all, (CPUS == 1) & (narrays == 160)),
##      axes = FALSE,
##      main = "Effect of number of genes (number of arrays = 160)",
##      lwd = 2,
##      col = "blue")
## axis(2)
## box()
## #par(las = 3)
## axis(1, at =  c(1000, 2000, 4000, 6000, 12000, 24000, 48000))

## ## points(time ~ ngenes,
## ##        type = "b", log = "y",
## ##        data = subset(fcms.all, (CPUS == 1)  & (narrays == 160)),
## ##        lwd = 2,
## ##        col = "blue")
## points(time ~ ngenes,
##        type = "b", log = "xy",
##        data = subset(fcms.all, (CPUS == 2)  & (narrays == 160)),
##        lwd = 2,
##        col = "green")
## points(time ~ ngenes,
##        type = "b", log = "xy",
##        data = subset(fcms.all, (CPUS == 10)  & (narrays == 160)),
##        lwd = 2,
##        col = "violet")
## points(time ~ ngenes,
##        type = "b", log = "xy",
##        data = subset(fcms.all, (CPUS == 20)  & (narrays == 160)),
##        lwd = 2,
##        col = "orange")
## points(time ~ ngenes,
##        type = "b", log = "xy",
##        data = subset(fcms.all, (CPUS == 60)  & (narrays == 160)),
##        lwd = 2,
##        col = "red")

## dev.off()
