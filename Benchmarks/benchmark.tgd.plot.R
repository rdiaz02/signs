
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
