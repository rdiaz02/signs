rm(list = ls())

system("lamhalt")
system("lamboot -v lamb-host-60slaves.def")

source("fcms.common.R")

fcms.paral.60cpu <- matrix(cbind(rep(NA, 8),
                                 narrays = c(rep(160, 4), 20, 40, 80, 100),
                                 ngenes = c(1000, 2000, 4000, 6000, rep(7399, 4))),
                           ncol = 3)
                         
colnames(fcms.paral.60cpu) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(fcms.paral.60cpu)) {
    fcms.paral.60cpu[i, 1] <- fParal("dlbcl", arrays = fcms.paral.60cpu[i, 2],
                                   genes = fcms.paral.60cpu[i, 3])
}

save(file = "fcms.paral.60cpu.RData", fcms.paral.60cpu)




