rm(list = ls())

system("lamhalt")
system("lamboot -v lamb-host-20slaves.def")

source("fcms.common.R")

fcms.paral.20cpu <- matrix(cbind(rep(NA, 8),
                                 narrays = c(rep(160, 4), 20, 40, 80, 100),
                                 ngenes = c(1000, 2000, 4000, 6000, rep(7399, 4))),
                           ncol = 3)
                         
colnames(fcms.paral.20cpu) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(fcms.paral.20cpu)) {
    cat(" \n\n arrays = ", fcms.paral.20cpu[i, 2],
        " genes = ", fcms.paral.20cpu[i, 3], "\n")

    fcms.paral.20cpu[i, 1] <- fParal("dlbcl", arrays = fcms.paral.20cpu[i, 2],
                                   genes = fcms.paral.20cpu[i, 3])
}

save(file = "fcms.paral.20cpu.RData", fcms.paral.20cpu)




