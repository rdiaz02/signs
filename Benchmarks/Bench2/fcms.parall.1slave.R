rm(list = ls())
system("lamhalt")
source("fcms.common.R")

fcms.paral.1cpu <- matrix(cbind(rep(NA, 8),
                                 narrays = c(rep(160, 4), 20, 40, 80, 100),
                                 ngenes = c(1000, 2000, 4000, 6000, rep(7399, 4))),
                           ncol = 3)
                         
colnames(fcms.paral.1cpu) <- c("time", "narrays", "ngenes")
                                          
for(i in 1:nrow(fcms.paral.1cpu)) {
    cat(" \n\n arrays = ", fcms.paral.1cpu[i, 2],
        " genes = ", fcms.paral.1cpu[i, 3], "\n")
    fcms.paral.1cpu[i, 1] <- fParal("dlbcl", arrays = fcms.paral.1cpu[i, 2],
                                   genes = fcms.paral.1cpu[i, 3])
}

save(file = "fcms.paral.1cpu.RData", fcms.paral.1cpu)




