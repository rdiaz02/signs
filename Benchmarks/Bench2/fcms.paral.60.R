ARGUMENT '60' __ignored__


R version 2.4.1 Patched (2007-02-28 r40809)
Copyright (C) 2007 The R Foundation for Statistical Computing
ISBN 3-900051-07-0

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> rm(list = ls())
> 
> ll <- length(commandArgs())
> nf <- commandArgs()[ll]
> numero <- as.numeric(commandArgs()[ll - 1])
Warning message:
NAs introduced by coercion 
> 
> cat("\n numero is ", numero, "\n")

 numero is  NA 
> 
> system("lamhalt")

LAM 7.1.2/MPI 2 C++/ROMIO - Indiana University

> 
> if(numero != 1)
+     system(paste("lamboot -v lamb-host-", numero, "slaves.def", sep = ""))
Error in if (numero != 1) system(paste("lamboot -v lamb-host-", numero,  : 
	missing value where TRUE/FALSE needed
Execution halted
