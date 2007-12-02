##############################################
##############################################
######                              ##########
######       Program starts         ##########
######                              ##########
##############################################
##############################################

## rm(list = ls())


cat("\nRunning\n", file = "Status.msg")

checkpoint.num <- scan("checkpoint.num", what = double(0), n = 1)




## From: http://ace.acadiau.ca/math/ACMMaC/Rmpi/sample.html
# In case R exits unexpectedly, have it automatically clean up 
# resources  taken up by Rmpi (slaves, memory, etc...)
## But does it really do it??
.Last <- function(){
    try(sink()) ## in case we are bailing out from within sink

    status <- file("Status.msg", "w")
    cat("Normal termination\n", file = status)
    flush(status)
    close(status)
    
    RterminatedOK <- file("RterminatedOK", "w")
    cat("\nNormal termination\n", file = RterminatedOK)
    flush(RterminatedOK)
    close(RterminatedOK)
    save.image()
    if (is.loaded("mpi_initialize")){ 
        if (mpi.comm.size(1) > 0){ 
        try(print("Please use mpi.close.Rslaves() to close slaves."), silent = TRUE)
        try(mpi.close.Rslaves() , silent = TRUE)
        } 
        try(print("Please use mpi.quit() to quit R"), silent = TRUE)
        cat("\n\n Normal termination\n")
        try(stopCluster(TheCluster), silent = TRUE)
        ##        .Call("mpi_finalize")
        try(system(paste("/http/mpi.log/killLAM.py", lamSESSION, "&")))
        cat(paste("\n Did the call to killLAM.py 1 ", date(), " \n"),
        file = "tmp.checks", append = TRUE)
        try(mpi.quit(save = "no"), silent = TRUE)
    }
    try(stopCluster(TheCluster), silent = TRUE)
    cat("\n\n Normal termination\n")
    ## In case the CGI is not called (user kills browser)
    ## have a way to stop lam
    try(system(paste("/http/mpi.log/killLAM.py", lamSESSION, "&")))
    cat(paste("\n Did the call to killLAM.py 2 ", date(), " \n"),
        file = "tmp.checks", append = TRUE)
    try(mpi.quit(save = "no"), silent = TRUE)
}

doCheckpoint <- function(num) {
    save.image()
    sink("checkpoint.num")
    cat(num)
    sink()
}




startExecTime <- format(Sys.time())

pid <- Sys.getpid()
write.table(file = "pid.txt", pid,
            row.names = FALSE,
            col.names = FALSE)


## attach pid to name in R.running.procs
hostn <- system("hostname", intern = TRUE)
new.name1 <- unlist(strsplit(getwd(), "\/"))
new.name1 <- paste(new.name1[length(new.name1)], "@", hostn, sep = "")
new.name <- paste("R.", new.name1, "%", pid, sep = "")
new.name1 <- paste("R.", new.name1, sep = "")
system(paste("mv ../../R.running.procs/", new.name1,
             " ../../R.running.procs/", new.name,
             sep = ""))

sink(file = "hostname")
cat(hostn)
sink()


### The above is no longer really needed. What follows is what buryThem2.py
##  uses
sink(file = "current_R_proc_info")
cat(hostn)
cat("   ")
cat(pid)
sink()





#library(CGIwithR)
library(Rmpi)
library(survival)
library(combinat)
library(MASS)
library(SignS2)
library(papply)
#library(snow)
library(R2HTML)
library(GDD)
#library(imagemap) ## FIXME: remove


imClose <- function (im) {
## prevent all the "Closing PNG file ..."
   dev.off(im$Device)
}
		


gdd.width <- png.width <- 480
gdd.height <- png.height <- 410
png.pointsize <- 10
png.family = "Helvetica"
graphDir <- paste(getwd(), "/", sep = "")



nfold <- 10


##############################################
##############################################
######                              ##########
######         Error checking       ##########
######                              ##########
##############################################
##############################################


caughtUserError <- function(message) {
    GDD("ErrorFigure.png", width = 600,
           height = 500, ps = png.pointsize)
    plot(x = c(0, 1), y = c(0, 1), 
         type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.7, "There was a PROBLEM with your data.")
    text(0.5, 0.5,
    "Please read carefully the error messages,")
    
    text(0.5, 0.3, "fix the problem, and try again.")
    dev.off()
    sink(file = "results.txt")
    cat(message)
    sink()
    sink(file = "exitStatus")
    cat("Error\n\n")
    cat(message)
    sink()
    quit(save = "no", status = 11, runLast = TRUE)
}


caughtOurError <- function(message) {
    GDD("ErrorFigure.png", width = 600,
           height = 500, ps = png.pointsize)
    plot(x = c(0, 1), y = c(0, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.7, "There was a PROBLEM with the code.")
    text(0.5, 0.5,
    "Please let us know (send us the URL),")
    
    text(0.5, 0.3, "so that we can fix it.")
    dev.off()
    sink(file = "results.txt")
    cat(message)
    sink()
    sink(file = "exitStatus")
    cat("Error\n\n")
    cat(message)
    sink()
    quit(save = "no", status = 11, runLast = TRUE)
}




methodSurv <- scan("methodSurv", what = "", n = 1)
if( methodSurv == "TGD") {
    epi <- scan("epi", what = double(0), n = 1)
    maxiter <- scan("maxiter", what = double(0), n = 1)
    tau <- scan("tau", what = double(0), n = 1)
} else if (methodSurv == "FCMS") {
    MaxSize <- scan("MaxSize", what = double(0), n = 1)
    MinSize <- scan("MinSize", what = double(0), n = 1)
    MinCor <- scan("MinCor", what = double(0), n = 1)
    Minp <- scan("Minp", what = double(0), n = 1)
} else if (methodSurv == "cforest") {
    ngenes <- scan("ngenes", what = double(0), n = 1)
} else { ## nothing else for now
    caughtUserError("This method is not yet implemented.")
}


## FIXME: add params for cforest



########################################################

########   Start MPI here to check if everything OK

#########################################################

MPI_MIN_UNIVERSE_SIZE <- 10

if (mpi.universe.size () < MPI_MIN_UNIVERSE_SIZE) {
    cat("\n\n mpi.universe.size () < MPI_MIN_UNIVERSE_SIZE \n\n")
    quit(save = "no", status = 11, runLast = FALSE)
}

if(methodSurv == "TGD") mpiSpawnAll()
if(methodSurv == "FCMS") mpiSpawnAll(10)
if(methodSurv == "cforest") mpiSpawnAll(10)


## if(methodSurv == "TGD") {
##     TheCluster <- makeCluster(mpi.universe.size(), "MPI")
##     } else if(methodSurv == "FCMS") {
##     mpiSpawnAll()
## }

sink(file = "mpiOK")
cat("MPI started OK\n")
sink()

trylam <- try(
              lamSESSION <- scan("lamSuffix", what = "character",
                                 sep = "\t", strip.white = TRUE))


#########################################################

########   HTML and other utility functions

#########################################################

cleanHTMLhead <- function(file, title = "", h1 = NULL, append = FALSE) {
    if(is.null(h1)) h1 <- title
    out <-
        paste("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">",
          "\n<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-15\">",
          "\n<title>", title, "</title> \n <body>\n<h1>", h1, "</h1>")
    write(out, file = file, append = append)
}
cleanHTMLtail <- function(file, append= TRUE) {
    write("</body></html>", file = file, append = append)
}






#########################################################

########   Continue reading data

#########################################################

if(checkpoint.num < 1) {

idtype <- "None"
organism <- "None"
idtype <- try(scan("idtype", what = "", n = 1))
organism <- try(scan("organism", what = "", n = 1))


useValidation <- try(scan("usevalidation", what = "", n = 1))

trytime <- try(
             Time <- scan("time", sep = "\t", strip.white = TRUE, nlines = 1))
if(class(trytime) == "try-error")
    caughtUserError("The time file is not of the appropriate format\n")

## to prevent problems with a space at end of classes
if(is.na(Time[length(Time)])) Time <- Time[-length(Time)]

tryevent <- 
    try(Event <- scan("event", sep = "\t", strip.white = TRUE, nlines = 1))
if(class(tryevent) == "try-error")
    caughtUserError("The status file is not of the appropriate format\n")

## to prevent problems with a space at end of classes
if(is.na(Event[length(Event)])) Event <- Event[-length(Event)]


if (! all(Event %in% c(0, 1))) {
  caughtUserError("Your status data is not valid; can only be 0 or 1\n")
}


if(length(Time) < 10) {
    caughtUserError(paste("Your data should contain at least ten (10) samples\n",
                       "but your data have only", length(Time), ".\n"))
}


## get rid of """ in data
tmp <- try(
       	   system("sed 's/\"//g' covariate > cvt; mv cvt covariate")
	   )

tmp <- try(
 	       system("sed 's/\"//g' validationcovariate > vcvt; mv vcvt validationcovariate")
       )   	   
	  
## zz: get substitution of ' ' too	   

num.cols.covariate <- count.fields("covariate", sep = "\t",
                                   quote = "",
                                   comment.char = "#")

if(length(unique(num.cols.covariate)) > 1) {
    message <-
    paste("The number of columns in your covariate file\n",
          "is not the same for all rows (genes).\n",
          "We find the following number of columns\n",
          paste(num.cols.covariate, collapse = ", "))
    caughtUserError(message)
}


tryxdata <- try(
                xdata <- read.table("covariate", header = FALSE, sep = "\t",
                                    strip.white = TRUE,
                                    comment.char = "#",
                                    quote = ""))
if(class(tryxdata) == "try-error")
    caughtUserError("The covariate file is not of the appropriate format\n")

geneNames <- xdata[, 1]
xdata <- xdata[, -1]

if(length(unique(geneNames)) < length(geneNames)) {
    dupnames <- which(duplicated(geneNames))
    message <- paste("Gene names are not unique.\n",
                      "Please change them so that they are unique.\n",
                      "The duplicated names are in rows", dupnames, "\n")
    caughtUserError(message)
}
    
rownames(xdata) <- geneNames

arrayNames <- scan("arrayNames", sep = "\t", what = "char", quote ="")

if(length(arrayNames) > 0) {
    arrayNames <- arrayNames[-1]
    if(length(unique(arrayNames)) < length(arrayNames)) {
        dupnames <- which(duplicated(arrayNames))
        message <- paste("Array names are not unique.\n",
                          "Please change them so that they are unique.\n",
                          "The duplicated names are ", dupnames, "\n")
        caughtUserError(message)
    }
    if(ncol(xdata) != length(arrayNames)) {
        emessage <- paste("We get that the number of columns in your data (", ncol(xdata), ")\n",
                          "is different from the number of column names (", length(arrayNames), ")\n",
                          "Check for things such as '#' or '#NULL!' in the middle of your data.\n")
        caughtUserError(emessage)
    }
    colnames(xdata) <- arrayNames
}
xdata <- t(xdata)


if(length(Time) != dim(xdata)[1]) {
    message <- paste("The survival time file and the covariate file\n",
                      "do not agree on the number of arrays: \n",
                      length(Time), " arrays according to the survival time file but \n",
                      dim(xdata)[1], " arrays according to the covariate data.\n",
                      "Please fix this problem and try again.\n")
    caughtUserError(message)  
}

if(length(Event) != dim(xdata)[1]) {
    message <- paste("The survival status file and the covariate file\n",
                      "do not agree on the number of arrays: \n",
                      length(Event), " arrays according to the survival time file but \n",
                      dim(xdata)[1], " arrays according to the covariate data.\n",
                      "Please fix this problem and try again.\n")
    caughtUserError(message)  
}

if(methodSurv == "cforest") {
    if(ngenes > ncol(xdata)) ngenes <- ncol(xdata)
}


if(methodSurv == "FCMS")
    if((ncol(xdata) < MinSize)) {
        caughtUserError(paste("You have requested an impossible situation:\n",
                              "your minimal cluster size ", MinSize," is larger than the \n",
                              "number of genes ", ncol(xdata), ".\n"))
    }

if(!(is.numeric(xdata))) {
    caughtUserError("Your covariate file contains non-numeric data. \n That is not allowed.\n")
}
if(any(is.na(xdata))) {
    caughtUserError("Your covariate file contains missing values. \n That is not allowed.\n")
}
if(!(is.numeric(Time))) {
    caughtUserError("Your survival time file contains non-numeric data. \n That is not allowed.\n")
}
if(any(is.na(Time))) {
    caughtUserError("Your survival status file contains missing values. \n That is not allowed.\n")
}
if(!(is.numeric(Event))) {
    caughtUserError("Your survival status file contains non-numeric data. \n That is not allowed.\n")
}
if(any(is.na(Event))) {
    caughtUserError("Your survival status file contains missing values. \n That is not allowed.\n")
}



if(useValidation == "yes") {
    trytime <- try(
                   validationTime <- scan("validationtime", sep = "\t", strip.white = TRUE, nlines = 1))
    if(class(trytime) == "try-error")
        caughtUserError("The validation time file is not of the appropriate format\n")
    
    ## to prevent problems with a space at end of classes
    if(is.na(validationTime[length(validationTime)])) validationTime <- validationTime[-length(validationTime)]
    
    tryevent <- 
        try(validationEvent <- scan("validationevent", sep = "\t", strip.white = TRUE, nlines = 1))
    if(class(tryevent) == "try-error")
        caughtUserError("The validation status file is not of the appropriate format\n")
    
    ## to prevent problems with a space at end of classes
    if(is.na(validationEvent[length(validationEvent)])) validationEvent <- validationEvent[-length(validationEvent)]
    
    
    num.cols.validationcovariate <- count.fields("validationcovariate", sep = "\t",
                                       quote = "",
                                       comment.char = "#")
    
    if(length(unique(num.cols.validationcovariate)) > 1) {
        message <-
            paste("The number of columns in your validation covariate file\n",
                  "is not the same for all rows (genes).\n",
                  "We find the following number of columns\n",
                  paste(num.cols.validationcovariate, collapse = ""))
        caughtUserError(message)
    }
    
    tryxdata <- try(
                    validationxdata <- read.table("validationcovariate", header = FALSE, sep = "\t",
                                        strip.white = TRUE,
                                        comment.char = "#",
                                        quote = ""))
    if(class(tryxdata) == "try-error")
        caughtUserError("The validation covariate file is not of the appropriate format\n")
    validationgeneNames <- validationxdata[, 1]
    validationxdata <- validationxdata[, -1]
    
    if(length(unique(validationgeneNames)) < length(validationgeneNames)) {
        dupnames <- which(duplicated(validationgeneNames))
        message <- paste("validation Gene names are not unique.\n",
                         "Please change them so that they are unique.\n",
                         "The duplicated names are in rows", dupnames, "\n")
        caughtUserError(message)
    }
    
    rownames(validationxdata) <- validationgeneNames
    
    validationarrayNames <- scan("validationarrayNames", sep = "\t", what = "char", quote = "")
    
    if(length(validationarrayNames) > 0) {
        validationarrayNames <- validationarrayNames[-1]
        if(length(unique(validationarrayNames)) < length(validationarrayNames)) {
            dupnames <- which(duplicated(validationarrayNames))
            message <- paste("validation Array names are not unique.\n",
                             "Please change them so that they are unique.\n",
                             "The duplicated names are ", dupnames, "\n")
            caughtUserError(message)
        }
        colnames(validationxdata) <- validationarrayNames
    }
    validationxdata <- t(validationxdata)
    
    
    if(length(validationTime) != dim(validationxdata)[1]) {
        message <- paste("The validation survival time file and the covariate file\n",
                         "do not agree on the number of arrays: \n",
                         length(validationTime), " arrays according to the survival time file but \n",
                         dim(validationxdata)[1], " arrays according to the covariate data.\n",
                         "Please fix this problem and try again.\n")
        caughtUserError(message)  
    }
    
    if(length(validationEvent) != dim(validationxdata)[1]) {
        message <- paste("The validation survival status file and the covariate file\n",
                         "do not agree on the number of arrays: \n",
                         length(validationEvent), " arrays according to the survival time file but \n",
                         dim(validationxdata)[1], " arrays according to the covariate data.\n",
                         "Please fix this problem and try again.\n")
        caughtUserError(message)  
    }
    
    if(!(is.numeric(validationxdata))) {
        caughtUserError("Your validation covariate file contains non-numeric data. \n That is not allowed.\n")
    }
    if(any(is.na(validationxdata))) {
        caughtUserError("Your validation covariate file contains missing values. \n That is not allowed.\n")
    }
    if(!(is.numeric(validationTime))) {
        caughtUserError("Your validation survival time file contains non-numeric data. \n That is not allowed.\n")
    }
    if(any(is.na(validationTime))) {
        caughtUserError("Your validation survival status file contains missing values. \n That is not allowed.\n")
    }
    if(!(is.numeric(validationEvent))) {
        caughtUserError("Your validation survival status file contains non-numeric data. \n That is not allowed.\n")
    }
    if(any(is.na(validationEvent))) {
        caughtUserError("Your validation survival status file contains missing values. \n That is not allowed.\n")
    }


    if(!(identical(colnames(xdata), colnames(validationxdata)))) {
        caughtUserError("Gene names for the validation and training data MUST be the same. ")
    }
}

doCheckpoint(1)

}

options(warn = -1)

if(methodSurv == "TGD") {#### Starting part for Threshold Gradient Descent

    u.threshold <- tau


if(checkpoint.num < 2) {
    xdata <- scale(xdata, center = TRUE, scale = TRUE)

    trycode <- try(
                   allDataRun <- tgdSingle(xdata, Time, Event,
                                            unique.thres =u.threshold,
                                            epi, maxiter,
                                            nfold = 10)
                   )

    if(class(trycode) == "try-error")
        caughtOurError(paste("Function tgdSingleP bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))
                       
    doCheckpoint(2)
}
if(checkpoint.num < 3) {
    trycode <- try(
                   cvTGDResults <- cvTGDP(xdata, Time, Event,
                                          thres = c(u.threshold, u.threshold),
                                          epi, thresGrid = 1, 
                                          maxiter, 
                                          nfold) 
                   )
 
    if(class(trycode) == "try-error")
        caughtOurError(paste("Function cvTGDP bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))


    doCheckpoint(3)
}
if(checkpoint.num < 4) {

    GDD(file = "kmplot-honest.png", width = png.width,
        height = png.height, ps = png.pointsize)
    KM.visualize(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good   ####  Fig 1
    dev.off()
    GDD(file = "kmplot-overfitt.png", width = png.width,
        height = png.height, ps = png.pointsize)
    KM.visualize(allDataRun$tgd.alldata$scores, Time,          ####  Fig 2
                 Event, ngroups = 2) ## Overfitt
    dev.off()

    pdf(file = "kmplot-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good   ####  Fig 1
    dev.off()
    pdf(file = "kmplot-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize(allDataRun$tgd.alldata$scores, Time,          ####  Fig 2
                 Event, ngroups = 2) ## Overfitt
    dev.off()




    GDD(file = "kmplot4-honest.png", width = png.width,
           height = png.height,ps = png.pointsize)
    KM.visualize4(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good   ####  Fig 1.4
    dev.off()
    GDD(file = "kmplot4-overfitt.png", width = png.width,
        height = png.height,ps = png.pointsize)
    KM.visualize4(allDataRun$tgd.alldata$scores, Time,          ####  Fig 2.4
                 Event, ngroups = 2) ## Overfitt
    dev.off()

    pdf(file = "kmplot4-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize4(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good   ####  Fig 1.4
    dev.off()
    pdf(file = "kmplot4-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize4(allDataRun$tgd.alldata$scores, Time,          ####  Fig 2.4
                 Event, ngroups = 2) ## Overfitt
    dev.off()


    pdf(file = "kmplot3-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize3(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    pdf(file = "kmplot3-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize3(allDataRun$tgd.alldata$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()
    GDD(file = "kmplot3-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(cvTGDResults$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    GDD(file = "kmplot3-overfitt.png", w=png.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(allDataRun$tgd.alldata$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()

    
    GDD(file = "cvpl.png", width = png.width,
           height = png.height,ps = png.pointsize)
    plot.cvpl.single(allDataRun$cvpl.mat, epi, tau)                       ####  Fig 3
    dev.off()


    pdf(file = "cvpl.pdf", width = png.width,
           height = png.height)
    plot.cvpl.single(allDataRun$cvpl.mat, epi, tau)                        ####  Fig 3
    dev.off()





    
    sink(file = "results.txt")
    cat("\n <h2>4. Model fitted to all data</h2>\n")
        

    trycode <- try(
                   outm <- summaryTGDrun(xdata, Time, Event, allDataRun,
                                         epi, thres = c(u.threshold, u.threshold),
                                         thresGrid = 1, plot = TRUE,
                                         genesOut = TRUE,
                                         outfile = "genes.all.out") )          #### Fig 4: fstdgrun


    if(class(trycode) == "try-error")
        caughtOurError(paste("Function summaryTGDrun bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))


    cat("\n <hr><h2>5. Cross-validation runs</h2>\n")
    trycode <- try(
                   summary.cvTGD(cvTGDResults, allDataRun, rownames(xdata))
                   )
    if(class(trycode) == "try-error")
        caughtOurError(paste("Function summary.cvTGD bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

    if(useValidation == "yes") {
        cat("\n <hr><h2>6. Validation data</h2>\n")
        
        valpred <- validationxdata %*% allDataRun$betas
       
        cat("\n\n Scores (linear predictor) for validation data\n")
        print(valpred)
        
    }      
    
    sink()

    doCheckpoint(4)
}
    if(useValidation == "yes") {
        pdf(file = "kmplot-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize(valpred, validationTime,
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        GDD(file = "kmplot-validation.png", width = png.width,
               height = png.height, ps = png.pointsize)
        KM.visualize(valpred, validationTime,                         
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()

        pdf(file = "kmplot4-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize4(valpred, validationTime,
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off() 
        GDD(file = "kmplot4-validation.png", width = png.width,
               height = png.height,ps = png.pointsize)
        KM.visualize4(valpred, validationTime,                         
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        pdf(file = "kmplot3-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize3(valpred, validationTime,
                      validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        GDD(file = "kmplot3-validation.png", w=gdd.width,
            h = gdd.height, ps = png.pointsize,
            type = "png")
        KM.visualize3(valpred, validationTime,                         
                      validationEvent, ngroups = 2, addmain = NULL)
        dev.off()    


    }


    
##    stopCluster(TheCluster)
    save.image()
    
    
} else if(methodSurv == "FCMS") {



    if(checkpoint.num < 2) {

        mpi.bcast.Robj2slave(idtype)
        mpi.bcast.Robj2slave(organism)
        
    
    MaxIterationsCox <- 200
    
    trycode <- try(
                   all.res1 <- dStep1.parallel(xdata, Time, Event,
                                               Minp, MaxIterationsCox)
                   )
    if(class(trycode) == "try-error")
        caughtOurError(paste("Function dStep1.parallel bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

#################    p-value tables      #################            

    p.values.original <- data.frame(Names = geneNames,
                                    p.value = all.res1[, 2],
                                    coeff = all.res1[, 1], 
                                    abs.coeff = abs(all.res1[, 1]),
                                    fdr = all.res1[, 6],
                                    Warning = all.res1[, 5])

    if (any(is.na(p.values.original))) {
        p.values.original[is.na(p.values.original)] <- 999.999
        cat("Oh, oh, some nas in p.values.originalll",
            file = "nas.in.p.values.WARN")
    }
    
    ## zz: we will need, either here or in python, to generate the links al IDClight
    write.table(file = "p.values.coeffs.txt",
                p.values.original, row.names = FALSE, 
                col.names = TRUE,
                quote = FALSE,
                sep = "\t")
    system(paste("/http/signs/cgi/order.html.py", idtype, organism)) ## call python to generate pre-sorted HTML tables
    


        doCheckpoint(2)
    }
        
#################    step AIC      #################                


    if(checkpoint.num < 3) {

    sink(file = "stepAIC.output.txt")
    trycode <- try(
    all.res3 <- fitDave.res1Given(xdata, Time, Event,              #### Fig 01:
                                  res1 = all.res1,             #### Fig 02
                                  Minp, MaxSize,                  ### ClusterNegativeCoeffs
                                  MinSize, MinCor,             ### ClusterPositiveCoeffs
                                  MaxIterationsCox,
                                  plot = TRUE,
                                  interactive = TRUE)
                   )
    if(class(trycode) == "try-error")
    caughtOurError(paste("Function fitDave.res1Given bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

    cat("\n\n Final model selected (beware: p-values are biased down!!)\n")

    print(all.res3$model)
    sink()

    doCheckpoint(3)
}

#################    run all the rest      #################


if(checkpoint.num < 4) {
    
    trycode <- try(
                   cvDaveRun <- cvDave.parallel3(x = xdata, time = Time,
                                                 event = Event,
                                                 p = Minp, maxSize = MaxSize,
                                                 minSize = MinSize,
                                                 minCor = MinCor,
                                                 MaxIterationsCox = MaxIterationsCox,
                                                 nfold = nfold)
                   )
    
    if(class(trycode) == "try-error")
        caughtOurError(paste("Function cvDaveRun bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

    sink(file = "results.txt")

    cat("<h2>4. Model fitted to all data</h2>\n")
    cat("<h3>4.1. Components, genes, coefficients</h3>\n")	
    trycode <- try(
                   all.res.out <- selectedSignatures(all.res3, colnames(xdata),
                                                     print = TRUE, out = TRUE,
                                                     html = TRUE)
                   )
    
    if(class(trycode) == "try-error")
        caughtOurError(paste("Function selectedSignatures bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

    sink()

    doCheckpoint(4)

}

if(checkpoint.num < 5) {   

    sink(file = "results.txt", append = TRUE)


    ### Output correlation matrix of clusters
    tmp.cor <- cor(all.res3$clusterResults$md)
    ## ugly hack to get bold for row names; yes, could use CSS
    rownames(tmp.cor) <- paste("<b>", rownames(tmp.cor), "</b>", sep ="")

    cleanHTMLhead(file = "correlationMatrixClusters.html",
                  title = "Correlation matrix of clusters or components")
    HTML(tmp.cor, file = "correlationMatrixClusters.html",
         digits = 4, align = "right")
    cleanHTMLtail(file = "correlationMatrixClusters.html")

    ## get decent spacing
    system("sed 's/<th>/<th width=50>/g' correlationMatrixClusters.html > cm; mv cm correlationMatrixClusters.html")

    cleanHTMLhead(file = "stepAIC.output.html",
                  title = "Steps from variable selection")
    write(paste("A \"+\" means that the variable was considered for addition and ",
                "a \"-\" means that the variable was considered for deletion. ",
                "The number shows the AIC after the addition/deletion of ",
                "each of the variables. Smaller AIC is better.\n"),
          file = "stepAIC.output.html",append = TRUE)
    write("<pre>", file = "stepAIC.output.html",append = TRUE)
    system("cat stepAIC.output.html stepAIC.output.txt > tmp.s; mv tmp.s stepAIC.output.html")
    write("</pre>", file = "stepAIC.output.html",append = TRUE)
    cleanHTMLtail(file = "stepAIC.output.html")

    
    cat("\n <hr align=\"left\" width=80>")
    cat("<h3>4.2. <a href=\"correlationMatrixClusters.html\" target=\"CorMatrix_window\">View</a> correlation matrix of clusters</h3>")
    cat("<h3>4.3. <a href=\"stepAIC.output.html\" target=\"stepAIC_window\">View</a> steps of variable selection</h3>")



    cat("\n<hr><h2>5. Cross-validation runs</h2>\n")
   
    trycode <- try(
                   summary.cvDave(cvDaveRun, all.res.out, rownames(xdata),
                                  colnames(xdata), html = TRUE)
                   )

    if(class(trycode) == "try-error")
        caughtOurError(paste("Function summary.cvDave bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))

    sink()

doCheckpoint(5)
}

   
    if(checkpoint.num < 6) {
        sink(file = "results.txt", append = TRUE)

    if(useValidation == "yes") {
        cat("\n <h2>6. Validation data</h2>\n")
        trycode <- try(
                       valpred <- dPredictNew(all.res3, validationxdata)
                       )
        if(class(trycode) == "try-error")
            caughtOurError(paste("Function dPredictNew bombed unexpectedly with error",
                                 trycode, ". \n Please let us know so we can fix the code."))


        cat("\n\n <h3>6.1. <a href=\"scores.validation.html\" target=\"vscores_window\">View</a>",
            "the scores (linear predictor) for validation data.</h3>\n")
        cleanHTMLhead(file = "scores.validation.html",
                      title = "Linear predictor scores for validation data")
        write(paste("<TABLE frame=\"box\">\n",
                    "<tr><th>Validation subject/array</th> <th>Linear score</th></tr>\n"),
              file = "scores.validation.html", append = TRUE)
        wout <- ""
        for(i in 1:length(valpred)) {
            wout <- paste(wout, "\n <tr align=right>",
            "<td>", rownames(valpred)[i], "</td><td>", valpred[i], "</td></tr>\n")
        }
        wout <- paste(wout, "</TABLE>")
        write(wout, file = "scores.validation.html", append = TRUE)
        cleanHTMLtail(file = "scores.validation.html")

       
    }      

    sink()

        doCheckpoint(6)
    }
    

#########################################################
########
########   Do the plots, and we are done
########
#########################################################

    if(checkpoint.num < 7) {
    gdd.width <- 480
    gdd.height <- 410

    pdf(file = "kmplot-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1
    dev.off()
    pdf(file = "kmplot-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2
    dev.off()
    GDD(file = "kmplot-honest.png", w=gdd.width,
           h = gdd.height, ps = png.pointsize,
           type = "png")
    KM.visualize(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1
    dev.off()
    GDD(file = "kmplot-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2
    dev.off()
    pdf(file = "kmplot4-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize4(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.4
    dev.off()
    pdf(file = "kmplot4-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize4(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.4
    dev.off()
    GDD(file = "kmplot4-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize4(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.4
    dev.off()
    GDD(file = "kmplot4-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize4(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.4
    dev.off()


    pdf(file = "kmplot3-honest.pdf", width = png.width,
        height = png.height)
    KM.visualize3(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    pdf(file = "kmplot3-overfitt.pdf", width = png.width,
        height = png.height)
    KM.visualize3(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()
    GDD(file = "kmplot3-honest.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(cvDaveRun$OOB.scores, Time,
                 Event, ngroups = 2, addmain = NULL) ## Good              #### Fig 1.3
    dev.off()
    GDD(file = "kmplot3-overfitt.png", w=gdd.width,
        h = gdd.height, ps = png.pointsize,
        type = "png")
    KM.visualize3(all.res3$scores, Time,                         
                 Event, ngroups = 2) ## Overfitt                   #### Fig 2.3
    dev.off()

    doCheckpoint(7)
}







    

    
    if(useValidation == "yes") {
        pdf(file = "kmplot-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize(valpred, validationTime,
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        GDD(file = "kmplot-validation.png", w=gdd.width,
               h = gdd.height, ps = png.pointsize,
               type = "png")
        KM.visualize(valpred, validationTime,                         
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()

        pdf(file = "kmplot4-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize4(valpred, validationTime,
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        GDD(file = "kmplot4-validation.png", w=gdd.width,
               h = gdd.height, ps = png.pointsize,
               type = "png")
        KM.visualize4(valpred, validationTime,                         
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()


        pdf(file = "kmplot3-validation.pdf", width = png.width,
            height = png.height)
        KM.visualize3(valpred, validationTime,
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
        GDD(file = "kmplot3-validation.png", w=gdd.width,
               h = gdd.height, ps = png.pointsize,
               type = "png")
        KM.visualize3(valpred, validationTime,                         
                     validationEvent, ngroups = 2, addmain = NULL)
        dev.off()
    }
    
    

    save.image()
##    try(mpi.close.Rslaves())
##    mpi.quit(save = "no")
} else if(methodSurv == "cforest") {
    if(checkpoint.num < 2) { ## Model for all data
        mpi.bcast.Robj2slave(idtype)
        mpi.bcast.Robj2slave(organism)
        MaxIterationsCox <- 200
        if(useValidation == "yes") {
            trycode <- try(
                            cf.all <- my.cforest(xdata, Time, Event, ngenes, validationxdata)
                           )
        } else {
            trycode <- try(
                           cf.all <- my.cforest(xdata, Time, Event, ngenes, NULL)
                           )
        }
        if(class(trycode) == "try-error")
            caughtOurError(paste("Function my.cforest bombed unexpectedly with error",
                                 trycode, ". \n Please let us know so we can fix the code."))
        doCheckpoint(2)
    }
    if(checkpoint.num < 3) { ## Cross-validation
        trycode <- try(
                       cf.cv.output <- my.cforest.cv(xdata, Time, Event, ngenes)
                       )
        if(class(trycode) == "try-error")
            caughtOurError(paste("Function my.forest.cv bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))
        doCheckpoint(3)
    }
    if(checkpoint.num < 4) { ## plots
        kmplots(cf.cv.output$OOB.scores, cf.all$overfit_predicted_surv_time,
                Time, Event)
        if(useValidation == "yes") {
            kmplots.validation(cf.all$predicted_surv_time, validationTime,
                               validationEvent)
        }
        doCheckpoint(4)
    }
    if(checkpoint.num < 5) { ## all text output
        print.selected.genes(cf.all, idtype, organism)
        print.cv.results(cf.cv.output, cf.all, rownames(xdata), html.level = 3)
        doCheckpoint(5)
    }
    if(checkpoint.num < 6) {
        if(useValidation == "yes") {
            print.validation.results(cf.all)
        }
        doCheckpoint(6)
    }
    save.image()
} else if(methodSurv == "glmboost") {
    if(checkpoint.num < 2) { ## Model for all data
        mpi.bcast.Robj2slave(idtype)
        mpi.bcast.Robj2slave(organism)
        if(useValidation == "yes") {
            trycode <- try(
                            glmb.all <- my.glmboost(xdata, Time, Event, validationxdata)
                           )
        } else {
            trycode <- try(
                           glmb.all <- my.glmboost(xdata, Time, Event, NULL)
                           )
        }
        if(class(trycode) == "try-error")
            caughtOurError(paste("Function my.glmboost bombed unexpectedly with error",
                                 trycode, ". \n Please let us know so we can fix the code."))
        doCheckpoint(2)
    }
    if(checkpoint.num < 3) { ## Cross-validation
        trycode <- try(
                       gb.cv.output <- my.glmboost.cv(xdata, Time, Event)
                       )
        if(class(trycode) == "try-error")
            caughtOurError(paste("Function my.glmboost.cv bombed unexpectedly with error",
                             trycode, ". \n Please let us know so we can fix the code."))
        doCheckpoint(3)
    }
    if(checkpoint.num < 4) { ## plots
        kmplots(gb.cv.output$OOB.scores, glmb.all$overfit_predicted_surv_time,
                Time, Event)
        if(useValidation == "yes") {
            kmplots.validation(glmb.all$predicted_surv_time, validationTime,
                               validationEvent)
        }
        doCheckpoint(4)
    }
    if(checkpoint.num < 5) { ## all text output
        print.selected.genes(glmb.all, idtype, organism)
        print.cv.results(gb.cv.output, glmb.all, rownames(xdata), html.level = 3)
        doCheckpoint(5)
    }
    if(checkpoint.num < 6) {
        if(useValidation == "yes") {
            print.validation.results(glmb.all)
        }
        doCheckpoint(6)
    }
    save.image()
}

    
    

##    try(mpi.close.Rslaves())
##    mpi.quit(save = "no")


## turn the html into txt with:
## html2text -width 200 -nobs -o Results.txt pre-results.html 
## but this is confussing because there exists a results.txt




