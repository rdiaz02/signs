#!/usr/bin/python

import glob
import socket
import sys
import os
import cgi 
import time
import shutil
import dircache
import random
import subprocess
from stat import ST_SIZE

# import cgitb
# cgitb.enable() 
# sys.stderr = sys.stdout

APP_NAME = "signs2"
sys.path.append("/asterias-web-apps/web-apps-common")
from web_apps_config import *
from web_apps_common_funcs import *

acceptedMethodSurvs = ('FCMS', 'TGD', 'cforest', 'glmboost')

def valueNumUpload(fieldName, testNumber = 'float', minValue = 0, APP_NAME = APP_NAME):
    """Upload and get the values and do some checking. For text and radio selections
    with positive numeric data.
    We assume there is an existing call to fs = cgi.FieldStorage()"""

    if not fs.has_key(fieldName):
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> ", fieldName, "value required </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if fs[fieldName].filename:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> ", fieldName, "should not be a file. </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    if type(fs[fieldName]) == type([]):
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> ", fieldName, "should be a single value.</p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    else:
        tmp = fs[fieldName].value

    ## Accept only numeric values that can be turned to floats or ints
    if testNumber == 'float':
        try:
            tmpn = float(tmp)
        except:
            commonOutput(APP_NAME)
            print "<h1> SignS INPUT ERROR </h1>"    
            print "<p> ", fieldName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()
    else:
        try:
            tmpn = int(tmp)
        except:
            commonOutput(APP_NAME)
            print "<h1> SignS INPUT ERROR </h1>"    
            print "<p> ", fieldName, "is not a valid numeric value.</p>"
            print "<p> Please fill up the required fields and try again</p>"
            print "</body></html>"
            sys.exit()

    if tmpn < minValue:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> ", fieldName, "smaller than smallest accepted value (", minValue, "). </p>"
        print "<p> Please fill up the required fields and try again</p>"
        print "</body></html>"
        sys.exit()
    # transferring files to final destination;
    fileInServer = tmpDir + "/" + fieldName
    srvfile = open(fileInServer, mode = 'w')
    srvfile.write(str(tmpn))
    srvfile.close()
    os.chmod(fileInServer, 0666)

    return tmpn


#########################################################
#########################################################

####          Execution starts here      ################

#########################################################
#########################################################



## Deleting tmp directories older than MAX_time
## NOT needed anymore; delete_old_dirs runs as cron job!
## YES, needed: cronjobs are a pain
currentTime = time.time()
currentTmp = dircache.listdir("/asterias-web-apps/signs2/www/tmp")
for directory in currentTmp:
    tmpS = "/asterias-web-apps/signs2/www/tmp/" + directory
    if (currentTime - os.path.getmtime(tmpS)) > MAX_time:
        shutil.rmtree(tmpS)


### Creating temporal directories
newDir = str(random.randint(1, 10000)) + str(os.getpid()) + str(random.randint(1, 100000)) + str(int(currentTime)) + str(random.randint(1, 10000))
redirectLoc = "/tmp/" + newDir
tmpDir = "/asterias-web-apps/signs2/www/tmp/" + newDir
os.mkdir(tmpDir)
os.chmod(tmpDir, 0700)


### File and parameter upload
fs = cgi.FieldStorage()

idtype = dummyUpload('idtype', 'None', tmpDir)
organism = dummyUpload('organism', 'None', tmpDir)


methodSurv = radioUpload('methodSurv', acceptedMethodSurvs, fs, tmpDir, APP_NAME)
if methodSurv == 'FCMS':
    maxsize = valueNumUpload('MaxSize', 'int', 2)
    minsize = valueNumUpload('MinSize', 'int', 2)
    if minsize >= maxsize:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> Max number of genes surely must be larger "
        print "than Min number of genes.</p>"
        print "<p> Please fill up the required fields and try again.</p>"
        print "minsize " + minsize + ";   maxsize" + maxsize
        print "</body></html>"
        sys.exit()
    mincor = valueNumUpload('MinCor', 'float', 0)
    if mincor > 1:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> Min. correlation should of course be between 0 and 1."
        print "<p> Please fill up the required fields and try again.</p>"
        print "mincor " + mincor
        print "</body></html>"
        sys.exit()
    Minp = valueNumUpload('Minp', 'float', 0)
    if Minp > 1:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> Minimal p for gene-wise Cox should of course be between 0 and 1."
        print "<p> Please fill up the required fields and try again.</p>"
        print "</body></html>"
        sys.exit()
if methodSurv == 'TGD':
    maxiter = valueNumUpload('maxiter', 'int', 100)
    epi = valueNumUpload('epi', 'float', 0)
    tau = valueNumUpload('tau', 'float', 0)
    if tau > 1:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"    
        print "<p> Tau should be between 0 and 1."
        print "<p> Please fill up the required fields and try again.</p>"
        print "</body></html>"
        sys.exit()
if methodSurv == 'cforest':
    ngenes = valueNumUpload('ngenes', 'int', 2)


##check if file coming from preP

# if(fs.getfirst("covariate2")!= None):
#     prep_tmpdir = fs.getfirst("covariate2")
#         ## an ugly hack, as prep not in this filesystem
#     os.system('wget http://prep.iib.uam.es/tmp/' + prep_tmpdir +
#               '/outdata.txt -O ' + tmpDir + '/covariate')
#     ## shutil.copy("/asterias-web-apps/prep/www/tmp/" + prep_tmpdir +"/outdata.txt",tmpDir + "/covariate")
# else:
## Uploading files and checking not abusively large
fileUpload('covariate', fs, tmpDir, APP_NAME)
if os.stat(tmpDir + '/covariate')[ST_SIZE] > MAX_covariate_size:
    shutil.rmtree(tmpDir)
    commonOutput(APP_NAME)
    print "<h1> SignS INPUT ERROR </h1>"
    print "<p> Covariate file way too large </p>"
    print "<p> Covariate files this size not allowed.</p>"
    print "</body></html>"
    sys.exit()

fileUpload('time', fs, tmpDir, APP_NAME)
if os.stat(tmpDir + '/time')[ST_SIZE] > MAX_time_size:
    shutil.rmtree(tmpDir)
    commonOutput(APP_NAME)
    print "<h1> SignS INPUT ERROR </h1>"
    print "<p> Survival time file way too large </p>"
    print "<p> This size is not allowed.</p>"
    print "</body></html>"
    sys.exit()

fileUpload('event', fs, tmpDir, APP_NAME)
if os.stat(tmpDir + '/event')[ST_SIZE] > MAX_time_size:
    shutil.rmtree(tmpDir)
    commonOutput(APP_NAME)
    print "<h1> SignS INPUT ERROR </h1>"
    print "<p> Survival status file way too large </p>"
    print "<p> This size is not allowed.</p>"
    print "</body></html>"
    sys.exit()

if fs.has_key('validation'):
    open(tmpDir + '/usevalidation', mode = 'w').write('yes\n')
    fileUpload('validationcovariate', fs, tmpDir, APP_NAME)
    if os.stat(tmpDir + '/validationcovariate')[ST_SIZE] > MAX_covariate_size:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"
        print "<p> Validation Covariate file way too large </p>"
        print "<p> Covariate files this size not allowed.</p>"
        print "</body></html>"
        sys.exit()
    fileUpload('validationtime', fs, tmpDir, APP_NAME)
    if os.stat(tmpDir + '/validationtime')[ST_SIZE] > MAX_time_size:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"
        print "<p> Validation Survival time file way too large </p>"
        print "<p> This size is not allowed.</p>"
        print "</body></html>"
        sys.exit()

    fileUpload('validationevent', fs, tmpDir, APP_NAME)
    if os.stat(tmpDir + '/validationevent')[ST_SIZE] > MAX_time_size:
        shutil.rmtree(tmpDir)
        commonOutput(APP_NAME)
        print "<h1> SignS INPUT ERROR </h1>"
        print "<p> Validation Survival status file way too large </p>"
        print "<p> This size is not allowed.</p>"
        print "</body></html>"
        sys.exit()
    
## Upload worked OK. We store the original names of the files in the
## browser for later report:
## We'll need to get this working for the validation data.zz
fileNamesBrowser = open(tmpDir + '/fileNamesBrowser', mode = 'w')
if(fs.getfirst("covariate2")== None):
   fileNamesBrowser.write(fs['covariate'].filename + '\n')
fileNamesBrowser.write(fs['time'].filename + '\n')
fileNamesBrowser.write(fs['event'].filename + '\n')
fileNamesBrowser.close()




## current number of processes > max number of processes?
## and we do it here, not before, so that we have the most
## current info about number of process right before we launch R.


## First, delete any R file left (e.g., from killing procs, etc).
RrunningFiles = dircache.listdir("/asterias-web-apps/signs2/www/R.running.procs")
for Rtouchfile in RrunningFiles:
    tmpS = "/asterias-web-apps/signs2/www/R.running.procs/" + Rtouchfile
    if (currentTime - os.path.getmtime(tmpS)) > R_MAX_time:
        os.remove(tmpS)

## Now, verify any processes left
numRsigns = len(glob.glob("/asterias-web-apps/signs2/www/R.running.procs/R.*@*%*"))
if numRsigns > MAX_signs:
    shutil.rmtree(tmpDir)
    commonOutput(APP_NAME)
    print "<h1> SignS problem: The servers are too busy </h1>"
    print "<p> Because of the popularity of the application "
    print " the maximum number of simultaneous runs of SignS has been reached.</p>"
    print "<p> Please try again later.</p>"
    print "<p> We apologize for the inconvenience.</p>"    
    print "</body></html>"
    sys.exit()
    

################        Launching R   ###############

# prepare the arrayNames file:

covarInServer = tmpDir + "/covariate"
arrayNames = tmpDir + "/arrayNames"
srvfile = open(covarInServer, mode = 'r')
arrayfile = open(arrayNames, mode = 'w')
num_name_lines = 0

while 1:
    line = srvfile.readline()
    if not line: break
    if (line.find("#name") == 0) or (line.find("#NAME") == 0) or (line.find("#Name") == 0) \
           or (line.find('"#name"') == 0) or (line.find('"#NAME"') == 0) or (line.find('"#Name"') == 0):
        num_name_lines = num_name_lines + 1
        if num_name_lines > 1:
            commonOutput(APP_NAME)
            print """ You have more than one line with #Name (or #NAME or #name), in the data matrix \
                   but only one is allowed."""
            sys.exit()
        arrayfile.write(line)
        arrayfile.write("\n\n")
    
srvfile.close()
arrayfile.close()   
os.chmod(arrayNames, 0600)



if fs.has_key('validation'):
    validationcovarInServer = tmpDir + "/validationcovariate"
    validationarrayNames = tmpDir + "/validationarrayNames"
    validationsrvfile = open(validationcovarInServer, mode = 'r')
    validationarrayfile = open(validationarrayNames, mode = 'w')
    num_name_lines = 0
    while 1:
        line = validationsrvfile.readline()
        if not line: break
        if (line.find("#name") == 0) or (line.find("#NAME") == 0) or (line.find("#Name") == 0) \
               or (line.find('"#name"') == 0) or (line.find('"#NAME"') == 0) or (line.find('"#Name"') == 0):
            num_name_lines = num_name_lines + 1
            if num_name_lines > 1:
                commonOutput(APP_NAME)
                print """ You have more than one line with #Name (or #NAME or #name), in the data matrix \
                but only one is allowed."""
                sys.exit()
            validationarrayfile.write(line)
            validationarrayfile.write("\n\n")
        
    validationsrvfile.close()
    validationarrayfile.close()   
    os.chmod(validationarrayNames, 0600)



# I do not understand why, but the next makes it wait eternally until apache timesout
# This is now in f1.R, but this is not the cause of the problem. It is spawning
# the other python process
# checkpoint = os.system("/bin/echo '0' > " + tmpDir + "/checkpoint.num")
# checkpointfile = open(tmpDir + "/checkpoint.num", mode = "w")
# checkpointfile.write('0')
# checkpointfile.close()
# os.chmod(tmpDir + "/checkpoint.num", 0666)

## touch Rout, o.w. checkdone can try to open a non-existing file
touchRout = os.system("/bin/touch " + tmpDir + "/f1.Rout") 
touchRrunning = os.system("/bin/touch /asterias-web-apps/signs2/www/R.running.procs/R." + newDir +
                          "@" + socket.gethostname())
shutil.copy("/asterias-web-apps/signs2/cgi/f1.R", tmpDir)
createResultsFile = os.system("/bin/touch " + tmpDir + "/results.txt")




###########   Creating a results.hmtl   ###############

## Copy to tmpDir a results.html that redirects to checkdone.cgi
## If communication gets broken, there is always a results.html
## that will do the right thing.
shutil.copy("/asterias-web-apps/signs2/cgi/results-pre.html", tmpDir)
os.system('echo copied_results-pre >> ' + tmpDir + '/checkdone2')
os.system("/bin/sed 's/sustituyeme/" + newDir + "/g' " +
          tmpDir + "/results-pre.html > " +
          tmpDir + "/results.html; rm " +
          tmpDir +"/results-pre.html")
# os.system('echo sed_results-pre >> ' + tmpDir + '/checkdone2')


## Launch the running and monitoring program
subprocess.Popen(['/asterias-web-apps/signs2/cgi/runAndCheck.py', tmpDir],
                 stdout = subprocess.PIPE, stdin = subprocess.PIPE, \
                 stderr = subprocess.PIPE)


##############    Return autorefresing results.hmtl    ##################
print "Location:"+ getQualifiedURL("/tmp/" + newDir + "/results.html")
print ""




