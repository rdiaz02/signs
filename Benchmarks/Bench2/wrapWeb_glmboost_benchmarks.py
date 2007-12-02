#!/usr/bin/env python


import os
import sys
import time

def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkSigns_glmboost.py Signs.' + test + '_glmboost')
        time.sleep(5)
    for uu in range(users):
        timef[uu] = float(t[uu].readlines()[-1].strip())

    return timef

def launchAll(test, lusers):
    outall = []
    for users in lusers:
        print 'users ' + str(users)
        outall = outall + launchUTests(test, users)
    return outall



def writeFile(testout, name):
    fout = open(name, mode = 'w')
    for result in testout:
        fout.write(str(result))
        fout.write('\t')
    fout.close()
    


## OK, this is a very ugly thing, but if the script crashes, no need to start from zero;
    ## we save partial work
    
breast_15_1 = [-9999999.11]
breast_10_1 = [-9999999.11]

dlbcl_10_1 = [-9999999.11]
dlbcl_15_1 = [-9999999.11]


try:
    breast_10_1 = launchUTests('breast', 10)
except:
    None
writeFile(breast_10_1, 'breast.web.glmboost.10_1.txt')
   
try:
    dlbcl_10_1 = launchUTests('dlbcl', 10)
except:
    None
writeFile(dlbcl_10_1, 'dlbcl.web.glmboost.10_1.txt')
   

try:
    breast_15_1 = launchUTests('breast', 15)
except:
    None
writeFile(breast_15_1, 'breast.web.glmboost.15_1.txt')
   
try:
    dlbcl_15_1 = launchUTests('dlbcl', 15)
except:
    None
writeFile(dlbcl_15_1, 'dlbcl.web.glmboost.15_1.txt')

breast_1_1 = launchUTests('breast', 1)
breast_1_2 = launchUTests('breast', 1)
breast_1_3 = launchUTests('breast', 1)
breast_1_4 = launchUTests('breast', 1)
breast_1_5 = launchUTests('breast', 1)
breast_5_1 = launchUTests('breast', 5)

breastA = breast_1_1 + breast_1_2 + breast_1_3 + breast_1_4 + breast_1_5 +\
          breast_5_1 
writeFile(breastA, 'breastA.web.glmboost.txt')

dlbcl_1_1 = launchUTests('dlbcl', 1)
dlbcl_1_2 = launchUTests('dlbcl', 1)
dlbcl_1_3 = launchUTests('dlbcl', 1)
dlbcl_1_4 = launchUTests('dlbcl', 1)
dlbcl_1_5 = launchUTests('dlbcl', 1)
dlbcl_5_1 = launchUTests('dlbcl', 5)

dlbclA = dlbcl_1_1 + dlbcl_1_2 + dlbcl_1_3 + dlbcl_1_4 + dlbcl_1_5 +\
         dlbcl_5_1 
writeFile(dlbclA, 'dlbclA.web.glmboost.txt')


# try:
#     dlbcl_20_1 = launchUTests('dlbcl', 20)
# except:
#     None
# writeFile(dlbcl_20_1, 'dlbcl.web.tests.reedition.20_1.txt')


# try:
#     breast_20_1 = launchUTests('breast', 20)
# except:
#      None
# writeFile(breast_20_1, 'breast.web.tests.reedition.20_1.txt')



# try:
#     breast_50_1 = launchUTests('breast', 50)
# except:
#     None
# writeFile(breast_50_1, 'breast.web.tests.reedition.50_1.txt')

# try:
#     dlbcl_50_1 = launchUTests('dlbcl', 50)
# except:
#     None
# writeFile(dlbcl_50_1, 'dlbcl.web.tests.reedition.50_1.txt')
