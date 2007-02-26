#!/usr/bin/env python


## To do: Use only 1, 2, 5, 10 simultaneous users.
## and save each one as run (and embed in try, and assign a default
## so we can save).


import os
import sys
import time

def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkSigns.py Signs.' + test)
        time.sleep(5)
    for uu in range(users):
        timef[uu] = float(t[uu].readlines()[1].strip())

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
    



breast_10_1 = [-9999999.11]
breast_20_1 = [-9999999.22]
breast_50_1 = [-9999999.22]

dlbcl_10_1 = [-9999999.11]
dlbcl_20_1 = [-9999999.22]
dlbcl_50_1 = [-9999999.22]


## running 50 completely stupid; 50 is more than the number we get over two weeks!

# try:
#     breast_50_1 = launchUTests('breast', 50)
# except:
#     None
# writeFile(breast_50_1, 'breast.web.tests.50_1.txt')

# try:
#     dlbcl_50_1 = launchUTests('dlbcl', 50)
# except:
#     None
# writeFile(dlbcl_50_1, 'dlbcl.web.tests.50_1.txt')



try:
    dlbcl_20_1 = launchUTests('dlbcl', 20)
except:
    None
writeFile(dlbcl_20_1, 'dlbcl.web.tests.20_1.txt')


try:
    breast_20_1 = launchUTests('breast', 20)
except:
     None
writeFile(breast_20_1, 'breast.web.tests.20_1.txt')

try:
    breast_10_1 = launchUTests('breast', 10)
except:
    None
writeFile(breast_10_1, 'breast.web.tests.10_1.txt')
   
try:
    dlbcl_10_1 = launchUTests('dlbcl', 10)
except:
    None
writeFile(dlbcl_10_1, 'dlbcl.web.tests.10_1.txt')
   




breast_1_1 = launchUTests('breast', 1)
breast_1_2 = launchUTests('breast', 1)
breast_1_3 = launchUTests('breast', 1)
breast_1_4 = launchUTests('breast', 1)
breast_1_5 = launchUTests('breast', 1)
breast_5_1 = launchUTests('breast', 5)

breastA = breast_1_1 + breast_1_2 + breast_1_3 + breast_1_4 + breast_1_5 +\
          breast_5_1 
writeFile(breastA, 'breastA.web.tests.txt')

dlbcl_1_1 = launchUTests('dlbcl', 1)
dlbcl_1_2 = launchUTests('dlbcl', 1)
dlbcl_1_3 = launchUTests('dlbcl', 1)
dlbcl_1_4 = launchUTests('dlbcl', 1)
dlbcl_1_5 = launchUTests('dlbcl', 1)
dlbcl_5_1 = launchUTests('dlbcl', 5)

dlbclA = dlbcl_1_1 + dlbcl_1_2 + dlbcl_1_3 + dlbcl_1_4 + dlbcl_1_5 +\
         dlbcl_5_1 
writeFile(dlbclA, 'dlbclA.web.tests.txt')

breast = breast_1_1 + breast_1_2 + breast_1_3 + breast_1_4 + breast_1_5 +\
         breast_5_1 + breast_10_1 + breast_20_1 ## + breast_50_1
writeFile(breast, 'breast.web.tests.txt')

dlbcl = dlbcl_1_1 + dlbcl_1_2 + dlbcl_1_3 + dlbcl_1_4 + dlbcl_1_5 +\
        dlbcl_5_1 + dlbcl_10_1 + dlbcl_20_1 ## + dlbcl_50_1
writeFile(dlbcl, 'dlbcl.web.tests.txt')