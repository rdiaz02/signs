#!/usr/bin/env python

import os
import sys

## outappend = sys.argv[1]

outappend = 'singlerun'

NUM_USERS = (1, 1, 1, 1, 1, 
##             2, 2, 2,
             5,
             10,
             20)


def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkSigns.py Signs.' + test)
    
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
    

# dlbcl = launchAll('dlbcl', NUM_USERS)
# writeFile(dlbcl, 'dlbcl.web.bnchmk.' + outappend + 'txt')


# breast = launchAll('breast', NUM_USERS)
# writeFile(breast, 'breast.web.bnchmk.' + outappend + 'txt')

## 2362.9028670799999, 2368.0365898599998: dlbcl with 1 sample

## OK, this is a very ugly thing, but I was getting random crashes
## that lost all work. Now I do this from an interactive python shell.

breast_1_1 = launchUTests('breast', 1)
breast_1_2 = launchUTests('breast', 1)
breast_1_3 = launchUTests('breast', 1)
breast_1_4 = launchUTests('breast', 1)
breast_1_5 = launchUTests('breast', 1)
breast_5_1 = launchUTests('breast', 5)
breast_10_1 = launchUTests('breast', 10)
breast_20_1 = launchUTests('breast', 20)

breast = []
breast = breast_1_1 + breast_1_2 + breast_1_3 + breast_1_4 + breast_1_5 +\
         breast_5_1 + breast_10_1 + breast_20_1

writeFile(breast, 'breast.web.tests.txt')

dlbcl_1_1 = launchUTests('dlbcl', 1)
dlbcl_1_2 = launchUTests('dlbcl', 1)
dlbcl_1_3 = launchUTests('dlbcl', 1)
dlbcl_1_4 = launchUTests('dlbcl', 1)
dlbcl_1_5 = launchUTests('dlbcl', 1)
dlbcl_5_1 = launchUTests('dlbcl', 5)
dlbcl_10_1 = launchUTests('dlbcl', 10)
dlbcl_20_1 = launchUTests('dlbcl', 20)

dlbcl = []
dlbcl = dlbcl_1_1 + dlbcl_1_2 + dlbcl_1_3 + dlbcl_1_4 + dlbcl_1_5 +\
         dlbcl_5_1 + dlbcl_10_1 + dlbcl_20_1

