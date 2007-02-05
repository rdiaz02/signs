#!/usr/bin/env python

import os
import sys

## outappend = sys.argv[1]

outappend = 'singlerun'

NUM_USERS = (1, 1, 1, 1, 1, 
##             2, 2, 2,
##             5,
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
    

dlbcl = launchAll('dlbcl', NUM_USERS)
writeFile(dlbcl, 'dlbcl.web.bnchmk.' + outappend + 'txt')


breast = launchAll('breast', NUM_USERS)
writeFile(breast, 'breast.web.bnchmk.' + outappend + 'txt')


## 2362.9028670799999, 2368.0365898599998: dlbcl with 1 sample
