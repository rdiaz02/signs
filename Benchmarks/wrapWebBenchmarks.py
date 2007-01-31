#!/usr/bin/env python

NUM_SAMPLES = 5

NUM_USERS = (1, 2, 10, 20)

TESTS = ('test1b')


def launchUTests(test, users):
    t = [-99999 for i in range(users)]
    timef = [-99999 for i in range(users)]
    for uu in range(users):
        iin, t[uu] = os.popen2('fl-run-test benchmarkSigns.py Signs.' + test)
    
    for uu in range(users):
        timef[uu] = float(t[uu].readlines()[1].strip())

    return timef

def launchTest(test, users, samples):
    out = []
    for samp in range(samples):
        out = out + launchUTests(test, users)
    return out


def launchAll(test, lusers, samples):
    outall = []
    for users in lusers:
        print 'users ' + str(users)
        print 'samples ' + str(samples)
        outall = outall + launchTest(test, users, samples)
    return outall




ta = launchTest('test1b', 2, 3)    

taa = launchAll('test1b', (1, 2), 3)
