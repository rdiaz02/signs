#!/usr/bin/python

## when ready, turn signs2 into signs
## or viceversa


import sys
import os

direction = sys.argv[1]

if direction=='signs2':
    os.system("sed 's/signs.bioinfo/signs2.bioinfo/g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/signs.bioinfo/signs2.bioinfo/g' signsR.cgi > tmp; mv tmp signsR.cgi")
    os.system("sed 's/signs.bioinfo/signs2.bioinfo/g' ../www/signs.html > tmp; mv tmp ../www/signs.html")
    os.system("sed 's/\/http\/signs\//\/http\/signs2\//g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/\/http\/signs\//\/http\/signs2\//g' signsR.cgi > tmp; mv tmp signsR.cgi")
    os.system("sed 's/\/http\/signs\//\/http\/signs2\//g' ../www/signs.html > tmp; mv tmp ../www/signs.html")
    
if direction=='signs':
    os.system("sed 's/signs2.bioinfo/signs.bioinfo/g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/signs2.bioinfo/signs.bioinfo/g' signsR.cgi > tmp; mv tmp signsR.cgi")
    os.system("sed 's/signs2.bioinfo/signs.bioinfo/g' ../www/signs.html > tmp; mv tmp ../www/signs.html")
    os.system("sed 's/\/http\/signs2\//\/http\/signs\//g' checkdone.cgi > tmp; mv tmp checkdone.cgi")
    os.system("sed 's/\/http\/signs2\//\/http\/signs\//g' signsR.cgi > tmp; mv tmp signsR.cgi")
    os.system("sed 's/\/http\/signs2\//\/http\/signs\//g' ../www/signs.html > tmp; mv tmp ../www/signs.html")

os.system('chmod u+x /asterias-web-apps/signs2/cgi/*.cgi')    
os.system('chown -R www-data /asterias-web-apps/signs2')
