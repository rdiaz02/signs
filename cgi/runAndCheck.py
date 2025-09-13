#!/usr/bin/python

## All this code is copyright Ramon Diaz-Uriarte.
## Released under the Affero GPL.


## This file is very similar to the one in ADaCGH, but with minor adaptations to signs.

import sys
import os
import cgi
import time
import shutil
import glob
import random
import socket
##import fcntl

import tarfile
# import cgitb; cgitb.enable()
# sys.stderr = sys.stdout

sys.path = sys.path + ['/home2/ramon/web-apps/web-apps-common']
import counterApplications
from web_apps_config import *

tmpDir = sys.argv[1]
ROOT_TMP_DIR = "/home2/ramon/web-apps/signs2/www/tmp/"
newDir = tmpDir.replace(ROOT_TMP_DIR, "")
runningProcs = tmpDir.split('/tmp/')[0] + '/R.running.procs/'


## I think we no longer check tmpDir is OK, because this is not launched
## by the user, byt by the signsR.cgi file.




## procTable = tmpDir.split('/tmp/')[0] + '/R.running.procs/procTable'

## Must ensure the procTable exists and has a valid value
## No longer used
# if not os.path.exists(procTable):
#     fo = open(procTable, mode = 'w')
#     fo.write('0')
#     fo.close()

# def set_defaults_lam(tmpDir):
#     """ Set defaults for lamboot and Rslaves and number procs
#     based on the size of the data file. This is all heuristic,
#     but works for us with 6 GB RAM per node. The key idea is to
#     prevent swapping. ncpu are the Rslaves spawned by lamd, or the cpu=ncpu
#     in the lamb-host file. max_num_procs is the maximum number of simultaneous
#     adacgh processes running at any time.
#     We return the tuple ncpu, max_num_procs"""
#     return(1, 5)
# ## In contrast to ADaCGH, we do not see swapping here, so no need to do
# ## size-dependent lam boots. FIXME: cforest and glmboost can use a lot of RAM ...
# ## Maybe enable again? Use (1, something); the 1, to launch each node
# ## only one process (though different runs could use the same node).

# ## All methods, except, TGD, use only 10 spawned procs. Maybe use this to
# ## optimize how to allocate procs. to nodes, etc?

# #     datsize1 = 0
# #     datsize2 = 0
# #     if os.path.exists(tmpDir + '/acghData'):
# #         datsize1 = int(os.popen('ls ' + tmpDir + '/acghData -sk').read().split()[0])
# #     if os.path.exists(tmpDir + '/acghAndPosition'):
# #         datsize2 = int(os.popen('ls ' + tmpDir + '/acghAndPosition -sk').read().split()[0])
# #     datsize = max(datsize2, datsize1)
# #     if datsize < 2000:
# #         return (2, 3)
# #     elif datsize < 6000:
# #         return (2, 2)
# #     elif datsize < 14000:
# #         return (2, 1)
# #     else:
# #         return (1, 1)


def collectZombies(k = 10):
    """ Make sure there are no zombies in the process tables.
    This is probably an overkill, but works.
    """
    for nk in range(k):
        try:
            tmp = os.waitpid(-1, os.WNOHANG)
        except:
            None


def issue_echo(fecho, tmpDir):
    """Silly function to output small tracking files"""
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S'))
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "' + fecho + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "    " >> ' + tmpDir + '/checkdone.echo')


def kill_pid_machine(pid):
    'as it says: to kill somehting somewhere'
    os.system("kill -s 9 " + pid )




# def clean_for_PaLS(file_in, file_out):
#     """ Make sure no file has two consecutive lines that start with '#',
#     so there are no lists without genes."""
#     f1 = open(file_in, mode = 'r').readlines()
#     f2 = open(file_out, mode = 'w')
#     maxi = len(f1) - 1
#     i = 0
#     tmp1 = f1[i]
#     while True:
#         if i == maxi:
#             break
#         tmp2 = f1[i + 1]
#         if not tmp1.startswith('#'):
#             f2.write(tmp1)
#         elif not tmp2.startswith('#'):
#             f2.write(tmp1)
#         tmp1 = tmp2
#         i += 1

#     ### make sure last one is written if not a "#"
#     if not tmp2.startswith('#'):
#         f2.write(tmp2)
#     f2.close()


# def extract_for_PaLS_from_Signs(file_in, file_out, all_runs = True):
#     """ We should be able to get the output from R directly but:
#     a) its much more of a mess, and we need to hack deep inside
#     functions; we do a lot with the fitted objects (the function to
#     look at is selectedSignatures, also when it is called from
#     summary.cvDave);
#     b) because of the above, a simple, clean, call to a function
#     outside is unlikely to be simple to write;
#     c) we need to mess around with the library, etc, which is a pain;
#     d) it is a lot harder to debug;
#     e) doing it this way, with pattern matchin, etc, is a lot more fun and
#     we are just picking up the final, digested, results. Should work.

#     That said, this function works on an object (Results.txt) which is the
#     output from a call of html2txt, based on results.html. So there are a
#     lof ot places where this can bomb if things change in the code.


#     The logic is to keep processing lines of outut and:
#     - look for where a set of results starts;
#     - keep adding the gene names until either that is done, or a new
#     component is added (mark the later when it happens).
#     - if the set is done, look again for where a set of results starts.

#     """

#     f1 = open(file_in, mode = 'r').readlines()
#     f2 = open(file_out, mode = 'w')

#     def can_make_float(x):
#         try:
#             float(x)
#         except:
#             return False
#         else:
#             return True

#     def is_start_run(x):
#         tmp = x.split()
#         if (len(tmp) == 4) and \
#            (tmp[0] == 'Component') and \
#            (tmp[1] == 'name') and \
#            (tmp[2] == 'Genes') and \
#            (tmp[3] == 'Coefficient'):
#             return True
#         else:
#             return False

#     def another_component(x):
#         tmp = x.split()
#         if (len(tmp) == 2) and \
#            ((tmp[0][0] == 'N') or (tmp[0][0] == 'P')) and \
#            can_make_float(tmp[1]):
#             return True
#         else:
#             return False

#     def end_run(x):
#         tmp = x.split()
#         if (len(tmp) == 1 and tmp[0].startswith("==================")):
#             return True
#         else:
#             return False


#     i = 0
#     maxi = len(f1) - 1

#     j = 1
#     k = 1
#     while True:
#         if is_start_run(f1[i]):
#             k = 1
#             f2.write("#Run." + str(j) + ".component." + str(k) + "\n")
#             i = i + 2
#             ltmp = f1[i]
#             while not end_run(ltmp):
#                 if not another_component(ltmp):
#                     f2.write(ltmp.lstrip())
#                 elif another_component(ltmp):
#                     k += 1
#                     f2.write("#Run." + str(j) + ".component." + str(k) + "\n")
#                 i += 1
#                 ltmp = f1[i]
#             if not all_runs: break
#             j += 1

#         i += 1
#         if i >= maxi:
#             f2.close()
#             break


# def printPalsURL(newDir,
#                  tmpDir,
#                  application_url = "http://signs2.iib.uam.es",
#                  f1 = "Selected.genes.txt",
#                  f2 = "Selected.and.CV.selected.txt",
#                  s1 = "genes selected in all components in main run",
#                  s2 = "genes selected in main run and in CV runs (all components in every run)"):
#     """ Based on Pomelo II's Send_to_Pals.cgi."""
#     f=open(tmpDir + "/idtype")
#     idtype = f.read().strip()
#     f.close()
#     f=open(tmpDir + "/organism")
#     organism = f.read().strip()
#     f.close()
#     if (idtype != "None" and organism != "None"):
#         url_org_id = "org=" + organism + "&idtype=" + idtype + "&"
#     else:
#         url_org_id = ""
#     gl_base = application_url + '/tmp/' + newDir + '/'
#     gl1 = gl_base + f1
#     gl2 = gl_base + f2

#     clean_for_PaLS(tmpDir + '/' + f1, tmpDir + '/' + f1)
#     clean_for_PaLS(tmpDir + '/' + f2, tmpDir + '/' + f2)

#     outstr0 = '<br /> <hr> ' + \
#               '<h3> Send results to <a href = "http://pals.iib.uam.es">' + \
#               '<IMG BORDER="0" SRC="../../palsfavicon40.png" align="middle"></a></h3>'
#     outstr = outstr0 + \
#              '<p> Send set of <a href="http://pals.iib.uam.es?' + \
#              url_org_id + 'datafile=' + gl1 + \
#              '">' + s1 + ' to PaLS</a></p>' + \
#              '<p> Send set of <a href="http://pals.iib.uam.es?' + \
#              url_org_id + 'datafile=' + gl2 + \
#              '">' + s2 + ' to PaLS</a></p>'
#     return(outstr)




def getQualifiedURL(uri = None):
    """ Return a full URL starting with schema, servername and port.

    *uri* -- append this server-rooted uri (must start with a slash)
    """
    schema, stdport = ('http', '80')
    host = os.environ.get('HTTP_HOST')
    if not host:
        host = os.environ.get('SERVER_NAME')
        port = os.environ.get('SERVER_PORT', '80')
        if port != stdport: host = host + ":" + port
    result = "%s://%s" % (schema, host)
    if uri: result = result + uri
    return result

def getScriptname():
    """ Return te scriptname part of the URL."""
    return os.environ.get('SCRIPT_NAME', '')

def getBaseURL():
    """ Return a fully qualified URL to this script. """
    return getQualifiedURL(getScriptname())

def commonOutput():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    <title>SignS results</title>
    </head>
    <body>
    """

def getScriptname():
    """ Return te scriptname part of the URL."""
    return os.environ.get('SCRIPT_NAME', '')


def relaunchCGI():
    issue_echo('inside relaunchCGI', tmpDir)

    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    """
    print '<meta http-equiv="Refresh"'
    print 'content="30; URL=' + getBaseURL() + '?newDir=' + newDir + '">'
    print '<title>SignS results</title>'
    print '</head> <body>'
    print '<p> This is an autorefreshing page; your results will eventually be displayed here.\n'
    print 'If your browser does not autorefresh, the results will be kept for five days at</p>'
    print '<p><a href="' + getBaseURL() + '?newDir=' + newDir + '">', 'http://signs2.iib.uam.es/tmp/'+ newDir + '/results.html</a>.'
    print '</p> </body> </html>'
    issue_echo('end of relaunchCGI', tmpDir)





## Output-generating functions

## Output-generating functions
def printErrorRun(errorfile):
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    errormsg = open(errorfile).read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>SignS results </title></head><body>\n")
    outf.write("<h1> ERROR: There was a problem with the R code </h1> \n")
    outf.write("<p>  This could be a bug on our code, or a problem  ")
    outf.write("with your data (that we hadn't tought of). Below is all the output from the execution ")
    outf.write("of the run. Unless it is obvious to you that this is a fault of your data ")
    outf.write("(and that there is no way we could have avoided the crash) ")
    outf.write("please let us know so we can fix the problem. ")
    outf.write("Please sed us this URL and the output below</p>")
    outf.write("<p> This is the results file:</p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("<p> And this the error messages file:</p>")
    outf.write("<pre>")
    outf.write(cgi.escape(errormsg))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



def printOKRun():
    issue_echo('starting printOKRun', tmpDir)
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
#    outf.write("Content-type: text/html\n\n")

    outf.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">')
    outf.write('\n<html><head>')
    outf.write('\n <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-15">')
    outf.write('\n <SCRIPT type="text/javascript" SRC="../../aqtree3clickable.js"></SCRIPT> ')
    outf.write('\n <LINK REL="stylesheet" HREF="../../aqtree3clickable.css"> ')
    outf.write('\n <LINK REL="stylesheet" HREF="../../style1.css"> ')

    outf.write("\n <title>SignS results </title></head><body>\n")

    if os.path.exists(tmpDir + "/ErrorFigure.png"):
        outf.write('<IMG BORDER="0" SRC="ErrorFigure.png">')
        outf.write("<br /><br /> <hr>")
        outf.write("<pre>")
        outf.write('<br /><br /><h2> Results <a href="http://signs2.iib.uam.es/help/signs-help.html#outputText">(help)</a></h2> \n')
        outf.write("<br /><br /> <hr>")
        outf.write(cgi.escape(resultsFile))
        outf.write("</pre>")
        outf.write("</body></html>")
        outf.close()
        Rresults.close()
        shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
    else:
        methodUsed = open(tmpDir + '/methodSurv').read()
##        listPNGS = glob.glob(tmpDir + "/*.png")


        if(methodUsed == 'cforest') or (methodUsed == 'cforest\n') or (methodUsed == 'glmboost') or (methodUsed == 'glmboost\n'):
            if(methodUsed == 'cforest') or (methodUsed == 'cforest\n'):
                outf.write("<h2> Results using Random forests (Hothorn et al., 2006a)</h2><br/ >\n")
            if (methodUsed == 'glmboost') or (methodUsed == 'glmboost\n'):
                outf.write("<h2> Results using boosting of component-wise Cox models (Hothorn et al., 2006b)</h2><br/ >\n")
            outf.write('<h2>1. Survival plots</h2>\n')
            outf.write('<h3>1.1. Survival plots using scores from final model <a href="http://signs2.iib.uam.es/help/signs-help.html#outKM">(help)</a></h3> \n')
            outf.write('<h4>Two groups</h4>')
            outf.write('<IMG BORDER="0" SRC="kmplot-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot-overfitt.png">')
            outf.write('<p>Please, DO NOT use the overfitt one, except for pedagogic purposes to show ')
            outf.write('consequencues of overfitting and not doing cross-validation.</p>')
            outf.write('<h4>Three groups</h4>\n')
            outf.write('<IMG BORDER="0" SRC="kmplot3-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot3-overfitt.png">')
            outf.write('<h4>Four groups</h4>')
            outf.write('<IMG BORDER="0" SRC="kmplot4-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot4-overfitt.png">')

            if os.path.exists(tmpDir + "/usevalidation"):
                outf.write('<h3>1.2. Survival curves for validation data <a href="http://signs2.iib.uam.es/help/signs-help.html#outKM">(help)</a></h3> \n')
                outf.write('<h4>Two groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot-validation.png">')
                outf.write('<h4>Three groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot3-validation.png">')
                outf.write('<h4>Four groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot4-validation.png">')

	    outf.write('<h2>2. Single gene for the chosen number of genes</h2>\n')
            outf.write('<ul class="aqtree3clickable">\n')
	    outf.write('<li><a>Single-gene Cox model p-values and statistics</a><ul>\n')
	    outf.write('<li><a href="p.v.sort.pv.a.html" target="pv_window">Sorted by p-value, ascending</a>')
	    outf.write('<li><a href="p.v.sort.pv.d.html" target="pv_window">Sorted by p-value, descending</a>')
	    outf.write('<li><a href="p.v.sort.coef.d.html" target="pv_window">Sorted by coefficient, descending</a>')
	    outf.write('<li><a href="p.v.sort.coef.a.html" target="pv_window">Sorted by coefficient, ascending</a>')
	    outf.write('<li><a href="p.v.sort.abscoef.a.html" target="pv_window">Sorted by absolute value of coefficient, ascending</a>')
	    outf.write('<li><a href="p.v.sort.abscoef.d.html" target="pv_window">Sorted by absolute value of coefficient, descending</a>')
	    outf.write('<li><a href="p.v.sort.name.d.html" target="pv_window">Sorted by name, descending</a>')
	    outf.write('<li><a href="p.v.sort.name.a.html" target="pv_window">Sorted by name, ascending</a></ul></ul></ul>')

            outf.write("<br /><br /> <hr>")
            outf.write(resultsFile)
            Rresults.close()
            if os.path.exists(tmpDir + '/results.txt'): os.remove(tmpDir + '/results.txt')
            allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
            os.chdir(tmpDir)
            os.system('html2text -width 200 -nobs -o Results.txt pre-results.html')
            lll = glob.glob('*')
            for flname in lll:
                allResults.add(flname)
            allResults.close()
            outf.write('<hr> <a href="all.results.tar.gz">Download</a> all figures and text results.')
            # outf.write(printPalsURL(newDir, tmpDir,
            #                         s1 = "Genes selected in run with all data",
            #                         s2 = "Genes selected in run with all data and in CV runs"))
            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")

        if (methodUsed == 'TGD') or (methodUsed == 'TGD\n'):
            outf.write("<h2> Results using the Threshold Gradient Descent method of Li and Gui</h2><br/ >\n")
            outf.write('<h2>1. Survival plots</h2>\n')
            outf.write('<h3>1.1. Survival plots using scores from final model <a href="http://signs2.iib.uam.es/help/signs-help.html#outKM">(help)</a></h3> \n')
            outf.write('<h4>Two groups</h4>')
            outf.write('<IMG BORDER="0" SRC="kmplot-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot-overfitt.png">')
            outf.write('<p>Please, DO NOT use the overfitt one, except for pedagogic purposes to show ')
            outf.write('consequencues of overfitting and not doing cross-validation.</p>')
            outf.write('<h4>Three groups</h4>\n')
            outf.write('<IMG BORDER="0" SRC="kmplot3-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot3-overfitt.png">')
            outf.write('<h4>Four groups</h4>')
            outf.write('<IMG BORDER="0" SRC="kmplot4-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot4-overfitt.png">')

            if os.path.exists(tmpDir + "/usevalidation"):
                outf.write('<h3>1.2. Survival curves for validation data <a href="http://signs2.iib.uam.es/help/signs-help.html#outKM">(help)</a></h3> \n')
                outf.write('<h4>Two groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot-validation.png">')
                outf.write('<h4>Three groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot3-validation.png">')
                outf.write('<h4>Four groups</h4>')
                outf.write('<IMG BORDER="0" SRC="kmplot4-validation.png">')

            outf.write('<br /> <br /><h2>2. Cross-validated partial likelihood')
            outf.write('<a href="http://signs2.iib.uam.es/help/signs-help.html#out.tgd">(help)</a></h2> \n')
            outf.write('<IMG BORDER="0" SRC="cvpl.png">')
            outf.write('<p>(Use this plot to asses if the chosen values of maximum iterations and &#8710;&#951; worked appropriately.)</p>')

            outf.write('<br /> <br /><h2>3. Coefficients of genes (at best cross-validated partial likelihood) ')
            outf.write('<a href="http://signs2.iib.uam.es/help/signs-help.html#out.tgd">(help)</a></h2> \n')
            outf.write('<IMG BORDER="0" SRC="fstdgrun.png">')

            outf.write("<br /><br /> <hr>")
            outf.write("<pre>")
#             outf.write('<br /><br /><h2> Results <a href="http://signs2.iib.uam.es/help/signs-help.html#outputText">(help)</a></h2> \n')
##            outf.write(cgi.escape(resultsFile))
            outf.write(resultsFile)
            outf.write("</pre>")

            allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
#            allResults.add(tmpDir + '/results.txt', 'results.txt')
            os.chdir(tmpDir)
            os.system('html2text -width 200 -nobs -o Results.txt pre-results.html')
            allResults.add(tmpDir + '/results.txt', 'Results.txt')

            if os.path.exists(tmpDir + "/kmplot-honest.png"): allResults.add(tmpDir + '/kmplot-honest.png', 'SurvivalPlot-honest.png')
            if os.path.exists(tmpDir + "/kmplot-overfitt.png"): allResults.add(tmpDir + '/kmplot-overfitt.png', 'SurvivalPlot-overfitt.png')
            if os.path.exists(tmpDir + "/fstdgrun.png"): allResults.add(tmpDir + '/fstdgrun.png', 'GeneCoefficients.png')
            if os.path.exists(tmpDir + "/cvpl.png"): allResults.add(tmpDir + '/cvpl.png', 'Cross.val.partial.likelihood.png')
            if os.path.exists(tmpDir + "/kmplot-honest.pdf"): allResults.add(tmpDir + '/kmplot-honest.pdf', 'SurvivalPlot-honest.pdf')
            if os.path.exists(tmpDir + "/kmplot-overfitt.pdf"): allResults.add(tmpDir + '/kmplot-overfitt.pdf', 'SurvivalPlot-overfitt.pdf')
            if os.path.exists(tmpDir + "/fstdgrun.pdf"): allResults.add(tmpDir + '/fstdgrun.pdf', 'GeneCoefficients.pdf')
            if os.path.exists(tmpDir + "/cvpl.pdf"): allResults.add(tmpDir + '/cvpl.pdf', 'Cross.val.partial.likelihood.pdf')
            if os.path.exists(tmpDir + "/genes.all.out"): allResults.add(tmpDir + '/genes.all.out', 'genes.all.out')

            if os.path.exists(tmpDir + "/kmplot4-honest.png"): allResults.add(tmpDir + '/kmplot4-honest.png', 'SurvivalPlot4-honest.png')
            if os.path.exists(tmpDir + "/kmplot4-overfitt.png"): allResults.add(tmpDir + '/kmplot4-overfitt.png', 'SurvivalPlot4-overfitt.png')
            if os.path.exists(tmpDir + "/kmplot4-honest.pdf"): allResults.add(tmpDir + '/kmplot4-honest.pdf', 'SurvivalPlot4-honest.pdf')
            if os.path.exists(tmpDir + "/kmplot4-overfitt.pdf"): allResults.add(tmpDir + '/kmplot4-overfitt.pdf', 'SurvivalPlot4-overfitt.pdf')

            if os.path.exists(tmpDir + "/kmplot-validation.png"): allResults.add(tmpDir + '/kmplot-validation.png', 'SurvivalPlot-validation.png')
            if os.path.exists(tmpDir + "/kmplot4-validation.png"): allResults.add(tmpDir + '/kmplot4-validation.png', 'SurvivalPlot4-validation.png')
            if os.path.exists(tmpDir + "/kmplot-validation.pdf"): allResults.add(tmpDir + '/kmplot-validation.pdf', 'SurvivalPlot-validation.pdf')
            if os.path.exists(tmpDir + "/kmplot4-validation.pdf"): allResults.add(tmpDir + '/kmplot4-validation.pdf', 'SurvivalPlot4-validation.pdf')

            allResults.close()
            outf.write('<hr> <a href="all.results.tar.gz">Download</a> all figures and text results.')
            # outf.write(printPalsURL(newDir, tmpDir,
            #                         s1 = "Genes selected in run with all data",
            #                         s2 = "Genes selected in run with all data and in CV runs"))
            outf.write("</body></html>")
            outf.close()
            Rresults.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")

        if (methodUsed == 'FCMS') or (methodUsed == 'FCMS\n'):
            outf.write("<h2> Results using FCMS (filter, cluster, and model selection, as in Dave et al.)</h2> <br /><hr><hr>\n")

	    outf.write('<h2>1. Single gene statistics and p-values</h2>\n')
            outf.write('<ul class="aqtree3clickable">\n')
	    outf.write('<li><a>Single-gene Cox model p-values and statistics</a><ul>\n')
	    outf.write('<li><a href="p.v.sort.pv.a.html" target="pv_window">Sorted by p-value, ascending</a>')
	    outf.write('<li><a href="p.v.sort.pv.d.html" target="pv_window">Sorted by p-value, descending</a>')
	    outf.write('<li><a href="p.v.sort.coef.d.html" target="pv_window">Sorted by coefficient, descending</a>')
	    outf.write('<li><a href="p.v.sort.coef.a.html" target="pv_window">Sorted by coefficient, ascending</a>')
	    outf.write('<li><a href="p.v.sort.abscoef.a.html" target="pv_window">Sorted by absolute value of coefficient, ascending</a>')
	    outf.write('<li><a href="p.v.sort.abscoef.d.html" target="pv_window">Sorted by absolute value of coefficient, descending</a>')
	    outf.write('<li><a href="p.v.sort.name.d.html" target="pv_window">Sorted by name, descending</a>')
	    outf.write('<li><a href="p.v.sort.name.a.html" target="pv_window">Sorted by name, ascending</a></ul></ul></ul>')

            outf.write('<hr><h2>2. Survival plots</h2>')
            outf.write('<h3>2.1. Survival plots using scores from final model</h3>\n')
            outf.write('<h4>Two groups</h4>\n')
            outf.write('<IMG BORDER="0" SRC="kmplot-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot-overfitt.png">')
            outf.write('<p>Please, DO NOT use the overfitt one, except for pedagogic purposes to show ')
            outf.write('consequencues of overfitting and not doing cross-validation.</p>')
            outf.write('<h4>Three groups</h4>\n')
            outf.write('<IMG BORDER="0" SRC="kmplot3-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot3-overfitt.png">')
            outf.write('<h4>Four groups</h4>\n')
            outf.write('<IMG BORDER="0" SRC="kmplot4-honest.png">')
            outf.write('<IMG BORDER="0" SRC="kmplot4-overfitt.png">')

            if os.path.exists(tmpDir + "/usevalidation"):
                outf.write('<h3>2.2. Survival plots for validation data</h3>\n')
                outf.write('<h4>Two groups</h4>\n')
                outf.write('<IMG BORDER="0" SRC="kmplot-validation.png">')
                outf.write('<h4>Three groups</h4>\n')
                outf.write('<IMG BORDER="0" SRC="kmplot3-validation.png">')
                outf.write('<h4>Four groups</h4>\n')
                outf.write('<IMG BORDER="0" SRC="kmplot4-validation.png">')

	    outf.write('<br /> <br /> <hr>')
	    outf.write('<h2>3. Dendrograms of gene clusters</h2>')
            outf.write('<p>(If needed, click on the "+" to expand the list. Then, clik on the figure you want to open. Once a dendrogram is opened, ' +
	    'click on the node name, or close to where the name would be placed, ' +
	    'to see additional information for each gene.)</p>')
	    outf.write('<ul class="aqtree3clickable">\n')
	    outf.write('<li><a>Dendrograms for genes with positive coefficients</a><ul>\n')
	    if os.path.exists(tmpDir + '/NoPositiveCluster'):
		outf.write('<p>There are no genes with positive coefficients that satisfy the p, ' +
		'minimum correlation and size restrictions.</p></ul>')
	    else:
		outf.write('<li><a>All genes named</a><ul>\n')
		outf.write('<li><a href="dend.P.factor0.5.alllabelsTRUE.html" target="dend_window">Half size</a>')
		outf.write('<li><a href="dend.P.factor1.alllabelsTRUE.html" target="dend_window">Normal (1200x800) size</a>')
		outf.write('<li><a href="dend.P.factor2.alllabelsTRUE.html" target="dend_window">Double size</a></ul>')
		outf.write('<li><a>Only cluster genes named</a><ul>\n')
		outf.write('<li><a href="dend.P.factor0.5.alllabelsFALSE.html" target="dend_window">Half size</a>')
		outf.write('<li><a href="dend.P.factor1.alllabelsFALSE.html" target="dend_window">Normal (1200x800) size</a>')
		outf.write('<li><a href="dend.P.factor2.alllabelsFALSE.html" target="dend_window">Double size</a></ul></ul>')
	    outf.write('<li><a>Dendrograms for genes with negative coefficients</a><ul>\n')
	    if os.path.exists(tmpDir + '/NoNegativeCluster'):
		outf.write('<p>There are no genes with negative coefficients that satisfy the p, ' +
		'minimum correlation and size restrictions.</p></ul></ul>')
	    else:
		outf.write('<li><a>All genes named</a><ul>\n')
		outf.write('<li><a href="dend.N.factor0.5.alllabelsTRUE.html" target="dend_window">Half size</a>')
		outf.write('<li><a href="dend.N.factor1.alllabelsTRUE.html" target="dend_window">Normal (1200x800) size</a>')
		outf.write('<li><a href="dend.N.factor2.alllabelsTRUE.html" target="dend_window">Double size</a></ul>')
		outf.write('<li><a>Only cluster genes named</a><ul>\n')
		outf.write('<li><a href="dend.N.factor0.5.alllabelsFALSE.html" target="dend_window">Half size</a>')
		outf.write('<li><a href="dend.N.factor1.alllabelsFALSE.html" target="dend_window">Normal (1200x800) size</a>')
		outf.write('<li><a href="dend.N.factor2.alllabelsFALSE.html" target="dend_window">Double size</a></ul></ul></ul>')

            outf.write("<br /><br /> <hr>")
            outf.write(resultsFile)

            #if os.path.exists(tmpDir + '/f1.R'): os.remove(tmpDir + '/f1.R')
            Rresults.close()
            if os.path.exists(tmpDir + '/results.txt'): os.remove(tmpDir + '/results.txt')

            allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
            os.chdir(tmpDir)
            if os.path.exists('correlationMatrixCluters.html'):
                os.system('html2text -width 200 -nobs -o correlationMatrixClusters.txt correlationMatrixCluters.html')
            if os.path.exists('scores.oob.html'):
                os.system('html2text -width 200 -nobs -o scores.oob.txt scores.oob.html')
            if os.path.exists('scores.validation.html'):
                os.system('html2text -width 200 -nobs  -o scores.validation.txt scores.validation.html')
            if os.path.exists('p.v.sort.name.a.html'):
                os.system('html2text -width 200 -nobs  -o single.gene.cox.coef.p.value.txt p.v.sort.name.a.html')
            outf.flush()
            os.system('html2text -width 200 -nobs -o Results.txt pre-results.html')

#             ll1 = glob.glob('*.log')
#             for dname in ll1:
#                 try: os.remove(dname)
#                 except: None

            lll = glob.glob('*')
            for flname in lll:
                allResults.add(flname)

            allResults.close()
            outf.write('<hr> <a href="all.results.tar.gz">Download</a> all figures and text results.')
            # extract_for_PaLS_from_Signs('Results.txt',
            #                             'Selected.genes.txt',
            #                             all_runs = False)
            # extract_for_PaLS_from_Signs('Results.txt',
            #                             'Selected.and.CV.selected.txt',
            #                             all_runs = True)
            # outf.write(printPalsURL(newDir, tmpDir))
            outf.write("</body></html>")
            outf.close()
            shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")
    issue_echo('end of printOKRun', tmpDir)



def printRKilled():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>SignS results </title></head><body>\n")
    outf.write("<h1> ERROR: R process killed </h1> \n")
    outf.write("<p>  The R process lasted longer than the maximum  allowed time, ")
    outf.write(str(R_MAX_time))
    outf.write(" seconds,  and was killed.")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


# def printMPIerror(tmpDir, numtries, application = 'SignS'):
#     if not os.path.exists('/home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog'):
#         os.system('touch /home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog')
#     outlog = open('/home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog', mode = 'a')
#     outlog.write('MPI fails on ' + time.ctime(time.time()) +
#                  ' Directory: ' + tmpDir + '\n')
#     outlog.close()
#     out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
#     out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
#     out1.write('MPI initialization error!!')
#     out2.write('MPI initialization error!!')
#     out1.close()
#     out2.close()
#     outf = open(tmpDir + "/pre-results.html", mode = "w")
#     outf.write("<html><head><title> MPI initialization problem.</title></head><body>\n")
#     outf.write("<h1> MPI initialization problem.</h1>")
#     outf.write("<p> After " + str(numtries) + " attempts we have been unable to ")
#     outf.write(" initialize MPI.</p>")
#     outf.write("<p> We will be notified of this error, but we would also ")
#     outf.write("appreciate if you can let us know of any circumstances or problems ")
#     outf.write("so we can diagnose the error.</p>")
#     outf.write("</body></html>")
#     outf.close()
#     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


# def printMPITooBusy(tmpDir, MAX_DURATION_TRY, application = 'SignS'):
#     if not os.path.exists('/home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog'):
#         os.system('touch /home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog')
#     outlog = open('/home2/ramon/web-apps/mpi.log/' + application + 'ErrorLog', mode = 'a')
#     outlog.write('Something fails on ' + time.ctime(time.time()) +
#                  ' Directory: ' + tmpDir + '\n')
#     outlog.close()
#     out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
#     out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
#     out1.write('Cannot start!!')
#     out2.write('Cannot start!!')
#     out1.close()
#     out2.close()
#     outf = open(tmpDir + "/pre-results.html", mode = "w")
#     outf.write("<html><head><title> Cannot start application.</title></head><body>\n")
#     outf.write("<h1> Cannot start application.</h1>")
#     outf.write("<p> After " + str(MAX_DURATION_TRY) + " seconds we have been unable to ")
#     outf.write(" start the application.</p>")
#     outf.write("<p> Most likely this means the servers are too busy and many ")
#     outf.write("are running ahead of yours. ")
#     outf.write("Please try again later. You can also get in touch with us ")
#     outf.write("if you think this is our error.</p>")
#     outf.write("</body></html>")
#     outf.close()
#     shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")

def printTooBusy(tmpDir, MAX_DURATION_TRY, application = 'SignS'):
    if not os.path.exists('/home2/ramon/web-apps/log/' + application + 'ErrorLog'):
        os.system('touch /home2/ramon/web-apps/log/' + application + 'ErrorLog')
    outlog = open('/home2/ramon/web-apps/log/' + application + 'ErrorLog', mode = 'a')
    outlog.write('Something fails on ' + time.ctime(time.time()) +
                 ' Directory: ' + tmpDir + '\n')
    outlog.close()
    out1 = open(tmpDir + "/natural.death.pid.txt", mode = "w")
    out2 = open(tmpDir + "/kill.pid.txt", mode = "w")
    out1.write('Cannot start!!')
    out2.write('Cannot start!!')
    out1.close()
    out2.close()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title> Cannot start application.</title></head><body>\n")
    outf.write("<h1> Cannot start application.</h1>")
    outf.write("<p> After " + str(MAX_DURATION_TRY) + " seconds we have been unable to ")
    outf.write(" start the application.</p>")
    outf.write("<p> Most likely this means the servers are too busy and many ")
    outf.write("are running ahead of yours. ")
    outf.write("Please try again later. You can also get in touch with us ")
    outf.write("if you think this is our error.</p>")
    outf.write("</body></html>")
    outf.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



# def lamboot(lamSuffix, ncpu, runningProcs = runningProcs):
#     'Boot a lam universe and leave a sentinel file behind'
#     issue_echo('before sentinel inside lamboot', tmpDir)
#     issue_echo('newDir is ' + newDir, tmpDir)
#     issue_echo('lamSuffix ' + lamSuffix, tmpDir)
#     issue_echo('runningProcs ' + runningProcs, tmpDir)
#     sentinel = os.open(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]),
#                        os.O_RDWR | os.O_CREAT | os.O_NDELAY)
#     issue_echo('before fullCommand inside lamboot', tmpDir)
#     fullCommand = 'export LAM_MPI_SESSION_SUFFIX="' + lamSuffix + \
#                   '"; /home2/ramon/web-apps/mpi.log/tryBootLAM2.py ' + lamSuffix + \
#                   ' ' + str(ncpu)
#     issue_echo('before os.system inside lamboot', tmpDir)
#     lboot = os.system(fullCommand)
#     issue_echo('after lboot ---os.system--- inside lamboot. Exiting lamboot', tmpDir)


# def check_tping(lamSuffix, tmpDir, tsleep = 15, nc = 2):
#     """ Use tping to verify LAM universe OK.
#     tsleep is how long we wait before checking output of tping.
#     Verify also using 'lamexec C hostname' """

#     tmp2 = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
#                      lamSuffix + '"; cd ' + tmpDir + \
#                      '; tping C N -c' + str(nc) + \
#                      ' > tping.out & ')
#     time.sleep(tsleep)
#     tmp = int(os.popen('cd ' + tmpDir + \
#                        '; wc tping.out').readline().split()[0])
#     os.system('rm ' + tmpDir + '/tping.out')
#     timeHuman = '##########   ' + \
#                 str(time.strftime('%d %b %Y %H:%M:%S'))
#     os.system('echo "' + timeHuman + \
#               '" >> ' + tmpDir + '/checkTping.out')
#     if tmp == 0:
#         os.system('echo "tping fails" >> ' + \
#                   tmpDir + '/checkTping.out')
#         return 0
#     elif tmp > 0:
#         os.system('echo "tping OK" >> ' + \
#                   tmpDir + '/checkTping.out')
#         lamexec = os.system('export LAM_MPI_SESSION_SUFFIX="' +\
#                             lamSuffix + '"; lamexec C hostname')
#         if lamexec == 0:
#             os.system('echo "lamexec OK" >> ' + \
#                       tmpDir + '/checkTping.out')
#             return 1
#         else:
#             os.system('echo "lamexec fails" >> ' + \
#                       tmpDir + '/checkTping.out')
#             return 0
#     else:
#         os.system('echo "tping weird ' + str(tmp) + '" >> ' + \
#                   tmpDir + '/checkTping.out')
#         return 0



# def lam_crash_log(tmpDir, value):
#     """ Write to the lam crash log, 'recoverFromLAMCrash.out' """
#     timeHuman = str(time.strftime('%d %b %Y %H:%M:%S'))
#     os.system('echo "' + value + '  at ' + timeHuman + \
#               '" >> ' + tmpDir + '/recoverFromLAMCrash.out')

def generic_crash_log(tmpDir, value):
    """ Write to the lam crash log, 'recoverFromLAMCrash.out' """
    timeHuman = str(time.strftime('%d %b %Y %H:%M:%S'))
    os.system('echo "' + value + '  at ' + timeHuman + \
              '" >> ' + tmpDir + '/recoverFromLAMCrash.out')

# def recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS, lamSuffix,
#                            runningProcs = runningProcs,
#                            machine_root = 'karl'):
#     """Check if lam crashed during R run. If it did, restart R
#     after possibly rebooting the lam universe.
#     Leave a trace of what happened."""

#     os.remove(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]))
#     del_mpi_logs(tmpDir, machine_root)
#     lam_crash_log(tmpDir, 'Crashed')
#     try:
#         os.system('mv ' + tmpDir + '/mpiOK ' + tmpDir + '/previous_mpiOK')
#     except:
#         None

#     check_room = my_queue(MAX_NUM_PROCS)
#     if check_room == 'Failed':
#         printMPITooBusy(tmpDir, MAX_DURATION_TRY = 5 * 3600)

#     lam_ok = check_tping(lamSuffix, tmpDir)
#     if lam_ok == 0:
#         lboot = lamboot(lamSuffix, NCPU)
#     Rrun(tmpDir, lamSuffix)
#     lam_crash_log(tmpDir, '..... recovering')


def Rrun(tmpDir, R_bin):
    """ Launch R, after setting the lam stuff."""
    issue_echo(' inside Rrun ', tmpDir)
    Rcommand = 'cd ' + tmpDir + \
               '; sleep 1; ' + \
               R_bin + ' --no-readline --no-save --slave <f1.R >>f1.Rout 2>> Status.msg &'
    issue_echo('the Rcommand is ' + Rcommand, tmpDir)
    Rtorun = os.system(Rcommand)




def status_run(tmpDir):
    """ Read Status.msg and return status."""
    status_r = open(tmpDir + '/Status.msg').read()
    if status_r.find('Normal termination\n') > -1:
        return('FinishedOK')
    if status_r.find('Execution halted\n') > -1:
        return('Halted')
    if status_r.find('Running\n') > -1:
        return('Running')
    if status_r.find('Rmpi error\n') > -1:
        return('Error_mpi')
    if status_r.find('Run out of time; killed\n') > -1:
        return('Out_of_time')


def did_R_crash_in_slaves(tmpDir, machine_root = 'karl'):
    """ Verify whether R crashed in any of the slaves by
    checking lam logs."""
    R_LAM_MSGS = 'Error:  Error in'
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    in_lam_logs = 0
    for lam_log in lam_logs:
        tmp1 = int(os.popen('grep "' + R_LAM_MSGS + '" ' + \
                            lam_log + ' | wc').readline().split()[0])
        if tmp1 > 0:
            in_lam_logs = 1
            break
    if in_lam_logs > 0:
        return True, lam_log
    else:
        return False, 'NA'



def did_lam_crash(tmpDir, machine_root = 'karl'):
    """ Verify whether LAM/MPI crashed by checking logs and f1.Rout
    for single universe lamboot."""
    issue_echo('          did_lam_crash ?', tmpDir)
    OTHER_LAM_MSGS = 'Call stack within LAM:'
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    in_error_msg = int(os.popen('grep MPI_Error_string ' + \
                                tmpDir + '/Status.msg | wc').readline().split()[0])
#     no_universe = int(os.popen('grep "Running serial version of papply" ' + \
#                                tmpDir + '/f1.Rout | wc').readline().split()[0])
## We do NOT want that, because sometimes a one node universe is legitimate!!!
    if in_error_msg > 0:
        for lam_log in lam_logs:
            os.system('rm ' + lam_log)
#     elif no_universe > 0:
#         os.system("sed -i 's/Running serial version of papply/already_seen:running serial version of papply/g'" + \
#                   tmpDir + "/f1.Rout")
    else: ## look in lam logs
        in_lam_logs = 0
        for lam_log in lam_logs:
            tmp1 = int(os.popen('grep "' + OTHER_LAM_MSGS + '" ' + \
                                lam_log + ' | wc').readline().split()[0])
            if tmp1 > 0:
                in_lam_logs = 1
                break
    if (in_error_msg > 0) or (in_lam_logs > 0):
        return True
    else:
        return False

def did_mpi_crash(tmpDir, machine_root = 'karl'):
    """ Either Rmpi or LAM crashed"""
    issue_echo('    did Rmpi or LAM crash (did_mpi_crash)?', tmpDir)
    if (status_run(tmpDir) == 'Error_mpi') or \
       did_lam_crash(tmpDir, machine_root):
        return True
    else:
        return False

def del_mpi_logs(tmpDir, machine_root = 'karl'):
    """ Delete logs from LAM/MPI."""
    lam_logs = glob.glob(tmpDir + '/' + machine_root + '*.*.*.log')
    try:
        os.system('rm ' + tmpDir + '/Status.msg')
    except:
        None
    try:
        for lam_log in lam_logs:
            os.system('rm ' + lam_log)
    except:
        None

def did_run_out_of_time(tmpDir, R_MAX_time):
    """ Did the process run longer than allowed?"""
    issue_echo('did we run out of time?', tmpDir)
    if not os.path.exists(tmpDir + "/pid.txt"):
        return False
    elif ((time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time) and \
       (status_run(tmpDir) == 'Running'):
        return True
    else:
        return False


def cleanups(tmpDir, newDir, newnamepid,
             runningProcs = runningProcs,
             appl = 'signs2'):
    """ Clean up actions; kill lam, delete running.procs files, clean process table."""
##    lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    rinfo = open(tmpDir + '/current_R_proc_info', mode = 'r').readline().split()

    try:
        kill_pid_machine(rinfo[1])
    except:
        None
    # try:
    #     os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
    #               '; lamhalt -H; lamwipe -H')
    # except:
    #     None
    try:
        os.system('rm /home2/ramon/web-apps/' + appl + '/www/R.running.procs/R.' + newDir + '*')
    except:
        None
    try:
        os.rename(tmpDir + '/pid.txt', tmpDir + '/' + newnamepid)
    except:
        None
    # try:
    #     os.remove(''.join([runningProcs, 'sentinel.lam.', newDir, '.', lamSuffix]))
    # except:
    #     None


def finished_ok(tmpDir):
    """ check ok termination and act accordingly."""
    if status_run(tmpDir) == 'FinishedOK':
        return True
    else:
        return False

def halted(tmpDir):
    """ check halted execution and act accordingly."""
    if status_run(tmpDir) == 'Halted':
        return True
    else:
        return False


def master_out_of_time(time_start):
    """If this process run longer than allowed, kill it and kill lam and R."""
    if (time.time () - time_start) > R_MAX_time:
        return True
    else:
        return False


# def add_to_proc_table(max_num_procs, add_procs = 1):
#     """Try to add add_procs to the process table. If it can
#     returns OK, otherwise (e.g., too many procs) return Failed.
#     Locking would be great ... but it does not work over NFS. """

#     fo = open(procTable, mode = 'r+')
#     fcntl.flock(fo.fileno(), fcntl.LOCK_EX)
#     currentProcs = int(fo.read())
#     if currentProcs >= max_num_procs:
#         fcntl.flock(fo.fileno(), fcntl.LOCK_UN)
#         fo.close()
#         return 'Failed'
#     else:
#         fo.seek(0)
#         fo.write(str(currentProcs + add_procs))
#         fcntl.flock(fo.fileno(), fcntl.LOCK_UN)
#         fo.close()
#         return 'OK'

# def add_to_proc_table(max_num_procs, add_procs = 1):
#     """Try to add add_procs to the process table. If it can
#     returns OK, otherwise (e.g., too many procs) return Failed."""
#     fo = open(procTable, mode = 'r+')
#     currentProcs = int(fo.read())
#     if currentProcs >= max_num_procs:
#         fo.close()
#         return 'Failed'
#     else:
#         fo.seek(0)
#         fo.write(str(currentProcs + add_procs))
#         fo.close()
#         return 'OK'




def del_from_proc_table(del_procs = 1):
    """Decrease count in the process table."""
    fo = open(procTable, mode = 'r+')
    currentProcs = int(fo.read())
    fo.seek(0)
    fo.write(str(currentProcs - del_procs))
    fo.close()
    return 'OK'


# def my_queue(MAX_NUM_PROCS,
#              runningProcs = runningProcs,
#              ADD_PROCS = 1,
#              CHECK_QUEUE = 23,
#              MAX_DURATION_TRY = 25 * 3600):
#     """ Wait here until the number of processes is smaller than
#     MAX_NUM_PROCS and number of slaves smaller than MAX_NUM_PROCS + ADD_PROCS
#     (so we allow for other apps. launching lamd).
#     But only wait for MAX_DURATION. Check
#     every CHECK_QUEUE seconds. If able to find an opening, return
#     OK, otherwise return Failed"""
#     out_value = 'OK'
#     startTime = time.time()
#     while True:
#         killedlamandr = os.system('/home2/ramon/web-apps/mpi.log/killOldLamAllMachines.py')
#         issue_echo('     inside my_queue ', tmpDir)
#         if (time.time() - startTime) > MAX_DURATION_TRY:
#             out_value = 'Failed'
#             break
#         num_lamd = int(os.popen('pgrep -u www lamd | wc').readline().split()[0])
#         num_sentinel = int(len(glob.glob(''.join([runningProcs, 'sentinel.lam.*']))))
#         if (num_lamd < (MAX_NUM_PROCS + ADD_PROCS)) and (num_sentinel < MAX_NUM_PROCS):
#             issue_echo('     OK; num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             break
#         else:
# 	    issue_echo('     wait:  num_lamd = ' + str(num_lamd) + \
#                        '; num_sentinel = ' + str(num_sentinel), tmpDir)
#             time.sleep(CHECK_QUEUE + random.uniform(0.1, 5))
#     return out_value

def my_queue(MAX_NUM_PROCS,
             runningProcs = runningProcs,
             ADD_PROCS = 1,
             CHECK_QUEUE = 23,
             MAX_DURATION_TRY = MAX_DURATION_TRY_Signs):
    """ Wait here until the number of processes is smaller than
    MAX_NUM_PROCS
    But only wait for MAX_DURATION. Check
    every CHECK_QUEUE seconds. If able to find an opening, return
    OK, otherwise return Failed"""
    out_value = 'OK'
    startTime = time.time()
    while True:
#        killedlamandr = os.system('/home2/ramon/web-apps/mpi.log/killOldLamAllMachines.py')
        issue_echo('     inside my_queue ', tmpDir)
        if (time.time() - startTime) > MAX_DURATION_TRY:
            out_value = 'Failed'
            break
 #       num_lamd = int(os.popen('pgrep -u www lamd | wc').readline().split()[0])
        num_sentinel = int(len(glob.glob(''.join([runningProcs, 'sentinel.lam.*']))))
        if (num_sentinel < MAX_NUM_PROCS):
            issue_echo(' num_sentinel = ' + str(num_sentinel), tmpDir)
            break
        else:
	    issue_echo('     wait:  ' + \
                       '; num_sentinel = ' + str(num_sentinel), tmpDir)
            time.sleep(CHECK_QUEUE + random.uniform(0.1, 5))
    return out_value


# def generate_lam_suffix(tmpDir):
#     """As it says. Generate and write it out"""
#     lamSuffix = str(int(time.time())) + \
#                 str(os.getpid()) + str(random.randint(10, 999999))
#     lamenvfile = open(tmpDir + '/lamSuffix', mode = 'w')
#     lamenvfile.write(lamSuffix)
#     lamenvfile.flush()
#     lamenvfile.close()
#     return lamSuffix



## Starting. First, the very first run.

issue_echo('starting', tmpDir)


## NCPU, MAX_NUM_PROCS = set_defaults_lam(tmpDir)

try:
    counterApplications.add_to_counter_log('Signs2', tmpDir, socket.gethostname())
except:
    None

issue_echo('at 2', tmpDir)

## lamSuffix = generate_lam_suffix(tmpDir)

issue_echo('at 3', tmpDir)

time.sleep(random.uniform(0.01, 1 )) ## Break ties if starting at identical times

check_room = my_queue(MAX_signs,
                      MAX_DURATION_TRY = MAX_DURATION_TRY_Signs)

issue_echo('after check_room', tmpDir)

if check_room == 'Failed':
    printTooBusy(tmpDir, MAX_DURATION_TRY_Signs)
    sys.exit()

issue_echo('before lamboot', tmpDir)


# lamboot(lamSuffix, NCPU)
# issue_echo('after lamboot', tmpDir)

# counterApplications.add_to_LAM_SUFFIX_LOG(lamSuffix, 'SignS2', tmpDir,
#                                           socket.gethostname())

issue_echo('before  Rrun', tmpDir)

Rrun(tmpDir, R_bin)

issue_echo('after Rrun', tmpDir)

time_start = time.time()
time.sleep(TIME_BETWEEN_CHECKS + random.uniform(0.1, 3))

count_mpi_crash = 0

while True:
    if did_run_out_of_time(tmpDir, R_MAX_time):
        issue_echo('run out of time', tmpDir)
        cleanups(tmpDir, newDir, 'killed.pid.txt')
        printRKilled()
        break
    elif finished_ok(tmpDir):
        issue_echo('finished OK', tmpDir)
        cleanups(tmpDir, newDir, 'natural.death.pid.txt')
        printOKRun()
        break
    elif halted(tmpDir):
        issue_echo('halted', tmpDir)
        cleanups(tmpDir, newDir, 'natural.death.pid.txt')
        printErrorRun(tmpDir + '/Status.msg')
        break
    # elif did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[0]:
    #     issue_echo('R crash in slaves', tmpDir)
    #     cleanups(tmpDir, newDir, 'natural.death.pid.txt')
    #     printErrorRun(did_R_crash_in_slaves(tmpDir, machine_root = 'karl')[1])
    #     break
    elif master_out_of_time(time_start):
        issue_echo('master out of time', tmpDir)
        cleanups(tmpDir, newDir, 'killed.pid.txt')
        printRKilled()
        break
    # elif did_mpi_crash(tmpDir, machine_root = 'karl'):
    #     count_mpi_crash += 1
    #     if count_mpi_crash > MAX_MPI_CRASHES:
    #         printMPIerror(tmpDir, MAX_MPI_CRASHES)
    #         cleanups(tmpDir, newDir, 'MPIerror.pid.txt')
    #         break
    #     else:
    #         recover_from_lam_crash(tmpDir, NCPU, MAX_NUM_PROCS,
    #                                lamSuffix,
    #                                machine_root = 'karl')
    else:
        generic_crash_log(tmpDir, 'NoCrash') ## if we get here, this much we know
    time.sleep(TIME_BETWEEN_CHECKS)



issue_echo('at the very end!', tmpDir)
