#!/usr/bin/python


import sys
import os
import cgi
import types
import time
import shutil
##import string
import signal
import re
import glob
import tarfile

import cgitb
cgitb.enable() ## zz: eliminar for real work?
## sys.stderr = sys.stdout ## eliminar?

R_MAX_time = 8 * 3600 ## max duration allowd for any process


def issue_echo(fecho, tmpDir):
    """Silly function to output small tracking files"""
    timeHuman = '##########   ' + \
                str(time.strftime('%d %b %Y %H:%M:%S'))
    os.system('echo "' + timeHuman + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "' + fecho + \
              '" >> ' + tmpDir + '/checkdone.echo')
    os.system('echo "    " >> ' + tmpDir + '/checkdone.echo')



def kill_lamcheck(pid, machine):
    'as it says: to kill lamcheck; actually, anything'
    os.system('ssh ' + machine + ' "kill -s 9 ' + pid + '"')

def clean_for_PaLS(file_in, file_out):
    """ Make sure no file has two consecutive lines that start with '#',
    so there are not lists without genes."""
    f1 = open(file_in, mode = 'r').readlines()
    f2 = open(file_out, mode = 'w')
    maxi = len(f1) - 1
    i = 0
    tmp1 = f1[i]
    while True:
        if i == maxi:
            break
        tmp2 = f1[i + 1]
        if not tmp1.startswith('#'):
            f2.write(tmp1)
        elif not tmp2.startswith('#'):
            f2.write(tmp1)
        tmp1 = tmp2
        i += 1

    ### make sure last one is written if not a "#"
    if not tmp2.startswith('#'):
        f2.write(tmp2)
    f2.close()


def extract_for_PaLS_from_Signs(file_in, file_out, all_runs = True):
    """ We should be able to get the output from R directly but:
    a) its much more of a mess, and we need to hack deep inside
    functions; we do a lot with the fitted objects (the function to
    look at is selectedSignatures, also when it is called from
    summary.cvDave);
    b) because of the above, a simple, clean, call to a function
    outside is unlikely to be simple to write;
    c) we need to mess around with the library, etc, which is a pain;
    d) it is a lot harder to debug;
    e) doing it this way, with pattern matchin, etc, is a lot more fun and
    we are just picking up the final, digested, results. Should work.

    That said, this function works on an object (Results.txt) which is the
    output from a call of html2txt, based on results.html. So there are a
    lof ot places where this can bomb if things change in the code.


    The logic is to keep processing lines of outut and:
    - look for where a set of results starts;
    - keep adding the gene names until either that is done, or a new
    component is added (mark the later when it happens).
    - if the set is done, look again for where a set of results starts.

    """

    f1 = open(file_in, mode = 'r').readlines()
    f2 = open(file_out, mode = 'w')

    def can_make_float(x):
        try:
            float(x)
        except:
            return False
        else:
            return True

    def is_start_run(x):
        tmp = x.split()
        if (len(tmp) == 4) and \
           (tmp[0] == 'Component') and \
           (tmp[1] == 'name') and \
           (tmp[2] == 'Genes') and \
           (tmp[3] == 'Coefficient'):
            return True
        else:
            return False

    def another_component(x):
        tmp = x.split()
        if (len(tmp) == 2) and \
           ((tmp[0][0] == 'N') or (tmp[0][0] == 'P')) and \
           can_make_float(tmp[1]):
            return True
        else:
            return False

    def end_run(x):
        tmp = x.split()
        if (len(tmp) == 1 and tmp[0].startswith("==================")):
            return True
        else:
            return False


    i = 0
    maxi = len(f1) - 1

    j = 1
    k = 1
    while True:
        if is_start_run(f1[i]):
            k = 1
            f2.write("#Run." + str(j) + ".component." + str(k) + "\n")
            i = i + 2
            ltmp = f1[i]
            while not end_run(ltmp):
                if not another_component(ltmp):
                    f2.write(ltmp.lstrip())
                elif another_component(ltmp):
                    k += 1
                    f2.write("#Run." + str(j) + ".component." + str(k) + "\n")
                i += 1
                ltmp = f1[i]
            if not all_runs: break
            j += 1

        i += 1
        if i >= maxi:
            f2.close()
            break


def printPalsURL(newDir,
                 tmpDir,
                 application_url = "http://signs2.iib.uam.es",
                 f1 = "Selected.genes.txt",
                 f2 = "Selected.and.CV.selected.txt",
                 s1 = "genes selected in all components in main run",
                 s2 = "genes selected in main run and in CV runs (all components in every run)"):
    """ Based on Pomelo II's Send_to_Pals.cgi."""
    f=open(tmpDir + "/idtype")
    idtype = f.read().strip()
    f.close()
    f=open(tmpDir + "/organism")
    organism = f.read().strip()
    f.close()
    if (idtype != "None" and organism != "None"):
        url_org_id = "org=" + organism + "&idtype=" + idtype + "&"
    else:
        url_org_id = ""
    gl_base = application_url + '/tmp/' + newDir + '/'
    gl1 = gl_base + f1
    gl2 = gl_base + f2

    clean_for_PaLS(tmpDir + '/' + f1, tmpDir + '/' + f1)
    clean_for_PaLS(tmpDir + '/' + f2, tmpDir + '/' + f2)

    outstr0 = '<br /> <hr> ' + \
              '<h3> Send results to <a href = "http://pals.iib.uam.es">' + \
              '<IMG BORDER="0" SRC="../../palsfavicon40.png" align="middle"></a></h3>'
    outstr = outstr0 + \
             '<p> Send set of <a href="http://pals.iib.uam.es?' + \
             url_org_id + 'datafile=' + gl1 + \
             '">' + s1 + ' to PaLS</a></p>' + \
             '<p> Send set of <a href="http://pals.iib.uam.es?' + \
             url_org_id + 'datafile=' + gl2 + \
             '">' + s2 + ' to PaLS</a></p>'
    return(outstr)




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


# def getPathinfo():
#     """ Return the remaining part of the URL. """
#     pathinfo = os.environ.get('PATH_INFO', '')
#     return pathinfo

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

## to keep executing myself:
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
def printErrorRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>SignS results </title></head><body>\n")
    outf.write("<h1> ERROR: There was a problem with the R code </h1> \n")
    outf.write("<p>  This could be a bug on our code, or a problem  ")
    outf.write("with your data (that we hadn't tought of). Below is all the output from the execution ")
    outf.write("of the run. Unless it is obvious to you that this is a fault of your data ")
    outf.write("(and that there is no way we could have avoided the crash) ")
    outf.write("please let us know so we can fix the problem. ")
    outf.write("Please sed us this URL and the output below</p>")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("<p> And this the error messages file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(error.msg))
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

        if (methodUsed == 'TGD') or (methodUsed == 'TGD\n'):
            outf.write("<h2> Results using the Threshold Gradient Descent method of Li and Gui</h2><br/ >\n")
            outf.write('<h2>1. Survival plots \n')
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

            ll1 = glob.glob('*.log')
            for dname in ll1:
                try: os.remove(dname)
                except: None

            lll = glob.glob('*')
            for flname in lll:
                allResults.add(flname)

            allResults.close()
            outf.write('<hr> <a href="all.results.tar.gz">Download</a> all figures and text results.')
            extract_for_PaLS_from_Signs('Results.txt',
                                        'Selected.genes.txt',
                                        all_runs = False)
            extract_for_PaLS_from_Signs('Results.txt',
                                        'Selected.and.CV.selected.txt',
                                        all_runs = True)
            outf.write(printPalsURL(newDir, tmpDir))
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
###     outf.write("<p> This is the output from the R run:<p>")
###     outf.write("<pre>")
###     outf.write(cgi.escape(soFar))
###     outf.write("</pre>")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")



## Changing to the appropriate directory

form = cgi.FieldStorage()
if form.has_key('newDir'):
   value=form['newDir']
   if type(value) is types.ListType:
       commonOutput()
       print "<h1> ERROR </h1>"
       print "<p> newDir should not be a list. </p>"
       print "<p> Anyone trying to mess with it?</p>"
       print "</body></html>"
       sys.exit()
   else:
       newDir = value.value
else:
    commonOutput()
    print "<h1> ERROR </h1>"
    print "<p> newDir is empty. </p>"
    print "</body></html>"
    sys.exit()

if re.search(r'[^0-9]', str(newDir)):
## newDir can ONLY contain digits.
    commonOutput()
    print "<h1> ERROR </h1>"
    print "<p> newDir does not have a valid format. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()

redirectLoc = "/tmp/" + newDir
tmpDir = "/asterias-web-apps/signs2/www/tmp/" + newDir

if not os.path.isdir(tmpDir):
    commonOutput()
    print "<h1> ERROR </h1>"
    print "<p> newDir is not a valid directory. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()



lam_check = open(tmpDir + '/lamCheckPID', mode = 'r'). readline().split()
lam_check_machine = lam_check[1]
lam_check_pid = lam_check[0]


## Were we already done in a previous execution?
## No need to reopen files or check anything else. Return url with results
## and bail out.
if os.path.exists(tmpDir + "/natural.death.pid.txt") or os.path.exists(tmpDir + "/killed.pid.txt"):
    print 'Location: http://signs2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'
    try:
        kill_lamcheck(lam_check_pid, lam_check_machine)
    except:
        None
    sys.exit()

## No, we were not done. Need to examine R output and possible crashes



Rrout = open(tmpDir + "/f1.Rout")
soFar = Rrout.read()
Rrout.close()
finishedOK = soFar.endswith("Normal termination\n")
errorRun = soFar.endswith("Execution halted\n")

issue_echo('right after Rrout', tmpDir)

if os.path.exists(tmpDir + '/RterminatedOK'):
    finishedOK = True

issue_echo('right after checking Rterminated', tmpDir)

## if os.path.exists(tmpDir + "/pid.txt"):
if (not finishedOK) and (not errorRun) and (os.path.exists(tmpDir + "/pid.txt")):
    ## do we need to kill an R process?
    issue_echo('did we run out of time?', tmpDir)
    if (time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
        try:
            os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                      '; lamhalt -H; lamwipe -H')
        except:
            None

        kill_lamcheck(lam_check_pid, lam_check_machine)
        printRKilled()
        os.rename(tmpDir + '/pid.txt', tmpDir + '/killed.pid.txt')
##        os.remove(tmpDir + '/f1.R')
        try:
            os.system("rm /asterias-web-apps/signs2/www/R.running.procs/R." + newDir + "*")
        except:
            None
        print 'Location: http://signs2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'
        sys.exit()

if errorRun > 0:
    issue_echo('errorRun is 1', tmpDir)
    kill_lamcheck(lam_check_pid, lam_check_machine)
    printErrorRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv + '; lamhalt -H; lamwipe -H')
    except:
        None
    try:
        os.system("rm /asterias-web-apps/signs2/www/R.running.procs/R." + newDir + "*")
    except:
        None
    print 'Location: http://signs2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'


elif finishedOK > 0:
    issue_echo('finishedOK is 1', tmpDir)

    kill_lamcheck(lam_check_pid, lam_check_machine)
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        lamkill = os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv + '; lamhalt -H; lamwipe -H')
    except:
        None
    printOKRun()
    try: os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    except: None
    try:
        os.system("rm /asterias-web-apps/signs2/www/R.running.procs/R." + newDir + "*")
    except:
        None
    print 'Location: http://signs2.iib.uam.es/tmp/'+ newDir + '/results.html \n\n'


else:
    issue_echo('relaunching', tmpDir)

    ## we only end up here if: we were not done in a previous run AND no process was overtime
    ## AND we did not just finish. So we must continue.
    relaunchCGI()


issue_echo('END of checkdone.cgi', tmpDir)
