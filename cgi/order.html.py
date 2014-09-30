#!/usr/bin/python

import sys

def permutation_indices(data):
    """ Title: list permutation order indices
    Submitter: Andrew Dalke
    # Based on a post by Peter Otten on comp.lang.python 2004/09/04.
    # This uses the 'key' option to sort and the 'sorted' function, both
    # new in Python 2.4.
    From python cookbook page """
    return sorted(range(len(data)), key = data.__getitem__)

def permutation_indices_reverse(data):
    """ Title: list permutation order indices
    Submitter: Andrew Dalke
    # Based on a post by Peter Otten on comp.lang.python 2004/09/04.
    # This uses the 'key' option to sort and the 'sorted' function, both
    # new in Python 2.4.
    From python cookbook page """
    return sorted(range(len(data)), key = data.__getitem__, reverse = True)

def for_print_p_value(x):
    if x < 1e-7:
        return '< 0.0000001'
    else:
        try:  return str(round(x, 7))
        except: return(x)

def print_warning(x):
    if x == 0:
        return ''
    elif x == 1:
        return 'Coefficient may be infinite'
    elif x == 2:
        return 'Did not converge'

    
def table_gen_sort(l1, l2, l3, l4, l5, l6, order, idtype,
                   organism, fileout = 'tabla.html'):
    fout = open(fileout, mode = 'w')
    fout.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">')
    fout.write('\n<html><head><meta http-equiv="Content-Type" content="text/html; charset=iso-8859-15">')
    fout.write('\n<title>Single gene Cox model p-values and coefficients </title>')
    fout.write('<body> <h1>Single gene Cox model p-values and coefficients </h1>\n')
    fout.write('<table frame="box">')
    outstring = ''.join(['<tr><td width=150>Gene Name</td><td width=150>p-value',
                         '</td><td width=100>Coefficient</td><td width=100>Absolute value of coefficient</td>',
                         '<td width=100>FDR-adjusted p-value</td><td width=150>Warnings?</td></tr>\n'])
    fout.write(outstring)
    for i in range(len(df11)-1):
	if (not(l2[order[i]] == 'NA')) and (l2[order[i]] > 98):
	    outstring = ''.join(['<tr><td>', linkGene(l1[order[i]], idtype, organism),
                            '</td><td>', 'NA',
                            '</td><td>', 'NA',
                            '</td><td>', 'NA',
                            '</td><td>', 'NA', 
                            '</td><td>',print_warning(l6[order[i]]),
                             '</td></tr>\n'])     
	else:
	    outstring = ''.join(['<tr><td>', linkGene(l1[order[i]], idtype, organism),
                                 '</td><td>',for_print_p_value(l2[order[i]]),
                                 '</td><td>',str(round(l3[order[i]], 4)),
                                 '</td><td>',str(round(l4[order[i]], 4)),
                                 '</td><td>',for_print_p_value(l5[order[i]]),
                                 '</td><td>',print_warning(l6[order[i]]),                             
                                 '</td></tr>\n'])
        fout.write(outstring)
    fout.write('</table></body></html>')


def linkGene(geneName, idtype, organism):
    if idtype == 'None' or organism == 'None':
        return geneName
    else:
        return ''.join(['<a href="http://idclight.iib.uam.es/IDClight.prog?idtype=',
                        idtype, '&id=', geneName, '&internal=0&org=',
                        organism,'" target="icl_window">',geneName,'</a>'])


idtype = sys.argv[1]
organism = sys.argv[2]

df1 = open("p.values.coeffs.txt", mode = 'r')
df11 = df1.read().splitlines()

l1=[]
l2=[]
l3=[]
l4=[]
l5=[]
l6=[]

def read_float_other(x):
    try:
        y = float(x)
    except:
        y = x
    return y
    

for i in range(1, len(df11)): ## skip first line
    splitted = df11[i].split('\t')
    l1.append(splitted[0])
    l2.append(read_float_other(splitted[1]))
    l3.append(read_float_other(splitted[2]))
    l4.append(read_float_other(splitted[3]))
    l5.append(read_float_other(splitted[4]))
    l6.append(read_float_other(splitted[5]))

df1.close()


nameAscending = permutation_indices(l1)
nameDescending = permutation_indices_reverse(l1)
pvAscending = permutation_indices(l2)
pvDescending = permutation_indices_reverse(l2)
coefAscending = permutation_indices(l3)
coefDescending = permutation_indices_reverse(l3)
acoeAscending = permutation_indices(l4)
acoeDescending = permutation_indices_reverse(l4)

table_gen_sort(l1, l2, l3, l4, l5, l6, nameAscending, idtype, organism, fileout = "p.v.sort.name.a.html")
table_gen_sort(l1, l2, l3, l4, l5, l6, nameDescending, idtype, organism, fileout ="p.v.sort.name.d.html")   
table_gen_sort(l1, l2, l3, l4, l5, l6, pvAscending, idtype, organism, fileout ="p.v.sort.pv.a.html")     
table_gen_sort(l1, l2, l3, l4, l5, l6, pvDescending, idtype, organism, fileout ="p.v.sort.pv.d.html")    
table_gen_sort(l1, l2, l3, l4, l5, l6, coefAscending, idtype, organism, fileout ="p.v.sort.coef.a.html")   
table_gen_sort(l1, l2, l3, l4, l5, l6, coefDescending, idtype, organism, fileout ="p.v.sort.coef.d.html")   
table_gen_sort(l1, l2, l3, l4, l5, l6, acoeAscending, idtype, organism, fileout ="p.v.sort.abscoef.a.html")
table_gen_sort(l1, l2, l3, l4, l5, l6, acoeDescending, idtype, organism, fileout ="p.v.sort.abscoef.d.html")


