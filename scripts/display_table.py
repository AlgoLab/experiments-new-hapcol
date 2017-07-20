usage = '''

usage: python display_table.py [-d] [-b] [-s] [-x] [file]

display a table annotated by study_table.py to stdout, slack, ..

positional arguments:
  file        input file (default: <stdin>)

optional arguments:
  -h, --help  show this help message and exit
  -d          whether we distinguish marked entries or not
  -b          output marked entries in blue for stdout
  -s          output marked entries in slack bold
  -x          output latex-formatted table (marked entries in bold)

note: only one display mode can be activated at a time, and if no
display mode is activated, it removes annotation from the table

'''

import sys

#
# format line for latex output
def latex_format(line, count) :

    line = line.replace('^','\\textbf{').replace('$','}')
    s = line.split()

    if not count : # must be at the header
        print('\\begin{tabular}{||c||'+'r'*(len(s)-1)+'||}')
        print('\\hline')
        s[0] = '\\emph{'+s[0].lstrip('#').replace('\\','\\textbackslash ')+'}'
        print(' & '.join(s), '\\\\')
        print('\\hline')
        
    elif not line.strip() : # must be b/w datasets
        print('\\hline')

    else : # o.w., a line
        print(' & '.join(s), '\\\\')
        
#
# PARSER
#----------------------------------------------------------------------

distinguished = False
blue = False
slack = False
latex = False
entree = sys.stdin

a = sys.argv[1:]
i = 0
while i < len(a) :

    if a[i].startswith('-') :
        if a[i] == '-d' :
            distinguished = True
        elif a[i] == '-b' :
            blue = True
        elif a[i] == '-s' :
            if blue :
                assert False, usage
            slack = True
        elif a[i] == '-x' :
            if blue or slack :
                assert False, usage
            latex = True
        else :
            assert False, usage
    else :
        entree = open(a[i],'r')

    i += 1

#
# MAIN
#----------------------------------------------------------------------

count = 0
for line in entree :

    if not distinguished :
        line = line.replace('^',' ').replace('$',' ')

    if blue :
        print(line.replace('^','  \033[94m').replace('$','\033[0m'), end = '')
    elif slack :
        print(line.replace('^','*').replace('$','*'), end = '')
    elif latex :
        latex_format(line, count)
    else :
        print(line.replace('^',' ').replace('$',' '), end = '')

    count += 1

if latex :
    print('\\end{tabular}')
