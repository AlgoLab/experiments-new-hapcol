usage = '''

usage: python display_table.py [-b] [-s] [file]

display a table annotated by study_table.py to stdout, slack, ..

positional arguments:
  file        input file (default: <stdin>)

optional arguments:
  -h, --help  show this help message and exit
  -b          output marked entries in blue for stdout
  -s          output marked entries in slack bold

note: only one display mode can be activated at a time, and if no
display mode is activated, it removes annotation from the table

'''

import sys

#
# PARSER
#----------------------------------------------------------------------

blue = False
slack = False
entree = sys.stdin

a = sys.argv[1:]
i = 0
while i < len(a) :

    if a[i].startswith('-') :
        if a[i] == '-b' :
            blue = True
            i += 1
        elif a[i] == '-s' :
            if blue :
                assert False, usage

            slack = True
            i += 1
        else :
            assert False, usage
    else :
        entree = open(a[i],'r')

    i += 1

#
# MAIN
#----------------------------------------------------------------------

for line in entree :
    
    if blue :
        print(line.replace('^','  \033[94m').replace('$','\033[0m'), end = '')
    elif slack :
        print(line.replace('^','*').replace('$','*'), end = '')
    else :
        print(line.replace('^',' ').replace('$',' '), end = '')
