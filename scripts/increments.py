usage = '''

usage: python increments.py [file.log]

given a log file file.log, the output of hapchat onto stdout, obtain
the sequence of increments at the sites where the k was incremented

'''

import sys

#
# Parser
#----------------------------------------------------------------------

entree = sys.stdin
a = sys.argv[1:]
i = 0
while i < len(a) : # bash-style argparse

    entree = open(a[i],'r')

    i += 1 # shift once by default

#
# Main
#----------------------------------------------------------------------

steps = {} # dictionary to store sequences of increments

for line in entree :

    if 'INCREMENT' in line : # .. STEP n INCREMENT k from a to b

        s = line.split()
        n, a, b = (int(x) for x in (s[-7], s[-3], s[-1])) # pull these values from line

        if n not in steps :
            steps[n] = set([])

        steps[n].add(a)
        steps[n].add(b)

# print dictionary at the end
print(*'step sequence'.split(), sep='\t')
for step in sorted(steps) :
    print(step, end='\t')
    print(*sorted(steps[step]))
