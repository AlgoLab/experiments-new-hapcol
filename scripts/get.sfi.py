usage = '''

usage: python get.sfi.py variants.var [reads.bam]

given a file variants.var containing a table of variants, output (to
stdout) for each line of file reads.bam (default: stdin) the
"snv/fragment info" (sfi) : snv positions -- along with corresponding
base call and quality -- that are covered by the (part of the) read
corresponding to this line, with respect to the variants, as well as
some other info, such as stats on indels

'''

import sys

#
# load set of variants (just the positions for now)
def load_variants(lines) :

    positions = set([])
    for line in lines :
        positions.add(int(line.split()[0]))

    return positions

#
# output a dictionary as a list of index:value pairs
def as_list(dic) :

    list = ''
    for i in sorted(dic) :
        list += ' ' + str(i) + ':' + str(dic[i])

    return list

#
# given a line of the bam file, return the corresponding line wrt the
# variants
def get_line(line, variants) :

    s = line.split()[:11]
    row = ' '.join([s[0],s[1],s[2],s[4],s[6],s[7],s[8]])
    row += ' $' # separate basic info from SNVs we will add

    pos = int(s[3])
    cigar = s[5]
    seq = s[9]
    qual = s[10]

    insertions = {} # for stats on indels
    deletions = {}

    ops = 'M I D N S H P = X'.split() # cigar string operations

    # should the cigar information not be available, it's as if there
    # are no variants, so we setup a bypass directly to return
    if cigar == '*' :
        cigar = ''
        seq = ''

    # eat up the cigar string, eating the sequence accordingly, while
    # advancing the position pointer along the reference.  The idea is
    # that we will always be progressing on at least one of the
    # sequence (i.e., ind) and the reference (i.e., pos).  Whenever we
    # hit a position that is in the set of SNVs, and the read is
    # active at that position (i.e., no clipping or deletion), we add
    # this base call (with its quality) to the row
    while cigar :

        # strip token off of cigar string
        p = min(filter(lambda x : x > 0, (cigar.find(o) for o in ops)))
        assert p > 0 # there should be such a p if cigar is not empty
        o = cigar[p]
        assert o in ops # sanity check
        a,b,c = cigar.partition(o) # leftmost such o, by def.
        cigar = c # strip this token off
        n = int(a) # should raise error if a is not purely digital

        if o in ['M', '=', 'X'] : # 'active' region
            for i in range(n) : # zip through range looking for SNVs
                if pos in variants : # add weighted SNV
                    row += ' ' + str(pos) + ',' + seq[0] + ':' + qual[0]

                # move forward in both seq and reference
                pos += 1
                seq = seq[1:]
                qual = qual[1:]

        elif o == 'I' : # insertion
            seq = seq[n:] # just eat n characters in seq
            qual = qual[n:]

            if n not in insertions : # for the stats
                insertions[n] = 0
            insertions[n] += 1

        elif o == 'D' : # deletion
            pos += n # just move n positions on the reference

            if n not in deletions : # for the stats
                deletions[n] = 0
            deletions[n] += 1

        elif o == 'N' : # skipped region
            pos += n # just move n positions on the reference

        elif o == 'S' : # soft clipping
            pos += n # just move everything forward
            seq = seq[n:]
            qual = qual[n:]

        else : # else it's hard clipping, in which we do nothing
            assert o in ['H', 'P']

    assert not seq # should be empty now, o.w., there's a problem

    row += ' #insertions' + as_list(insertions) # add indel stats
    row += ' #deletions' + as_list(deletions)

    return row

#
# PARSER
#

variants = None
entree = sys.stdin

a = sys.argv[1:]
assert len(a), usage
i = 0
while i < len(a) : # bash-style argparse

    if not variants :
        variants = load_variants(open(a[i],'r'))
    else :
        entree = open(a[i],'r')

    i += 1 # shift once by default in any case

#
# MAIN
#

for line in entree : # move through lines of bam file getting the lines
    print(get_line(line, variants))
