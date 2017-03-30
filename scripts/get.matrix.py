usage = '''

usage: python get.matrix.py [-s] [-e] [-w] variants.var [file.sfi]

given an sfi file file.sfi (default: stdin) and a set of variants
variants.var output (to stdout) the resulting SNV/fragment matrix
(optional [-w] to output matrix in the wif format of whatshap),
optional [-s] to restrict it down to SNPs, i.e., outputs SNP/fragment
matrix.  Optional [-e] to build 'extended' matrix, with +'s for those
positions that are not gaps but are multiallelic (and would be set to
gaps in the [-s] context).  Note that [-e] is overridden by [-s] and
has undefined behaviour in the wif [-w] context

'''

import sys

#
# load the variants as dic[pos] = [ref,alt]
def load_variants(lines) :

    variants = {}
    for line in lines :
        s = line.split()
        pos = int(s[0])

        assert pos not in variants
        variants[pos] = s[1:]

    return variants

#
# build index into position for O(1)-time access
def index_positions(positions) :

    index = {}
    for a,b in enumerate(positions) :
        index[b] = a

    return index

#
# from a line of an sfi file, return the index (line number in sfi
# file) and a dictionary dic[pos] = base:qual of the SNVs for the
# corresponding read ... restricted down to those SNVs that are also
# SNPs wrt variants (REF or ALT) if snps is True
def get_read(line, variants, snps) :

    name = line.split()[0]

    left, right = line.find(' $'), line.find('#insertions')
    assert right > left > 0
    segment = line[left:right]

    read = {}
    for token in segment.split()[1:] :
        a,b,c = token.partition(',')
        pos = int(a)
        allele = c.split(':')[0]

        # if not REF or ALT
        if allele not in variants[pos] :
            if snps : continue # it is not a SNP
            if extended : c = '+' # nor is it a gap

        read[pos] = c

    return name, read

#
# given a read: its set of SNVs, positions, and an index into the
# positions, return the run of SNVs consecutive wrt to the positions
def get_run(read, positions, pos_index) :

    first, last = min(read), max(read)
    i,j = pos_index[first], pos_index[last]
    run = positions[i:j+1]
    assert set(read) <= set(run)    

    return run

#
# given a read and a run on this read, return the row of the
# SNV/fragment matrix wrt the positions
def get_row(read, run) :

    if not read : return []

    row = []
    for pos in run :
        if pos in read :
            row.append(read[pos])
        else :
            row.append('-')

    return row

#
# given a read and a run on this read, return the corresponding line
# of a (whatshap) wif file
def get_wif(read, run, variants) :

    if not read : return ''

    wif = []
    for pos in run :
        ref,alt = variants[pos]
        nuc, allele, phred = ref, '0', '0'

        if pos in read :
            nuc,qual = read[pos].split(':')
            phred = str(ord(qual))
            if nuc == alt : allele = '1'

        wif.append(' '.join([str(pos), nuc, allele, phred, ':']))

    wif.append('# X : X') # dummy endpiece (for now)

    return ' '.join(wif)

#
# PARSER
#

snps = False
extended = False
wif = False
variants = None
entree = sys.stdin

a = sys.argv[1:]
assert len(a), usage
i = 0
while i < len(a) : # bash-style argparse

    if a[i].startswith('-') :
        if a[i] == '-s' :
            snps = True
        elif a[i] == '-e' :
            extended = True
        elif a[i] == '-w' :
            wif = True
        else :
            assert False, usage
    elif not variants :
        variants = load_variants(open(a[i],'r'))
    else :
        entree = open(a[i],'r')

    i += 1 # shift once by default in any case

#
# MAIN
#

positions = sorted(variants)
pos_index = index_positions(positions)

for line in entree :

    name, read = get_read(line, variants, snps)
    if not read : continue # don't consider empty reads (for whatever reason)

    run = get_run(read, positions, pos_index)
    assert len(run) # should be, since read is not empty

    if wif :
        print(get_wif(read, run, variants))

    else :
        pos = pos_index[run[0]]
        row = get_row(read, run)
        print(name + '\t' + str(pos+1) + '\t' + ' '.join(row))
