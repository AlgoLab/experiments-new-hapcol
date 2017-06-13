import sys

# setup
#----------------------------------------------------------------------
columns = 'Tool Data Technology Individual Chr MeanCov RawRealigned WhDowns Merge MergeE MergeM MergeT MergeN MergeTimeSec MergeMaxMemMB RndDowns RndDownsSeed RndDownsMaxCov RndDownsTimeSec RndDownsMaxMemMB FurtherMerging Epsilon Alpha BalThr BalRatio SwErrRatePerc HamDistPerc MecScore PhasTimeSec PhasMaxMemMB'.split()

tools = 'HapChat WhatsHap'.split()
data = 'ashk sim'.split()
chrs = 'chr1 chr21'.split()
meancovs = '15 20 25 30 40 50 60'.split()
modes = 'raw realigned'.split()
whdowns = 'N 15 20'.split()
merges = ['no NA NA NA NA',
          'yes 0.15 0.25 1e+6 1e+3',
          'yes 0.15 0.25 1e+17 1e+3']
rnddowns = ['no NA NA', 'yes 1 15', 'yes 1 20']
alphas = 'NA 0.1 0.01 0.001 0.0001 0.00001'.split()
betas = ['NA NA', 'N 0']

table = {}

# parse table
#----------------------------------------------------------------------
entree = sys.stdin
a = sys.argv[1:]
if a :
    entree = open(a[0],'r')

count = 0
for line in entree :

    table[count] = {}
    s = line.split(',')
    assert len(s) == len(columns), 'discordant number of columns'

    row = table[count]
    for i,column in enumerate(columns) :
        row[column] = s[i]

    count += 1

# test if table is correct and complete
hccount = 0
whcount = 0
for entry in table :
    if entry == 0 :
        continue

    row = table[entry]

    # correctness
    assert row['Tool'] in tools, 'unknown Tool: '+row['Tool']
    assert row['Data'] in data, 'unknown Data: '+row['Data']
    assert row['Chr'] in chrs, 'unknown Chr: '+row['Chr']
    assert row['MeanCov'] in meancovs, 'unknown MeanCov: '+row['MeanCov']
    assert row['RawRealigned'] in modes, 'unknown RawRealigned'+row['RawRealigned']
    assert row['WhDowns'] in whdowns, 'unknown WhDowns: '+row['WhDowns']
    assert ' '.join([row['Merge'], row['MergeE'], row['MergeM'],
                     row['MergeT'], row['MergeN']]) in merges
    assert ' '.join([row['RndDowns'], row['RndDownsSeed'],
                     row['RndDownsMaxCov']]) in rnddowns
    assert row['Alpha'] in alphas, 'unkown Alpha: '+row['Alpha']
    assert ' '.join([row['BalThr'], row['BalRatio']]) in betas

    # completeness
    if row['Tool'] == 'HapChat' :
        hccount += 1
    if row['Tool'] == 'WhatsHap' :
        whcount += 1

common = len(data) * len(chrs) * len(meancovs) * len(modes)
wh = common * (len(whdowns) - 1)
hc = common * (len(whdowns) - 1) * (len(merges) - 1) * (len(rnddowns) - 1) * (len(alphas) - 1)
assert whcount == wh, 'discordant number of WhatsHap records'
assert hccount == hc, 'discordant number of HapChat records'

# add your own stuff here ...
#----------------------------------------------------------------------
