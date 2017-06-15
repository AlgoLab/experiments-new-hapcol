import sys

# setup
#----------------------------------------------------------------------
columns = 'Tool Data Technology Individual Chr MeanCov RawRealigned WhDowns WhDownsTimeSec Merge MergeE MergeM MergeT MergeN MergeTimeSec MergeMaxMemMB RndDowns RndDownsSeed RndDownsMaxCov RndDownsTimeSec RndDownsMaxMemMB FurtherMerging Epsilon Alpha BalThr BalRatio SwErrRatePerc HamDistPerc MecScore PhasTimeSec PhasMaxMemMB'.split()
ncols = len(columns)
print('number of columns:', ncols, file = sys.stderr)

tools = 'HapChat WhatsHap'.split()
data = 'ashk sim'.split()
chrs = 'chr1 chr21'.split()
meancovs = '15 20 25 30 40 50 60'.split()
modes = 'raw realigned'.split()
whdowns = 'N 15 20'.split()
no_merging = 'no NA NA NA NA'
mergings = [no_merging,
          'yes 0.15 0.25 1e+6 1e+3',
          'yes 0.15 0.25 1e+17 1e+3']
no_downs = 'no NA NA'
rnddowns = [no_downs, 'yes 1 15', 'yes 1 20']
alphas = 'NA 0.1 0.01 0.001 0.0001 0.00001'.split()
no_beta = 'NA NA'
betas = [no_beta, 'N 0']

# initialize table
table = {}
for tool in tools :

    table[tool] = {}
    t_tool = table[tool]
    for datum in data :

        t_tool[datum] = {}
        t_datum = t_tool[datum]
        for chr in chrs :

            t_datum[chr] = {}
            t_chr = t_datum[chr]
            for cov in meancovs :

                t_chr[cov] = {}
                t_cov = t_chr[cov]
                for mode in modes :

                    t_cov[mode] = {}
                    t_mode = t_cov[mode]
                    for whdown in whdowns :

                        t_mode[whdown] = {}
                        t_whdown = t_mode[whdown]
                        for merging in mergings :

                            t_whdown[merging] = {}
                            t_merging = t_whdown[merging]
                            for rnddown in rnddowns :

                                t_merging[rnddown] = {}
                                t_rnddown = t_merging[rnddown]
                                for alpha in alphas :

                                    t_rnddown[alpha] = {}
                                    t_alpha = t_rnddown[alpha]
                                    for beta in betas :

                                        t_alpha[beta] = {}

# parse csv file
#----------------------------------------------------------------------
entree = sys.stdin
a = sys.argv[1:]
if a :
    entree = open(a[0],'r')

count = 0
for line in entree :

    s = line.strip().split(',')
    if not count :
        assert s == columns, 'discordant header: '+str(s)
        count += 1
        continue

    assert len(s) == len(columns), 'discordant number of columns: '+len(s)

    # index and fill row
    row = {}
    for i,column in enumerate(columns) :
        row[column] = s[i]

    # annotate parameters of row
    tool = row['Tool']
    datum = row['Data']
    chr = row['Chr']
    cov = row['MeanCov']
    mode = row['RawRealigned']
    whdown = row['WhDowns']
    merging = ' '.join([row['Merge'], row['MergeE'], row['MergeM'],
                        row['MergeT'], row['MergeN']])
    rnddown = ' '.join([row['RndDowns'], row['RndDownsSeed'],
                       row['RndDownsMaxCov']])
    alpha = row['Alpha']
    beta = ' '.join([row['BalThr'], row['BalRatio']])

    # test if row is correct
    assert tool in tools, 'unknown Tool: '+tool
    assert datum in data, 'unknown Data: '+datum
    assert chr in chrs, 'unknown Chr: '+chr
    assert cov in meancovs, 'unknown MeanCov: '+cov
    assert mode in modes, 'unknown RawRealigned'+mode
    assert whdown in whdowns, 'unknown WhDowns: '+whdown
    assert merging in mergings, 'unknown merging: '+merging
    assert rnddown in rnddowns, 'unknown random downsampling: '+rnddown
    assert alpha in alphas, 'unknown Alpha: '+alpha
    assert beta in betas, 'unknown beta'+beta

    # add row to the table
    table[tool][datum][chr][cov][mode][whdown][merging][rnddown][alpha][beta] = row

    count += 1

print('read {} entries'.format(count - 1), file = sys.stderr)

# test if table is complete
hccount = 0
whcount = 0
for datum in data :
    for chr in chrs :
        for cov in meancovs :
            for mode in modes :
                dataset = '{}.{}.cov{}.{}'.format(datum,chr,cov,mode)
                t_hc = table['HapChat'][datum][chr][cov][mode]

                for whdown in whdowns[1:] :
                    t_wh = table['WhatsHap'][datum][chr][cov][mode][whdown]

                    assert t_wh[no_merging][no_downs]['NA'][no_beta], 'no record for WhatsHap, dataset: {}, WhDowns: {}'.format(dataset, whdown)
                    whcount += 1

                    for alpha in alphas[1:] :
                        for beta in betas[1:] :
                            assert t_hc[whdown][no_merging][no_downs][alpha][beta], 'no record for HapChat, dataset: {}, WhDowns: {}, Alpha: {}'.format(dataset, whdown, alpha)
                            hccount += 1

                for merging in mergings :
                    for rnddown in rnddowns[1:] :
                        t_m = t_hc['N'][merging][rnddown]

                        for alpha in alphas[1:] :
                            for beta in betas[1:] :
                                assert t_m[alpha][beta], 'no record for HapChat, dataset: {}, merge: {}, random downsample: {}, Alpha: {}'.format(dataset, merging, rnddown, alpha)
                                hccount += 1

print('number of WhatsHap records:', whcount, file = sys.stderr)
print('number of HapChat records:', hccount, file = sys.stderr)

# aux functions
#----------------------------------------------------------------------
def msg(line) :
    print(line, file = sys.stderr)

def head() :
    msg('')

def tail() :
    msg(' '+70*'-')
    msg('')

def modemax(mode, maxcov) :
    msg(' '+mode)
    msg(' final cov = {}'.format(maxcov))

def emph(string) :
    return '^{}$'.format(string)

def emph_winners(values) :
    row = []
    winning_value = min([float(value) for value in values])

    for value in values :
        if float(value) == winning_value :
            row.append(emph(value))
        else :
            row.append(value)

    return ' '.join(row)
            
# some shortcuts in the table, etc.
#----------------------------------------------------------------------
pipelines = 'whdowns merge-t6 merge-t17 rnddowns'.split()
    
pipeline_name = {'whdowns' : 'whatshap downsampling',
                  'merge-t6' : 'merging threshold 6',
                  'merge-t17' : 'merging threshold 17',
                  'rnddowns' : 'random downsampling'}

def pipeline_record(pipeline, t_data, maxcov, alpha) :

    if pipeline == 'whdowns' :
        return t_data[str(maxcov)][no_merging][no_downs][alpha]['N 0']
    elif pipeline == 'merge-t6' :
        return t_data['N'][mergings[1]]['yes 1 {}'.format(maxcov)][alpha]['N 0']
    elif pipeline == 'merge-t17' :
        return t_data['N'][mergings[2]]['yes 1 {}'.format(maxcov)][alpha]['N 0']
    elif pipeline == 'rnddowns' :
        return t_data['N'][no_merging]['yes 1 {}'.format(maxcov)][alpha]['N 0']
    else :
        assert False, 'unknown pipeline: '+pipeline

def whatshap_record(datum, chr, cov, mode, maxcov) :
    return table['WhatsHap'][datum][chr][cov][mode][str(maxcov)][no_merging][no_downs]['NA'][no_beta]

# compare the different pipelines for hapchat
#----------------------------------------------------------------------
def compare_pipelines(mode, maxcov, alpha) :

    head()
    msg(' HapChat -- switch error for each pipeline')
    modemax(mode, maxcov)
    msg(' alpha = {}'.format(alpha))
    tail()

    print('#dataset', ' '.join(pipelines))
    print()
    for datum in data :
        for chr in chrs :
            for cov in meancovs :

                t_data = table['HapChat'][datum][chr][cov][mode]
                row = []
                for pipeline in pipelines :
                    record = pipeline_record(pipeline, t_data, maxcov, alpha)
                    row.append(record['SwErrRatePerc'])

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how does a hapchat pipeline vary with alpha
#----------------------------------------------------------------------
def vary_alpha(pipeline, mode, maxcov) :

    head()
    msg(' HapChat -- switch error as a function of alpha')
    modemax(mode, maxcov)
    msg(' pipeline = {}'.format(pipeline_name[pipeline]))
    tail()

    print('#dataset', ' '.join(alphas[1:]))
    print()
    for datum in data :
        for chr in chrs :
            for cov in meancovs :

                t_data = table['HapChat'][datum][chr][cov][mode]
                row = []
                for alpha in alphas[1:] :
                    record = pipeline_record(pipeline, t_data, maxcov, alpha)
                    row.append(record['SwErrRatePerc'])

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how does hapchat compare to whatshap
#----------------------------------------------------------------------
def hapchat_whatshap(pipeline, mode, maxcov, alpha) :

    head()
    msg(' HapChat vs. WhatsHap')
    modemax(mode, maxcov)
    msg(' HapChat pipeline = {}'.format(pipeline_name[pipeline]))
    msg(' HapChat alpha = {}'.format(alpha))
    tail()

    print('#dataset HapChat WhatsHap')
    print()
    for datum in data :
        for chr in chrs :
            for cov in meancovs :

                t_data = table['HapChat'][datum][chr][cov][mode]
                record = pipeline_record(pipeline, t_data, maxcov, alpha)
                wh_record = whatshap_record(datum,chr,cov,mode,maxcov)
                winline = emph_winners([record['SwErrRatePerc'],
                                        wh_record['SwErrRatePerc']])

                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# add your own stuff here ...
#----------------------------------------------------------------------

mode = 'realigned'
maxcov = 20
alpha = '0.1'

#compare_pipelines(mode, maxcov, alpha)
#vary_alpha(pipelines[0], mode, maxcov)
hapchat_whatshap('merge-t17', mode, maxcov, alpha)
