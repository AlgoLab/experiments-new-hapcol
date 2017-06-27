import sys

# setup
#----------------------------------------------------------------------
columns = 'Tool Data Technology Individual Chr MeanCov RawRealigned WhDowns WhDownsTimeSec Merge MergeE MergeM MergeT MergeN MergeTimeSec MergeMaxMemMB RndDowns RndDownsSeed RndDownsMaxCov RndDownsTimeSec RndDownsMaxMemMB FurtherMerging Epsilon Alpha BalThr BalRatio SwErrRatePerc HamDistPerc MecScore PhasTimeSec PhasMaxMemMB'.split()
ncols = len(columns)
print('number of columns:', ncols, file = sys.stderr)

tools = 'HapChat WhatsHap HapCol'.split()
data = 'ashk sim'.split()
chrs = 'chr1 chr21'.split()
meancovs = '15 20 25 30 40 50 60'.split()
modes = 'raw realigned'.split()
whdowns = 'N 15 20 25 30'.split()
no_merging = 'no NA NA NA NA'
mergings = [no_merging,
            'yes 0.15 0.25 1e+6 1e+3',
            'yes 0.15 0.25 1e+17 1e+3']
no_downs = 'no NA NA'
rnddowns = [no_downs] + ['yes 1 {}'.format(maxcov) for maxcov in whdowns[1:]]
alphas = 'NA 0.1 0.01 0.001 0.0001 0.00001'.split()
no_beta = 'NA NA'
null_beta = 'N 0'
betas = [no_beta, null_beta]

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
hlcount = 0
for datum in data :
    for chr in chrs :
        for cov in meancovs :
            for mode in modes :
                dataset = '{}.{}.cov{}.{}'.format(datum,chr,cov,mode)
                t_hc = table['HapChat'][datum][chr][cov][mode]

                for whdown in whdowns[1:] :
                    for alpha in alphas[1:] :
                        if whdown in whdowns[-2:] : # high covs for alpha=0.1, realigned only
                            if mode != 'realigned' :
                                continue
                            if alpha != '0.1' :
                                continue

                        for beta in betas[1:] :
                            assert t_hc[whdown][no_merging][no_downs][alpha][beta], 'no record for HapChat, dataset: {}, WhDowns: {}, Alpha: {}'.format(dataset, whdown, alpha)
                            hccount += 1

                    if whdown in whdowns[-1] : # no highest cov for HapCol
                        continue

                    t_hl = table['HapCol'][datum][chr][cov]['realigned'][whdown]
                    assert t_hl[no_merging][no_downs]['NA'][no_beta], 'no record for HapCol, dataset: {}, WhDowns: {}'.format(dataset, whdown)
                    if mode == 'realigned' :
                        hlcount += 1

                    if whdown in whdowns[-2:] : # no high covs for whatsHap
                        continue

                    t_wh = table['WhatsHap'][datum][chr][cov][mode][whdown]
                    assert t_wh[no_merging][no_downs]['NA'][no_beta], 'no record for WhatsHap, dataset: {}, WhDowns: {}'.format(dataset, whdown)
                    whcount += 1

                for merging in mergings :
                    for rnddown in rnddowns[1:] :
                        t_m = t_hc['N'][merging][rnddown]

                        for alpha in alphas[1:] :
                            if rnddown in rnddowns[-2:] : # high covs for alpha=0.1, realigned only
                                if mode != 'realigned' :
                                    continue
                                if alpha != '0.1' :
                                    continue

                            for beta in betas[1:] :
                                assert t_m[alpha][beta], 'no record for HapChat, dataset: {}, merge: {}, random downsample: {}, Alpha: {}'.format(dataset, merging, rnddown, alpha)
                                hccount += 1

print('number of WhatsHap records:', whcount, file = sys.stderr)
print('number of HapCol records:', hlcount, file = sys.stderr)
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
    msg(' realignment mode = {}'.format(mode))
    msg(' final coverage = {}'.format(maxcov))

def emph(string) :
    return '^{}$'.format(string)

def emph_winners(values) :

    comparables = []
    for value in values :
        if value not in nonvalues :
            comparables.append(value)

    winning_value = None
    if comparables :
        winning_value = min([float(value) for value in comparables])

    row = []
    for value in values :
        if value in comparables and float(value) == winning_value :
            row.append(emph(value))
        else :
            row.append(value)

    return ' '.join(row)

# some preamble for the display of the data
#----------------------------------------------------------------------
pipelines = 'whdowns merge-t6 merge-t17 rnddowns'.split()

pipeline_name = {'whdowns' : 'whatshap downsampling',
                 'merge-t6' : 'merging threshold 6',
                 'merge-t17' : 'merging threshold 17',
                 'rnddowns' : 'random downsampling'}

measures = 'swerr hamming time mem'.split()

measure_name = {'swerr' : 'switch error percentage',
                'hamming' : 'hamming distance percentage',
                'time' : 'time in seconds',
                'mem' : 'memory used in Mb'}

steps = 'prep phasing total'.split()

def step_name(tool, measure, step) :
    if tool == 'WhatsHap' :
        return ''
    if measure in ['swerr', 'hamming'] :
        return ''
    if step == 'prep' :
        return 'preprocessing '
    return step+' '

nonvalues = '? _'.split()

variants = 'mode alpha maxcov'.split()

variant_name = {'mode' : 'realignment mode',
                'alpha' : 'alpha',
                'maxcov' : 'final coverage'}

short = {'WhatsHap' : 'WH',
         'HapCol' : 'HC',
         'HapChat' : 'HX'}

def variant_vals(tool) :

    if tool == 'WhatsHap' :
        return {'mode' : modes,
                'alpha' : alphas[:1],
                'maxcov' : whdowns[1:-2]}

    elif tool == 'HapCol' :
        return {'mode' : ['realigned'],
                'alpha' : alphas[:1],
                'maxcov' : whdowns[1:-1]}

    elif tool == 'HapChat' :
        return {'mode' : modes,
                'alpha' : alphas[1:],
                'maxcov' : whdowns[1:]}

def apply_variant(variant, mode, maxcov, alpha, value) :

    if variant == 'mode' :
        return value, maxcov, alpha
    elif variant == 'maxcov' :
        return mode, value, alpha
    elif variant == 'alpha' :
        return mode, maxcov, value
    else :
        assert False, 'unknown variant: '+variant

def display_invariant(tool, variant, mode, maxcov, alpha) :

    mode_ = ' realignment mode = {}'.format(mode)
    maxcov_ = ' final coverage = {}'.format(maxcov)
    alpha_ = ' alpha = {}'.format(alpha)

    if tool == 'HapChat' :
        if variant == 'mode' :
            msg(maxcov_)
            msg(alpha_)
        elif variant == 'maxcov' :
            msg(mode_)
            msg(alpha_)
        elif variant == 'alpha' :
            msg(mode_)
            msg(maxcov_)
        else :
            assert False, 'unknown variant: '+variant
    else :
        assert tool in ['WhatsHap', 'HapCol'], 'unknown tool: '+tool

        if variant == 'mode' :
            assert tool == 'WhatsHap', 'mode an illegal variant for '+tool
            msg(maxcov_)
        elif variant == 'maxcov' :
            msg(mode_)
        elif variant == 'alpha' :
            assert False, 'alpha an illegal variant for '+tool
        else :
            assert False, 'unknown variant: '+variant

# some shortcuts in the table, etc.
#----------------------------------------------------------------------
def pipeline_record(pipeline, t_data, maxcov, alpha) :

    if pipeline == 'whdowns' :
        return t_data[str(maxcov)][no_merging][no_downs][alpha][null_beta]
    elif pipeline == 'merge-t6' :
        return t_data['N'][mergings[1]]['yes 1 {}'.format(maxcov)][alpha][null_beta]
    elif pipeline == 'merge-t17' :
        return t_data['N'][mergings[2]]['yes 1 {}'.format(maxcov)][alpha][null_beta]
    elif pipeline == 'rnddowns' :
        return t_data['N'][no_merging]['yes 1 {}'.format(maxcov)][alpha][null_beta]
    else :
        assert False, 'unknown pipeline: '+pipeline

def pipeline_time(pipeline, step, t_data, maxcov, alpha) :

    record = pipeline_record(pipeline, t_data, maxcov, alpha)
    phasetime = record['PhasTimeSec']
    times = []

    if pipeline == 'whdowns' :
        times = [record['WhDownsTimeSec']]
    elif pipeline in 'merge-t6 merge-t17'.split() :
        times = [record['MergeTimeSec'], record['RndDownsTimeSec']]
    elif pipeline == 'rnddowns' :
        times = [record['RndDownsTimeSec']]
    else :
        assert False, 'unkown pipeline: '+pipeline

    if step == 'phasing' :
        times = [phasetime]
    elif step == 'total' :
        times += [phasetime]
    else :
        assert step == 'prep', 'unknown step: '+step

    return '{:.3f}'.format(sum([float(time) for time in times]))

def pipeline_mem(pipeline, step, t_data, maxcov, alpha) :

    record = pipeline_record(pipeline, t_data, maxcov, alpha)
    phasemem = record['PhasMaxMemMB']
    mems = []

    if pipeline == 'whdowns' :
        mems = ['?'] # unsupported at the moment
    elif pipeline in 'merge-t6 merge-t17'.split() :
        mems = [record['MergeMaxMemMB'], record['RndDownsMaxMemMB']]
    elif pipeline == 'rnddowns' :
        mems = [record['RndDownsMaxMemMB']]
    else :
        assert False, 'unkown pipeline: '+pipeline

    if step == 'phasing' :
        mems = [phasemem]
    elif step == 'total' :
        mems += [phasemem]
    else :
        assert step == 'prep', 'unknown step: '+step

    if '?' in mems :
        return '?'
    return '{:.3f}'.format(max([float(mem) for mem in mems]))

def pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha) :

    if measure == 'swerr' :
        return pipeline_record(pipeline, t_data, maxcov, alpha)['SwErrRatePerc']
    elif measure == 'hamming' :
        return pipeline_record(pipeline, t_data, maxcov, alpha)['HamDistPerc']
    elif measure == 'time' :
        return pipeline_time(pipeline, step, t_data, maxcov, alpha)
    elif measure == 'mem' :
        return pipeline_mem(pipeline, step, t_data, maxcov, alpha)
    else :
        assert False, 'unknown measure: '+measure

def whatshap_record(datum, chr, cov, mode, maxcov) :
    return table['WhatsHap'][datum][chr][cov][mode][str(maxcov)][no_merging][no_downs]['NA'][no_beta]

def whatshap_time(datum, chr, cov, mode, maxcov) :
    record = whatshap_record(datum, chr, cov, mode, maxcov)
    times = [record['WhDownsTimeSec'], record['PhasTimeSec']]

    return '{:.3f}'.format(sum([float(time) for time in times]))

def whatshap_measure(measure, datum, chr, cov, mode, maxcov) :

    if measure == 'swerr' :
        return whatshap_record(datum, chr, cov, mode, maxcov)['SwErrRatePerc']
    elif measure == 'hamming' :
        return whatshap_record(datum, chr, cov, mode, maxcov)['HamDistPerc']
    elif measure == 'time' :
        return whatshap_time(datum, chr, cov, mode, maxcov)
    elif measure == 'mem' :
        return '?' # unsupported at the moment
    else :
        assert False, 'unknown measure: '+measure

def hapcol_record(datum, chr, cov, maxcov) :
    return table['HapCol'][datum][chr][cov]['realigned'][str(maxcov)][no_merging][no_downs]['NA'][no_beta]

def hapcol_time(step, datum, chr, cov, maxcov) :
    record = hapcol_record(datum, chr, cov, maxcov)
    times = [record['WhDownsTimeSec']]
    phasetime = record['PhasTimeSec']

    if step == 'phasing' :
        times = [phasetime]
    elif step == 'total' :
        times += [phasetime]
    else :
        assert step == 'prep', 'unknown step: '+step

    return '{:.3f}'.format(sum([float(time) for time in times]))

def hapcol_mem(step, datum, chr, cov, maxcov) :
    record = hapcol_record(datum, chr, cov, maxcov)

    if step == 'phasing' :
        return '{:.3f}'.format(float(record['PhasMaxMemMB']))
    elif step in ['prep', 'total'] :
        return '?'
    else :
        assert False, 'unknown step: '+step

def hapcol_measure(measure, step, datum, chr, cov, maxcov) :

    if measure == 'swerr' :
        return hapcol_record(datum, chr, cov, maxcov)['SwErrRatePerc']
    elif measure == 'hamming' :
        return hapcol_record(datum, chr, cov, maxcov)['HamDistPerc']
    elif measure == 'time' :
        return hapcol_time(step, datum, chr, cov, maxcov)
    elif measure == 'mem' :
        return hapcol_mem(step, datum, chr, cov, maxcov)

def tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha) :

    if tool == 'HapChat' :
        t_data = table['HapChat'][datum][chr][cov][mode]
        return pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha)
    elif tool == 'WhatsHap' :
        return whatshap_measure(measure, datum, chr, cov, mode, maxcov)
    else :
        assert tool == 'HapCol', 'unknown tool: '+tool
        return hapcol_measure(measure, step, datum, chr, cov, maxcov)

# compare the different pipelines for hapchat
#----------------------------------------------------------------------
def compare_pipelines(measure, step, mode, maxcov, alpha) :

    head()
    msg(' HapChat -- {}{} for each pipeline'.format(step_name('HapChat',measure,step), measure_name[measure]))
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
                    row.append(pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha))

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how does a tool (+ pipeline) vary with some parameter
#----------------------------------------------------------------------
def vary_param(tool, variant, measure, pipeline, step, mode, maxcov, alpha) :

    head()
    msg(' {} -- {}{} as a function of {}'.format(tool, step_name(tool,measure,step), measure_name[measure], variant_name[variant]))
    display_invariant(tool, variant, mode, maxcov, alpha)
    if tool == 'HapChat' :
        msg(' pipeline = {}'.format(pipeline_name[pipeline]))
    tail()

    print('#dataset\{}'.format(variant), ' '.join(variant_vals(tool)[variant]))
    print()
    for datum in data :
        for chr in chrs :
            for cov in meancovs :

                row = []
                for value in variant_vals(tool)[variant] :
                    mode, maxcov, alpha = apply_variant(variant, mode, maxcov, alpha, value)
                    col = tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha)
                    row.append(col)

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how do tools compare to each other over a set of coverages
#----------------------------------------------------------------------
def compare_tools(tools, maxcovs, measure, pipeline, mode, alpha) :

    head()
    tools_str = ' vs. '.join(tools)
    msg(' {} in terms of {}'.format(tools_str, measure_name[measure]))
    msg(' realignment mode = {}'.format(mode))
    if 'HapChat' in tools :
        msg(' HapChat alpha = {}'.format(alpha))
        msg(' HapChat pipeline = {}'.format(pipeline_name[pipeline]))
    tail()

    table_header = ' '.join(['{},cov={}'.format(short[tool], maxcov)
                             for tool in tools
                             for maxcov in maxcovs[tool]])

    print('#dataset', table_header)
    print()
    for datum in data :
        for chr in chrs :
            for cov in meancovs :

                row = []
                for tool in tools :
                    for maxcov in maxcovs[tool] :
                        col = tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha)
                        row.append(col)

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# add your own stuff here ...
#----------------------------------------------------------------------

tool = 'HapChat' # WhatsHap, HapCol, HapChat
variant = 'maxcov' # mode, alpha, maxcov
measure = 'time' # swerr, hamming, time, mem
pipeline = 'rnddowns' # whdowns, merge-t6, merge-t17, rnddowns
step = 'total' # prep, phasing, total
mode = 'realigned' # raw, realigned
maxcov = 20 # 15, 20, 25, 30
alpha = '0.1' # 0.1, 0.01, 0.001, 0.0001, 0.00001

tools = ['HapChat', 'HapCol', 'WhatsHap']
maxcovs = {'WhatsHap' : [15], # 15, 20
           'HapCol' : [15], # 15, 20, 25
           'HapChat' : [20]} # 15, 20, 25, 30

# tables
#vary_param(tool, variant, measure, pipeline, step, mode, maxcov, alpha)
#compare_pipelines(measure, step, mode, maxcov, alpha)
#compare_tools(tools, maxcovs, measure, pipeline, mode, alpha)
