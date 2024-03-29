import sys

# to change as more data arrives
hxmaxs = '15 20 25 30'.split()
hxas = '0.1 0.01 0.001'.split()
whmaxs = '15 20'.split()
hcmaxs = '15 20 25 30'.split()

# setup
#----------------------------------------------------------------------
columns = 'Tool Data Technology Individual Chr MeanCov RawRealigned WhDowns WhDownsTimeSec Merge MergeE MergeM MergeT MergeN MergeTimeSec MergeMaxMemMB RndDowns RndDownsSeed RndDownsMaxCov RndDownsTimeSec RndDownsMaxMemMB FurtherMerging Epsilon Alpha BalThr BalRatio IndelMode SwErrRatePerc HamDistPerc MecScore PhasTimeSec PhasMaxMemMB CleanFinish FeasibleSoln'.split()
ncols = len(columns)
print('number of columns:', ncols, file = sys.stderr)

tools = 'HapChat WhatsHap HapCol HapCUT2'.split()
data = 'ashk sim'.split()
chrs = 'chr1 chr21'.split()
meancovs = ['{}'.format(c) for c in range(15, 65, 5)]
chr_covs = { 'chr1' : meancovs, 'chr21' : meancovs[:-2] }
modes = 'raw realigned'.split()
maxcovs = '15 20 25 30'.split()
whdowns = maxcovs + ['N']
no_merging = 'no NA NA NA NA'
mergings = [no_merging,
            'yes 0.15 0.25 1e+6 1e+3',
            'yes 0.15 0.25 1e+17 1e+3']
no_downs = 'no NA NA'
rnddowns = [no_downs] + ['yes 1 {}'.format(m) for m in maxcovs]
alphas = 'NA 0.1 0.01 0.001 0.0001 0.00001'.split()
no_beta = 'NA NA'
null_beta = 'N 0'
betas = [no_beta, null_beta]
indelmodes = 'NA no yes'.split()

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
                                        t_beta = t_alpha[beta]
                                        for indelmode in indelmodes :

                                            t_beta[indelmode] = {}

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
    indelmode = row['IndelMode']

    # test if row is correct
    assert tool in tools, 'unknown Tool: '+tool
    assert datum in data, 'unknown Data: '+datum
    assert chr in chrs, 'unknown Chr: '+chr
    assert cov in chr_covs[chr], 'unknown MeanCov {} for chr {}'.format(cov, chr)
    assert mode in modes, 'unknown RawRealigned '+mode
    assert whdown in whdowns, 'unknown WhDowns: '+whdown
    assert merging in mergings, 'unknown merging: '+merging
    assert rnddown in rnddowns, 'unknown random downsampling: '+rnddown
    assert alpha in alphas, 'unknown Alpha: '+alpha
    assert beta in betas, 'unknown beta '+beta
    assert indelmode in indelmodes, 'unknown IndelMode '+indelmode

    # add row to the table
    table[tool][datum][chr][cov][mode][whdown][merging][rnddown][alpha][beta][indelmode] = row

    count += 1

print('read {} entries'.format(count - 1), file = sys.stderr)

# test if table is complete
hxcount = 0
whcount = 0
hccount = 0
hucount = 0
for datum in data :
    for chr in chrs :
        for cov in chr_covs[chr] :

            # hapcut2
            dataset = '{}.{}.cov{}'.format(datum,chr,cov)
            for indelmode in indelmodes[1:] :
                t_hu = table['HapCUT2'][datum][chr][cov]['raw']['N'][no_merging][no_downs]
                assert t_hu['NA'][no_beta][indelmode], 'no record for HapCUT2, dataset: {}, IndelMode: {}'.format(dataset, indelmode)
                hucount += 1

            # hapcol
            for whdown in hcmaxs :
                for mode in modes :
                    t_hc = table['HapCol'][datum][chr][cov][mode][whdown]
                    assert t_hc[no_merging][no_downs]['NA'][no_beta]['NA'], 'no record for HapCol, dataset: {}, RawRealigend: {}, WhDowns: {}'.format(dataset, mode, whdown)
                    hccount += 1

            # whatshap
            for whdown in whmaxs :
                t_wh = table['WhatsHap'][datum][chr][cov]['realigned'][whdown]
                assert t_wh[no_merging][no_downs]['NA'][no_beta]['NA'], 'no record for WhatsHap in realignment mode, dataset: {}, WhDowns: {}'.format(dataset, whdown)
                whcount += 1

            # hapchat
            t_hx = table['HapChat'][datum][chr][cov]['realigned']
            for whdown in hxmaxs :
                for alpha in hxas :
                    assert t_hx[whdown][no_merging][no_downs][alpha][null_beta]['NA'], 'no record for HapChat in realignment mode, dataset: {}, WhDowns: {}, Alpha: {}'.format(dataset, whdown, alpha)
                    hxcount += 1

            for merging in mergings[:-1] :
                for rnddown in ['yes 1 {}'.format(m) for m in hxmaxs] :
                    t_m = t_hx['N'][merging][rnddown]

                    for alpha in hxas :
                        assert t_m[alpha][null_beta]['NA'], 'no record for HapChat, dataset: {}, merge: {}, random downsample: {}, Alpha: {}'.format(dataset, merging, rnddown, alpha)
                        hxcount += 1

print('number of WhatsHap records:', whcount, file = sys.stderr)
print('number of HapCol records:', hccount, file = sys.stderr)
print('number of HapChat records:', hxcount, file = sys.stderr)
print('number of HapCUT2 records:', hucount, file = sys.stderr)

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
pipelines = 'whdowns merge-t6 rnddowns'.split()

pipeline_name = {'whdowns' : 'whatshap downsampling',
                 'merge-t6' : 'merging threshold 6',
                 'rnddowns' : 'random downsampling'}

measures = 'swerr hamming time mem'.split()

measure_name = {'swerr' : 'switch error percentage',
                'hamming' : 'hamming distance percentage',
                'time' : 'time in seconds',
                'mem' : 'memory used in Mb'}

steps = 'prep phasing total'.split()

def step_name(tool, measure, step) :
    if tool in ['WhatsHap', 'HapCUT2'] :
        return ''
    if measure in ['swerr', 'hamming'] :
        return ''
    if step == 'prep' :
        return 'preprocessing '
    return step+' '

nonvalues = '? -'.split()

variants = 'mode alpha maxcov indelmode'.split()

variant_name = {'mode' : 'realignment mode',
                'alpha' : 'alpha',
                'maxcov' : 'final coverage',
                'indelmode' : 'indel mode'}

short = {'WhatsHap' : 'WH',
         'HapCol' : 'HC',
         'HapChat' : 'HX',
         'HapCUT2' : 'HU'}

def variant_vals(tool) :

    if tool == 'WhatsHap' :
        return {'mode' : ['realigned'],
                'alpha' : [],
                'maxcov' : whmaxs,
                'indelmode' : []}

    elif tool == 'HapCol' :
        return {'mode' : modes,
                'alpha' : ['0.01'],
                'maxcov' : hcmaxs,
                'indelmode' : []}

    elif tool == 'HapChat' :
        return {'mode' : ['realigned'],
                'alpha' : hxas,
                'maxcov' : hxmaxs,
                'indelmode' : []}

    elif tool == 'HapCUT2' :
        return {'mode' : [],
                'alpha' : [],
                'maxcov' : [],
                'indelmode' : indelmodes[1:]}

def apply_variant(variant, mode, maxcov, alpha, indelmode, value) :

    if variant == 'mode' :
        return value, maxcov, alpha, indelmode
    elif variant == 'maxcov' :
        return mode, value, alpha, indelmode
    elif variant == 'alpha' :
        return mode, maxcov, value, indelmode
    elif variant == 'indelmode' :
        return mode, maxcov, alpha, value
    else :
        assert False, 'unknown variant: '+variant

def display_invariant(tool, variant, mode, maxcov, alpha, indelmode) :

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
        elif variant == 'indelmode' :
            assert False, 'indelmode an illegal variant for HapChat'
        else :
            assert False, 'unknown variant: '+variant

    elif tool == 'HapCUT2' :
        assert variant == 'indelmode', variant+' an illegal variant for HapCUT2'

    else :
        assert tool in ['WhatsHap', 'HapCol'], 'unknown tool: '+tool

        if variant == 'mode' :
            msg(maxcov_)
        elif variant == 'maxcov' :
            msg(mode_)
        elif variant in ['alpha', 'indelmode'] :
            assert False, variant+' an illegal variant for '+tool
        else :
            assert False, 'unknown variant: '+variant

# some shortcuts in the table, etc.
#----------------------------------------------------------------------
def form(entry) :
    if entry in nonvalues :
        return entry
    return '{:.{prec}f}'.format(float(entry), prec = precision)

def pipeline_record(pipeline, t_data, maxcov, alpha) :

    if pipeline == 'whdowns' :
        return t_data[str(maxcov)][no_merging][no_downs][alpha][null_beta]['NA']
    elif pipeline == 'merge-t6' :
        return t_data['N'][mergings[1]]['yes 1 {}'.format(maxcov)][alpha][null_beta]['NA']
    elif pipeline == 'rnddowns' :
        return t_data['N'][no_merging]['yes 1 {}'.format(maxcov)][alpha][null_beta]['NA']
    else :
        assert False, 'unknown pipeline: '+pipeline

def pipeline_time(pipeline, step, t_data, maxcov, alpha) :

    record = pipeline_record(pipeline, t_data, maxcov, alpha)
    phasetime = record['PhasTimeSec']
    times = []

    if pipeline == 'whdowns' :
        times = [record['WhDownsTimeSec']]
    elif pipeline in 'merge-t6'.split() :
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

    if '-' in times :
        return '-'
    return form(sum([float(time) for time in times]))

def pipeline_mem(pipeline, step, t_data, maxcov, alpha) :

    record = pipeline_record(pipeline, t_data, maxcov, alpha)
    phasemem = record['PhasMaxMemMB']
    mems = []

    if pipeline == 'whdowns' :
        mems = ['?'] # unsupported at the moment
    elif pipeline in 'merge-t6'.split() :
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

    if '-' in mems :
        return '-'
    if '?' in mems :
        return '?'
    return form(max([float(mem) for mem in mems]))

def pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha) :

    if measure == 'swerr' :
        return form(pipeline_record(pipeline, t_data, maxcov, alpha)['SwErrRatePerc'])
    elif measure == 'hamming' :
        return form(pipeline_record(pipeline, t_data, maxcov, alpha)['HamDistPerc'])
    elif measure == 'time' :
        return pipeline_time(pipeline, step, t_data, maxcov, alpha)
    elif measure == 'mem' :
        return pipeline_mem(pipeline, step, t_data, maxcov, alpha)
    else :
        assert False, 'unknown measure: '+measure

def whatshap_record(datum, chr, cov, mode, maxcov) :
    return table['WhatsHap'][datum][chr][cov][mode][str(maxcov)][no_merging][no_downs]['NA'][no_beta]['NA']

def whatshap_time(datum, chr, cov, mode, maxcov) :
    record = whatshap_record(datum, chr, cov, mode, maxcov)
    times = [record['WhDownsTimeSec'], record['PhasTimeSec']]

    if '-' in times :
        return '-'
    return form(sum([float(time) for time in times]))

def whatshap_measure(measure, datum, chr, cov, mode, maxcov) :

    if measure == 'swerr' :
        return form(whatshap_record(datum, chr, cov, mode, maxcov)['SwErrRatePerc'])
    elif measure == 'hamming' :
        return form(whatshap_record(datum, chr, cov, mode, maxcov)['HamDistPerc'])
    elif measure == 'time' :
        return whatshap_time(datum, chr, cov, mode, maxcov)
    elif measure == 'mem' :
        return '?' # unsupported at the moment
    else :
        assert False, 'unknown measure: '+measure

def hapcol_record(datum, chr, cov, mode, maxcov) :
    return table['HapCol'][datum][chr][cov][mode][str(maxcov)][no_merging][no_downs]['NA'][no_beta]['NA']

def hapcol_time(step, datum, chr, cov, mode, maxcov) :
    record = hapcol_record(datum, chr, cov, mode, maxcov)
    times = [record['WhDownsTimeSec']]
    phasetime = record['PhasTimeSec']

    if step == 'phasing' :
        times = [phasetime]
    elif step == 'total' :
        times += [phasetime]
    else :
        assert step == 'prep', 'unknown step: '+step

    if '-' in times :
        return '-'
    return form(sum([float(time) for time in times]))

def hapcol_mem(step, datum, chr, cov, mode, maxcov) :
    record = hapcol_record(datum, chr, cov, mode, maxcov)

    if step == 'phasing' :
        return form(float(record['PhasMaxMemMB']))
    elif step in ['prep', 'total'] :
        return '?'
    else :
        assert False, 'unknown step: '+step

def hapcol_measure(measure, step, datum, chr, cov, mode, maxcov) :

    if measure == 'swerr' :
        return form(hapcol_record(datum, chr, cov, mode, maxcov)['SwErrRatePerc'])
    elif measure == 'hamming' :
        return form(hapcol_record(datum, chr, cov, mode, maxcov)['HamDistPerc'])
    elif measure == 'time' :
        return hapcol_time(step, datum, chr, cov, mode, maxcov)
    elif measure == 'mem' :
        return hapcol_mem(step, datum, chr, cov, mode, maxcov)

def hapcut2_record(datum, chr, cov, indelmode) :
    return table['HapCUT2'][datum][chr][cov]['raw']['N'][no_merging][no_downs]['NA'][no_beta][indelmode]

def hapcut2_measure(measure, datum, chr, cov, indelmode) :
    record = hapcut2_record(datum, chr, cov, indelmode)

    if measure == 'swerr' :
        return form(record['SwErrRatePerc'])
    elif measure == 'hamming' :
        return form(record['HamDistPerc'])
    elif measure == 'time' :
        return form(float(record['PhasTimeSec']))
    elif measure == 'mem' :
        return form(record['PhasMaxMemMB'])
    else :
        assert False, 'unknown measure '+measure

def tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha, indelmode) :

    if tool == 'HapChat' :
        t_data = table['HapChat'][datum][chr][cov]['realigned']
        return pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha)
    elif tool == 'WhatsHap' :
        return whatshap_measure(measure, datum, chr, cov, 'realigned', maxcov)
    elif tool == 'HapCol' :
        return hapcol_measure(measure, step, datum, chr, cov, mode, maxcov)
    else :
        assert tool == 'HapCUT2', 'unknown tool: '+tool
        return hapcut2_measure(measure, datum, chr, cov, indelmode)

# clear some fields in the table, based on other fields
def clear_fields(clearinfeasible, clearunfinished = True) :

    non_parameters = 'WhDownsTimeSec MergeTimeSec MergeMaxMemMB RndDownsTimeSec RndDownsMaxMemMB SwErrRatePerc HamDistPerc MecScore PhasTimeSec PhasMaxMemMB'.split()

    for tool in table :
        t_tool = table[tool]
        for datum in t_tool :
            t_datum = t_tool[datum]
            for chr in t_datum :
                t_chr = t_datum[chr]
                for cov in t_chr :
                    t_cov = t_chr[cov]
                    for mode in t_cov :
                        t_mode = t_cov[mode]
                        for whdown in t_mode :
                            t_whdown = t_mode[whdown]
                            for merging in t_whdown :
                                t_merging = t_whdown[merging]
                                for rnddown in t_merging :
                                    t_rnddown = t_merging[rnddown]
                                    for alpha in t_rnddown :
                                        t_alpha = t_rnddown[alpha]
                                        for beta in t_alpha :
                                            t_beta = t_alpha[beta]
                                            for indelmode in indelmodes :
                                                t_indelmode = t_beta[indelmode]
                                                if(t_indelmode) :

                                                    if(clearunfinished) :
                                                        if t_indelmode['CleanFinish'] == 'no' :
                                                            for field in non_parameters :
                                                                t_indelmode[field] = '-'

                                                    if(clearinfeasible) :
                                                        if t_indelmode['FeasibleSoln'] == 'no' :
                                                            for field in non_parameters :
                                                                t_indelmode[field] = '-'

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
            for cov in chr_covs[chr] :

                t_data = table['HapChat'][datum][chr][cov][mode]
                row = []
                for pipeline in pipelines :
                    row.append(pipeline_measure(measure, pipeline, step, t_data, maxcov, alpha))

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how does a tool (+ pipeline) vary with some parameter
#----------------------------------------------------------------------
def vary_param(tool, variant, measure, pipeline, step, mode, maxcov, alpha, indelmode) :

    head()
    msg(' {} -- {}{} as a function of {}'.format(tool, step_name(tool,measure,step), measure_name[measure], variant_name[variant]))
    display_invariant(tool, variant, mode, maxcov, alpha, indelmode)
    if tool == 'HapChat' :
        msg(' pipeline = {}'.format(pipeline_name[pipeline]))
    tail()

    print('#dataset\{}'.format(variant), ' '.join(variant_vals(tool)[variant]))
    print()
    for datum in data :
        for chr in chrs :
            for cov in chr_covs[chr] :

                row = []
                for value in variant_vals(tool)[variant] :
                    mode, maxcov, alpha, indelmode = apply_variant(variant, mode, maxcov, alpha, indelmode, value)
                    col = tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha, indelmode)
                    row.append(col)

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# how do tools compare to each other over a set of coverages
#----------------------------------------------------------------------
def compare_tools(tools, maxcovs, measure, pipeline, mode, alpha, indelmode) :

    head()
    tools_str = ' vs. '.join(['{} ({})'.format(tool, short[tool]) for tool in tools])
    msg(' {} in terms of {}'.format(tools_str, measure_name[measure]))
    if 'HapCol' in tools :
        msg(' HapCol realignment mode = {}'.format(mode))
    if 'HapChat' in tools :
        msg(' HapChat alpha = {}'.format(alpha))
        msg(' HapChat pipeline = {}'.format(pipeline_name[pipeline]))
    if 'HapCUT2' in tools :
        msg(' HapCUT2 indelmode = {}'.format(indelmode))
    tail()

    table_header = ' '.join(['{}-{}'.format(short[tool], maxcov)
                             for tool in tools
                             for maxcov in maxcovs[tool]])

    print('#Data\Tool-MaxC', table_header)
    print()
    for datum in data :
        for chr in chrs :
            for cov in chr_covs[chr] :

                row = []
                for tool in tools :
                    for maxcov in maxcovs[tool] :
                        col = tool_measure(tool, measure, pipeline, step, datum, chr, cov, mode, maxcov, alpha, indelmode)
                        row.append(col)

                winline = emph_winners(row)
                print('{}.{}.cov{}'.format(datum,chr,cov), winline)
            print()

# add your own stuff here ...
#----------------------------------------------------------------------

precision = 3 # 2, 3, 4, ..
tool = 'HapChat' # WhatsHap, HapCol, HapChat, HapCUT2
variant = 'alpha' # mode, alpha, maxcov, indelmode
measure = 'swerr' # swerr, hamming, time, mem
pipeline = 'merge-t6' # whdowns, merge-t6, rnddowns
step = 'total' # prep, phasing, total
mode = 'raw' # raw, realigned
maxcov = 30 # 15, 20, 25, ..
alpha = '0.01' # 0.1, 0.01, 0.001, 0.0001, ..
indelmode = 'no' # yes, no

tools = ['HapChat', 'HapCol', 'WhatsHap', 'HapCUT2']
maxcovs = {'WhatsHap' : whmaxs,
           'HapCol' : [25,30],
           'HapChat' : [25,30],
           'HapCUT2' : ['N']}

clearinfeasible = False # True, False
clear_fields(clearinfeasible)

# redefine in order to display a subset
meancovs = ['{}'.format(c) for c in range(25, 65, 5)]
chr_covs = { 'chr1' : meancovs, 'chr21' : meancovs[:-2] }

# tables
#vary_param(tool, variant, measure, pipeline, step, mode, maxcov, alpha, indelmode)
#compare_pipelines(measure, step, mode, maxcov, alpha)
compare_tools(tools, maxcovs, measure, pipeline, mode, alpha, indelmode)
