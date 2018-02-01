#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
include : 'setup.snake'

time = '/usr/bin/time'
timeout = '/usr/bin/timeout'
compare = 'programs/whatshap/venv/bin/whatshap compare'
hapcut2vcf = 'programs/whatshap/venv/bin/whatshap hapcut2vcf'

# programs
_corewh_ = 'programs/core_whatshap/build/dp'
_hapchat_ = 'programs/balancing-hapcol/build/hapcol'
_hapcol_ = 'programs/HapCol/build/hapcol'
_hapcut2_ = 'programs/HapCUT2/build/HAPCUT2'
_probhap_ = 'programs/ProbHap/probhap.py'

# limits on memory usage and runtime
memlimit = 64 * 1024 * 1024 # 64GB limit (in KB)
timelimit = '24h' # 24 hour time limit

# methods
sih_methods = ['refhap', 'fasthare']
hapcut_methods = ['hapcut', 'hapcut2'] + sih_methods + ['probhap']
hap_methods = ['core_wh', 'hapcol', 'hapchat']
all_methods = hapcut_methods + hap_methods + ['whatshap']

def list_regex(a) :
	return '|'.join([x for x in a])

sih_pattern = '{method,(' + list_regex(sih_methods) + ')}'
hapcut_pattern = '{method,(' + list_regex(hapcut_methods) + ')}'
hap_pattern = '{method,(' + list_regex(hap_methods) + ')}'
methods_pattern = '{method,(' + list_regex(all_methods) + ')}'

def sih_method(wildcards) :

	sih = 'programs/refhap/SingleIndividualHaplotyper/SIH.jar'
	main = 'mpg.molgen.sih.main.SIH'
	opt = '-a FastHare' if wildcards.method == 'fasthare' else ''

        return 'java -cp {} {} {}'.format(sih, main, opt)

# epislon / alpha pairs for hapcol, and variants
ea_vals = ['05_1', '05_01', '05_001', '05_0001', '05_00001']
ea_two = ['01_1', '01_01', '01_001', '01_0001', '1_1', '1_01', '1_001', '1_0001']

# pattern (taking into account all input/output types)
full_pattern = post_pattern + '{ea,(|.[0-9]+_[0-9]+)}{balancing,(|.b([0-9]+|N)_[0-9]+)}{indelmode,(|.indels|.noindels)}'
output_pattern = methods_pattern + '/' + full_pattern

#
# useful list-defining functions (and lists)
#----------------------------------------------------------------------

# datasets processed by whatshap to a specified list of max cov.
def whatshap(datasets_, modes_, maxs_) :
	return ['{}.{}.h{}.no_merging.no_downs.no_merging'.format(dataset_, mode_, h_)
		for dataset_ in datasets_
		for mode_ in modes_
		for h_ in maxs_]

# partial paramerization of a merging (threshold and neg. threshold)
def merging(thrs_, negthrs_) :
	return ['merged_e{}_m{}_t{}_n{}'.format(err_, max_, thresh_, neg_)
		for err_ in error_rates
		for max_ in max_errs
		for thresh_ in thrs_
		for neg_ in negthrs_]

# downsampling to a specified list of max coverages
def downs(maxs_) :
	return ['downs_s{}_m{}'.format(seed_, max_)
		for seed_ in seeds
		for max_ in maxs_]

# datasets postprocessed to a specified list of max cov
def postproc(datasets_, modes_, thrs_, negthrs_, maxs_, rnddowns = False) :
	only_rnddowns = ['.no_merging'] if rnddowns else []
	return ['{}.{}.hN.{}.{}.no_merging'.format(dataset_, mode_, merging_, downsampling_)
		for dataset_ in datasets_
		for mode_ in modes_
		for merging_ in merging(thrs_, negthrs_) + only_rnddowns
		for downsampling_ in downs(maxs_)]

# datasets processed by a hairs method (just a shortcut, really)
def hairs(datasets_, modes_, indelmodes_) :
	return ['{}.{}.hN.no_merging.no_downs.no_merging.{}'.format(dataset_, mode_, indelmode_)
		for dataset_ in datasets_
		for mode_ in modes_
		for indelmode_ in indelmodes_]

# datasets both processed by whatshap and postprocessed to list of max cov.
def sliceof(datasets_, modes_, thrs_, negthrs_, maxs_) :
	return whatshap(datasets_, modes_, maxs_) + postproc(datasets_, modes_, thrs_, negthrs_, maxs_, True)

# define a subset of the datasets in terms of chromosomes and coverages
def datasubset(chr_covs_) :
	return ['{}.pacbio.child.chr{}.cov{}'.format(data_, chromosome_, coverage_)
		for data_ in data
		for chromosome_ in chr_covs_
		for coverage_ in chr_covs_[chromosome_]]

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('output/whatshap/{pattern}.sum',
			pattern = whatshap(datasets, ['realigned'],
				[15, 20, 25])),

		expand('output/hapcol/{pattern}.sum',
			pattern = whatshap(datasets, ['raw'],
				[15, 20, 25, 30])),

		expand('output/hapchat/{pattern}.05_01.bN_0.{ext}',
			pattern = postproc(datasets, ['realigned'], [6], [3],
				[15, 20, 25, 30, 35, 40]),
			ext = ['sum', 'inc']),

		expand('output/{method}/{pattern}.out',
			method = sih_methods + ['probhap'],
			pattern = hairs(datasets, ['raw'], indelmodes))

# coming up ..
rule next :
	input :
		expand('output/hapcut2/{pattern}.out',
			pattern = hairs(datasets, modes, indelmodes)),

		expand('output/core_wh/{pattern}.sum',
			pattern = whatshap(
				datasubset(
					{21:[5,10]}),
				['realigned'],
				[15, 20])),

		expand('output/hapchat/{pattern}.05_00001.{bal}.sum',
			pattern = sliceof(
				datasubset(
                                        {1:[5,10,15,20],21:[5,10,15,20]}),
				['realigned'], [6], [3],
				[25]),
			bal = ['bN_0', 'b20_45'])

#
# run whatshap on a bam / vcf pair
#----------------------------------------------------------------------
rule run_whatshap :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf',
		ref = 'reference/human_g1k_v37.fasta'

	params :
		realignment = lambda wildcards, input :
                        '--reference '+input.ref if wildcards.realignment == 'realigned' else '',
		h = lambda wildcards :
			'1000' if wildcards.h == 'N' else wildcards.h

	output : 'output/whatshap/' + post_pattern + '.phased.vcf'

	log :
		log = 'output/whatshap/' + post_pattern + '.log',
		time = 'output/whatshap/' + post_pattern + '.time'

	message : '''

   running whatshap on bam/vcf pair :

   {input.bam} / {input.vcf}

   selecting coverage down to {wildcards.h}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      whatshap phase -o {output} {params.realignment} -H {params.h} \
         {input.vcf} {input.bam} > {log.log} 2>&1 || true
   touch {output} '''

#
# run the core whatshap dp on a wif file
#----------------------------------------------------------------------
rule run_core_whatshap :
	input : 'wif/' + post_pattern + '.wif'
	output : 'output/core_wh/' + post_pattern + '.hap'

	log :
		log = 'output/core_wh/' + post_pattern + '.log',
		time = 'output/core_wh/' + post_pattern + '.time',
		wif = 'output/core_wh/' + post_pattern + '.wif'

	message : '''

   running core whatshap dp on :

   {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_corewh_} -h {output} -a {input} \
         > {log.wif} 2> {log.log} || true
   touch {output} '''

#
# run hapchat with increase-k and balancing on a wif file
#----------------------------------------------------------------------
rule run_hapchat :
	input : 'wif/' + post_pattern + '.wif'

	params :
		epsilon = lambda wildcards :
			'0' + wildcards.ea.split('_')[0],
		alpha = lambda wildcards :
			'0.' + wildcards.ea.split('_')[1],
		balance_cov = lambda wildcards :
			'1000' if wildcards.balancing.split('_')[0].split('b')[1] == 'N' else wildcards.balancing.split('_')[0].split('b')[1],
		balance_ratio = lambda wildcards :
			'0.' + wildcards.balancing.split('_')[1]

	output :
		hap = 'output/hapchat/' + full_pattern + '.hap',
		log = 'output/hapchat/' + full_pattern + '.log',

	log :
		time = 'output/hapchat/' + full_pattern + '.time'

	message : '''

   running hapchat on :

   {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_hapchat_} -i {input} -o {output.hap} -A \
         -e {params.epsilon} -a {params.alpha} \
         -b {params.balance_cov} -r {params.balance_ratio} \
            > {output.log} 2>&1 || true
   touch {output} '''

#
# run hapcol on a wif file (from within a script that adapts alpha)
#----------------------------------------------------------------------
rule run_hapcol :
	input :
		script = 'scripts/run.hapcol.bash',
		wif = 'wif/' + post_pattern + '.wif'

	output : 'output/hapcol/' + post_pattern + '.hap'

	log :
		log = 'output/hapcol/' + post_pattern + '.log',
		time = 'output/hapcol/' + post_pattern + '.time'

	message : '''

   running hapcol on :

   {input.wif}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      bash {input.script} {_hapcol_} {input.wif} {output} \
         > {log.log} 2>&1 || true
   touch {output} '''

#
# run hapcut2 on a hairs file
#----------------------------------------------------------------------
rule run_hapcut2 :
	input :
		hairs = 'hairs2/' + hairs_pattern + '.hairs',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output :
		'output/hapcut2/' + full_pattern + '.out'

	log :
		log = 'output/hapcut2/' + full_pattern + '.log',
		time = 'output/hapcut2/' + full_pattern + '.time'

	message : '''

   running hapcut2 on hairs file: {input.hairs}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_hapcut2_} --fragments {input.hairs} --VCF {input.vcf} \
         --output {output} > {log.log} 2>&1 || true
   touch {output} '''

#
# run probhap on a hairs file
#----------------------------------------------------------------------
rule run_probhap :
	input :
		'hairs/' + hairs_pattern + '.hairs'

	output :
		'output/probhap/' + full_pattern + '.out'

	log :
		reads = 'output/probhap/' + full_pattern + '.reads',
		assignments = 'output/probhap/' + full_pattern + '.assignments',
		log = 'output/probhap/' + full_pattern + '.log',
		time = 'output/probhap/' + full_pattern + '.time'

	message : '''

   running probhap on hairs file: {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      python2 {_probhap_} \
         --reads {input} \
         --parsed-reads {log.reads} \
         --phase {output} \
         --assignments {log.assignments} > {log.log} 2>&1 || true
   touch {output} '''

#
# run a SIH method (refhap, fasthare, ..)
#----------------------------------------------------------------------
rule run_sih_method :
	input :
		'hairs/' + hairs_pattern + '.hairs',

	params :
		method = sih_method

	output :
		'output/' + sih_pattern + '/' + full_pattern + '.out'

	log :
		log = 'output/{method}/' + full_pattern + '.log',
		time = 'output/{method}/' + full_pattern + '.time'

	message : '''

   running {wildcards.method} on hairs file: {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {params.method} {input} {output} > {log.log} 2>&1 || true
   touch {output} '''

#
# rule that asks for all the different results that we want for a run
#----------------------------------------------------------------------
rule gather_summary :
	input :
		diff = 'output/' + output_pattern + '.diff',
		mec = 'output/' + output_pattern + '.mec',
		sites = 'output/' + output_pattern + '.sites'

	output : 'output/' + output_pattern + '.sum'

	message : '''

   gather summary from:

   {input.diff}

   and:

   {input.mec} '''

	shell : '''

   grep -m 1 "switch error rate:" {input.diff} > {output} || true
   grep -m 1 "Block-wise Hamming distance " {input.diff} >> {output} || true
   cat {input.mec} | \
      awk '{{ print "                            "$1" "$2"  "$3 }}' \
         >> {output} '''

#
# compare phased vcfs to true phasing
#----------------------------------------------------------------------
rule compare_vcfs :
	input :
		true = 'vcf/' + vcf_pattern + '.phased.vcf',
		vcf = 'output/' + output_pattern + '.phased.vcf'

	output :
		diff = 'output/' + output_pattern + '.diff',
		bed = 'output/' + output_pattern + '.bed'

	log : 'output/' + output_pattern + '.diff.log'

	message : '''

   comparing inferred phasing:

   {input.vcf}

   to true phasing: {input.true} '''

	shell : '''

   {compare} --switch-error-bed {output.bed} \
      {input.true} {input.vcf} \
         > {output.diff} 2> {log} || true
   touch {output} '''

# convert old-skool hap format to a phased vcf
rule phase_vcf :
	input :
		script = 'scripts/subvcf.py',
		hap = 'output/' + hap_pattern + '/' + full_pattern + '.hap',
		blocks = 'wif/' + post_pattern + '.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'output/' + hap_pattern + '/' + full_pattern + '.phased.vcf'

	log : 'output/' + hap_pattern + '/' + full_pattern + '.phased.vcf.log'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} 2> {log} '''

# convert hapcut output format to (phased) vcf
rule hapcut_to_vcf :
	input :
		out = 'output/' + hapcut_pattern + '/' + full_pattern + '.out',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'output/' + hapcut_pattern + '/' + full_pattern + '.phased.vcf'

	log : 'output/' + hapcut_pattern + '/' + full_pattern + '.phased.vcf.log'

	message : 'converting hapcut {input} to {output}'

	shell : '''

   {hapcut2vcf} {input.vcf} {input.txt} > {output} 2> {log} '''

#
# compute MEC score of phased vcf wrt instance, as a wif file
#----------------------------------------------------------------------
rule mec_score :
	input :
		script = 'scripts/wiftools.py',
		vcf = 'output/' + output_pattern + '.phased.vcf',
		wif = 'wif/' + post_pattern + '.wif'

	output : 'output/' + output_pattern + '.mec'

	log : 'output/' + output_pattern + '.mec.log'

	message : '''

   infer mec score of:

   {input.vcf}

   with respect to: {input.wif} '''

	shell : '''

   python {input.script} -v {input.vcf} {input.wif} \
      > {output} 2> {log} '''

#
# get sitewise details of input/run, e.g., coverage, switch error, etc.
#----------------------------------------------------------------------
rule sitewise_details :
	input :
		script = 'scripts/sitesinfo.py',
		sites = 'wif/' + post_pattern + '.wif.info_/sites_',
		swerrs = 'output/' + output_pattern + '.bed'

	output : 'output/' + output_pattern + '.sites'

	log : 'output/' + output_pattern + '.sites.log'

	message : 'obtain sitewise details {output}'

	shell : '''

   python {input.script} -s {input.swerrs} {input.sites} \
      > {output} 2> {log} '''

#
# get details on increasing k from a hapchat log file
#----------------------------------------------------------------------
rule increments :
	input :
		script = 'scripts/increments.py',
		log = 'output/hapchat/' + full_pattern + '.log'

	output : 'output/hapchat/' + full_pattern + '.inc'

	log : 'output/hapchat/' + full_pattern + '.inc.log'

	message : 'obtain details on increasing k from {input.log}'

	shell : '''

   printf "%s Cov. %s: " {wildcards.dataset} {wildcards.coverage} > {output}
   python {input.script} -r {input.log} >> {output} 2> {log} '''
