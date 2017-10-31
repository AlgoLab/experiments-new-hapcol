#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
include : 'setup.snake'
time = '/usr/bin/time'
timeout = '/usr/bin/timeout'
compare = 'programs/whatshap/venv/bin/whatshap compare'
hapcut2vcf = 'programs/whatshap/venv/bin/whatshap hapcut2vcf'

# softwares
corewh = 'programs/core_whatshap/build/dp'
hapchat = 'programs/balancing-hapcol/build/hapcol'
hapcol = 'programs/HapCol/build/hapcol'
hapcut2 = 'programs/HapCUT2/build/HAPCUT2'

# limits on memory usage and runtime
memlimit = 64 * 1024 * 1024 # 64GB limit (in KB)
timelimit = '24h' # 24 hour time limit

# pattern (taking into account all input/output types)
full_pattern = post_pattern + '{ea,(|.[0-9]+_[0-9]+)}{balancing,(|.b([0-9]+|N)_[0-9]+)}{indelmode,(|.indels|.noindels)}'

# output directory pattern (for all methods that output .hap files)
outdir_pattern = '{dir,(output/hapchat|output/hapcol|output/core_wh)}'

# epislon / alpha pairs for hapcol, and variants
ea_vals = ['05_1', '05_01', '05_001', '05_0001', '05_00001']
ea_two = ['01_1', '01_01', '01_001', '01_0001', '1_1', '1_01', '1_001', '1_0001']

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
	return ['.merged_e{}_m{}_t{}_n{}'.format(err_, max_, thresh_, neg_)
		for err_ in error_rates
		for max_ in max_errs
		for thresh_ in thrs_
		for neg_ in negthrs_]

# downsampling to a specified list of max coverages
def downs(maxs_) :
	return ['.downs_s{}_m{}'.format(seed_, max_)
		for seed_ in seeds
		for max_ in maxs_]

# datasets postprocessed to a specified list of max cov
def postproc(datasets_, modes_, thrs_, negthrs_, maxs_, rnddowns = False) :
	only_rnddowns = ['.no_merging'] if rnddowns else []
	return ['{}.{}.hN{}{}.no_merging'.format(dataset_, mode_, merging_, downsampling_)
		for dataset_ in datasets_
		for mode_ in modes_
		for merging_ in merging(thrs_, negthrs_) + only_rnddowns
		for downsampling_ in downs(maxs_)]

# datasets processed by hapcut2 (just a shortcut, really)
def hapcut2(datasets_, modes_) :
	return ['{}.raw.hN.no_merging.no_downs.no_merging.{}'.format(dataset_, mode_)
		for dataset_ in datasets_
		for mode_ in modes_]

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
				[15, 20, 25, 30, 35]),
			ext = ['sum', 'inc']),

		expand('output/hapcut2/{pattern}.sum',
			pattern = hapcut2(datasets, indelmodes))

# coming up ..
rule next :
	input :
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

	output : 'output/whatshap/' + whatshap_pattern + '.no_merging.no_downs.no_merging.phased.vcf'

	log :
		log = 'output/whatshap/' + whatshap_pattern + '.no_merging.no_downs.no_merging.log',
		time = 'output/whatshap/' + whatshap_pattern + '.no_merging.no_downs.no_merging.time'

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
      {corewh} -h {output} -a {input} \
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
      {hapchat} -i {input} -o {output.hap} -A \
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
      bash {input.script} {hapcol} {input.wif} {output} \
         > {log.log} 2>&1 || true
   touch {output} '''

#
# run hapcut2 on a hairs file
#----------------------------------------------------------------------
rule run_hapcut2 :
	input :
		txt = 'hairs/' + hapcut_pattern + '.txt',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output :
		'output/hapcut2/' + dataset_pattern + '.raw.hN.no_merging.no_downs.no_merging.{indelsornot,(indels|noindesl)}.txt'

	log :
		log = 'output/hapcut2/' + hapcut_pattern + '.raw.hN.no_merging.no_downs.no_merging.{indelsornot}.log',
		time = 'output/hapcut2/' + hapcut_pattern + '.raw.hN.no_merging.no_downs.no_merging.{indelsornot}.time'

	message : '''

   running hapcut2 on hairs file: {input.txt}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {hapcut2} --fragments {input.txt} --VCF {input.vcf} \
         --output {output} &> {log.log} '''

#
# rule that asks for all the different results that we want for a run
#----------------------------------------------------------------------
rule gather_summary :
	input :
		diff = '{dir}/' + full_pattern + '.diff',
		mec = '{dir}/' + full_pattern + '.mec',
		sites = '{dir}/' + full_pattern + '.sites'

	output : '{dir}/' + full_pattern + '.sum'

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
		vcf = '{dir}/' + full_pattern + '.phased.vcf'

	output :
		diff = '{dir}/' + full_pattern + '.diff',
		bed = '{dir}/' + full_pattern + '.bed'

	log : '{dir}/' + full_pattern + '.diff.log'

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
		hap = outdir_pattern + '/' + full_pattern + '.hap',
		blocks = 'wif/' + post_pattern + '.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : outdir_pattern + '/' + full_pattern + '.phased.vcf'

	log : outdir_pattern + '/' + full_pattern + '.phased.vcf.log'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} 2> {log} '''

# convert hapcut2 output to (phased) vcf
rule hapcut2_to_vcf :
	input :
		txt = 'output/hapcut2/' + full_pattern + '.txt',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'output/hapcut2/' + full_pattern + '.phased.vcf'

	log : 'output/hapcut2/' + full_pattern + '.phased.vcf.log'

	message : 'converting hapcut {input} to {output}'

	shell : '''

   {hapcut2vcf} {input.vcf} {input.txt} > {output} 2> {log} '''

#
# compute MEC score of phased vcf wrt instance, as a wif file
#----------------------------------------------------------------------
rule mec_score :
	input :
		script = 'scripts/wiftools.py',
		vcf = '{dir}/' + full_pattern + '.phased.vcf',
		wif = 'wif/' + post_pattern + '.wif'

	output : '{dir}/' + full_pattern + '.mec'

	log : '{dir}/' + full_pattern + '.mec.log'

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
		swerrs = '{dir}/' + full_pattern + '.bed'

	output : '{dir}/' + full_pattern + '.sites'

	log : '{dir}/' + full_pattern + '.sites.log'

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

	shell : 'python {input.script} {input.log} > {output} 2> {log}'
