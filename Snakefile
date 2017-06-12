#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
include : 'setup.snake'
time = '/usr/bin/time'
compare = 'programs/whatshap/venv/bin/whatshap compare'

# softwares
corewh = 'programs/core_whatshap/build/dp'
hapchat = 'programs/balancing-hapcol/build/hapcol'

# pattern (taking into account hapcol)
full_pattern = post_pattern + '{ea,(|.[0-9]+_[0-9]+)}{balancing,(|.b([0-9]+|N)_[0-9]+)}'

# epislon / alpha pairs for hapcol, and variants
ea_vals = ['05_1', '05_01', '05_001', '05_0001', '05_00001']
ea_two = ['01_1', '01_01', '01_001', '01_0001', '1_1', '1_01', '1_001', '1_0001']

# extensions desired
exts = ['diff', 'mec', 'sites']

#
# useful list-defining functions
#----------------------------------------------------------------------

# downsampling to a specified list of max coverages
def downs(maxs_) :
	return ['.downs_s{}_m{}'.format(seed_, max_)
		for seed_ in seeds
		for max_ in maxs_]

# datasets processed by whatshap to a specified list of max cov.
def whatshap(datasets_, maxs_) :
	return ['{}.h{}.no_merging.no_downs.no_merging'.format(dataset_, h_)
		for dataset_ in datasets_
		for h_ in maxs_]

# partial paramerization of a merging (threshold and neg. threshold)
def merging(thrs_, negthrs_) :
	return ['.merged_e{}_m{}_t{}_n{}'.format(err_, max_, thresh_, neg_)
		for err_ in error_rates
		for max_ in max_errs
		for thresh_ in thrs_
		for neg_ in negthrs_] + ['.no_merging']

# datasets postprocessed to a specified list of max cov
def postproc(datasets_, thrs_, negthrs_, maxs_) :
	return ['{}.hN{}{}.no_merging'.format(dataset_, merging_, downsampling_)
		for dataset_ in datasets_
		for merging_ in merging(thrs_, negthrs_)
		for downsampling_ in downs(maxs_)]

# datasets both processed by whatshap and postprocessed to list of max cov.
def sliceof(datasets_, thrs_, negthrs_, maxs_) :
	return whatshap(datasets_, maxs_) + postproc(datasets_, thrs_, negthrs_, maxs_)

# define a subset of the datasets in terms of chromosomes and coverages
def datasubset(chromosomes_, coverages_, modes_) :
	return ['{}.pacbio.child.chr{}.cov{}.{}'.format(data_, chromosome_, coverage_, mode_)
		for data_ in data
		for chromosome_ in chromosomes_
		for coverage_ in coverages_
		for mode_ in modes_]

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('output/whatshap/{pattern}.{ext}',
			pattern = whatshap(datasets, [15, 20]),
			ext = exts),

		expand('output/hapchat/{pattern}.{ea}.bN_0.{ext}',
			pattern = sliceof(datasets, [6, 17], [3], [15, 20]),
			ea = ea_vals,
			ext = exts),

# coming up ..
rule next :
	input :
		expand('output/core_wh/{pattern}.{ext}',
			pattern = whatshap(
				datasubset([21], [5, 10], realignment),
				[15, 20]),
			ext = exts),

		expand('output/hapchat/{pattern}.05_00001.{bal}.{ext}',
			pattern = sliceof(
				datasubset(
                                        chromosomes,
					[5, 10, 15, 20],
					realignment),
				[6], [3], [25]),
			bal = ['bN_0', 'b20_45'],
			ext = exts)

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

	output : 'output/whatshap/' + whatshap_pattern + '.phased.vcf'

	log :
		log = 'output/whatshap/' + whatshap_pattern + '.log',
		time = 'output/whatshap/' + whatshap_pattern + '.time'

	message : '''

   running whatshap on bam/vcf pair :

   {input.bam} / {input.vcf}

   selecting coverage down to {wildcards.h} '''

	shell : '''

   {time} -v -o {log.time} \
      whatshap phase -o {output} {params.realignment} -H {params.h} \
         {input.vcf} {input.bam} > {log.log} 2>&1 '''

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

	message : 'running core whatshap dp on {input}'

	shell : '''

   {time} -v -o {log.time} \
      {corewh} -h {output} -a {input} > {log.wif} 2> {log.log} '''

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

	output : 'output/hapchat/' + full_pattern + '.hap'

	log :
		log = 'output/hapchat/' + full_pattern + '.log',
		time = 'output/hapchat/' + full_pattern + '.time'

	message : 'running hapchat on {input}'

	shell : '''

   {time} -v -o {log.time} \
      {hapchat} -i {input} -o {output} -A \
         -e {params.epsilon} -a {params.alpha} \
         -b {params.balance_cov} -r {params.balance_ratio} \
            > {log.log} 2>&1 '''

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

	message : '''

   comparing inferred phasing:

   {input.vcf}

   to true phasing: {input.true} '''

	shell : '''

   {compare} --switch-error-bed {output.bed} \
      {input.true} {input.vcf} > {output.diff} '''

# convert old-skool hap format to a phased vcf
rule phase_vcf :
	input :
		script = 'scripts/subvcf.py',
		hap = '{dir}/' + full_pattern + '.hap',
		blocks = 'wif/' + post_pattern + '.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : '{dir}/' + full_pattern + '.phased.vcf'

	log : '{dir}/' + full_pattern + '.phased.vcf.log'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} 2> {log} '''

#
# compute MEC score of phased vcf wrt instance, as a wif file
#----------------------------------------------------------------------
rule mec_score :
	input :
		script = 'scripts/wiftools.py',
		vcf = '{dir}/' + full_pattern + '.phased.vcf',
		wif = 'wif/' + post_pattern + '.wif'

	output : '{dir}/' + full_pattern + '.mec'

	message : '''

   infer mec score of:

   {input.vcf}

   with respect to: {input.wif} '''

	shell : '''

   python {input.script} -v {input.vcf} {input.wif} > {output} '''

#
# get sitewise details of input/run, e.g., coverage, switch error, etc.
#----------------------------------------------------------------------
rule sitewise_details :
	input :
		script = 'scripts/sitesinfo.py',
		sites = 'wif/' + post_pattern + '.wif.info_/sites_',
		swerrs = '{dir}/' + full_pattern + '.bed'

	output : '{dir}/' + full_pattern + '.sites'

	message : 'obtain sitewise details {output}'

	shell : '''

   python {input.script} -s {input.swerrs} {input.sites} > {output} '''
