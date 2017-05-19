#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
include : 'setup.snake'
time = '/usr/bin/time'

# softwares
corewh = 'programs/core_whatshap/build/dp'
hapchat = 'programs/increase-k-hapcol/build/hapcol'

# pattern (taking into account hapcol)
full_pattern = post_pattern + '{ea,(|.[0-9]+_[0-9]+)}{balancing,(|.balanced)}'

# epislon / alpha pairs for hapcol, and variants
ea_vals = ['05_1', '05_01', '05_001', '05_0001', '05_00001']
ea_two = ['01_1', '01_01', '01_001', '01_0001', '1_1', '1_01', '1_001', '1_0001']

# downsamplings to a max cov. of 20
downsamplings_20 = ['.downs_s{}_m{}'.format(seed, max)
	for seed in seeds
	for max in [15, 20]]

#
# everything of max coverage 20
#----------------------------------------------------------------------
whatshap_max20 = ['{}.h{}'.format(dataset, h)
	for dataset in datasets
        for h in [15, 20]]
post_max20 = ['{}.hN{}{}'.format(dataset, merging, downsampling)
	for dataset in datasets
	for merging in [''] + mergings
	for downsampling in downsamplings_20]
slice_max20 = whatshap_max20 + post_max20

#
# subset of datasets for chr21 and some of the smaller avg. coverages
#----------------------------------------------------------------------
subset_one = ['{}.pacbio.child.chr{}.cov{}.{}'.format(data, chromosome, coverage, mode)
        for data in data
        for chromosome in [21]
	for coverage in [5, 10, 15]
	for mode in realignment]

# subsets of max coverage 20
whatshap_s1_max20 = ['{}.h{}'.format(dataset, h)
	for dataset in subset_one
        for h in [15, 20]]
post_s1_max20 = ['{}.hN{}{}'.format(dataset, merging, downsampling)
	for dataset in subset_one
	for merging in [''] + mergings
    	for downsampling in downsamplings_20]
slice_s1_max20 = whatshap_s1_max20 + post_s1_max20

# subsets of max coverage 25
whatshap_s1_max25 = ['{}.h{}'.format(dataset, h)
	for dataset in subset_one
        for h in [15, 20, 25]]
post_s1_max25 = ['{}.hN{}{}'.format(dataset, merging, downsampling)
	for dataset in subset_one
	for merging in [''] + mergings
    	for downsampling in downsamplings_20]
slice_s1_max25 = whatshap_s1_max25 + post_s1_max25

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('output/whatshap/{pattern}.diff',
			pattern = whatshap_max20),

		expand('output/core_wh/{pattern}.diff',
			pattern = slice_s1_max20),

		expand('output/hapchat/{pattern}.{ea}.diff',
			pattern = slice_max20,
			ea = ea_vals)

#
# run whatshap on a bam / vcf pair
#----------------------------------------------------------------------
rule run_whatshap :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	params :
		realignment = '', # TODO: add this function
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
# run hapchat with increase-k on a wif file
#----------------------------------------------------------------------
rule run_hapchat :
	input : 'wif/' + post_pattern + '.wif'

	params :
		e = lambda wildcards : '0' + wildcards.ea.split('_')[0],
		a = lambda wildcards : '0.' + wildcards.ea.split('_')[1]

	output : 'output/hapchat/' + full_pattern + '.hap'

	log :
		log = 'output/hapchat/' + full_pattern + '.log',
		time = 'output/hapchat/' + full_pattern + '.time'

	message : 'running hapchat on {input}'

	shell : '''

   {time} -v -o {log.time} \
      {hapchat} -i {input} -o {output} -e {params.e} -a {params.a} -A \
         > {log.log} 2>&1 '''

#
# compare phased vcfs to true phasing
#----------------------------------------------------------------------
rule compare_vcfs :
	input :
		true = 'vcf/' + vcf_pattern + '.phased.vcf',
		vcf = '{dir}/' + full_pattern + '.phased.vcf'

	output : '{dir}/' + full_pattern + '.diff'

	message : '''

   comparing inferred phasing:

   {input.vcf}

   to true phasing: {input.true} '''

	shell : 'whatshap compare {input.true} {input.vcf} > {output}'

# convert old-skool hap format to a phased vcf
rule phase_vcf :
	input :
		script = 'scripts/subvcf.py',
		hap = '{dir}/' + full_pattern + '.hap',
		blocks = 'wif/' + post_pattern + '.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : '{dir}/' + full_pattern + '.phased.vcf'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} '''

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
