#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
include : 'setup.snake'
time = '/usr/bin/time'

# softwares
corewh = 'programs/core_whatshap/build/dp'

# everything of max coverage 20
whatshap_one = ['{}.raw.h{}'.format(dataset, h)
	for dataset in datasets
        for h in [15, 20]]
post_one = ['{}.raw.hN{}sh1-max{}'.format(dataset, merge, max)
	for dataset in datasets
	for merge in ['.', '.merged.']
    	for max in [15, 20]]
slice_one = whatshap_one + post_one

# datasets for chr21 and some of the smaller average coverages
datasubset_one = ['{}.pacbio.child.chr{}.cov{}'.format(data, chromosome, coverage)
        for data in data
        for chromosome in [21]
	for coverage in [5, 10, 15, 20]]

# subsets of max coverage 20
whatshap_subset_one = ['{}.raw.h{}'.format(dataset, h)
	for dataset in datasubset_one
        for h in [15, 20]]
post_subset_one = ['{}.raw.hN{}sh1-max{}'.format(dataset, merge, max)
	for dataset in datasubset_one
	for merge in ['.', '.merged.']
    	for max in [15, 20]]
slice_subset_one = whatshap_subset_one + post_subset_one

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('output/whatshap/{pattern}.diff',
			pattern = whatshap_one),

		expand('output/core_wh/{pattern}.diff',
			pattern = slice_subset_one)

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
	input : 'wif/' + whatshap_pattern + '{pattern}.wif'
	output : 'output/core_wh/' + whatshap_pattern + '{pattern}.hap'

	log :
		log = 'output/core_wh/' + whatshap_pattern + '{pattern}.log',
		time = 'output/core_wh/' + whatshap_pattern + '{pattern}.time',
		wif = 'output/core_wh/' + whatshap_pattern + '{pattern}.wif'

	message : 'running core whatshap dp on {input}'

	shell : '''

   {time} -v -o {log.time} \
      {corewh} -h {output} -a {input} > {log.wif} 2> {log.log} '''

#
# compare phased vcfs to true phasing
#----------------------------------------------------------------------
rule compare_vcfs :
	input :
		true = 'vcf/' + vcf_pattern + '.phased.vcf',
		vcf = '{dir}/' + whatshap_pattern + '{pattern}.phased.vcf'

	output : '{dir}/' + whatshap_pattern + '{pattern}.diff'

	message : '''

   comparing inferred phasing:

   {input.vcf}

   to true phasing: {input.true} '''

	shell : 'whatshap compare {input.true} {input.vcf} > {output}'

# convert old-skool hap format to a phased vcf
rule phase_vcf :
	input :
		script = 'scripts/subvcf.py',
		hap = '{dir}/' + whatshap_pattern + '{pattern}.hap',
		blocks = 'wif/' + whatshap_pattern + '{pattern}.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : '{dir}/' + whatshap_pattern + '{pattern}.phased.vcf'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} '''
