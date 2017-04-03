#
# for testing and comparing the new hapcol
#----------------------------------------------------------------------

# settings
data_dir = '/data/haplotyping_data/phasing-comparison-experiments'

datasets = ['ashk', 'sim']
individuals = ['child'] # mother, father, ..
coverages = [10, 20, 40, 'all']
chromosomes = [1, 21]

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}'

# versions (commits) of hapcol
hapcol_versions = {
        'original' : '68a9f3fbce84020e2faef054fd07dfb1bd86052f',
	'balanced_only' : '8651eb1782c32f77f638317dad7095d831b864af',
	'increase_k_only' : '0a5415edc697304d0eb9a80b5832ad0a572c514b',
        'increase_k_and_balancing' : '7a9cea5b63c5e0d48a16ad15a373aa294b73f119' }

# master rule
rule master :
	input :
		expand('input_wif/{dataset}.pacbio.{individual}.chr{chromosome}.cov{coverage}.wif',
			dataset = datasets,
			individual = individuals,
			chromosome = chromosomes,
			coverage = coverages),
		expand('hapcol_builds/{version}/hapcol',
			version = hapcol_versions)

#
# link to files from phasing comparison experiments directory
#----------------------------------------------------------------------
rule link_vcf :
	input : data_dir + '/vcf/' + vcf_pattern + '.{phase,(phased|unphased)}.vcf'
	output : 'vcf/' + vcf_pattern + '.{phase,(phased|unphased)}.vcf'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_bam_bai :
	input : data_dir + '/bam/' + dataset_pattern + '.{ext,(bam|bai)}'
	output : 'bam/' + dataset_pattern + '.{ext,(bam|bai)}'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

#
# obtain a wif file from a bam / vcf pair
#----------------------------------------------------------------------
rule link_wif :
        input : 'wif/' + dataset_pattern + '.wif'
	output : 'input_wif/' + dataset_pattern + '.wif'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

# take only reads with flag 0 or 16 and non-empyt set of SNVs, number
# the lines and convert this to wif, using the var file
rule get_wif :
	input :
		sfi = 'wif/' + dataset_pattern + '.sfi',
		var = 'vcf/' + vcf_pattern + '.var'

	output : 'wif/' + dataset_pattern + '.wif'
	log :
		log = 'wif/' + dataset_pattern + '.wif.log',
		time = 'wif/' + dataset_pattern + '.wif.time'

	message : '''

   obtaining wif file {output} from {input.sfi} / {input.var} pair '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      awk '(($2==0 || $2==16) && ($9 != "#insertions")) {{print NR,$0}}' \
         {input.sfi} | python scripts/get.matrix.py -s -w {input.var} | \
            sort -nk1,1 > {output} 2> {log.log} '''

# get a snv/fragment info (snv) file from a bam / var pair
rule get_sfi :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		var = 'vcf/' + vcf_pattern + '.var'

	output : 'wif/' + dataset_pattern + '.sfi'
	log :
		log = 'wif/' + dataset_pattern + '.sfi.log',
		time = 'wif/' + dataset_pattern + '.sfi.time'

	message : '''

   obtaining SNV/fragment info {output} from {input.bam} / {input.var} pair '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      samtools view {input.bam} | python scripts/get.sfi.py {input.var} \
         > {output} 2> {log.log} '''

# get a variants (SNVs) file from a vcf file
rule get_var :
	input : 'vcf/' + vcf_pattern + '.unphased.vcf'
	output : 'vcf/' + vcf_pattern + '.var'
	log :
		log = 'vcf/' + vcf_pattern + '.var.log',
		time = 'vcf/' + vcf_pattern + '.var.time'

	message : 'obtaining SNVs file {output} from {input}'
	shell : '''

   /usr/bin/time -v -o {log.time} python scripts/get.variants.py \
      tmp_{wildcards.dataset}_{wildcards.individual} {input} > {log.log} 2>&1
   mv tmp_{wildcards.dataset}_{wildcards.individual}_{wildcards.chromosome}.var \
      {output} '''

#
# build a version of hapcol
#----------------------------------------------------------------------
rule build_hapcol :
        input : 'hapcol_builds/.hapcol_obtained'
	output : 'hapcol_builds/{version}/hapcol'
	params : version = lambda wildcards : hapcol_versions[wildcards.version]
	message : 'building a version of hapcol called: {wildcards.version}'
	shell : '''

   mkdir -p hapcol_builds/{wildcards.version}
   mkdir hapcol_builds/HapCol/build
   cd hapcol_builds/HapCol/build
   git checkout {params.version}
   cmake ../src
   make -j 16
   mv hapcol ../../{wildcards.version}
   rm -rf ../build '''

# obtain hapcol (from its git repo)
rule get_hapcol :
        output : 'hapcol_builds/.hapcol_obtained'
	message : 'cloning hapcol from its git repository'
	shell : '''

   mkdir -p hapcol_builds
   cd hapcol_builds
   touch .hapcol_obtained
   git clone git@github.com:AlgoLab/HapCol.git '''
