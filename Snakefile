#
# for testing and comparing the new hapcol
#----------------------------------------------------------------------
data_dir = '/data/haplotyping_data/phasing-comparison-experiments'
hap_dir = '/home/prj_rnabwt/haplotyping'

# datasets
datasets = ['ashk', 'sim']
individuals = ['child'] # mother, father, ..
coverages = [10, 20, 40, 'all']
chromosomes = [1, 21]

# scripts
scripts = ['get.matrix.py', 'get.sfi.py', 'get.variants.py', 'wiftools.py']
scregex = '('+'|'.join([s for s in scripts])+')'

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}'

datasets_exp = ['{}.pacbio.{}.chr{}.cov{}'.format(dataset, individual, chromosome, coverage)
	for dataset in datasets
	for individual in individuals
	for chromosome in chromosomes
	for coverage in coverages]

# versions (commits) of hapcol
hapcol_versions = {
        'original' : '68a9f3fbce84020e2faef054fd07dfb1bd86052f',
	'balanced_only' : '8651eb1782c32f77f638317dad7095d831b864af',
	'increase_k_only' : '0a5415edc697304d0eb9a80b5832ad0a572c514b',
        'increase_k_and_balancing' : '7a9cea5b63c5e0d48a16ad15a373aa294b73f119' }

# test modes (parameter sets)
test_modes = {
        'none' : '',
        'with_beta' : '-b 0' }

# master rule
rule master :
	input :
		expand('input_wif/{paths}.wif', paths = datasets_exp),
		expand('merged_wif/{paths}.merged.wif', paths = datasets_exp),
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
# link to a script in the haplotyping/scripts directory, etc.
#----------------------------------------------------------------------
rule link_script :
        input : hap_dir + '/scripts/{script}'
	output : 'scripts/{script,'+scregex+'}'
	message : 'linking script {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_red_blue :
        input : hap_dir + '/red_blue/build-red-blue-graph.py'
	output : 'scripts/build-red-blue-graph.py'
	message : 'linking script {input} to {output}'
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
		scr = 'scripts/get.matrix.py',
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
         {input.sfi} | python {input.scr} -s -w {input.var} | \
            sort -nk1,1 > {output} 2> {log.log} '''

# get a snv/fragment info (snv) file from a bam / var pair
rule get_sfi :
	input :
		scr = 'scripts/get.sfi.py',
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
      samtools view {input.bam} | python {input.scr} {input.var} \
         > {output} 2> {log.log} '''

# get a variants (SNVs) file from a vcf file
rule get_var :
	input :
		scr = 'scripts/get.variants.py',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'vcf/' + vcf_pattern + '.var'
	log :
		log = 'vcf/' + vcf_pattern + '.var.log',
		time = 'vcf/' + vcf_pattern + '.var.time'

	message : 'obtaining SNVs file {output} from {input.vcf}'
	shell : '''

   /usr/bin/time -v -o {log.time} python {input.scr} \
      tmp_{wildcards.dataset}_{wildcards.individual} {input.vcf} \
         > {log.log} 2>&1
   mv tmp_{wildcards.dataset}_{wildcards.individual}_{wildcards.chromosome}.var \
      {output} '''

#
# obtain a (red-blue-) merged wif from a wif
#----------------------------------------------------------------------
rule link_merged_wif :
	input : 'wif/' + dataset_pattern + '.merged.wif'
	output : 'merged_wif/' + dataset_pattern + '.merged.wif'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

# merge a wif according to connected components
rule merge_wif :
	input :
		scr = 'scripts/wiftools.py',
		wif = 'wif/' + dataset_pattern + '.wif',
		ccs = 'wif/' + dataset_pattern + '.ccs'

	output : 'wif/' + dataset_pattern + '.merged.wif'
	log :
		log = 'wif/' + dataset_pattern + '.merged.wif.log',
		time = 'wif/' + dataset_pattern + '.merged.wif.time'

	message : '''

   merge reads of {input.wif} according to {input.ccs}, producing {output} '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.scr} -c {input.ccs} {input.wif} \
         > {output} 2> {log.log} '''

# build red-blue graph, and obtain connected components (one per line)
rule get_redblue_ccs :
	input :
		scr = 'scripts/build-red-blue-graph.py',
		mat = 'wif/' + dataset_pattern + '.mat'

	output : 'wif/' + dataset_pattern + '.ccs'
	log :
		log = 'wif/' + dataset_pattern + '.ccs.log',
		time = 'wif/' + dataset_pattern + '.ccs.time'

	message : '''

   obtain connected components {output} from {input.mat} using red-blue graph '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.scr} {input.mat} > {output} 2> {log.log} '''

# convert wif file to a (zygosity) matrix
rule get_zygosity_matrix :
	input :
		scr = 'scripts/wiftools.py',
		wif = 'wif/' + dataset_pattern + '.wif'

	output : 'wif/' + dataset_pattern + '.mat'
	log :
		log = 'wif/' + dataset_pattern + '.mat.log',
		time = 'wif/' + dataset_pattern + '.mat.time'

	message : 'converting {input.wif} to (zygosity) matrix: {output}'
	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.scr} -z {input.wif} | awk '{{ \
         for(i=3; i<= NF; ++i) {{ \
            if($i=="0"){{$i="G"}} ; if($i=="1"){{$i="C"}} }} ; \
         print }}' > {output} 2> {log.log} '''

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
   git pull
   git checkout {params.version}
   cmake ../src
   make -j 16
   mv hapcol ../../{wildcards.version}
   git checkout master
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

#
# test a version of hapcol
#----------------------------------------------------------------------
rule test_hapcol_version :
	input : 'wif/sim.pacbio.child.chr21.cov10.wif'
	output : 'test_output/hapcol-{version}.mode-{mode}.hap'
	params : lambda wildcards : test_modes[wildcards.mode]
	log :
		log = 'test_output/hapcol-{version}.mode-{mode}.log',
		time = 'test_output/hapcol-{version}.mode-{mode}.time'

	message : '''

   testing hapcol version {wildcards.version} with parameters {params}
   on {input} '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      ./hapcol_builds/{wildcards.version}/hapcol -i {input} -o {output} -A \
         {params} > {log.log} 2>&1 '''
