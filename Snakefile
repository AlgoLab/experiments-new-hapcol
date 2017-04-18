#
# for testing and comparing the new hapcol
#----------------------------------------------------------------------
data_dir = '/data/phasing-comparison-experiments'
hap_dir = '/home/prj_rnabwt/haplotyping'

# datasets
datasets = ['ashk', 'sim']
individuals = ['child'] # mother, father, ..
coverages = [5, 10, 15, 20, 25, 30, 'all']
chromosomes = [1, 21]

# scripts
scripts = ['get.matrix.py', 'subsam.py', 'get.variants.py', 'wiftools.py']
scregex = '('+'|'.join([s for s in scripts])+')'

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}'
pattern_ext = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.shuf{seed,[0-9]+}.max{max,[0-9]+}'

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
		expand('wif/{paths}{merge}.wif.info_/blocks_',
			paths = datasets_exp,
			merge = ['', '.merged']),

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

# take only reads with flag 0 or 16 and non-empty set of SNVs, number
# the lines and convert this to wif, using the var file
rule get_wif :
	input :
		script = 'scripts/get.matrix.py',
		sfi = 'wif/' + dataset_pattern + '.sfi',
		var = 'vcf/' + vcf_pattern + '.var'

	output :
		wif = 'wif/' + dataset_pattern + '.wif',
		sub = 'wif/' + dataset_pattern + '.subset'

	log :
		log = 'wif/' + dataset_pattern + '.wif.log',
		time = 'wif/' + dataset_pattern + '.wif.time'

	message : '''

   obtaining wif file {output.wif} from {input.sfi} / {input.var} pair '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      awk '(($2==0 || $2==16) && ($9 != "#insertions")) {{print NR,$0}}' \
         {input.sfi} | python {input.script} -s -w -t {output.sub} {input.var} | \
            sort -nk1,1 > {output.wif} 2> {log.log} '''

# get a snv/fragment info (snv) file from a bam / var pair
rule get_sfi :
	input :
		script = 'scripts/subsam.py',
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
      samtools view {input.bam} | python {input.script} -v {input.var} \
         > {output} 2> {log.log} '''

# get a variants (SNVs) file from a vcf file
rule get_var :
	input :
		script = 'scripts/get.variants.py',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'vcf/' + vcf_pattern + '.var'
	log :
		log = 'vcf/' + vcf_pattern + '.var.log',
		time = 'vcf/' + vcf_pattern + '.var.time'

	message : 'obtaining SNVs file {output} from {input.vcf}'
	shell : '''

   /usr/bin/time -v -o {log.time} python {input.script} \
      tmp_{wildcards.dataset}_{wildcards.individual} {input.vcf} \
         > {log.log} 2>&1
   mv tmp_{wildcards.dataset}_{wildcards.individual}_{wildcards.chromosome}.var \
      {output} '''

#
# downsample a wif file to a specified max coverage
#----------------------------------------------------------------------
rule link_downsampled_wif :
	input : 'wif/' + pattern_ext + '.wif'
	output : 'input_wif/' + pattern_ext + '.wif'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

# extract from wif (or a set of reads) file (the lines of) the sample
rule extract_sample :
	input :
		source = 'wif/' + pattern_ext + '.{ext}',
		sample = 'wif/' + pattern_ext + '.wif.sample'

	output : 'wif/' + pattern_ext + '.{ext,(wif|subset)}'
	message : 'extract lines {input.sample} from {input.source}'
	shell : '''

   awk '{{printf "%.20d %s\\n", NR, $0}}' {input.source} | join - \
      <(awk '{{printf "%.20d\\n", $1}}' {input.sample} | sort) | \
         sed 's/^[0-9]* //' > {output} '''

# greedily downsample wif to a coverage according to a shuffle
rule downsample :
	input :
		script = 'scripts/wiftools.py',
		wif = 'wif/' + dataset_pattern + '.wif',
		shuf = 'wif/' + dataset_pattern + '.wif.lines.shuf{seed}'

	output : 'wif/' + pattern_ext + '.wif.sample'
	log :
		log = 'wif/' + pattern_ext + '.wif.sample.log',
		time = 'wif/' + pattern_ext + '.wif.sample.time'

	message : '''

   downsampling {input.wif} to coverage {wildcards.cov}
   according to {input.shuf} '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} -s {wildcards.cov} {input.shuf} {input.wif} \
         > {output} 2> {log.log} '''

# seeded pseudorandom shuffle of lines of a file (cf. gnu.org)
rule permute_lines :
	input : '{path}.lines'
	output : '{path}.lines.shuf{seed,[0-9]+}'
	message : 'pseudorandom shuffle of {input} with seed {wildcards.seed}'
	shell : '''

   shuf {input} --random-source=<(openssl enc -aes-256-ctr \
      -pass pass:"{wildcards.seed}" -nosalt </dev/zero 2>/dev/null) > {output} '''

# get lines (numbers) from a file
rule get_lines :
	input : '{path}'
	output : '{path}.lines'
	message : 'obtain lines (numbers) from {input}'
	shell : ''' awk '{{print NR}}' {input} > {output} '''

#
# select reads from bam file (1:1) corresponding to downsampled wif
#----------------------------------------------------------------------
rule select_reads :
	input :
		script = 'scripts/subsam.py',
		bam = 'bam/' + dataset_pattern + '.bam',
		sub = 'wif/' + pattern_ext + '.subset'

	output : 'bam/' + pattern_ext + '.bam'
	log :
		log = 'bam/' + pattern_ext + '.bam.log',
		time = 'bam/' + pattern_ext + '.bam.time'

	message : 'selecting reads from {input.bam} according to {input.sub}'
	shell : '''

   /usr/bin/time -v -o {log.time} \
      samtools view -h {input.bam} | python {input.script} -s {input.sub} | \
         samtools view -hb - > {output} 2> {log.log} '''

#
# obtain a (red-blue-) merged wif from a wif
#----------------------------------------------------------------------
rule link_merged_wif :
	input : 'wif/{pattern}.merged.wif'
	output : 'merged_wif/{pattern}.merged.wif'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

# merge a wif according to connected components
rule merge_wif :
	input :
		script = 'scripts/wiftools.py',
		wif = 'wif/{pattern}.wif',
		ccs = 'wif/{pattern}.ccs'

	output : 'wif/{pattern}.merged.wif'
	log :
		log = 'wif/{pattern}.merged.wif.log',
		time = 'wif/{pattern}.merged.wif.time'

	message : '''

   merge reads of {input.wif} according to {input.ccs}, producing {output} '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} -c {input.ccs} {input.wif} \
         > {output} 2> {log.log} '''

# build red-blue graph, and obtain connected components (one per line)
rule get_redblue_ccs :
	input :
		script = 'scripts/build-red-blue-graph.py',
		mat = 'wif/{pattern}.mat'

	output : 'wif/{pattern}.ccs'
	log :
		log = 'wif/{pattern}.ccs.log',
		time = 'wif/{pattern}.ccs.time'

	message : '''

   obtain connected components {output} from {input.mat} using red-blue graph '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} {input.mat} > {output} 2> {log.log} '''

# convert wif file to a (zygosity) matrix
rule get_zygosity_matrix :
	input :
		script = 'scripts/wiftools.py',
		wif = 'wif/{pattern}.wif'

	output : 'wif/{pattern}.mat'
	log :
		log = 'wif/{pattern}.mat.log',
		time = 'wif/{pattern}.mat.time'

	message : 'converting {input.wif} to (zygosity) matrix: {output}'
	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} -z {input.wif} | awk '{{ \
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

#
# obtain properties of a wif file
#----------------------------------------------------------------------
rule wif_info :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif'
	output : '{path}.wif.info_/blocks_'
	message : 'obtaining info for {input.wif}'
	shell : 'python {input.script} -i {input.wif}'
