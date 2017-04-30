#
# for generating all of the data needed for the new hapcol
#----------------------------------------------------------------------
#
data_dir = '/data/phasing-comparison-experiments'
hap_dir = '/home/prj_rnabwt/haplotyping'

# datasets
data = ['ashk', 'sim']
individuals = ['child'] # mother, father, ..
chromosomes = [1, 21]
states = ['full', 'hetero']
coverages = [5, 10, 15, 20, 25, 30, 'all']
seeds = [1] # 2, 3, .. for downsampling
max_covs = [20, 25, 30]

# scripts
scripts = ['subsam.py', 'subvcf.py', 'subvert.py', 'wiftools.py']
scripts_regex = '('+'|'.join([s for s in scripts])+')'

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.{state,(hetero|full)}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.{state,(hetero|full)}.cov{coverage,(all|[0-9]+)}'
pattern_ext = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.{state,(hetero|full)}.cov{coverage,(all|[0-9]+)}.shuf{seed,[0-9]+}.max{max,[0-9]+}'

# lists (the datasets in the form of a list)
datasets = ['{}.pacbio.{}.chr{}.{}.cov{}'.format(data, individual, chromosome, state, coverage)
	for data in data
	for individual in individuals
	for chromosome in chromosomes
	for state in states
	for coverage in coverages]

datasets_ext = ['{}.shuf{}.max{}'.format(dataset, seed, max)
	for dataset in datasets
	for seed in seeds
	for max in max_covs]

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('input_wif/{pattern}.wif',
			pattern = datasets + datasets_ext),
		expand('merged_wif/{pattern}.merged.wif',
			pattern = datasets + datasets_ext),
		expand('bam/{pattern}.bam',
			pattern = datasets_ext),

		expand('wif/{pattern}.{ext}.info_/blocks_',
			pattern = datasets + datasets_ext,
			ext = ['wif', 'merged.wif']),
		expand('vcf/{data}.child.chr{chr}.phased.vcf',
			data = data,
			chr = chromosomes)

#
# link to files from phasing comparison experiments directory
#----------------------------------------------------------------------
rule link_vcf :
	input : data_dir + '/vcf/{dataset}.{individual}.chr{chromosome}.{phase}.vcf'
	output : 'vcf/{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.{phase,(phased|unphased)}.vcf'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_bam_bai :
	input : data_dir + '/bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.{ext}'
	output : 'bam/{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.{ext,(bam|bai)}'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

#
# link to a script in the haplotyping/scripts directory, etc.
#----------------------------------------------------------------------
rule link_script :
        input : hap_dir + '/scripts/{script}'
	output : 'scripts/{script,' + scripts_regex + '}'
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
		script = 'scripts/subvert.py',
		var = 'vcf/' + vcf_pattern + '.var',
		sfi = 'wif/' + dataset_pattern + '.sfi'

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
      awk '(($2==0 || $2==16) && ($9 != "#insertions")) {{print NR,$0}}' {input.sfi} | \
         python {input.script} -w -t {output.sub} {input.var} > {output.wif} '''

# get a snv/fragment info (snv) file from a bam / var pair
rule get_sfi :
	input :
		script = 'scripts/subsam.py',
		bam = 'bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bam',
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
		script = 'scripts/subvcf.py',
		vcf = 'vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf'

	params : lambda wildcards :
		'-H' if wildcards.state == 'hetero' else '-v'

	output : 'vcf/' + vcf_pattern + '.var'
	log :
		log = 'vcf/' + vcf_pattern + '.var.log',
		time = 'vcf/' + vcf_pattern + '.var.time'

	message : 'obtaining SNVs file {output} from {input.vcf}'
	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} {params} {input.vcf} > {output} 2> {log.log} '''

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
		source = 'wif/' + dataset_pattern + '.{ext}',
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

   downsampling {input.wif} to coverage {wildcards.max}
   according to {input.shuf} '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} -s {wildcards.max} {input.shuf} {input.wif} \
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
		bam = 'bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bam',
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

# perform a sanity check to ensure that we downsample the bam correctly
rule sanity_check :
	input :
                wif = 'wif/' + pattern_ext + '.wif',
		check = 'wif/' + pattern_ext + '.sfi.wif'

	output : 'wif/' + pattern_ext + '.sfi.wif.diff'
	message : 'check {input.wif} against {input.check}'
	shell : 'diff {input.wif} {input.check} > {output}'

# extract the wif from the selected sfi for a sanity check against the
# downsampled wif
rule sanity_wif :
	input :
		script = 'scripts/subvert.py',
		var = 'vcf/' + vcf_pattern + '.var',
		sfi = 'wif/' + pattern_ext + '.sfi'

	output : 'wif/' + pattern_ext + '.sfi.wif'
	message : 'obtain wif from {input.sfi} for sanity check'
	shell : '''

   python {input.script} -w {input.var} {input.sfi} > {output} '''

# select reads from an sfi file in the same way we would for the bam
# file (e.g., for the purposes of a sanity check)
rule select_lines :
	input :
		sfi = 'wif/' + dataset_pattern + '.sfi',
		sub = 'wif/' + pattern_ext + '.subset'

	output : 'wif/' + pattern_ext + '.sfi'
	message : 'select lines {input.sub} from {input.sub}'
	shell : '''

   awk '{{printf "%.20d %s\\n", NR, $0}}' {input.sfi} | join - \
      <(awk '{{printf "%.20d\\n", $1}}' {input.sub} | sort) | \
         sed 's/^[0-9]* //' > {output} '''

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
# obtain properties of a wif file
#----------------------------------------------------------------------
rule wif_info :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif'
	output : '{path}.wif.info_/blocks_'
	message : 'obtaining info for {input.wif}'
	shell : 'python {input.script} -i {input.wif}'
