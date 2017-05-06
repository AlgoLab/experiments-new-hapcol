#
# for generating all of the data needed for the experiments
#----------------------------------------------------------------------
#
data_dir = '/data/phasing-comparison-experiments'
hap_dir = '/home/prj_rnabwt/haplotyping'

# scripts and programs
whatshap = 'programs/whatshap/venv/bin/whatshap'
scripts = ['wiftools.py']
scripts_regex = '('+'|'.join([s for s in scripts])+')'

# datasets
data = ['ashk', 'sim']
platforms = ['pacbio']
individuals = ['child'] # mother, father, ..
chromosomes = [1, 21]
coverages = [5, 10, 15, 20, 25, 30, 'all']

# whatshap processing
realignment = ['raw', 'realigned'] # whatshap realignment
hs = [15, 20, 25, 30, 'N'] # whatshap read selection

# remaining processing
seeds = [1] # 2, 3, .. for (pseudo-) random downsampling
max_covs = [15, 20, 25, 30]

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}'
whatshap_pattern = dataset_pattern + '.{realignment,(raw|realigned)}.h{h,([0-9]+|N)}'

# common lists
datasets = ['{}.pacbio.child.chr{}.cov{}'.format(data, chromosome, coverage)
	for data in data
	for chromosome in chromosomes
	for coverage in coverages]

whatshap_downsample = ['{}.raw.h{}'.format(dataset, h)
	for dataset in datasets
	for h in hs]

outside_whatshap = ['{}.raw.hN{}sh1-max{}'.format(dataset, merge, max)
	for dataset in datasets
	for merge in ['.', '.merged.']
        for max in max_covs]

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('wif/{pattern}.wif.info_/blocks_',
			pattern = whatshap_downsample + outside_whatshap),

		expand('vcf/{data}.child.chr{chr}.phased.vcf',
			data = data,
			chr = chromosomes)

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
	output : 'scripts/{script,' + scripts_regex + '}'
	message : 'linking script {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_red_blue :
        input : hap_dir + '/red_blue/build-red-blue-graph.py'
	output : 'scripts/build-red-blue-graph.py'
	message : 'linking script {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

#
# obtain a wif file from a bam / vcf pair using whatshap
#----------------------------------------------------------------------
rule get_wif :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	params :
                realignment = '', # TODO: add this function
		h = lambda wildcards :
			'1000' if wildcards.h == 'N' else wildcards.h

	output : 'wif/' + whatshap_pattern + '.wif'

	log :
		transcript = 'wif/' + whatshap_pattern + '.wif.transcript',
		log = 'wif/' + whatshap_pattern + '.wif.log',
		time = 'wif/' + whatshap_pattern + '.wif.time'

	message : '''

   obtaining wif file {output} from {input.bam} / {input.vcf} pair '''

	shell : '''
   
   /usr/bin/time -v -o {log.time} \
      {whatshap} phase -o /dev/null {params.realignment} \
         --output-wif {output} -H {params.h} \
         {input.vcf} {input.bam} > {log.transcript} 2> {log.log} '''

#
# downsample a wif file to a specified max coverage
#----------------------------------------------------------------------
rule extract_sample :
	input :
		source = '{path}.wif',
		sample = '{path}.wif.sample.sh{seed}-max{max}'

	output : '{path}.sh{seed,[0-9]+}-max{max,[0-9]+}.wif'
	message : 'extract lines {input.sample} from {input.source}'

	shell : '''

   awk '{{printf "%.20d %s\\n", NR, $0}}' {input.source} | join - \
      <(awk '{{printf "%.20d\\n", $1}}' {input.sample} | sort) | \
         sed 's/^[0-9]* //' > {output} '''

# greedily downsample wif to a coverage according to a shuffle
rule downsample :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif',
		shuf = '{path}.wif.lines.sh{seed}'

	output : '{path}.wif.sample.sh{seed,[0-9]+}-max{max,[0-9]+}'

	log :
		log = '{path}.wif.sample.sh{seed}-max{max}.log',
		time = '{path}.wif.sample.sh{seed}-max{max}.time'

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
	output : '{path}.lines.sh{seed,[0-9]+}'
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
# obtain a (red-blue-) merged wif from a wif
#----------------------------------------------------------------------
rule merge_wif :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif',
		ccs = '{path}.ccs'

	output : '{path}.merged.wif'

	log :
		log = '{path}.merged.wif.log',
		time = '{path}.merged.wif.time'

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
		mat = '{path}.mat'

	output : '{path}.ccs'

	log :
		log = '{path}.ccs.log',
		time = '{path}.ccs.time'

	message : '''

   obtain connected components {output} from {input.mat} using red-blue graph '''

	shell : '''

   /usr/bin/time -v -o {log.time} \
      python {input.script} {input.mat} > {output} 2> {log.log} '''

# convert wif file to a (zygosity) matrix
rule get_zygosity_matrix :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif'

	output : '{path}.mat'

	log :
		log = '{path}.mat.log',
		time = '{path}.mat.time'

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
