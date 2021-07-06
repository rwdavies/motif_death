if GENOTYPER == "HaplotypeCaller":
    GENOTYPING_MULTITHREAD_FLAG="-nct"
elif GENOTYPER == "UnifiedGenotyper":
    GENOTYPING_MULTITHREAD_FLAG="-nt"

CHR_CHUNKS = range(1, 2) ## disable - too complicated as crashes if out of rang - run 1 chunk only
TREEMIX_MIGRANT_RANGE = list(range(0, 10))

IBAMS=""
for species in SPECIES_LIST:
    IBAMS = IBAMS + " -I mapping/" + species + "/" + species + "." + BAM_SUFFIX

MERGE_VCF_GATK_INPUT=""
for piece in CHR_CHUNKS:
    for chr in CHR_LIST_ONLY_AUTOS:
    	MERGE_VCF_GATK_INPUT = MERGE_VCF_GATK_INPUT + " -V vcf/" + VCF_PREFIX + ".chr" + str(chr) + ".filtered.piece" + str(piece) + ".vcf.gz"
	
rule downstream_all:
    input:
        [
	expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX),
	expand("coverage/coverage.{species}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS, species = SPECIES_LIST),
	expand("treemix/{treemix_prefix}.treemix.migrants.{migrants}.out.treeout.gz", treemix_prefix = TREEMIX_PREFIX, migrants = TREEMIX_MIGRANT_RANGE)
	]

# rule call_and_tree:
#     input:
#         [expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX), expand("treemix/{treemix_prefix}.treemix.migrants.{migrants}.out.treeout.gz", treemix_prefix = TREEMIX_PREFIX, migrants = TREEMIX_MIGRANT_RANGE)]

# rule all_calling:
#     input:
#         expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX)

# rule all_calling_chr:
#     input:
#         expand("vcf/{vcf_prefix}.chr{chr}.piece{piece}.vcf.gz", vcf_prefix = VCF_PREFIX, chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS)

# rule treemix_input:
#     input:
#         expand("treemix/{vcf_prefix}.treemix.frq.gz", vcf_prefix = VCF_PREFIX)

rule treemix:
    input:
        expand("treemix/{treemix_prefix}.treemix.migrants.{migrants}.out.treeout.gz", treemix_prefix = TREEMIX_PREFIX, migrants = TREEMIX_MIGRANT_RANGE)

# rule all_doc:
#     input:
#         expand("coverage/coverage.{species}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS, species = SPECIES_LIST)


rule call_chr:
    input:
        ref = REF_DIR + REFNAME + ".fa",
        bams = expand("mapping/{species}/{species}.{bam_suffix}", species = SPECIES_LIST, bam_suffix = BAM_SUFFIX),
        bais = expand("mapping/{species}/{species}.{bam_suffix}.bai", species = SPECIES_LIST, bam_suffix = BAM_SUFFIX)
    output:
        vcf = expand("vcf/{vcf_prefix}.chr{{chr}}.piece{{piece}}.vcf.gz", vcf_prefix = VCF_PREFIX),
        tbi = expand("vcf/{vcf_prefix}.chr{{chr}}.piece{{piece}}.vcf.gz.tbi", vcf_prefix = VCF_PREFIX)
    params:
        N='call_chr',
        threads = GENOTYPING_THREADS,
        queue = GENOTYPING_QUEUE
    wildcard_constraints:
        chr = WILDCARD_CHR_CONSTRAINT,
        piece='\d{1,3}'	
    shell:
        'mkdir -p vcf/logs && '
        'if [ \"{OPERATE_GATK_PER_CHR}\" == \"FALSE" ]; then '
        'minus_L=\" \" ; else '
        'minus_L=\"-L {GATK_CHR_PREFIX}{wildcards.chr}\" ; fi && '
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T {GENOTYPER} '
        '-R {input.ref} '
        '{GENOTYPING_MULTITHREAD_FLAG} {params.threads} '
        '{IBAMS} '
        '${{minus_L}} '
        '-o vcf/{VCF_PREFIX}.chr{wildcards.chr}.piece{wildcards.piece}.vcf.gz '
        '2> '
        'vcf/logs/{VCF_PREFIX}.chr{wildcards.chr}.piece{wildcards.piece}.vcf.gz.log '

# TODO Robbie Temp - delete
##        'REGION_START=$(({wildcards.piece} * 10000000 + 1)) && '
##        'REGION_END=$((({wildcards.piece} + 1) * 10000000)) && '
##        '-L {GATK_CHR_PREFIX}{wildcards.chr}:${{REGION_START}}-${{REGION_END}} '


rule filter_chr:
    input:
        input_vcf = expand("vcf/{vcf_prefix}.chr{{chr}}.piece{{piece}}.vcf.gz", vcf_prefix = VCF_PREFIX),
        input_tbi = expand("vcf/{vcf_prefix}.chr{{chr}}.piece{{piece}}.vcf.gz.tbi", vcf_prefix = VCF_PREFIX),
        ref = REF_DIR + REFNAME + ".fa"
    output:
        filtered_vcf = expand("vcf/{vcf_prefix}.chr{{chr}}.filtered.piece{{piece}}.vcf.gz", vcf_prefix = VCF_PREFIX),
        filtered_tbi = expand("vcf/{vcf_prefix}.chr{{chr}}.filtered.piece{{piece}}.vcf.gz.tbi", vcf_prefix = VCF_PREFIX)
    params: N='filter_chr', threads=1, queue = "short.qc"
    wildcard_constraints:
        chr = WILDCARD_CHR_CONSTRAINT,
        piece='\d{1,3}'
    shell:
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T VariantFiltration '
        '-R {input.ref} '
        '-V {input.input_vcf} '
        '--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" '
        '--filterName "HardFilter" '
        '-o {output.filtered_vcf} 2> '
        'vcf/logs/{VCF_PREFIX}.chr{chr}.filtered.piece{piece}.vcf.gz.log && '
        'rm {input.input_vcf} && rm {input.input_tbi} '


rule merge_chr:
    input:
        vcf = expand("vcf/{vcf_prefix}.chr{chr}.filtered.piece{piece}.vcf.gz", vcf_prefix = VCF_PREFIX, chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        tbi = expand("vcf/{vcf_prefix}.chr{chr}.filtered.piece{piece}.vcf.gz.tbi", vcf_prefix = VCF_PREFIX, chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        ref = REF_DIR + REFNAME + ".fa"	
    output:
        merged_vcf = expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX),
        merged_tbi = expand("vcf/{vcf_prefix}.filtered.vcf.gz.tbi", vcf_prefix = VCF_PREFIX)
    params: N='merge_chr', threads=1, queue = "short.qc"
    wildcard_constraints:
        chr = WILDCARD_CHR_CONSTRAINT
    shell:
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T CombineVariants '
        '-R {input.ref} '
        '{MERGE_VCF_GATK_INPUT} '
        '--assumeIdenticalSamples '
        '-o {output.merged_vcf} 2> '
        '{output.merged_vcf}.log && '
        'rm {input.vcf} && rm {input.tbi}'

rule calculate_doc:
    input:
        ref = REF_DIR + REFNAME + ".fa",
        bam = expand("mapping/{{species}}/{{species}}.{bam_suffix}", bam_suffix = BAM_SUFFIX),
        bais = expand("mapping/{{species}}/{{species}}.{bam_suffix}.bai", bam_suffix = BAM_SUFFIX)
    output:
        outfile = "coverage/coverage.{species}.chr{chr}.txt.gz"
    params:
        N='get_depth_of_coverage',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
        species='\w{1,40}'	
    shell:
        'mkdir -p coverage && '
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T DepthOfCoverage '
        '-R {input.ref} '
        '-L {GATK_CHR_PREFIX}{wildcards.chr} '
        '--minMappingQuality 17 --minBaseQuality 0 '
        '-I {input.bam} '
        '--downsample_to_coverage 1000 '
        '-o coverage/coverage.{wildcards.species}.chr{wildcards.chr}.txt '
        '> {output.outfile}.log && '
        'gzip -1 coverage/coverage.{wildcards.species}.chr{wildcards.chr}.txt'


rule prepare_reference:
    input:
        amb = REF_DIR + REFNAME + ".fa.amb",
        ann = REF_DIR + REFNAME + ".fa.ann"
    output:
        summary = REF_DIR + REFNAME + ".fa.summary.txt"    
    params:
        N='prepare_reference',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
    shell:
        'R -f {R_GET_GENOME_STATS} --args {REF_DIR} {REFNAME}.fa'


rule get_callable_regions:
    input:
        infile = expand("coverage/coverage.{{species}}.chr{chr}.txt.gz", chr = CHR_LIST_ONLY_AUTOS),
        summary = REF_DIR + REFNAME + ".fa.summary.txt"
    output:
        beds = expand("coverage/coverage.{{species}}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS)
    params:
        N='get_callable_regions',
        threads = 2,
        queue = "short.qc"
    wildcard_constraints:
        species='\w{1,40}'
    shell:
        'echo start && '
        'R -f {R_GET_PER_SAMPLE_AVERAGE_COV} --args {wildcards.species} {REF_DIR} {REFNAME}.fa && '
        'echo done'


rule prepare_treemix:
    input:
        merged_vcf = expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX),
        merged_tbi = expand("vcf/{vcf_prefix}.filtered.vcf.gz.tbi", vcf_prefix = VCF_PREFIX)	
    output:
        treemix_input = expand("treemix/{treemix_prefix}.treemix.frq.gz", treemix_prefix = TREEMIX_PREFIX)
    params:
        N='prepare_treemix',
        threads=1,
        queue = "short.qc"
    shell:
        'mkdir -p treemix && '
        'echo recode as treemix format && date && '
        '${{PYTHON_352}} {VCF2TREEMIX} {input.merged_vcf} {output.treemix_input} > {output.treemix_input}.log && echo done && date'
        
        
rule run_treemix:
    input:
        treemix_input = expand("treemix/{treemix_prefix}.treemix.frq.gz", treemix_prefix = TREEMIX_PREFIX)
    output:
        treemix_output = expand("treemix/{treemix_prefix}.treemix.migrants.{{migrants}}.out.treeout.gz", treemix_prefix = TREEMIX_PREFIX)
    params:
        N='run_treemix',
        threads=TREEMIX_THREADS,
        queue = "short.qc"
    wildcard_constraints:
        migrants='\d{1,2}'
    shell:
        'echo begin treemix && date &&'
        '${{TREEMIX}} -i {input} -m {wildcards.migrants} -noss -se -k {TREEMIX_K} '
        '-root {TREEMIX_OUTGROUP} '
        '-n_warn 10 '
        '-o treemix/{TREEMIX_PREFIX}.treemix.migrants.{wildcards.migrants}.out &&'
        'echo done treemix && date && '
        'R -f {TREEMIX_R_PLOT} '
        '--args "./" {TREEMIX_R_PLOT_HELPER} '
        'treemix/{TREEMIX_PREFIX}.treemix.migrants.{wildcards.migrants}.out '
        'treemix/{TREEMIX_PREFIX}.treemix.frq.gz  && echo done && date'


# TODO: obsolete?
# config = {
#   "simple_repeat" : {
#      "yes" : HATBAG_SIMPLE_REPEAT,
#      "no" : "NULL"
#    },
#   "callable_bed" : {
#      "yes" : HATBAG_CALLABLE_BED,
#      "no" : "NULL"
#    }
# }

