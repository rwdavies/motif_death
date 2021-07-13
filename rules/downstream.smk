rule downstream_all:
    input:
        [
	f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz",
	expand("coverage/coverage.{species}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS, species = SPECIES_LIST),
	expand(f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz", migrants = TREEMIX_MIGRANT_RANGE)
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
        expand(f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz", migrants = TREEMIX_MIGRANT_RANGE)

# rule all_doc:
#     input:
#         expand("coverage/coverage.{species}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS, species = SPECIES_LIST)


rule call_chr:
    input:
        ref = REF_DIR + REF_NAME + ".fa",
        bams = expand("mapping/{species}/{species}.{bam_suffix}", species = SPECIES_LIST, bam_suffix = BAM_SUFFIX),
        bais = expand("mapping/{species}/{species}.{bam_suffix}.bai", species = SPECIES_LIST, bam_suffix = BAM_SUFFIX)
    output:
        vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz"),
        tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz.tbi")
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
        '-o {output.vcf} '
        # '-o vcf/{SPECIES_ORDER}/{RUN_ID}/chr{wildcards.chr}.piece{wildcards.piece}.vcf.gz '
        '2> '
        # 'vcf/{SPECIES_ORDER}/{RUN_ID}/chr{wildcards.chr}.piece{wildcards.piece}.vcf.gz.log '
        '{output.vcf}.log'

# TODO Robbie Temp - delete
##        'REGION_START=$(({wildcards.piece} * 10000000 + 1)) && '
##        'REGION_END=$((({wildcards.piece} + 1) * 10000000)) && '
##        '-L {GATK_CHR_PREFIX}{wildcards.chr}:${{REGION_START}}-${{REGION_END}} '


rule filter_chr:
    input:
        input_vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz"),
        input_tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz.tbi"),
        ref = REF_DIR + REF_NAME + ".fa"
    output:
        filtered_vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.filtered.piece{{{{piece}}}}.vcf.gz"),
        filtered_tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.filtered.piece{{{{piece}}}}.vcf.gz.tbi")
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
        # 'vcf/logs/{VCF_PREFIX}.chr{chr}.filtered.piece{piece}.vcf.gz.log && '
        '{output.filtered_vcf}.log && '
        'rm {input.input_vcf} && rm {input.input_tbi} '


rule merge_chr:
    input:
        vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.filtered.piece{{piece}}.vcf.gz", chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.filtered.piece{{piece}}.vcf.gz.tbi", chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        ref = REF_DIR + REF_NAME + ".fa"	
    output:
        merged_vcf = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz",
        merged_tbi = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz.tbi"
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
        ref = REF_DIR + REF_NAME + ".fa",
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
        amb = REF_DIR + REF_NAME + ".fa.amb",
        ann = REF_DIR + REF_NAME + ".fa.ann"
    output:
        summary = REF_DIR + REF_NAME + ".fa.summary.txt"    
    params:
        N='prepare_reference',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
    shell:
        'R -f {R_GET_GENOME_STATS} --args {REF_DIR} {REF_NAME}.fa'


rule get_callable_regions:
    input:
        infile = expand("coverage/coverage.{{species}}.chr{chr}.txt.gz", chr = CHR_LIST_ONLY_AUTOS),
        summary = REF_DIR + REF_NAME + ".fa.summary.txt"
    output:
        beds = expand("coverage/coverage.{{species}}.chr{chr}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS),
	merged_bed = "coverage/coverage.{species}.callableOnly.bed"
    params:
        N='get_callable_regions',
        threads = 2,
        queue = "short.qc"
    wildcard_constraints:
        species='\w{1,40}'
    shell:
        'echo start && '
        'R -f {R_GET_PER_SAMPLE_AVERAGE_COV} --args {wildcards.species} {REF_DIR} {REF_NAME}.fa '
        '&& echo done'

CALLABLE_MIN_N=CALLABLE_MIN_FRACTION * len(SPECIES_LIST)

rule get_single_callable_regions:
    input:
        beds = expand("coverage/coverage.{species}.callableOnly.bed", species = SPECIES_LIST)
    output:
        bed = f"coverage/coverage.{SPECIES_ORDER}.all.callableOnly.bed"
    params:
        N='get_single_callable_regions',
        threads = 1,
        queue = "short.qc"
    shell:
        """
        multiIntersectBed -i {input.beds} | awk '{{if($4 >= ({CALLABLE_MIN_N})) {{print $1"\t"$2"\t"$3}}}}' > {output.bed}
        """



rule prepare_treemix:
    input:
        merged_vcf = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz",
        merged_tbi = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz.tbi"	
    output:
        treemix_input = f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.frq.gz"
    params:
        N='prepare_treemix',
        threads=1,
        queue = "short.qc@@short.hge" # TODO: change back to all nodes if can figure out 'bash illegal instruction...' error on `${{PYTHON_352}} ...`
    shell:
        'mkdir -p treemix && '
        'echo recode as treemix format && date && '
        '${{PYTHON_352}} {VCF2TREEMIX} {input.merged_vcf} {output.treemix_input} > {output.treemix_input}.log && echo done && date'
        
        
rule run_treemix:
    input:
        treemix_input = f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.frq.gz"
    output:
        treemix_output = f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz"
    params:
        N='run_treemix',
        threads=TREEMIX_THREADS,
        queue = "short.qc@@short.hge" # TODO: change back to all nodes if can figure out 'bash illegal instruction...' error on `${{TREEMIX}} ...`
    wildcard_constraints:
        migrants='\d{1,2}'
    shell:
        'echo begin treemix && date && '
        '${{TREEMIX}} -i {input} -m {wildcards.migrants} -noss -se -k {TREEMIX_K} '
        '-root {TREEMIX_OUTGROUP} '
        '-n_warn 10 '
        '-o treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{wildcards.migrants}.out &&'
        'echo done treemix && date && '
        'R -f {TREEMIX_R_PLOT} '
        '--args "./" {TREEMIX_R_PLOT_HELPER} '
        'treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{wildcards.migrants}.out '
        'treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.frq.gz  && echo done && date'


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

