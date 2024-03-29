rule downstream_all:
    input:
        f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz",
        f"coverage/{SPECIES_ORDER}/coverage.{SPECIES_ORDER}.{RUN_ID}.all.callableOnly.bed",
        expand(f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz", migrants = TREEMIX_MIGRANT_RANGE)

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
        ref = f"{REF_DIR}/{REF_NAME}.fa",
        bams = expand(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}", species = SPECIES_LIST),
        bais = expand(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}.bai", species = SPECIES_LIST)
    output:
        vcf = f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.piece{{piece}}.vcf.gz",
        tbi = f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.piece{{piece}}.vcf.gz.tbi"
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
        '2> '
        '{output.vcf}.log'

# TODO Robbie Temp - delete
##        'REGION_START=$(({wildcards.piece} * 10000000 + 1)) && '
##        'REGION_END=$((({wildcards.piece} + 1) * 10000000)) && '
##        '-L {GATK_CHR_PREFIX}{wildcards.chr}:${{REGION_START}}-${{REGION_END}} '


rule filter_chr:
    input:
        input_vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz"),
        input_tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.piece{{{{piece}}}}.vcf.gz.tbi"),
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        filtered_vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.filtered.piece{{{{piece}}}}.vcf.gz"),
        filtered_tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{{{chr}}}}.filtered.piece{{{{piece}}}}.vcf.gz.tbi")
    params: N='filter_chr', threads=1, queue = "short"
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
        '{output.filtered_vcf}.log && '
        'rm {input.input_vcf} && rm {input.input_tbi} '


rule merge_chr:
    input:
        vcf = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.filtered.piece{{piece}}.vcf.gz", chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        tbi = expand(f"vcf/{SPECIES_ORDER}/{RUN_ID}/chr{{chr}}.filtered.piece{{piece}}.vcf.gz.tbi", chr = CHR_LIST_ONLY_AUTOS, piece = CHR_CHUNKS),
        ref = f"{REF_DIR}/{REF_NAME}.fa"	
    output:
        merged_vcf = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz",
        merged_tbi = f"vcf/{SPECIES_ORDER}/{RUN_ID}/filtered.vcf.gz.tbi"
    params: N='merge_chr', threads=1, queue = "short"
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
        ref = f"{REF_DIR}/{REF_NAME}.fa",
        bam = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}",
        bais = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}.bai"
    output:
        outfile = f"coverage/{SPECIES_ORDER}/coverage.{{species}}.chr{{chr}}.txt.gz"
    params:
        N='get_depth_of_coverage',
        threads=1,
        queue = "short"
    wildcard_constraints:
        species='\w{1,40}'	
    shell:
        'mkdir -p coverage && '
        '${{JAVA}} -Xmx14g -jar ${{GATK}} '
        '-T DepthOfCoverage '
        '-R {input.ref} '
        '-L {GATK_CHR_PREFIX}{wildcards.chr} '
        '--minMappingQuality 17 --minBaseQuality 0 '
        '-I {input.bam} '
        '--downsample_to_coverage 1000 '
        '-o coverage/{SPECIES_ORDER}/coverage.{wildcards.species}.chr{wildcards.chr}.txt '
        '> {output.outfile}.log && '
        'gzip -1 coverage/{SPECIES_ORDER}/coverage.{wildcards.species}.chr{wildcards.chr}.txt'


rule prepare_reference:
    input:
        amb = f"{REF_DIR}/{REF_NAME}.fa.amb",
        ann = f"{REF_DIR}/{REF_NAME}.fa.ann"
    output:
        summary = f"{REF_DIR}/{REF_NAME}.fa.summary.txt"
    params:
        N='prepare_reference',
        threads=1,
        queue = "short"
    wildcard_constraints:
    shell:
        """
        {R_GET_GENOME_STATS} --ref_dir={REF_DIR} --ref={REF_NAME}.fa --chr_prefix={GATK_CHR_PREFIX} {CHR_LIST_ONLY_AUTOS}
        """


rule get_callable_regions:
    input:
        infile = expand(f"coverage/{SPECIES_ORDER}/coverage.{{{{species}}}}.chr{{chr}}.txt.gz", chr = CHR_LIST_ONLY_AUTOS),
        summary = f"{REF_DIR}/{REF_NAME}.fa.summary.txt"
    output:
        beds = expand(f"coverage/{SPECIES_ORDER}/coverage.{{{{species}}}}.chr{{chr}}.callableOnly.bed", chr = CHR_LIST_ONLY_AUTOS),
        merged_bed = f"coverage/{SPECIES_ORDER}/coverage.{{species}}.callableOnly.bed",
        average = f"coverage/{SPECIES_ORDER}/average.{{species}}.txt"
    params:
        N='get_callable_regions',
        threads = 3,
        queue = "short"
    wildcard_constraints:
        species='\w{1,40}'
    shell:
        """
        echo start get_callable_regions
        {R_GET_PER_SAMPLE_AVERAGE_COV} --order={SPECIES_ORDER} --species={wildcards.species} --ref_dir={REF_DIR} --ref={REF_NAME}.fa --chr_prefix={GATK_CHR_PREFIX} {CHR_LIST_ONLY_AUTOS}
        echo done get_callable_regions
        """

CALLABLE_MIN_N=CALLABLE_MIN_FRACTION * len(SPECIES_LIST)

EASY_SPECIES_LIST=""
for species in SPECIES_LIST:
    EASY_SPECIES_LIST = EASY_SPECIES_LIST + species + " "


rule get_single_callable_regions:
    input:
        beds = expand(f"coverage/{SPECIES_ORDER}/coverage.{{species}}.callableOnly.bed", species = SPECIES_LIST),
        avs = expand(f"coverage/{SPECIES_ORDER}/average.{{species}}.txt", species = SPECIES_LIST)
    output:
        bed = f"coverage/{SPECIES_ORDER}/coverage.{SPECIES_ORDER}.{RUN_ID}.all.callableOnly.bed",
        av = f"coverage/{SPECIES_ORDER}/coverage.{SPECIES_ORDER}.{RUN_ID}.all.average.txt",
        checkfile = f"coverage/{SPECIES_ORDER}/coverage.{SPECIES_ORDER}.{RUN_ID}.all.callableOnly.bed.clean.txt"
    params:
        N='get_single_callable_regions',
        threads = 1,
        queue = "short"
    shell:
        """
        multiIntersectBed -i {input.beds} | awk '{{if($4 >= ({CALLABLE_MIN_N})) {{print $1"\t"$2"\t"$3}}}}' > {output.bed}
        if [ ! -s {output.bed} ] # raise error if output is empty
        then
            echo "Error: file {output.bed} is empty"
            exit 1
        fi
	rm -f {output.av}
	for sps in {EASY_SPECIES_LIST}
	do
            cov=`cat coverage/{SPECIES_ORDER}/average.${{sps}}.txt`
	    echo coverage is ${{cov}}
	    echo -e "${{sps}}\t${{cov}}" >> {output.av}
        done
	echo remove temp files callable regions
        rm -f coverage/{SPECIES_ORDER}/coverage.*.chr*.txt*
	rm -f coverage/{SPECIES_ORDER}/coverage*.chr*.callableOnly.bed
	for sps in {EASY_SPECIES_LIST}
	do
            rm -f coverage/{SPECIES_ORDER}/coverage.${{sps}}.chr*.callableOnly.bed
            rm -f coverage/{SPECIES_ORDER}/coverage.${{sps}}.callableOnly.bed	    
        done
	touch {output.checkfile}
	echo done cleanup
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
        queue = "short"
    shell:
        'mkdir -p treemix && '
        'echo recode as treemix format && date && '
        'echo ${{PYTHON_352}} {VCF2TREEMIX} {input.merged_vcf} {output.treemix_input} && '
        '${{PYTHON_352}} {VCF2TREEMIX} {input.merged_vcf} {output.treemix_input} > {output.treemix_input}.log && echo done && date'
        
        
rule run_treemix:
    input:
        treemix_input = f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.frq.gz"
    output:
        treemix_output = f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz"
    params:
        N='run_treemix',
        threads=TREEMIX_THREADS,
        queue = "short"
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

