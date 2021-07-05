SPECIES = config["species"]
UNITS = config["units"].keys()
N_MAPPING_PIECES = config["n_mapping_pieces"]
N_MAPPING_PIECES_SE = N_MAPPING_PIECES * 10
MAPPING_PIECES = range(1, int(N_MAPPING_PIECES + 1), 1)
MAPPING_PIECES_SE = range(1, int(N_MAPPING_PIECES_SE + 1), 1)
MAPPING_QUEUE = config["mapping_queue"]
FASTQ_SUFFIX = "fastq.gz"
if SPECIES == "marmoset":
    FASTQ_SUFFIX = "standardized_phred.fastq.gz"

if OPERATE_GATK_PER_CHR == "FALSE":
    INDEL_REALIGNMENT_QUEUE="long.qc"
else:
    INDEL_REALIGNMENT_QUEUE="short.qc"

rule all:
    input:
        ["mapping/" + SPECIES + "/" + SPECIES + ".realigned.rmdup.bam"]

## make head this for faster version "| head -n 400"
rule extract_and_map_fastq_pieces:
    input:
        pen1 = expand("mapping/{species}/{{units}}_1.{fastq_suffix}", species = SPECIES, fastq_suffix = FASTQ_SUFFIX),
        pen2 = expand("mapping/{species}/{{units}}_2.{fastq_suffix}", species = SPECIES, fastq_suffix = FASTQ_SUFFIX)
    output:
        bam = expand("mapping/{species}/{{units}}_piece{{piece}}.bam", species = SPECIES),
        bai = expand("mapping/{species}/{{units}}_piece{{piece}}.bam.bai", species = SPECIES)
    params:
        N='map',
        threads=1,
        n_mapping_pieces=N_MAPPING_PIECES,
        queue = MAPPING_QUEUE,
        lb = lambda wildcards: config["units"][wildcards.units]["lb"],
        lb_insert_size = lambda wildcards: config["units"][wildcards.units]["lb_insert_size"],	
        flowcell_barcode = lambda wildcards: config["units"][wildcards.units]["flowcell_barcode"],
        flowcell_lane = lambda wildcards: config["units"][wildcards.units]["flowcell_lane"],
        head=" "
    wildcard_constraints:
        units='[A-Za-z0-9]+',
        piece='\d{1,3}'
    shell:
        'set +o pipefail && '
        'gunzip -c {input.pen1} {params.head} | awk \'{{ if( ((int((NR - 1) / 4) ) % {params.n_mapping_pieces}) == ({wildcards.piece} - 1)) {{print $0 }} }} \' | gzip -1 > {input.pen1}.temp.{wildcards.piece}.gz && '
        'gunzip -c {input.pen2} {params.head} | awk \'{{ if( ((int((NR - 1) / 4) ) % {params.n_mapping_pieces}) == ({wildcards.piece} - 1)) {{print $0 }} }} \' | gzip -1 > {input.pen2}.temp.{wildcards.piece}.gz && '
        'set -o pipefail && '
        '{BWA} mem -R "@RG\\tID:{params.flowcell_barcode}.{params.flowcell_lane}\\tSM:{SPECIES}\\tPL:ILLUMINA\\tPI:{params.lb_insert_size}\\tPU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES}\\tLB:{params.lb}" '
        ' {REF_DIR}/{REFNAME}.fa -t{params.threads} '
        '{input.pen1}.temp.{wildcards.piece}.gz '
        '{input.pen2}.temp.{wildcards.piece}.gz | '
        '{SAMTOOLS} view -Sb - > {output.bam}.temp1.bam && '
        '{PYTHON_278} {STAMPY} '
        '--readgroup='
        'ID:{params.flowcell_barcode}.{params.flowcell_lane},'
        'SM:{SPECIES},'
        'PL:ILLUMINA,'
        'PI:{params.lb_insert_size},'
        'PU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES},'
        'LB:{params.lb}'
        '-t{params.threads} -g {REF_DIR}/{REFNAME} -h {REF_DIR}/{REFNAME} --overwrite --bamkeepgoodreads '
        '-M {output.bam}.temp1.bam | '
        '{SAMTOOLS} view -Sb - > {output.bam}.temp2.bam && '
        'echo STAMPY DONE && '
        '{SAMTOOLS} sort -T mapping/{SPECIES}/ -o {output.bam} {output.bam}.temp2.bam && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {output.bam}.temp1.bam && '	
        'rm {output.bam}.temp2.bam && '
        'rm {input.pen1}.temp.{wildcards.piece}.gz {input.pen2}.temp.{wildcards.piece}.gz'


## after STAMPY command maybe --sensitive
## make head this for faster version "| head -n 400"

rule merge_mapped_pieces:
    input:
        bams = expand("mapping/{species}/{{units}}_piece{piece}.bam", species = SPECIES, piece = MAPPING_PIECES),
        bais = expand("mapping/{species}/{{units}}_piece{piece}.bam.bai", species = SPECIES, piece = MAPPING_PIECES)
    output:
        bam = expand("mapping/{species}/{{units}}.bam", species = SPECIES),
        bai = expand("mapping/{species}/{{units}}.bam.bai", species = SPECIES)
    params:
        N='merge_pieces',
        threads=1,
        queue = "long.qc"
    wildcard_constraints:
        units='[A-Za-z0-9]+'
    shell:
        '{SAMTOOLS} merge {output.bam} {input.bams} && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {input.bams} && '
        'rm {input.bais}'

 
rule merge_units:
    input:
        bams = expand("mapping/{species}/{units}.bam", species = SPECIES, units = UNITS),
        bais = expand("mapping/{species}/{units}.bam.bai", species = SPECIES, units = UNITS)
    output:
        bam = expand("mapping/{species}/{species}.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.bam.bai", species = SPECIES)
    params:
        N='merge_units',
        threads=1,
        queue = "long.qc"
    wildcard_constraints:
        units='[A-Za-z0-9]+'
    shell:
        '{SAMTOOLS} merge {output.bam} {input.bams} && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {input.bams} && '
        'rm {input.bais}'


# can be done at the end due to library usage
rule mark_duplicates:
    input:
        bam = expand("mapping/{species}/{species}.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.bam.bai", species = SPECIES)
    output:
        bam = expand("mapping/{species}/{species}.rmdup.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.rmdup.bam.bai", species = SPECIES)
    params:
        N='rmdup',
        threads=3,
        queue = "long.qc"
    shell:
        'ulimit -n 4096 && {JAVA} -Xmx16G -jar {PICARD} MarkDuplicates '
        'ASSUME_SORTED=false '
        'READ_NAME_REGEX=null '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        'METRICS_FILE={output.bam}.metrics.txt && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {input.bam} && '
        'rm {input.bai}'


rule identify_indels:
    input:
        bam = expand("mapping/{species}/{species}.rmdup.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.rmdup.bam.bai", species = SPECIES),
        ref = REF_DIR + REFNAME + ".fa"	
    output:
        outfile = expand("mapping/{species}/intervals/{species}.chr{{chr}}.intervals.list", species = SPECIES)
    params:
        N='indel_identify',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
    shell:
        'if [ \"{OPERATE_GATK_PER_CHR}\" == \"FALSE" ]; then '
        'minus_L=\" \" ; else '
        'minus_L=\"-L {GATK_CHR_PREFIX}{wildcards.chr}\" ; fi && '
        'mkdir -p mapping/{SPECIES}/intervals && '
        '{JAVA} -Xmx8g -jar {GATK} '
        '-T RealignerTargetCreator '
        '-R {input.ref} '
        '${{minus_L}} '
        '-I {input.bam} '
        '-o {output.outfile} '
        '2> {output.outfile}.log'


rule realign_around_indels:
    input:
        bam = expand("mapping/{species}/{species}.rmdup.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.rmdup.bam.bai", species = SPECIES),
        ref = REF_DIR + REFNAME + ".fa",
        intervals = expand("mapping/{species}/intervals/{species}.chr{{chr}}.intervals.list", species = SPECIES)	
    output:
        bam = expand("mapping/{species}/{species}.chr{{chr}}.realigned.rmdup.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.chr{{chr}}.realigned.rmdup.bai", species = SPECIES)
    params:
        N='apply_indel',
        threads=1,
        queue = INDEL_REALIGNMENT_QUEUE
    shell:
        'if [ \"{OPERATE_GATK_PER_CHR}\" == \"FALSE" ]; then '
        'minus_L=\" \" ; else '
        'minus_L=\"-L {GATK_CHR_PREFIX}{wildcards.chr}\" ; fi && '
        '{JAVA} -Xmx8g -jar {GATK} '
        '-T IndelRealigner '
        '-R {input.ref} '
        '-targetIntervals {input.intervals} '
        '-I {input.bam} '
        '${{minus_L}} '
        '-o {output.bam} '
        '2> {input.intervals}.realignment.log '

rule merge_realigned_bams:
    input:
        bams = expand("mapping/{species}/{species}.chr{chr}.realigned.rmdup.bam", species = SPECIES, chr = CHR_LIST),
        bais = expand("mapping/{species}/{species}.chr{chr}.realigned.rmdup.bai", species = SPECIES, chr = CHR_LIST)
    output:
        bam = expand("mapping/{species}/{species}.realigned.rmdup.bam", species = SPECIES),
        bai = expand("mapping/{species}/{species}.realigned.rmdup.bam.bai", species = SPECIES)
    params:
        N='merge_units',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
    shell:
        '{SAMTOOLS} merge {output.bam} {input.bams} && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {input.bams} && '
        'rm {input.bais} && '
        'rm mapping/{SPECIES}/{SPECIES}.rmdup.bam && '
        'rm mapping/{SPECIES}/{SPECIES}.rmdup.bam.bai '	








##
## single end stuff below
##
rule all_se:
    input:
        [expand("mapping/{species}/{units}.se.bam", species = SPECIES, units = UNITS)]

rule extract_and_map_single_fastq_pieces:
    input:
        pen = expand("mapping/{species}/{{units}}.{fastq_suffix}", species = SPECIES, fastq_suffix = FASTQ_SUFFIX)
    output:
        bam = expand("mapping/{species}/{{units}}.se_piece{{piece}}.bam", species = SPECIES),
        bai = expand("mapping/{species}/{{units}}.se_piece{{piece}}.bam.bai", species = SPECIES)
    params:
        N='map',
        threads=1,
        n_mapping_pieces=N_MAPPING_PIECES_SE,
        queue = MAPPING_QUEUE,
        lb = lambda wildcards: config["units"][wildcards.units]["lb"],
        lb_insert_size = lambda wildcards: config["units"][wildcards.units]["lb_insert_size"],	
        flowcell_barcode = lambda wildcards: config["units"][wildcards.units]["flowcell_barcode"],
        flowcell_lane = lambda wildcards: config["units"][wildcards.units]["flowcell_lane"],
        head=" "
    wildcard_constraints:
        units='[A-Za-z0-9]+',
        piece='\d{1,3}'
    shell:
        'set +o pipefail && '
        'gunzip -c {input.pen} {params.head} | awk \'{{ if( ((int((NR - 1) / 4) ) % {params.n_mapping_pieces}) == ({wildcards.piece} - 1)) {{print $0 }} }} \' | gzip -1 > {input.pen}.se.temp.{wildcards.piece}.gz && '
        'set -o pipefail && '
        'echo ====BWA START==== && '
        '{BWA} mem -R "@RG\\tID:{params.flowcell_barcode}.{params.flowcell_lane}\\tSM:{SPECIES}\\tPL:ILLUMINA\\tPI:{params.lb_insert_size}\\tPU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES}\\t@LB{params.lb}" '
        ' {REF_DIR}/{REFNAME}.fa -t{params.threads} '
        '{input.pen}.se.temp.{wildcards.piece}.gz | '
        '{SAMTOOLS} view -Sb - > {output.bam}.temp1.bam && '
        'rm {input.pen}.se.temp.{wildcards.piece}.gz && '
        'echo ======STAMPY START====== && '	
        '{PYTHON_278} {STAMPY} '
        '--readgroup='
        'ID:{params.flowcell_barcode}.{params.flowcell_lane},'
        'SM:{SPECIES},'
        'PL:ILLUMINA,'
        'PI:{params.lb_insert_size},'
        'PU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES},'
        'LB:{params.lb}'
        '-t{params.threads} -g {REF_DIR}/{REFNAME} -h {REF_DIR}/{REFNAME} --overwrite --bamkeepgoodreads '
        '-M {output.bam}.temp1.bam | '
        '{SAMTOOLS} view -Sb - > {output.bam}.temp2.bam && '
        'rm {output.bam}.temp1.bam && '	
        'echo =====STAMPY DONE===== && '
        '{SAMTOOLS} sort -T mapping/{SPECIES}/ -o {output.bam} {output.bam}.temp2.bam && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {output.bam}.temp2.bam '

rule merge_se_pieces:
    input:
        bams = expand("mapping/{species}/{{units}}.se_piece{piece}.bam", species = SPECIES, piece = MAPPING_PIECES_SE),
        bais = expand("mapping/{species}/{{units}}.se_piece{piece}.bam.bai", species = SPECIES, piece = MAPPING_PIECES_SE)
    output:
        bam = expand("mapping/{species}/{{units}}.se.bam", species = SPECIES),
        bai = expand("mapping/{species}/{{units}}.se.bam.bai", species = SPECIES)
    params:
        N='merge_se_pieces',
        threads=1,
        queue = "long.qc"
    wildcard_constraints:
        units='[A-Za-z0-9]+'
    shell:
        '{SAMTOOLS} merge {output.bam} {input.bams} && '
        '{SAMTOOLS} index {output.bam} && '
        'rm {input.bams} && '
        'rm {input.bais}'
