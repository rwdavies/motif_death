# TODO: figure out if below needs to be set as species level?
# if SPECIES == "marmoset":
#     FASTQ_SUFFIX = "standardized_phred.fastq.gz"

# TODO: replace with one function
def get_bam_pieces(wildcards):
    unit_dir = checkpoints.chunk_fastq.get(**wildcards).output.unit_dir
    pieces = glob_wildcards(f"{unit_dir}/1.{FASTQ_SUFFIX}.temp.{{piece}}.gz").piece
    filenames = expand(rules.map_fastq_pieces.output.bam, **wildcards, piece=pieces)
    return filenames

def get_bai_pieces(wildcards):
    unit_dir = checkpoints.chunk_fastq.get(**wildcards).output.unit_dir
    pieces = glob_wildcards(f"{unit_dir}/1.{FASTQ_SUFFIX}.temp.{{piece}}.gz").piece
    filenames = expand(rules.map_fastq_pieces.output.bai, **wildcards, piece=pieces)
    return filenames

# TODO: replace with one function
def get_bam_units(wildcards):
    units = list(order_df.loc[wildcards.species, "units"])
    filenames = [f"mapping/{SPECIES_ORDER}/{wildcards.species}/{wildcards.species}.unit{u}.bam" for u in units]
    return filenames

def get_bai_units(wildcards):
    units = list(order_df.loc[wildcards.species, "units"])
    filenames = [f"mapping/{SPECIES_ORDER}/{wildcards.species}/{wildcards.species}.unit{u}.bam.bai" for u in units]
    return filenames


rule map_all:
    input:
        expand(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.realigned.rmdup.bam", species = order_df["species"])




checkpoint chunk_fastq:
    input:
        pen1 = f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_1.{FASTQ_SUFFIX}",
        pen2 = f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_2.{FASTQ_SUFFIX}",
    output:
        unit_dir = directory(f"mapping/{SPECIES_ORDER}/{{species}}/temp/{{units}}")
    params:
        N='chunk_fastq',
        threads=1,
        queue = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "mapping_queue"],
        n_mapping_pieces= lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "n_mapping_pieces"]
    wildcard_constraints:
        units=WILDCARD_UNIT_CONSTRAINT
    shell:
        """
        mkdir -p {output.unit_dir}
        # Note: bash math below. division returns integer only. must be multiple of 4 as fasta is in 4 line pieces.
        lines_per_chunk=$(($(gunzip -c {input.pen1} | wc -l) / {params.n_mapping_pieces}))
        lines_per_chunk=$(($lines_per_chunk - $(($lines_per_chunk % 4))))
        zcat {input.pen1} | split -l ${{lines_per_chunk}} --filter=\'gzip > $FILE.gz\' - {output.unit_dir}/1.{FASTQ_SUFFIX}.temp.
        zcat {input.pen2} | split -l ${{lines_per_chunk}} --filter=\'gzip > $FILE.gz\' - {output.unit_dir}/2.{FASTQ_SUFFIX}.temp.
        """


## make head this for faster version "| head -n 400"
rule map_fastq_pieces:
    input:
        pen1_chunk = f"mapping/{SPECIES_ORDER}/{{species}}/temp/{{units}}/1.{FASTQ_SUFFIX}.temp.{{piece}}.gz",
        pen2_chunk = f"mapping/{SPECIES_ORDER}/{{species}}/temp/{{units}}/2.{FASTQ_SUFFIX}.temp.{{piece}}.gz",
        ref = f"{REF_DIR}/{REF_NAME}.fa",
        ref_sa = f"{REF_DIR}/{REF_NAME}.fa.sa",
        ref_fai = f"{REF_DIR}/{REF_NAME}.fa.fai",
        ref_dict = f"{REF_DIR}/{REF_NAME}.dict",
        ref_stidx = f"{REF_DIR}/{REF_NAME}.stidx"
    output:
        bam = temp(f"mapping/{SPECIES_ORDER}/{{species}}/temp/{{units}}/piece.{{piece}}.bam"),
        bai = temp(f"mapping/{SPECIES_ORDER}/{{species}}/temp/{{units}}/piece.{{piece}}.bam.bai")
    params:
        N='map',
        threads=1,
        queue = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "mapping_queue"],
        lb = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "lb"],
        lb_insert_size = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "lb_insert_size"],
        flowcell_barcode = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "flowcell_barcode"],
        flowcell_lane = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "flowcell_lane"]
    wildcard_constraints:
        units=WILDCARD_UNIT_CONSTRAINT,
        piece='\D{1,3}'
    shell:
        'bwa mem -R "@RG\\tID:{params.flowcell_barcode}.{params.flowcell_lane}\\tSM:{wildcards.species}\\tPL:ILLUMINA\\tPI:{params.lb_insert_size}\\tPU:{params.flowcell_barcode}.{params.flowcell_lane}.{wildcards.species}\\tLB:{params.lb}" '
        ' {input.ref} -t{params.threads} '
        '{input.pen1_chunk} '
        '{input.pen2_chunk} | '
        'samtools view -Sb - > {output.bam}.temp1.bam && '
        '${{PYTHON_278}} ${{STAMPY}} '
        '--readgroup='
        'ID:{params.flowcell_barcode}.{params.flowcell_lane},'
        'SM:{wildcards.species},'
        'PL:ILLUMINA,'
        'PI:{params.lb_insert_size},'
        'PU:{params.flowcell_barcode}.{params.flowcell_lane}.{wildcards.species},'
        'LB:{params.lb}'
        '-t{params.threads} -g {REF_DIR}/{REF_NAME} -h {REF_DIR}/{REF_NAME} --overwrite --bamkeepgoodreads '
        '-M {output.bam}.temp1.bam | '
        'samtools view -Sb - > {output.bam}.temp2.bam && '
	'rm {output.bam}.temp1.bam && '
        'echo STAMPY DONE && '
        'samtools sort -T mapping/{SPECIES_ORDER}/{wildcards.species}/temp/{wildcards.units}/ -o {output.bam} {output.bam}.temp2.bam && '
        'samtools index {output.bam} && '
        'rm {output.bam}.temp2.bam && '
        'rm {input.pen1_chunk} && '
        'rm {input.pen2_chunk} '

## after STAMPY command maybe --sensitive
## make head this for faster version "| head -n 400"
rule merge_mapped_pieces:
    input:
        bams = get_bam_pieces,
        bais = get_bai_pieces
    output:
        bam = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.unit{{units}}.bam"),
        bai = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.unit{{units}}.bam.bai")
    params:
        N='merge_pieces',
        threads=1,
        queue = "short.qc@@short.hge"
    wildcard_constraints:
        units=WILDCARD_UNIT_CONSTRAINT
    shell:
        """
        samtools merge {output.bam} {input.bams}
        samtools index {output.bam}
        rm {input.bams}
        """


rule merge_units:
    input:
        bams = get_bam_units,
        bais = get_bai_units
    output:
        bam = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.bam"),
        bai = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.bam.bai")
    params:
        N='merge_units',
        threads=1,
        queue = "short.qc@@short.hge"
    wildcard_constraints:
        units=WILDCARD_UNIT_CONSTRAINT
    shell:
        """
        samtools merge {output.bam} {input.bams}
        samtools index {output.bam}
        # TODO: below is hack as snakemake temp(directory(...)) not working
        rm -rf mapping/{SPECIES_ORDER}/{wildcards.species}/temp
        """


# can be done at the end due to library usage
rule mark_duplicates:
    input:
        bam = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.bam",
        bai = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.bam.bai"
    output:
        bam = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam"),
        bai = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam.bai")
    params:
        N='rmdup',
        threads=3,
        queue = "short.qc@@short.hge"
    shell:
        'ulimit -n 4096 && ${{JAVA}} -Xmx16G -jar ${{PICARD}} MarkDuplicates '
        'ASSUME_SORTED=false '
        'READ_NAME_REGEX=null '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
        'METRICS_FILE={output.bam}.metrics.txt && '
        'samtools index {output.bam}'


rule identify_indels:
    input:
        bam = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam",
        bai = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam.bai",
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        outfile = f"mapping/{SPECIES_ORDER}/{{species}}/intervals/{{species}}.chr{{chr}}.intervals.list"
    params:
        N='indel_identify',
        threads=1,
        queue = "short.qc"
    wildcard_constraints:
    shell:
        'if [ \"{OPERATE_GATK_PER_CHR}\" == \"FALSE" ]; then '
        'minus_L=\" \" ; else '
        'minus_L=\"-L {GATK_CHR_PREFIX}{wildcards.chr}\" ; fi && '
        'mkdir -p mapping/{SPECIES_ORDER}/{wildcards.species}/intervals && '
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T RealignerTargetCreator '
        '-R {input.ref} '
        '${{minus_L}} '
        '-I {input.bam} '
        '-o {output.outfile} '
        '2> {output.outfile}.log'


rule realign_around_indels:
    input:
        bam = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam",
        bai = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.rmdup.bam.bai",
        ref = f"{REF_DIR}/{REF_NAME}.fa",
        intervals = f"mapping/{SPECIES_ORDER}/{{species}}/intervals/{{species}}.chr{{chr}}.intervals.list"
    output:
        bam = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.chr{{chr}}.realigned.rmdup.bam"),
        bai = temp(f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.chr{{chr}}.realigned.rmdup.bai")
    params:
        N='apply_indel',
        threads=1,
        queue = INDEL_REALIGNMENT_QUEUE
    shell:
        'if [ \"{OPERATE_GATK_PER_CHR}\" == \"FALSE" ]; then '
        'minus_L=\" \" ; else '
        'minus_L=\"-L {GATK_CHR_PREFIX}{wildcards.chr}\" ; fi && '
        '${{JAVA}} -Xmx8g -jar ${{GATK}} '
        '-T IndelRealigner '
        '-R {input.ref} '
        '-targetIntervals {input.intervals} '
        '-I {input.bam} '
        '${{minus_L}} '
        '-o {output.bam} '
        '2> {input.intervals}.realignment.log '

rule merge_realigned_bams:
    input:
        bams = expand(f"mapping/{SPECIES_ORDER}/{{{{species}}}}/{{{{species}}}}.chr{{chr}}.realigned.rmdup.bam", chr = CHR_LIST),
        bais = expand(f"mapping/{SPECIES_ORDER}/{{{{species}}}}/{{{{species}}}}.chr{{chr}}.realigned.rmdup.bai", chr = CHR_LIST)
    output:
        bam = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}",
        bai = f"mapping/{SPECIES_ORDER}/{{species}}/{{species}}.{BAM_SUFFIX}.bai"
    params:
        N='merge_bams',
        threads=1,
        queue = "short.qc@@short.hge"
    wildcard_constraints:
    shell:
        """
        samtools merge {output.bam} {input.bams}
        samtools index {output.bam}
        rm mapping/{SPECIES_ORDER}/{wildcards.species}/*{FASTQ_SUFFIX}
        """








##
## single end stuff below
##
# rule all_se:
#     input:
#         [expand("mapping/{species}/{units}.se.bam", species = SPECIES, units = UNITS)]

# rule extract_and_map_single_fastq_pieces:
#     input:
#         pen = expand("mapping/{species}/{{units}}.{fastq_suffix}", species = SPECIES, fastq_suffix = FASTQ_SUFFIX)
#     output:
#         bam = expand("mapping/{species}/{{units}}.se_piece{{piece}}.bam", species = SPECIES),
#         bai = expand("mapping/{species}/{{units}}.se_piece{{piece}}.bam.bai", species = SPECIES)
#     params:
#         N='map',
#         threads=1,
#         n_mapping_pieces=N_MAPPING_PIECES_SE,
#         queue = MAPPING_QUEUE,
#         lb = lambda wildcards: config["units"][wildcards.units]["lb"],
#         lb_insert_size = lambda wildcards: config["units"][wildcards.units]["lb_insert_size"],	
#         flowcell_barcode = lambda wildcards: config["units"][wildcards.units]["flowcell_barcode"],
#         flowcell_lane = lambda wildcards: config["units"][wildcards.units]["flowcell_lane"],
#         head=" "
#     wildcard_constraints:
#         units='[A-Za-z0-9]+',
#         piece='\d{1,3}'
#     shell:
#         'set +o pipefail && '
#         'gunzip -c {input.pen} {params.head} | awk \'{{ if( ((int((NR - 1) / 4) ) % {params.n_mapping_pieces}) == ({wildcards.piece} - 1)) {{print $0 }} }} \' | gzip -1 > {input.pen}.se.temp.{wildcards.piece}.gz && '
#         'set -o pipefail && '
#         'echo ====BWA START==== && '
#         '{BWA} mem -R "@RG\\tID:{params.flowcell_barcode}.{params.flowcell_lane}\\tSM:{SPECIES}\\tPL:ILLUMINA\\tPI:{params.lb_insert_size}\\tPU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES}\\t@LB{params.lb}" '
#         ' {REF_DIR}/{REF_NAME}.fa -t{params.threads} '
#         '{input.pen}.se.temp.{wildcards.piece}.gz | '
#         '{SAMTOOLS} view -Sb - > {output.bam}.temp1.bam && '
#         'rm {input.pen}.se.temp.{wildcards.piece}.gz && '
#         'echo ======STAMPY START====== && '	
#         '{PYTHON_278} {STAMPY} '
#         '--readgroup='
#         'ID:{params.flowcell_barcode}.{params.flowcell_lane},'
#         'SM:{SPECIES},'
#         'PL:ILLUMINA,'
#         'PI:{params.lb_insert_size},'
#         'PU:{params.flowcell_barcode}.{params.flowcell_lane}.{SPECIES},'
#         'LB:{params.lb}'
#         '-t{params.threads} -g {REF_DIR}/{REF_NAME} -h {REF_DIR}/{REF_NAME} --overwrite --bamkeepgoodreads '
#         '-M {output.bam}.temp1.bam | '
#         '{SAMTOOLS} view -Sb - > {output.bam}.temp2.bam && '
#         'rm {output.bam}.temp1.bam && '	
#         'echo =====STAMPY DONE===== && '
#         '{SAMTOOLS} sort -T mapping/{SPECIES}/ -o {output.bam} {output.bam}.temp2.bam && '
#         '{SAMTOOLS} index {output.bam} && '
#         'rm {output.bam}.temp2.bam '

# rule merge_se_pieces:
#     input:
#         bams = expand("mapping/{species}/{{units}}.se_piece{piece}.bam", species = SPECIES, piece = MAPPING_PIECES_SE),
#         bais = expand("mapping/{species}/{{units}}.se_piece{piece}.bam.bai", species = SPECIES, piece = MAPPING_PIECES_SE)
#     output:
#         bam = expand("mapping/{species}/{{units}}.se.bam", species = SPECIES),
#         bai = expand("mapping/{species}/{{units}}.se.bam.bai", species = SPECIES)
#     params:
#         N='merge_se_pieces',
#         threads=1,
#         queue = "long.qc"
#     wildcard_constraints:
#         units='[A-Za-z0-9]+'
#     shell:
#         '{SAMTOOLS} merge {output.bam} {input.bams} && '
#         '{SAMTOOLS} index {output.bam} && '
#         'rm {input.bams} && '
#         'rm {input.bais}'
