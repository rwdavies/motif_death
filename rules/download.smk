rule download_all:
    input:
        expand(f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_1.fastq.gz", zip, species = order_df["species"], units = order_df["units"]),
        expand(f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_2.fastq.gz", zip, species = order_df["species"], units = order_df["units"]),
        f"{EXTERNAL_DIR}/{REF_NAME}.simpleRepeat.gz",
        f"{EXTERNAL_DIR}/{REF_NAME}.rmsk.gz",
        ref = f"{REF_DIR}/{REF_NAME}.fa.gz"
        
# TODO obsolete?
# rule all_standardized_phred:
#     input:
#         expand("mapping/{species}/{units}_{pen}.standardized_phred.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)

# TODO: combine these into one job with variable [1, 2] !
rule download_fastq_1:
    output:
        f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_1.fastq.gz"
    params:
        N='download_fastq',
        threads=1,
        path = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "1"]
    wildcard_constraints:
        units='\D{1,8}\d{0,9}',
        pen='\d',
        queue="short.qc"
    shell:
        """
        mkdir -p mapping/{SPECIES_ORDER}
        cd mapping/{SPECIES_ORDER}
        mkdir -p {wildcards.species}
        cd {wildcards.species}
        wget {params.path}
        """

rule download_fastq_2:
    output:
        f"mapping/{SPECIES_ORDER}/{{species}}/{{units}}_2.fastq.gz"
    params:
        N='download_fastq',
        threads=1,
        path = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "2"]
    wildcard_constraints:
        units='\D{1,8}\d{0,9}',
        pen='\d',
        queue="short.qc"
    shell:
        """
        mkdir -p mapping/{SPECIES_ORDER}
        cd mapping/{SPECIES_ORDER}
        mkdir -p {wildcards.species}
        cd {wildcards.species}
        wget {params.path}
        """

rule download_ref:
    input:
    output:
        ref = f"{REF_DIR}/{REF_NAME}.fa.gz"
    params: N='make_ref', threads=1, queue = "short.qc"
    shell:
        """
        mkdir -p {REF_DIR}
        cd {REF_DIR}
        wget -O {REF_NAME}.fa.gz {REF_URL}
        """

rule download_rmask:
    # Note: Adds header if not present
    input:
    output:
        rmask = f"{EXTERNAL_DIR}/{REF_NAME}.rmsk.gz"
    params: N='download_rmask', threads=1, queue = "short.qc"
    shell:
        """
        mkdir -p {EXTERNAL_DIR}
        wget -O {output.rmask} {rmask_URL}
        if [ $(gunzip -c {output.rmask} | head -1 | cut -d$'\t' -f1) != '#bin' ]
        then
            echo Adding header to rmask file
            gunzip -c {output.rmask} > {EXTERNAL_DIR}/temp.rmsk
            echo -e {RMASK_HEADER} | cat - {EXTERNAL_DIR}/temp.rmsk > {EXTERNAL_DIR}/{REF_NAME}.rmsk
            gzip -f {EXTERNAL_DIR}/{REF_NAME}.rmsk
        fi
        rm {EXTERNAL_DIR}/temp.rmsk
        """

rule download_simple_repeat:
    # Note: Adds header if not present
    input:
    output:
        simpleRepeat = f"{EXTERNAL_DIR}/{REF_NAME}.simpleRepeat.gz"
    params:
        N='download_simple_repeat',
        threads=1,
        queue = "short.qc"
    shell:
        """
        mkdir -p {EXTERNAL_DIR} 
        wget -O {output.simpleRepeat} {simpleRepeat_URL} 
        if [ $(gunzip -c {output.simpleRepeat} | head -1 | cut -d$'\t' -f1) != '#bin' ]
        then
            echo Adding header to simpleRepeat file
            gunzip -c {output.simpleRepeat} > {EXTERNAL_DIR}/temp.simpleRepeat
            echo -e {SIMPLE_REPEAT_HEADER} | cat - {EXTERNAL_DIR}/temp.simpleRepeat > {EXTERNAL_DIR}/{REF_NAME}.simpleRepeat
            gzip -f {EXTERNAL_DIR}/{REF_NAME}.simpleRepeat
        fi
        rm {EXTERNAL_DIR}/temp.simpleRepeat
        """

# TODO Obsolete?
# rule standardize_fastq:
#     input:
#         fastq = expand("mapping/{species}/{{units}}_{{pen}}.fastq.gz", species = SPECIES)
#     output:
#         fastq = expand("mapping/{species}/{{units}}_{{pen}}.standardized_phred.fastq.gz", species = SPECIES)
#     params:
#         N='standardize_fastq',
#         threads=1,
#         path=lambda wildcards: config["units"][wildcards.units][wildcards.pen],
#         queue='short.qc'	
#     wildcard_constraints:
#         units='\D{1,3}\d{1,9}',
#         pen='\d'
#     shell:
#         'mkdir -p mapping && cd mapping && '
#         'mkdir -p {SPECIES} && cd {SPECIES} && '
#         'gzip -d {wildcards.units}_{wildcards.pen}.fastq.gz && '
#         '{FASTQSOL2PHRED} '
#         '{wildcards.units}_{wildcards.pen}.fastq '
#         '{wildcards.units}_{wildcards.pen}.standardized_phred.fastq && '
#         'gzip -1 {wildcards.units}_{wildcards.pen}.fastq && '
#         'gzip -1 {wildcards.units}_{wildcards.pen}.standardized_phred.fastq'

# TODO Obsolete?
# rule expand_fastq:
#     input:
#         expand("mapping/{species}/{units}_{pen}.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)
#     output:
#         expand("mapping/{species}/{{units}}_{{pen}}.fastq.gz", species = SPECIES)
#     params:
#         N='expand_fastq',
#         threads=1
#     wildcard_constraints:
#         units='\D{1,3}\d{1,9}',
#         pen='\d',
#         queue="short.qc"
#     shell:
#         'cd mapping && cd {SPECIES} && '
# 	'{PYTHON_352} SCRIPT ARGUMENTS'
