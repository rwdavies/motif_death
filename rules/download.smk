rule download_all:
    input:
        expand("mapping/{species}/{units}_1.fastq.gz", zip, species = order_df["species"], units = order_df["units"]),
        expand("mapping/{species}/{units}_2.fastq.gz", zip, species = order_df["species"], units = order_df["units"])
        ref = REF_DIR + REF_NAME + ".fa.gz"
        
# TODO obsolete?
# rule all_standardized_phred:
#     input:
#         expand("mapping/{species}/{units}_{pen}.standardized_phred.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)

# TODO: combine these into one job with variable [1, 2] !
rule download_fastq_1:
    output:
        expand("mapping/{{species}}/{{units}}_1.fastq.gz")
    params:
        N='download_fastq',
        threads=1,
        # path=lambda wildcards: config["units"][wildcards.units][wildcards.pen],
        path = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "1"]
    wildcard_constraints:
        units='\D{1,8}\d{0,9}',
        pen='\d',
        queue="short.qc"
    shell:
        'mkdir -p mapping && cd mapping && '
        'mkdir -p {wildcards.species} && cd {wildcards.species} && '
        'wget {params.path}'

rule download_fastq_2:
    output:
        expand("mapping/{{species}}/{{units}}_2.fastq.gz")
    params:
        N='download_fastq',
        threads=1,
        # path=lambda wildcards: config["units"][wildcards.units][wildcards.pen],
        path = lambda wildcards: order_df.loc[(wildcards.species, wildcards.units), "2"]
    wildcard_constraints:
        units='\D{1,8}\d{0,9}',
        pen='\d',
        queue="short.qc"
    shell:
        'mkdir -p mapping && cd mapping && '
        'mkdir -p {wildcards.species} && cd {wildcards.species} && '
        'wget {params.path}'

rule download_ref:
    input:
    output:
        ref = REF_DIR + REF_NAME + ".fa.gz"
    params: N='make_ref', threads=1, queue = "short.qc"
    shell:
        'mkdir -p {REF_DIR} && cd {REF_DIR} && '
        'wget {REF_URL}'
        
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
