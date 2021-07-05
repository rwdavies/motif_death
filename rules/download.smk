SPECIES = config["species"]
UNITS = config["units"].keys()
N_MAPPING_PIECES = config["n_mapping_pieces"]
MAPPING_PIECES = range(1, int(N_MAPPING_PIECES + 1), 1)
MAPPING_QUEUE = config["mapping_queue"]

rule all:
    input:
        expand("mapping/{species}/{units}_{pen}.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)
        
# TODO obsolete?
rule all_standardized_phred:
    input:
        expand("mapping/{species}/{units}_{pen}.standardized_phred.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)

rule download_fastq:
    output:
        expand("mapping/{species}/{{units}}_{{pen}}.fastq.gz", species = SPECIES)
    params:
        N='download_fastq',
        threads=1,
        path=lambda wildcards: config["units"][wildcards.units][wildcards.pen],
    wildcard_constraints:
        units='\D{1,3}\d{1,9}',
        pen='\d',
        queue="short.qc"
    shell:
        'mkdir -p mapping && cd mapping && '
        'mkdir -p {SPECIES} && cd {SPECIES} && '
        'wget {params.path}'

# TODO Obsolete?
rule standardize_fastq:
    input:
        fastq = expand("mapping/{species}/{{units}}_{{pen}}.fastq.gz", species = SPECIES)
    output:
        fastq = expand("mapping/{species}/{{units}}_{{pen}}.standardized_phred.fastq.gz", species = SPECIES)
    params:
        N='standardize_fastq',
        threads=1,
        path=lambda wildcards: config["units"][wildcards.units][wildcards.pen],
        queue='short.qc'	
    wildcard_constraints:
        units='\D{1,3}\d{1,9}',
        pen='\d'
    shell:
        'mkdir -p mapping && cd mapping && '
        'mkdir -p {SPECIES} && cd {SPECIES} && '
        'gzip -d {wildcards.units}_{wildcards.pen}.fastq.gz && '
        '{FASTQSOL2PHRED} '
        '{wildcards.units}_{wildcards.pen}.fastq '
        '{wildcards.units}_{wildcards.pen}.standardized_phred.fastq && '
        'gzip -1 {wildcards.units}_{wildcards.pen}.fastq && '
        'gzip -1 {wildcards.units}_{wildcards.pen}.standardized_phred.fastq'

# TODO Obsolete?
rule expand_fastq:
    input:
        expand("mapping/{species}/{units}_{pen}.fastq.gz", species = SPECIES, pen = [1, 2], units = UNITS)
    output:
        expand("mapping/{species}/{{units}}_{{pen}}.fastq.gz", species = SPECIES)
    params:
        N='expand_fastq',
        threads=1
    wildcard_constraints:
        units='\D{1,3}\d{1,9}',
        pen='\d',
        queue="short.qc"
    shell:
        'cd mapping && cd {SPECIES} && '
	'{PYTHON_352} SCRIPT ARGUMENTS'