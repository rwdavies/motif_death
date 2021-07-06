VCF_PREFIX = config["VCF_PREFIX"]
HATBAG_OUTPUT_DIR = config["HATBAG_OUTPUT_DIR"]
HATBAG_OUTPUT_DATE=config["HATBAG_OUTPUT_DATE"]
R_DIR=config["R_DIR"]
SPECIES_ORDER=config["SPECIES_ORDER"]
HATBAG_DIR=config["HATBAG_DIR"]

##
## current hack
##
rule HATBAG_HACK:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/F_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)


rule HATBAG_HACK_A:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/A_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
rule HATBAG_HACK_B:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/B_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
rule HATBAG_HACK_C:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/C_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
rule HATBAG_HACK_D:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/D_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
rule HATBAG_HACK_E:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/E_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
rule HATBAG_HACK_F:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/F_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)



##
## note, not sure if these should be changed
## if the cluster is busy, this will take a while
##
def get_hatbag_n_threads(wildcards):
    run = wildcards.run
    if run == "A":
        return(8)
    if run == "B":
        return(8)
    if run == "C":
        return(8)
    if run == "D":
        return(4)
    if run == "E":
        return(6)
    if run == "F":
        return(1)


def get_input_bam(wildcards):
    return(scenarios[wildcards.scenario]["impute_bam_list"][wildcards.sample_name])

def get_hatbag_input(wildcards):
    run = wildcards.run
    if run == "A":
        return()
    if run == "B":
        need = "A"
    if run == "C":
        need = "B"
    if run == "D":
        need = "C"
    if run == "E":
        need = "D"
    if run == "F":
        need = "E"
    return(expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/" + need + "_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE))


rule HATBAG_HACK_FUNCTION:
    input:
        get_hatbag_input
    output:
        decoy = expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/{{run}}_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
    params:
        N='hatbag_test',
        threads=get_hatbag_n_threads,
        queue="short.qc@@short.hge" # TODO: change back if can figure out problem on other nodes
    wildcard_constraints:
        run='[A-Z]{1,6}'
    shell:
        'R -f ' + R_DIR + 'run_all.R '
	'--args '
	' ' + R_DIR + ' '
        ' ' + SPECIES_ORDER + ' '
	' {wildcards.run} '
	' ' + HATBAG_OUTPUT_DATE + ' '
	' {params.threads} '
	' vcf/' + VCF_PREFIX + '.filtered.vcf.gz '
        ' ' + HATBAG_DIR + ' '
	' && touch {output.decoy} '








##
##
## UNCLEAR OPERATIONAL STATUS BELOW
##
## 

rule all_hatbag:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)


rule test_hatbag:
    input:
        expand("hatbag/{hatbag_output_dir}--{simple_repeat}_simple_repeat--{use_one_sided_pvalue}_use_one_sided_pvalue--mrle_{mrle}--ndge_{ndge}--callable_bed_{callable_bed}/{hatbag_output_date}/complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE, simple_repeat = ["yes", "no"], use_one_sided_pvalue = ["TRUE", "FALSE"], mrle = [6, 10], ndge = [3, 0], callable_bed = ["yes", "no"])


## simple - basic options
rule HATBAG_test:
    input:
        ref = "ref/" + REFNAME + ".fa",
        vcf = expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX),
        tbi = expand("vcf/{vcf_prefix}.filtered.vcf.gz.tbi", vcf_prefix = VCF_PREFIX)
    output:
        decoy = expand("hatbag/{hatbag_output_dir}--{{simple_repeat}}_simple_repeat--{{use_one_sided_pvalue}}_use_one_sided_pvalue--mrle_{{mrle}}--ndge_{{ndge}}--callable_bed_{{callable_bed}}/{hatbag_output_date}/complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
    params:
        N='hatbag_test',
        threads=8,
        queue = "short.qc",
        repeat_filename = lambda wildcards: config["simple_repeat"][wildcards.simple_repeat],
        callable_bed_filename = lambda wildcards: config["callable_bed"][wildcards.callable_bed]
    wildcard_constraints:
        piece='\d{1,3}',
        simple_repeat='[a-z]{1,5}',
        use_one_sided_pvalue='[A-Za-z]{1,5}',
	callable_bed='[a-z]{1,5}',	
        mrle='\d{1,2}',
        ndge='\d{1,1}'
    shell:
        'export PATH=/apps/well/htslib/1.4.1/bin/:${{PATH}} && '
        'which tabix && '
        'hostname && '
        'echo SGE_TASK_ID=${{SGE_TASK_ID}} && '
        'output_dir=hatbag/{HATBAG_OUTPUT_DIR}--{wildcards.simple_repeat}_simple_repeat--{wildcards.use_one_sided_pvalue}_use_one_sided_pvalue--mrle_{wildcards.mrle}--ndge_{wildcards.ndge}--callable_bed_{wildcards.callable_bed}/ && '
        'mkdir -p ${{output_dir}} && '
        '{R_HATBAG} '
        '--species={HATBAG_OUTPUT_DIR} '
        '--run=ABCDEF '
        '--nCores={params.threads} '
        '--outputDate={HATBAG_OUTPUT_DATE} '
        '--outputDir=${{output_dir}} '
        '--vcf_file={input.vcf} '
        '--reference={input.ref} '
        '--rmask_file={HATBAG_RMASK} '
        '--Klist=10 '
        '--chrlist="{HATBAG_CHRS}" '
        '--genomeSize={HATBAG_GENOME_SIZE} '
        '--lineages=\'{HATBAG_LINEAGES}\' '
        '--lineages_to_build=\'{HATBAG_LINEAGES_TO_BUILD}\' '	
        '--ancestral_lineage=\'{HATBAG_ANCESTRAL_LINEAGE}\' '
        '--outgroups=\'{HATBAG_OUTGROUPS}\' '
        '--simpleRepeat_file={params.repeat_filename} '
        '--use_one_sided_pvalue={wildcards.use_one_sided_pvalue} '
        '--callable_bed={params.callable_bed_filename} '
        '--mrle={wildcards.mrle} '
        '--ndge={wildcards.ndge} '	
        '&> '
        '{output.decoy}.log && '
        'touch {output.decoy}'


## A B C D 
rule HATBAG:
    input:
        ref = "ref/" + REFNAME + ".fa",
        vcf = expand("vcf/{vcf_prefix}.filtered.vcf.gz", vcf_prefix = VCF_PREFIX),
        tbi = expand("vcf/{vcf_prefix}.filtered.vcf.gz.tbi", vcf_prefix = VCF_PREFIX),
    output:
        decoy = expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE)
    params:
        N='hatbag',
        threads=8,
        queue = "short.qc"
    shell:
        'mkdir -p /tmp/{HATBAG_OUTPUT_DIR} && '
        'for RUN in E F; do '
        '(export PATH=/apps/well/htslib/1.4.1/bin/:${{PATH}} && '
        'which tabix && '
        'hostname && '
        'output_dir=hatbag/{HATBAG_OUTPUT_DIR}/ && '
        'mkdir -p ${{output_dir}} && '
        '{R_HATBAG} '
        '--species={HATBAG_OUTPUT_DIR} '
        '--run=${{RUN}} '
        '--nCores={params.threads} '
        '--outputDate={HATBAG_OUTPUT_DATE} '
        '--outputDir=${{output_dir}} '
        '--vcf_file={input.vcf} '
        '--reference={input.ref} '
        '--rmask_file={HATBAG_RMASK} '
        '--Klist=10 '
        '--chrlist=\'{HATBAG_CHRS}\' '
        '--genomeSize={HATBAG_GENOME_SIZE} '
        '--lineages=\'{HATBAG_LINEAGES}\' '
        '--lineages_to_build=\'{HATBAG_LINEAGES_TO_BUILD}\' '	
        '--ancestral_lineage=\'{HATBAG_ANCESTRAL_LINEAGE}\' '
        '--outgroups=\'{HATBAG_OUTGROUPS}\' '
        '--simpleRepeat_file={HATBAG_SIMPLE_REPEAT} '
        '--use_one_sided_pvalue=TRUE '
        '--mrle=6 '
        '--ndge=3 '
        ')&> '
        '{output.decoy}.${{RUN}}.log ; '
        'done && '
        'touch {output.decoy}'


