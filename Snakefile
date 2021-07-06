import pandas as pd

VCF_PREFIX = config["VCF_PREFIX"]
HATBAG_OUTPUT_DIR = config["HATBAG_OUTPUT_DIR"]
HATBAG_OUTPUT_DATE=config["HATBAG_OUTPUT_DATE"]
R_DIR=config["R_DIR"]
SPECIES_ORDER=config["SPECIES_ORDER"]
HATBAG_DIR=config["HATBAG_DIR"]
REF_DIR = config["REF_DIR"]
REFNAME = config["REFNAME"]
REF_URL = config["REF_URL"]
BAM_SUFFIX = config["BAM_SUFFIX"]
WILDCARD_CHR_CONSTRAINT = config["WILDCARD_CHR_CONSTRAINT"]
## WILDCARD_CHR_CONSTRAINT = '[A-Z0-9-_.]+'
GENOTYPING_QUEUE=config["GENOTYPING_QUEUE"]
GENOTYPING_THREADS=config["GENOTYPING_THREADS"] ## might get lucky!
GENOTYPER=config["GENOTYPER"]
SPECIES_LIST=config["SPECIES_LIST"]
CHR_LIST_ONLY_AUTOS=config["CHR_LIST_ONLY_AUTOS"]
VCF_PREFIX=config["VCF_PREFIX"]
TREEMIX_PREFIX=config["TREEMIX_PREFIX"]
TREEMIX_K=config["TREEMIX_K"]
TREEMIX_OUTGROUP=config["TREEMIX_OUTGROUP"]
TREEMIX_THREADS=config["TREEMIX_THREADS"]
GATK_CHR_PREFIX = config["GATK_CHR_PREFIX"]
OPERATE_GATK_PER_CHR = config["OPERATE_GATK_PER_CHR"]
REF_DIR = config["REF_DIR"]
REFNAME = config["REFNAME"]
REF_URL = config["REF_URL"]
FASTQ_SUFFIX = config["FASTQ_SUFFIX"]
CHR_LIST = config["CHR_LIST"]

R_GET_GENOME_STATS=R_DIR + config["R_GET_GENOME_STATS"]
R_GET_PER_SAMPLE_AVERAGE_COV=R_DIR + config["R_GET_PER_SAMPLE_AVERAGE_COV"]
VCF2TREEMIX=config["PYTHON_DIR"] + config["VCF2TREEMIX"]
TREEMIX_R_PLOT=R_DIR + config["TREEMIX_R_PLOT"]
TREEMIX_R_PLOT_HELPER=R_DIR + config["TREEMIX_R_PLOT_HELPER"]

order_df = pd.read_csv('/users/davies/zri347/proj/motif_death/test.csv').set_index(["species", "units"], drop=False)
# order_df.index = order_df.index.set_levels([i.astype(str) for i in order_df.index.levels])  # enforce str in index

rule all:
    input:
        expand("hatbag/{hatbag_output_dir}/{hatbag_output_date}/F_complete", hatbag_output_dir = HATBAG_OUTPUT_DIR, hatbag_output_date = HATBAG_OUTPUT_DATE),
        expand("treemix/{treemix_prefix}.treemix.migrants.{migrants}.out.treeout.gz", treemix_prefix = TREEMIX_PREFIX, migrants = TREEMIX_MIGRANT_RANGE)

include: "rules/download.smk"
include: "rules/prep_reference.smk"
include: "rules/mapping.smk"
include: "rules/downstream.smk"
include: "rules/HATBAG.smk"
