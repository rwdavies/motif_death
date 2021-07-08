import os
import pandas as pd

# From Activate or run.sh
R_DIR=os.environ["R_DIR"]
HATBAG_DIR=os.environ["HATBAG_DIR"]
PYTHON_DIR=os.environ["PYTHON_DIR"]
ORDER_CSV=os.environ["ORDER_CSV"]

# From config/{order}.json
VCF_PREFIX = config["VCF_PREFIX"]
HATBAG_OUTPUT_DIR = config["HATBAG_OUTPUT_DIR"]
HATBAG_OUTPUT_DATE=config["HATBAG_OUTPUT_DATE"]
HATBAG_PARAMS=config["HATBAG_PARAMS"]
SPECIES_ORDER=config["SPECIES_ORDER"]
REF_DIR = config["REF_DIR"]
REF_NAME = config["REF_NAME"]
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
FASTQ_SUFFIX = config["FASTQ_SUFFIX"]
CHR_LIST = config["CHR_LIST"]

# From config/filenames.json
R_GET_GENOME_STATS=R_DIR + config["R_GET_GENOME_STATS"]
R_GET_PER_SAMPLE_AVERAGE_COV=R_DIR + config["R_GET_PER_SAMPLE_AVERAGE_COV"]
VCF2TREEMIX=PYTHON_DIR + config["VCF2TREEMIX"]
TREEMIX_R_PLOT=R_DIR + config["TREEMIX_R_PLOT"]
TREEMIX_R_PLOT_HELPER=R_DIR + config["TREEMIX_R_PLOT_HELPER"]

# Downstream stuff
CHR_CHUNKS = range(1, 2) ## disable - too complicated as crashes if out of rang - run 1 chunk only
TREEMIX_MIGRANT_RANGE = list(range(0, 10))

if GENOTYPER == "HaplotypeCaller":
    GENOTYPING_MULTITHREAD_FLAG="-nct"
elif GENOTYPER == "UnifiedGenotyper":
    GENOTYPING_MULTITHREAD_FLAG="-nt"

IBAMS=""
for species in SPECIES_LIST:
    IBAMS = IBAMS + " -I mapping/" + species + "/" + species + "." + BAM_SUFFIX

MERGE_VCF_GATK_INPUT=""
for piece in CHR_CHUNKS:
    for chr in CHR_LIST_ONLY_AUTOS:
    	MERGE_VCF_GATK_INPUT = MERGE_VCF_GATK_INPUT + " -V vcf/" + VCF_PREFIX + ".chr" + str(chr) + ".filtered.piece" + str(piece) + ".vcf.gz"

if OPERATE_GATK_PER_CHR == "FALSE":
    INDEL_REALIGNMENT_QUEUE="long.qc"
else:
    INDEL_REALIGNMENT_QUEUE="short.qc"

order_df = pd.read_csv(ORDER_CSV).set_index(["species", "units"], drop=False)
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
