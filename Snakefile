import os
import pandas as pd

# From Activate or run.sh
R_DIR=os.environ["R_DIR"]
HATBAG_DIR=os.environ["HATBAG_DIR"]
PYTHON_DIR=os.environ["PYTHON_DIR"]
ORDER_CSV=os.environ["ORDER_CSV"]
ORDER_CONFIG=os.environ["ORDER_CONFIG"]

# From config/{order}_{run_id}.json
RUN_ID=config["RUN_ID"]
HATBAG_OUTPUT_DIR=config["HATBAG_OUTPUT_DIR"]
SPECIES_ORDER=config["SPECIES_ORDER"]
REF_DIR = config["REF_DIR"]
EXTERNAL_DIR = config["EXTERNAL_DIR"]
REF_NAME = config["REF_NAME"]
REF_URL = config["REF_URL"]
BAM_SUFFIX = config["BAM_SUFFIX"]
SPECIES_LIST = config["SPECIES_LIST"].keys()
CHR_LIST = config["CHR_LIST"]
CHR_LIST_ONLY_AUTOS = [c for c in CHR_LIST if c not in config["SEX_CHR_LIST"]]
TREEMIX_OUTGROUP=config["TREEMIX_OUTGROUP"]
GATK_CHR_PREFIX = config["GATK_CHR_PREFIX"]
OPERATE_GATK_PER_CHR = config["OPERATE_GATK_PER_CHR"]
FASTQ_SUFFIX = config["FASTQ_SUFFIX"]
WILDCARD_CHR_CONSTRAINT = config["WILDCARD_CHR_CONSTRAINT"]
HATBAG_PARAMS=config["HATBAG_PARAMS"]

TREEMIX_THREADS=config["DEFAULTS"]["TREEMIX_THREADS"]
CALLABLE_MIN_FRACTION=config["DEFAULTS"]["CALLABLE_MIN_FRACTION"]
TREEMIX_K=config["DEFAULTS"]["TREEMIX_K"]
GENOTYPING_QUEUE=config["DEFAULTS"]["GENOTYPING_QUEUE"]
GENOTYPING_THREADS=config["DEFAULTS"]["GENOTYPING_THREADS"] ## might get lucky!
GENOTYPER=config["DEFAULTS"]["GENOTYPER"]
WILDCARD_UNIT_CONSTRAINT='[A-Za-z0-9]+' # Note: cannot include _ in here, otherwise considers {unit}_1 as unit?

simpleRepeat_URL = config["simpleRepeat_URL"]
rmask_URL = config["rmask_URL"]
# Note: need escaped single quotes on below to evaluate correctly
SIMPLE_REPEAT_HEADER = "\'#bin\tchrom\tchromStart\tchromEnd\tname\tperiod\tcopyNum\tconsensusSize\tperMatch\tperIndel\tscore\tA\tC\tG\tT\tentropy\tsequence\'"
RMASK_HEADER = "\'#bin\tswScore\tmilliDiv\tmilliDel\tmilliIns\tgenoName\tgenoStart\tgenoEnd\tgenoLeft\tstrand\trepName\trepClass\trepFamily\trepStart\trepEnd\trepLeft\tid\'"

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
    IBAMS = IBAMS + f" -I mapping/{SPECIES_ORDER}/{species}/{species}.{BAM_SUFFIX}"

MERGE_VCF_GATK_INPUT=""
for piece in CHR_CHUNKS:
    for chr in CHR_LIST_ONLY_AUTOS:
    	MERGE_VCF_GATK_INPUT = MERGE_VCF_GATK_INPUT + f" -V vcf/{SPECIES_ORDER}/{RUN_ID}/" + "chr" + str(chr) + ".filtered.piece" + str(piece) + ".vcf.gz"

if OPERATE_GATK_PER_CHR == "FALSE":
    INDEL_REALIGNMENT_QUEUE="long"
else:
    INDEL_REALIGNMENT_QUEUE="short"

order_df = pd.read_csv(ORDER_CSV).set_index(["species", "units"], drop=False)
# order_df.index = order_df.index.set_levels([i.astype(str) for i in order_df.index.levels])  # enforce str in index

rule all:
    input:
        f"hatbag/{SPECIES_ORDER}/{RUN_ID}/{HATBAG_OUTPUT_DIR}/F_complete",
        expand(f"treemix/{SPECIES_ORDER}/{RUN_ID}/treemix.migrants.{{migrants}}.out.treeout.gz", migrants = TREEMIX_MIGRANT_RANGE)

include: "rules/download.smk"
include: "rules/prep_reference.smk"
include: "rules/mapping.smk"
include: "rules/downstream.smk"
include: "rules/HATBAG.smk"
