#!/usr/bin/env bash
## e.g. run.sh config/test.json all local --dryrun

set -e

if [ -z "$ANALYSIS_DIR" ]
then
      echo "\$ANALYSIS_DIR is undefined, are you sure you ran '. activate'?"
      exit 1
fi

order_config=$1
what=$2
where=$3
other=${@:4} # all other CLI arguments passed to snakemake

if [ "$where" == "cluster" ]
then
    jobs=1000
elif [ "$where" == "local" ]
then
     jobs=8
fi

# if $other is integer, use as number of cores
if [ "$other" -eq "$other" ] 2> /dev/null 
then
    jobs=$other
    other="" # so not passed to Snakemake
fi

echo "Running with ${jobs} cores"

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

export ORDER_CONFIG="${SCRIPTPATH}/${order_config}"

# From jq docs on how to load variables
# SPECIES_LIST=( $(jq -r '.SPECIES_LIST | keys[]' ${order_config}) )
eval "$(jq -r '@sh "SPECIES_ORDER=\(.SPECIES_ORDER)"' ${order_config})"
eval "$(jq -r '@sh "SPECIES_MAP_DIR_NAME=\(.SPECIES_MAP_DIR_NAME)"' config/filenames.json)"

echo "Motif Death output in ${ANALYSIS_DIR}"

export ORDER_CSV="${SCRIPTPATH}/${SPECIES_MAP_DIR_NAME}/${SPECIES_ORDER}.csv"

if [ -f $ORDER_CSV ]
then
    rm $ORDER_CSV
fi

R -f R/create_order_csv.R --args $ORDER_CONFIG $ORDER_CSV

LOG_DIR=${ANALYSIS_DIR}logs/
mkdir -p ${ANALYSIS_DIR}
mkdir -p ${LOG_DIR}

SNAKEFILE="${SCRIPTPATH}/Snakefile"

cd ${ANALYSIS_DIR}

if [ $where == "cluster" ]
then
    ##--verbose 
    ##--printshellcmds
    ${SNAKEMAKE} \
        --snakefile ${SNAKEFILE} \
        -w 30 \
	 --max-status-checks-per-second 0.1 \
        --cluster "qsub -cwd -V -N {params.N} -pe shmem {params.threads} -q {params.queue} -P davies.prjc -j Y -o ${LOG_DIR}" --jobs ${jobs} \
         ${other} ${what} \
        --configfiles ${ORDER_CONFIG} "${SCRIPTPATH}/config/filenames.json"
    ## qsub -V -N {params.N} -j oe -o ${LOG_DIR} -l nodes=1:ppn={params.threads}
elif [ $where == "local" ]
then
    ${SNAKEMAKE} \
        --snakefile ${SNAKEFILE} \
        --jobs ${jobs} \
        ${other} ${what} \
        --configfiles ${ORDER_CONFIG} "${SCRIPTPATH}/config/filenames.json"
echo done
else
    echo bad where: ${where}
    exit 1
fi
