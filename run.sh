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
other=("${@:4}") # all other CLI arguments passed to snakemake. if is_test, should be first item.

if [ "$where" == "cluster" ]
then
    jobs=1000
elif [ "$where" == "local" ]
then
     jobs=8
fi

# if $other first element is_test, run in test mode
is_test=false
if [ "${other[0]}" == "is_test" ]
then
    is_test=true
    other=("${other[@]:1}")
fi
# if $other next element is integer, use as number of cores
if [ "${other[0]}" -eq "${other[0]}" ] 2> /dev/null 
then
    jobs=${other[0]}
    other=("${other[@]:1}")
fi
echo "Running with ${jobs} cores"

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

export ORDER_CONFIG="${SCRIPTPATH}/${order_config}"

# Note: From jq docs on how to load variables
eval "$(jq -r '@sh "SPECIES_ORDER=\(.SPECIES_ORDER)"' ${order_config})"
eval "$(jq -r '@sh "SPECIES_MAP_DIR_NAME=\(.SPECIES_MAP_DIR_NAME)"' config/filenames.json)"

echo "Motif Death output in ${ANALYSIS_DIR}"

if [ "$is_test" = false ]
then
    export ORDER_CSV="${SCRIPTPATH}/${SPECIES_MAP_DIR_NAME}/${SPECIES_ORDER}.csv"
    ## only update order CSV if newer than order config
    if [ "${ORDER_CONFIG}" -nt "${ORDER_CSV}" ]
    then
	echo remake order CSV
	R -f R/create_order_csv.R --args $ORDER_CONFIG $ORDER_CSV
    fi
fi

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
         ${other[*]} ${what} \
        --configfiles ${ORDER_CONFIG} "${SCRIPTPATH}/config/filenames.json"
    ## qsub -V -N {params.N} -j oe -o ${LOG_DIR} -l nodes=1:ppn={params.threads}
elif [ $where == "local" ]
then
    ${SNAKEMAKE} \
        --snakefile ${SNAKEFILE} \
        --jobs ${jobs} \
        ${other[*]} ${what} \
        --configfiles ${ORDER_CONFIG} "${SCRIPTPATH}/config/filenames.json"
echo done
else
    echo bad where: ${where}
    exit 1
fi
