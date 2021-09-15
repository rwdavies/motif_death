#!/usr/bin/env bash

set -e

if [ -z "$ANALYSIS_DIR" ]
then
      echo "\$ANALYSIS_DIR is undefined, are you sure you ran '. activate'?"
      exit 1
fi

other=$1 # n_cores if integer, or can be `--dryrun`, `--debug-dag`, `--reason`

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

SIMULATE_DIR="${SCRIPTPATH}/../"simulate_input/ # in top level of motif_death
export ANALYSIS_DIR="${SCRIPTPATH}/../"simulate_results/ # in top level of motif_death
TEST_CONFIG_PATH="${SCRIPTPATH}/../config/test_run1.json"
eval "$(jq -r '@sh "SPECIES_ORDER=\(.SPECIES_ORDER) RUN_ID=\(.RUN_ID) HATBAG_OUTPUT_DIR=\(.HATBAG_OUTPUT_DIR)"' ${TEST_CONFIG_PATH})"
eval "$(jq -r '@sh "REF_DIR=\(.REF_DIR) EXTERNAL_DIR=\(.EXTERNAL_DIR) SPECIES_MAP_DIR=\(.SPECIES_MAP_DIR_NAME)"' ${SCRIPTPATH}/../config/filenames.json)"
FULL_HATBAG_OUTPUT_DIR="${ANALYSIS_DIR}hatbag/${SPECIES_ORDER}/${RUN_ID}/${HATBAG_OUTPUT_DIR}"

export ORDER_CSV="${SCRIPTPATH}/../${SPECIES_MAP_DIR}/${SPECIES_ORDER}.csv"

if [ -f $ORDER_CSV ]
then
    rm $ORDER_CSV
fi

echo "species,units,fastq_ftp,n_mapping_pieces,mapping_queue,lb,lb_insert_size,flowcell_barcode,flowcell_lane" > $ORDER_CSV

if [ -d "${SIMULATE_DIR}" ]
then
    rm -r "${SIMULATE_DIR}"
fi
if [ -d "${ANALYSIS_DIR}" ]
then
    rm -r "${ANALYSIS_DIR}"
fi

# to ensure sufficient test coverage, split outgroup into two fastq units
# as the code does something different if there are one or two
n_entries=`gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | wc -l`
half=`$((echo ${n_entries} / 2))`
gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | head -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupA_1.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | tail -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupB_1.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_2.fastq.gz" | head -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupA_2.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_2.fastq.gz" | tail -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupB_2.fastq.gz"


$HATBAG_DIR/HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir="$SIMULATE_DIR" --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

cd "$SIMULATE_DIR"

TEST_NAMES=$(find *.fastq.gz | cut -d _ -f 1 | uniq)

# TODO: specify all analysis subfolders up front in config file
for name in $TEST_NAMES; do
    # copy fastqs
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/${SPECIES_ORDER}/test_${name}"
    mkdir -p $SPECIES_DIR
    rsync -a "${name}"* $SPECIES_DIR

    echo "\"test_${name}\",\"${name}\",\"\",20,\"short.qc\",\"dummy_lb\",100,\"X1\",1" >> $ORDER_CSV
done

## Copy reference files
TEST_REF_DIR="${ANALYSIS_DIR}/${REF_DIR}/"
TEST_EXTERNAL_DIR="${ANALYSIS_DIR}/${EXTERNAL_DIR}/"
mkdir -p $TEST_REF_DIR $TEST_EXTERNAL_DIR
rsync -a ref.fa.gz ${TEST_REF_DIR}
rsync -a rmask.gz ${TEST_EXTERNAL_DIR}ref.rmsk.gz
rsync -a simpleRepeat.gz ${TEST_EXTERNAL_DIR}ref.simpleRepeat.gz

cd "${SCRIPTPATH}/../"

./run.sh config/test_run1.json all local is_test $other

# Check test here

# Test test_file exists
test_file=${FULL_HATBAG_OUTPUT_DIR}/F_document/qq.losslin.all.K6.png
if test -f "${test_file}"; then
    echo "test passed; ${test_file} exists"
else
    echo "test failed; ${test_file} does not exist"
    exit 1
fi

# Test significant k-mers
./R/test_HATBAG.R ${FULL_HATBAG_OUTPUT_DIR}
