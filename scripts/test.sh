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
eval "$(jq -r '@sh "REF_DIR=\(.REF_DIR) EXTERNAL_DIR=\(.EXTERNAL_DIR) SPECIES_MAPPING_INFO_DIR_NAME=\(.SPECIES_MAP_DIR_NAME)"' ${SCRIPTPATH}/../config/filenames.json)"
FULL_HATBAG_OUTPUT_DIR="${ANALYSIS_DIR}hatbag/${SPECIES_ORDER}/${RUN_ID}/${HATBAG_OUTPUT_DIR}"

if [ -d "${SIMULATE_DIR}" ]
then
    rm -r "${SIMULATE_DIR}"
fi
if [ -d "${ANALYSIS_DIR}" ]
then
    rm -r "${ANALYSIS_DIR}"
fi

$HATBAG_DIR/HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir="$SIMULATE_DIR" --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

cd "$SIMULATE_DIR"

TEST_NAMES=$(find *.fastq.gz | cut -d _ -f 1 | uniq)

# TODO: specify all analysis subfolders up front in config file
for name in $TEST_NAMES; do
    # echo $name
    # copy fastqs
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/${SPECIES_ORDER}/test_${name}"
    mkdir -p $SPECIES_DIR
    rsync -a "${name}"* $SPECIES_DIR
    # create json
    jq --null-input --arg species "test_$name" --arg unit "$name" \
        '{"species": $species, "n_mapping_pieces": 20, "platform": "Illumina", "mapping_queue": "short.qc", "units": {($unit): {"1": "dummy_URL", "2": "dummy_URL", "lb": "dummyout", "lb_insert_size": "100", "flowcell_barcode": "X1", "flowcell_lane": "1"}}}' > \
        "${SCRIPTPATH}/../${SPECIES_MAP_DIR_NAME}/test_${name}.json"
done

## Copy reference files
TEST_REF_DIR="${ANALYSIS_DIR}/${REF_DIR}/"
TEST_EXTERNAL_DIR="${ANALYSIS_DIR}/${EXTERNAL_DIR}/"
mkdir -p $TEST_REF_DIR $TEST_EXTERNAL_DIR
rsync -a ref.fa.gz ${TEST_REF_DIR}
rsync -a rmask.gz ${TEST_EXTERNAL_DIR}ref.rmsk.gz
rsync -a simpleRepeat.gz ${TEST_EXTERNAL_DIR}ref.simpleRepeat.gz

cd "${SCRIPTPATH}/../"

./run.sh config/test_run1.json all local $other

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
