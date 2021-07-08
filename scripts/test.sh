#!/usr/bin/env bash

set -e

if [ -z "$ANALYSIS_DIR" ]
then
      echo "\$ANALYSIS_DIR is undefined, are you sure you ran '. activate'?"
      exit 1
fi

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

SIMULATE_DIR="${SCRIPTPATH}/../"simulate_input/ # in top level of motif_death
export ANALYSIS_DIR="${SCRIPTPATH}/../"simulate_results/ # in top level of motif_death
TEST_CONFIG_PATH="${SCRIPTPATH}/../config/test.json"
eval "$(jq -r '@sh "HATBAG_OUTPUT_DIR=\(.HATBAG_OUTPUT_DIR) HATBAG_OUTPUT_DATE=\(.HATBAG_OUTPUT_DATE)"' ${TEST_CONFIG_PATH})"
HATBAG_TEST_DIR="${ANALYSIS_DIR}hatbag/${HATBAG_OUTPUT_DIR}/${HATBAG_OUTPUT_DATE}"

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
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/test_${name}"
    mkdir -p $SPECIES_DIR
    rsync -a "${name}"* $SPECIES_DIR
    # create json
    jq --null-input --arg species "test_$name" --arg unit "$name" \
        '{"species": $species, "n_mapping_pieces": 20, "platform": "Illumina", "mapping_queue": "short.qc", "units": {($unit): {"1": "dummy_URL", "2": "dummy_URL", "lb": "dummyout", "lb_insert_size": "100", "flowcell_barcode": "X1", "flowcell_lane": "1"}}}' > \
        "${SCRIPTPATH}/../species_mapping_info/test_${name}.json" # TODO: change back to SPECIES_MAP_DIR
done

# # Copy reference files
REF_DIR="${ANALYSIS_DIR}/ref/"
mkdir -p $REF_DIR
rsync -a . $REF_DIR --exclude=*.fastq.gz

cd "${SCRIPTPATH}/../"

## TODO: be able to set number of cores when running
./run.sh config/test.json all local

# Check test here

# Test test_file exists
test_file=${HATBAG_TEST_DIR}/F_document/qq.losslin.all.K6.png
if test -f "${test_file}"; then
    echo "test passed; ${test_file} exists"
else
    echo "test failed; ${test_file} does not exist"
    exit 1
fi

# Test significant k-mers
./R/test_HATBAG.R ${HATBAG_TEST_DIR}
