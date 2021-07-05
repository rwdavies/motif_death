#!/usr/bin/env bash

set -e
. activate

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

echo ${SCRIPTPATH}

SIMULATE_DIR="${SCRIPTPATH}/../"simulate_input/ # in top level of motif_death

## HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir="$SIMULATE_DIR" --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

cd "$SIMULATE_DIR"
TEST_NAMES=$(find *.fastq.gz | cut -d _ -f 1 | uniq)

# TODO: specify all analysis subfolders up front in config file
for name in $TEST_NAMES; do
    echo $name
    # copy fastqs
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/${name}"
    mkdir -p $SPECIES_DIR
    cp "${name}"* $SPECIES_DIR
    # create json
    jq --null-input --arg species "test_$name" --arg unit "$name" \
        '{"species": $species, "n_mapping_pieces": 20, "platform": "Illumina", "mapping_queue": "short.qc", "units": {($unit): {"lb": "dummyout", "lb_insert_size": "100", "flowcell_barcode": "X1", "flowcell_lane": "1"}}}' > \
        "${SCRIPTPATH}/../${SPECIES_MAP_DIR_NAME}/test_${name}.json"
done

# Copy reference files
REF_DIR="${ANALYSIS_DIR}/ref/"
mkdir -p $REF_DIR
rsync -av . $REF_DIR --exclude=*.fastq.gz

# ? Run test here?

# Check test here
test_file="F_document/qq.losslin.all.K6.png"

export TEST_HATBAG_DIR=$1

# Test test_file exists
if test -f "${TEST_HATBAG_DIR}/${test_file}"; then
    echo "test passed; ${test_file} exists"
else
    echo "test failed; ${test_file} does not exist"
    exit 1
fi

# Test significant k-mers
./R/test_HATBAG.R
