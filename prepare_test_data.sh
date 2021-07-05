#!/usr/bin/env bash

set -e
. activate

SIMULATE_DIR=./simulate_input/ # with reference to HATBAG_DIR
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

cd $HATBAG_DIR
./HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir=$SIMULATE_DIR --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

cd $SIMULATE_DIR
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
        "${SCRIPTPATH}/${SPECIES_MAP_DIR_NAME}/test_${name}.json"
done

# Copy reference files
REF_DIR="${ANALYSIS_DIR}/ref/"
mkdir -p $REF_DIR
rsync -av . $REF_DIR --exclude=*.fastq.gz
