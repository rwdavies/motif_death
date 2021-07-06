#!/usr/bin/env bash

set -e
. activate


SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

SIMULATE_DIR="${SCRIPTPATH}/../"simulate_input/ # in top level of motif_death
export ANALYSIS_DIR="${SCRIPTPATH}/../"simulate_results/ # in top level of motif_death

if [ -d "${SIMULATE_DIR}" ]
then
    rm -r "${SIMULATE_DIR}"
fi
if [ -d "${ANALYSIS_DIR}" ]
then
    rm -r "${ANALYSIS_DIR}"
fi

HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir="$SIMULATE_DIR" --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

cd "$SIMULATE_DIR"

TEST_NAMES=$(find *.fastq.gz | cut -d _ -f 1 | uniq)

# TODO: specify all analysis subfolders up front in config file
for name in $TEST_NAMES; do
    echo $name
    # copy fastqs
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/test_${name}"
    mkdir -p $SPECIES_DIR
    rsync -av "${name}"* $SPECIES_DIR
    # create json
    jq --null-input --arg species "test_$name" --arg unit "$name" \
        '{"species": $species, "n_mapping_pieces": 20, "platform": "Illumina", "mapping_queue": "short.qc", "units": {($unit): {"lb": "dummyout", "lb_insert_size": "100", "flowcell_barcode": "X1", "flowcell_lane": "1"}}}' > \
        "${SCRIPTPATH}/../${SPECIES_MAP_DIR_NAME}/test_${name}.json"
done



# Copy reference files
REF_DIR="${ANALYSIS_DIR}/ref/"
mkdir -p $REF_DIR
rsync -av . $REF_DIR --exclude=*.fastq.gz

cd "${SCRIPTPATH}/../"


## TODO: be able to set number of cores when running
## ./run.sh preprocess all local test # not needed
./run.sh prep_reference all local test_pop1
./run.sh mapping all local test_pop1
./run.sh mapping all local test_pop2
./run.sh mapping all local test_outgroup
./run.sh downstream all local test_pop1
./run.sh HATBAG HATBAG_HACK local test_pop1


# Check test here
# path hardcoded in Snakefile_reference_test
test_file="${ANALYSIS_DIR}hatbag/test/dummy_date/F_document/qq.losslin.all.K6.png"

# Test test_file exists
if test -f "${test_file}"; then
    echo "test passed; ${test_file} exists"
else
    echo "test failed; ${test_file} does not exist"
    exit 1
fi

# Test significant k-mers
Rscript ./R/test_HATBAG.R
