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

if [ -d "${SIMULATE_DIR}" ]
then
    rm -r "${SIMULATE_DIR}"
fi
if [ -d "${ANALYSIS_DIR}" ]
then
    rm -r "${ANALYSIS_DIR}"
fi

$HATBAG_DIR/HATBAG_simulate_input.R --chr_lengths='c(2e5,2e5)' --seed=145 --verbose=TRUE --outputDir="$SIMULATE_DIR" --model=1 --enriched_in_tails=FALSE --simulate_reads=TRUE

# to ensure sufficient test coverage, split outgroup into two fastq units
# as the code does something different if there are one or two
n_entries=`gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | wc -l`
half=$((${n_entries}/2))
gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | head -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupA_1.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_1.fastq.gz" | tail -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupB_1.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_2.fastq.gz" | head -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupA_2.fastq.gz"
gunzip -c "${SIMULATE_DIR}outgroup_2.fastq.gz" | tail -n ${half} | gzip -1 > "${SIMULATE_DIR}outgroupB_2.fastq.gz"
rm "${SIMULATE_DIR}outgroup_1.fastq.gz"
rm "${SIMULATE_DIR}outgroup_2.fastq.gz"

cd "$SIMULATE_DIR"


TEST_UNIT_NAMES=$(find *.fastq.gz | cut -d _ -f 1 | uniq)

echo "species,units,fastq_ftp,n_mapping_pieces,mapping_queue,lb,lb_insert_size,flowcell_barcode,flowcell_lane" > $ORDER_CSV
# TODO: specify all analysis subfolders up front in config file
for unit_name in $TEST_UNIT_NAMES; do
    species="test_${unit_name}"
    if [ "${unit_name}" == "outgroupA" ]
    then
	species="test_outgroup"
    elif [ "${unit_name}" == "outgroupB" ]
    then
	species="test_outgroup"
    fi
    echo "\"${species}\",\"${unit_name}\",\"\",20,\"short\",\"dummy_lb\",100,\"X1\",1" >> $ORDER_CSV
    # copy fastqs
    SPECIES_DIR="${ANALYSIS_DIR}/mapping/${SPECIES_ORDER}/${species}"
    mkdir -p $SPECIES_DIR
    rsync -a "${unit_name}"* $SPECIES_DIR
done

## Copy reference files
TEST_REF_DIR="${ANALYSIS_DIR}/${REF_DIR}/"
TEST_EXTERNAL_DIR="${ANALYSIS_DIR}/${EXTERNAL_DIR}/"
mkdir -p $TEST_REF_DIR $TEST_EXTERNAL_DIR
rsync -a ref.fa.gz ${TEST_REF_DIR}
rsync -a rmask.gz ${TEST_EXTERNAL_DIR}ref.rmsk.gz
rsync -a simpleRepeat.gz ${TEST_EXTERNAL_DIR}ref.simpleRepeat.gz

cd "${SCRIPTPATH}/../"

## split into a few bits, make sure it is OK to re-run partway through
WHERE="local"
WHERE="cluster"
./run.sh config/test_run1.json mapping/test/test_pop1/test_pop1.unitpop1.bam ${WHERE} is_test $other
./run.sh config/test_run1.json downstream_all ${WHERE} is_test $other
./run.sh config/test_run1.json all ${WHERE} is_test $other

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
