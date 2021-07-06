#!/usr/bin/env bash

set -e
## . activate # if load activate again, need to account for test ANALYSIS_DIR

## e.g. run.sh config/test.json all local --dryrun

# component=$1 ## either mapping, preprocess, or 
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

# put this all in config
## if [ "${species}" == "human" ] || [ "${species}" == "chimp" ] || [ "${species}" == "gorilla" ] || [ "${species}" == "orangutan" ] || [ "${species}" == "baboon" ] || [ "${species}" == "marmoset" ] ||  [ "${species}" == "baboon" ] || [ "${species}" == "macaque" ] || [ "${species}" == "snubnosedmonkey" ] || [ "${species}" == "vervet" ] || [ "${species}" == "squirrelmonkey" ] || [ "${species}" == "neanderthal" ]  || [ "${species}" == "denisovan" ]

# if [ "${species}" == "human" ] || [ "${species}" == "chimp" ] || [ "${species}" == "gorilla" ] || [ "${species}" == "orangutan" ] || [ "${species}" == "neanderthal" ]  || [ "${species}" == "denisovan" ]
# then
#     SPECIES_ORDER="hominoidea"
# elif [ "${species}" == "baboon" ] || [ "${species}" == "marmoset" ] ||  [ "${species}" == "baboon" ] || [ "${species}" == "macaque" ]
# then
#     SPECIES_ORDER="cercopithecoidea"
# elif [ "${species}" == "cow" ] || [ "${species}" == "goat" ] || [ "${species}" == "giraffe" ] || [ "${species}" == "okapi" ] || [ "${species}" == "whitetaileddeer" ]  || [ "${species}" == "reddeer" ] || [ "${species}" == "buffalo" ]
# then
#     SPECIES_ORDER="artiodactyla"
# elif [ "${species}" == "FAM" ] || [ "${species}" == "caroli" ] || [ "${species}" == "CAST_EiJ" ] || [ "${species}" == "WSB_EiJ" ] || [ "${species}" == "SPRET_EiJ" ] || [ "${species}" == "PWK_PhJ" ]
# then
#     SPECIES_ORDER="mice"
# elif [ "${species}" == "felidae" ]
# then
#     SPECIES_ORDER="felidae"
# elif [ "${species}" == "zebra_finch" ] ||  [ "${species}" == "rifleman" ] ||  [ "${species}" == "medium_ground_finch" ] ||  [ "${species}" == "long_tailed_finch" ] || [ "${species}" == "kea" ] || [ "${species}" == "golden_collared_manakin" ] ||  [ "${species}" == "double_barrelled_finch" ] || [ "${species}" == "american_crow" ]
# then
#     SPECIES_ORDER="avian"
# elif [ "${species}" == "arctic_char" ] ||  [ "${species}" == "atlantic_salmon" ] ||  [ "${species}" == "chinook" ] ||  [ "${species}" == "coho" ] || [ "${species}" == "rainbow_trout" ]
# then
#     SPECIES_ORDER="salmon"
# elif [ "${species}" == "big_brown_bat" ] ||  [ "${species}" == "brandts_bat" ] ||  [ "${species}" == "davids_bat" ] ||  [ "${species}" == "red_bat" ]
# then
#     SPECIES_ORDER="bat"
# elif [ "${species}" == "lizards" ] || [ "${species}" == "A_angusticeps" ] || [ "${species}" == "A_chlorocyanus" ] || [ "${species}" == "A_cristatellus" ] || [ "${species}" == "A_cybotes" ] || [ "${species}" == "A_evermanni" ] || [ "${species}" == "A_grahami" ] || [ "${species}" == "A_insolitus" ] || [ "${species}" == "A_lineatopus" ] || [ "${species}" == "A_occultus" ] || [ "${species}" == "A_sagrei" ] || [ "${species}" == "A_valencienni" ]
# then
#     SPECIES_ORDER="lizards"
# elif [ "${species}" == "canidae" ] || [ "${species}" == "african_golden_wolf" ] || [ "${species}" == "african_hunting_dog" ] || [ "${species}" == "andean_fox" ] || [ "${species}" == "coyote" ] || [ "${species}" == "dhole" ] || [ "${species}" == "dog" ] || [ "${species}" == "golden_jackal" ] || [ "${species}" == "grey_fox" ] || [ "${species}" == "island_fox" ] || [ "${species}" == "red_fox" ]
# then
#     SPECIES_ORDER="canidae"
# elif [ "${species}" == "ursidae" ]
# then
#     SPECIES_ORDER="ursidae"    
# elif [ "${species}" == "test_pop1" ] || [ "${species}" == "test_pop2" ] || [ "${species}" == "test_outgroup" ]
# then
#     SPECIES_ORDER="test"     
# else
#     echo Cannot determine order
#     exit 1
# fi

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

# From jq docs on how to load variables
eval "$(jq -r '@sh "SPECIES_LIST=( \([.SPECIES_LIST[]]) ) SPECIES_ORDER=\(.SPECIES_ORDER)"' ${order_config})"

# if [ $other = "test" ]
# then
#     # use ANALYSIS_DIR defined in scripts/test.sh
#     other="" # so not passed to snakemake

# fi
echo "Output in ${ANALYSIS_DIR}"

# eval "$(jq -r '@sh "SPECIES_MAP_DIR_NAME=\(.SPECIES_MAP_DIR_NAME)"' config/filenames.json)"

# json_list=""
# for name in ${SPECIES_LIST[@]}; do
#     # echo $name
#     json_list+=" ${SPECIES_MAP_DIR_NAME}/${name}.json"
# done

# TODO: change this to ouput csv. 
# jq -rs 'reduce .[] as $item ({}; . * $item) ' $json_list > "${SPECIES_ORDER}.json"

LOG_DIR="${ANALYSIS_DIR}/logs/"
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
        --cluster "qsub -cwd -V -N {params.N} -pe shmem {params.threads} -q {params.queue} -P davies.prjc -j Y -o "${ANALYSIS_DIR}"logs/" --jobs ${jobs} \
         ${other} ${what} \
        --configfiles "${SCRIPTPATH}/${order_config}" "${SCRIPTPATH}/config/filenames.json"
    ## "qsub -V -N {params.N} -j oe -o ${LOG_DIR} -l nodes=1:ppn={params.threads}
elif [ $where == "local" ]
then
    ${SNAKEMAKE} \
        --snakefile ${SNAKEFILE} \
        --jobs ${jobs} \
        ${other} ${what} \
        --configfiles "${SCRIPTPATH}/${order_config}" "${SCRIPTPATH}/config/filenames.json"
echo done
else
    echo bad where: ${where}
    exit 1
fi
