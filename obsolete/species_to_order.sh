# if [ "${species}" == "human" ] || [ "${species}" == "chimp" ] || [ "${species}" == "gorilla" ] || [ "${species}" == "orangutan" ] || [ "${species}" == "neanderthal" ]  || [ "${species}" == "denisovan" ]
# then
#     SPECIES_ORDER="hominoidea"
# elif [ "${species}" == "baboon" ] || [ "${species}" == "marmoset" ] ||  [ "${species}" == "baboon" ] || [ "${species}" == "macaque" ]
# then
#     SPECIES_ORDER="cercopithecoidea"
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
# else
#     echo Cannot determine order
#     exit 1
# fi