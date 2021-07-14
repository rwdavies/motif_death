get_params <- function(outputDate) {
    ## "default"
    num_non_missing_outgroups_required <- NULL
    num_missing_lineages_allowed <- NULL
    similar_kmer_criterion <- c("cpg", "at")
    mrle <- 6
    mncdnle <- 100 ## make irrelevant
    use_one_sided_pvalue <- FALSE    
    if (outputDate == "2018_04_27") {
        ## everything as default above
    } 
    if ((outputDate == "2018_05_30") | (outputDate == "2018_07_18") | (outputDate == "2018_11_21")) {
        num_non_missing_outgroups_required <- 1 ## only require 1 outgroup
        num_missing_lineages_allowed <- 1 ## allow up to 1 missing lineage
        similar_kmer_criterion <- c("cpg", "at", "mrle") ## add in mrle as something to control for
        mrle <- 5 ## shrink max homo to 5
        mncdnle <- 3 ## maximum number of consecutive dinucleotides
        use_one_sided_pvalue <- TRUE ## I think this is OK
    }
    return(
        list(
            num_non_missing_outgroups_required =num_non_missing_outgroups_required,
            num_missing_lineages_allowed =num_missing_lineages_allowed,
            similar_kmer_criterion = similar_kmer_criterion,
            mrle = mrle,
            mncdnle = mncdnle,
            use_one_sided_pvalue = use_one_sided_pvalue
        )
    )
}



convert_dict_to_chrlist <- function(dict_file, min_size = 1e6) {
    dict <- read.table(dict_file, header = FALSE, skip = 1)
    dict$LN <- as.numeric(substr(dict[, 3], 4, 100))
    chrlist <- substr(dict[dict[, "LN"] >= min_size, 2], 4, 100)
    return(chrlist)
}

# TODO: species should be order?
get_per_species_params <- function(species) {
    # if (species == "primates") {
    #     nCores <- 10
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.simpleRepeat.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/primates.ndhcgobmsvms.GATKug.filtered.vcf.gz"
    #     reference = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.fa.gz"
    #     chrlist = paste0("chr", c(1:22))
    #     genomeSize = 3137144693
    #     rmask_file = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.rmsk.gz"
    #     lineages <- list(
    #         "neanderthal" = "neanderthal",
    #         "human" = c("human"),
    #         "AHN" = c("human", "neanderthal"),
    #         "chimp" = c("chimp"),
    #         "AHNC" = c("AHN", "chimp"),
    #         "gorilla" = c("gorilla"),
    #         "AHNCG" = c("AHNC", "gorilla"),
    #         "orangutan" = "orangutan",
    #         "AHNCGO" = c("AHNCG", "orangutan"),
    #         "baboon" = c("baboon"),
    #         "macaque" = c("macaque"),
    #         "ABM" = c("baboon", "macaque"),
    #         "vervet" = c("vervet"),
    #         "ABMV" = c("ABM", "vervet"),
    #         "snubnosedmonkey" = c("snubnosedmonkey"),
    #         "ABMVS" = c("ABMV", "snubnosedmonkey")
    #     )
    #     ancestral_lineage = list(
    #         "AHNCGOBMVS" = c("ABMVS", "AHNCGO")
    #     )
    #     outgroups = c("marmoset", "squirrelmonkey")
    # }
    # if (species == "hominoidea") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.simpleRepeat.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/primates.ndhcgobmsv.GATKug.filtered.vcf.gz"
    #     reference = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.fa.gz"
    #     chrlist = paste0("chr", c(1:22))
    #     genomeSize = 3137144693
    #     rmask_file = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.rmsk.gz"
    #     lineages <- list(
    #         "neanderthal" = "neanderthal",
    #         "human" = c("human"),
    #         "AHN" = c("human", "neanderthal"),
    #         "chimp" = c("chimp"),
    #         "AHNC" = c("AHN", "chimp"),
    #         "gorilla" = c("gorilla"),
    #         "AHNCG" = c("AHNC", "gorilla"),
    #         "orangutan" = "orangutan"
    #     )
    #     ancestral_lineage = list(
    #         "AHNCGO" = c("AHNCG", "orangutan")
    #     )
    #     outgroups = c("baboon", "macaque", "vervet", "snubnosedmonkey")
    # }
    # if (species == "cercopithecidae") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.simpleRepeat.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/primates.ndhcgobmsv.GATKug.filtered.vcf.gz"
    #     reference = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.fa.gz"
    #     chrlist = paste0("chr", c(1:22))
    #     genomeSize = 3137144693
    #     rmask_file = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.rmsk.gz"
    #     lineages <- list(
    #         "baboon" = c("baboon"),
    #         "macaque" = c("macaque"),
    #         "ABM" = c("baboon", "macaque"),
    #         "vervet" = c("vervet"),
    #         "ABMV" = c("ABM", "vervet"),
    #         "snubnosedmonkey" = c("snubnosedmonkey")
    #     )
    #     ancestral_lineage = list(
    #         "ABMVS" = c("ABMV", "snubnosedmonkey")
    #     )
    #     outgroups = c("human", "chimp", "gorilla", "orangutan")
    # }
    # if ((species == "primates_nean") | (species == "primates_deni")) {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.simpleRepeat.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/primates/primates.ndhcgo.GATKug.filtered.vcf.gz"
    #     reference = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.fa.gz"
    #     chrlist = paste0("chr", c(1:22))
    #     genomeSize = 3137144693
    #     rmask_file = "/data/smew1/rdavies/motifLossAnalysis/primates/hg38.rmsk.gz"
    #     if (species == "primates_nean") {
    #         lineages = list(
    #             "neanderthal" = "neanderthal",
    #             "human" = c("human"),
    #             "AHN" = c("human", "neanderthal"),
    #             "chimp" = c("chimp"),
    #             "AHNC" = c("AHN", "chimp"),
    #             "gorilla" = c("gorilla")
    #         )
    #         ancestral_lineage = list(
    #             "AHNCG" = c("AHNC", "gorilla")
    #         )
    #         outgroups = c("orangutan")
    #     } else if (species == "primates_deni") {
    #         lineages = list(
    #             "denisovan" = "denisovan",
    #             "human" = c("human"),
    #             "AHD" = c("human", "denisovan"),
    #             "chimp" = c("chimp"),
    #             "AHDC" = c("AHD", "chimp"),
    #             "gorilla" = c("gorilla")
    #         )
    #         ancestral_lineage = list(
    #             "AHDCG" = c("AHDC", "gorilla")
    #         )
    #         outgroups = c("orangutan")
    #     }
    # }
    # if (species == "mice") {
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/mice/mm9.simpleRepeat.no_chr.gz"
    #     nCores <- 16
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/mice/mice.cfswpc.filtered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/wildmice/ref/NCBIM37_um.fa"
    #     chrlist = c(1:19)
    #     genomeSize =  2654895218
    #     rmask_file = "/data/smew1/rdavies/wildmice/downloads/mm9.rmask.nochr.gz"
    #     lineages = list(
    #         "WSBEiJ" = c("WSBEiJ"),
    #         "PWKPhJ" = c("PWKPhJ"),
    #         "CASTEiJ" = c("CASTEiJ"),
    #         "WSBEiJ.PWKPhJ" = c("WSBEiJ", "PWKPhJ"),
    #         "WSBEiJ.CASTEiJ" = c("WSBEiJ", "CASTEiJ"),
    #         "PWKPhJ.CASTEiJ" = c("PWKPhJ", "CASTEiJ"),
    #         "AM" = c("WSBEiJ", "PWKPhJ", "CASTEiJ", "WSBEiJ.PWKPhJ", "WSBEiJ.CASTEiJ", "PWKPhJ.CASTEiJ"),
    #         "Spretus" = "Spretus",
    #         "AMS" = c("AM", "Spretus"),
    #         "FAM" = "FAM"
    #     )
    #     lineages_to_build = list(
    #         "WSBEiJ.PWKPhJ" = c("WSBEiJ", "PWKPhJ"),
    #         "WSBEiJ.CASTEiJ" = c("WSBEiJ", "CASTEiJ"),
    #         "PWKPhJ.CASTEiJ" = c("PWKPhJ", "CASTEiJ")
    #     )
    #     ancestral_lineage = list(
    #         "AMSF" = c("FAM", "AMS")
    #     )
    #     outgroups = c("Caroli")
    # }
    if (species == "ruminantia" | species == "artiodactyla") {
        return(
            list(
                ## have manually rsynced these over
                reference = file.path("ref/bosTau8.fa.gz"),
                simpleRepeat_file = file.path("external/bosTau8.simpleRepeat.gz"),
                rmask_file = file.path("external/bosTau8.rmsk.gz"),
                # vcf_file = NA, ##  now specified elsewhere
                ## /data/smew1/rdavies/motifLossAnalysis/ruminantia
                ## reference <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/ref/bosTau8.fa.gz"
                ## vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bovidae.cbgogwd.GATKug.filtered.vcf.gz"
                # nCores = NA, ## now specified elsewhere
                ## simpleRepeat_file <- "/bosTau8.simpleRepeat.gz"
                ## reference <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/ref/bosTau8.fa.gz"
                ## vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bovidae.cbgogwd.GATKug.filtered.vcf.gz"
                ## rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bosTau8.rmsk.gz"
                chrlist = paste0("chr", 1:29),
                genomeSize = 2670044500,
                lineages_to_build = NULL,
                lineages = list(
                    "cow" = "cow",
                    "buffalo" = "buffalo",
                    "ACB" = c("cow", "buffalo"),
                    "goat" = "goat",
                    "ACBG" = c("ACB", "goat"),
                    "giraffe" = "giraffe",
                    "okapi" = "okapi",
                    "AOGi" = c("giraffe", "okapi")
                ),
                ancestral_lineage = list(
                    "AOGiACBG" = c("AOGi", "ACBG")
                ),
                # outgroups = c("reddeer", "whitetaileddeer")
                ## if no reddeer (yet)
                outgroups = c("whitetaileddeer")
            )
        )
    }
    if (species == "whippomorpha") {
        return(
            list(
                ## have manually rsynced these over
                reference = file.path("ref/balAcu1.fa.gz"),
                simpleRepeat_file = file.path("external/balAcu1.simpleRepeat.gz"),
                rmask_file = file.path("external/balAcu1.rmsk.gz"),
                # vcf_file = NA, ##  now specified elsewhere
                ## /data/smew1/rdavies/motifLossAnalysis/ruminantia
                ## reference <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/ref/bosTau8.fa.gz"
                ## vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bovidae.cbgogwd.GATKug.filtered.vcf.gz"
                # nCores = NA, ## now specified elsewhere
                ## simpleRepeat_file <- "/bosTau8.simpleRepeat.gz"
                ## reference <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/ref/bosTau8.fa.gz"
                ## vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bovidae.cbgogwd.GATKug.filtered.vcf.gz"
                ## rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/ruminantia/bosTau8.rmsk.gz"
                chrlist = paste0("chr", 1:21),
                genomeSize = 2670044500,
                lineages_to_build = NULL,
                lineages = list(
                    "cow" = "cow",
                    "buffalo" = "buffalo",
                    "ACB" = c("cow", "buffalo"),
                    "goat" = "goat",
                    "ACBG" = c("ACB", "goat"),
                    "giraffe" = "giraffe",
                    "okapi" = "okapi",
                    "AOGi" = c("giraffe", "okapi")
                ),
                ancestral_lineage = list(
                    "AOGiACBG" = c("AOGi", "ACBG")
                ),
                # outgroups = c("reddeer", "whitetaileddeer")
                ## if no reddeer (yet)
                outgroups = c("whitetaileddeer")
            )
        )
    }
    if (species == "test") {
        return(
            list(
                lineages_to_build = NULL,
                ## have manually rsynced these over
                reference = file.path("ref/ref.fa.gz"),
                simpleRepeat_file = file.path("external/ref.simpleRepeat.gz"),
                rmask_file = file.path("external/ref.rmsk.gz"),
                chrlist = paste0(1:2),
                genomeSize = 4e5,
                lineages = list(
                    "test_pop1" = "test_pop1",
                    "test_pop2" = "test_pop2"
                ),
                ancestral_lineage = list(
                    "P1P2Anc" = c("test_pop1", "test_pop2")
                ),
                outgroups = c("test_outgroup"),
                vcf_load_split_num_files = 1,
                Klist = 6,
                cgte = 0,
                rgte = 0,
                n_extra_random_starts = 1, ## since slow
                max_iters_atToGC = 10, ## also since slow, who cares what parameters are! 
                use_gradient_and_hessian_for_ATGC_model_fitting = FALSE, ## disable as uses nlminb which is slow with these
                ancestral_map_window_size = 1e4,
                n_initial_atToGC_fitting_reps = 3,
                skip_at_to_gc_ci_fitting = TRUE
            )
        )
    }
    # if (species == "felidae") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/felidae/felis_catus_6.2.simpleRepeat.no_chr.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/felidae/from_Yuki_2018_01_24/Cats.HardFiltered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/felidae/Felis_catus.Felis_catus_6.2.dna.toplevel.fa.gz"
    #     chrlist = c("A1", "A2", "A3", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "F1", "F2")
    #     ## no chrX
    #     genomeSize =  2365745914
    #     rmask_file = "/data/smew1/rdavies/motifLossAnalysis/felidae/felis_catus_6.2.rmsk.nochr.gz"
    #     lineages = list(
    #         "ChineseShorttailed_FelisCatus6.2" = "ChineseShorttailed_FelisCatus6.2",
    #         "Sand_FelisCatus6.2" = "Sand_FelisCatus6.2",
    #         "ACS" = c("ChineseShorttailed_FelisCatus6.2", "Sand_FelisCatus6.2"),
    #         "BlackFooted_FelisCatus6.2" = "BlackFooted_FelisCatus6.2",
    #         "ACSB" = c("ACS", "BlackFooted_FelisCatus6.2"),
    #         "Jungle_FelisCatus6.2" = "Jungle_FelisCatus6.2",
    #         "ACSBJ" = c("ACSB", "Jungle_FelisCatus6.2"),
    #         "Asian_leo_FelisCatus6.2" = "Asian_leo_FelisCatus6.2",
    #         "ACSBJA" = c("ACSBJ", "Asian_leo_FelisCatus6.2"),
    #         "Cheetah_FelisCatus6.2" = "Cheetah_FelisCatus6.2",
    #         "ACSBJAC" = c("ACSBJA", "Cheetah_FelisCatus6.2"),
    #         "EurasianLynx_FelisCatus6.2" = "EurasianLynx_FelisCatus6.2",
    #         "ACSBJACE" = c("ACSBJAC", "EurasianLynx_FelisCatus6.2"),
    #         "Serval_FelisCatus6.2" = "Serval_FelisCatus6.2"
    #     )
    #     ancestral_lineage = list(
    #         "ACSBJACES" = c("ACSBJACE", "Serval_FelisCatus6.2")
    #     )
    #     outgroups = c("Leopard_FelisCatus6.2", "Tiger_FelisCatus6.2", "SnowLeo_FelisCatus6.2", "Jaguar_FelisCatus6.2")
    # }
    # if (species == "avian") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/avian/taeGut2.simpleRepeat.gz"
    #     rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/avian/taeGut2.rmask.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/avian/birds.dzlmargk.GATKug.filtered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/avian/taeGut2.fa.gz"
    #     chrlist <- paste0("chr", c('1', '1A', '1B', '2', '3', '4', '4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28'))
    #     genomeSize <- 1232135591
    #     lineages = list(
    #         "zebra_finch" = "zebra_finch",
    #         "long_tailed_finch" = "long_tailed_finch",
    #         "AZL" = c("zebra_finch", "long_tailed_finch"),
    #         "double_barrelled_finch" = "double_barrelled_finch",
    #         "AZLD" = c("AZL", "double_barrelled_finch"),
    #         "medium_ground_finch" = "medium_ground_finch",
    #         "AZLDM" = c("AZLD", "medium_ground_finch"),
    #         "american_crow" = "american_crow",
    #         "AZLDMA" = c("AZLDM", "american_crow"),
    #         "golden_collared_manakin" = "golden_collared_manakin"
    #     )
    #     ##  ['double_barrelled_finch', 'zebra_finch', 'long_tailed_finch', 'medium_ground_finch', 'american_crow', 'rifleman', 'golden_collared_manakin', 'kea']
    #     ancestral_lineage = list(
    #         "AZLDMAG" = c("AZLDMA", "golden_collared_manakin")
    #     )
    #     outgroups = c("rifleman", "kea")
    # }
    # if (species == "salmon") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/salmon/Otsh_v1.0.simpleRepeat.txt.gz"
    #     rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/salmon/Otsh_v1.0.rmask2.txt.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/salmon/salmon.aracc.GATKug.filtered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/salmon/Otsh_v1.0.fa.gz"
    #     chrlist <- c('NC_037097.1', 'NC_037098.1', 'NC_037099.1', 'NC_037100.1', 'NC_037101.1', 'NC_037102.1', 'NC_037103.1', 'NC_037104.1', 'NC_037105.1', 'NC_037106.1', 'NC_037107.1', 'NC_037108.1', 'NC_037109.1', 'NC_037110.1', 'NC_037111.1', 'NC_037112.1', 'NC_037113.1', 'NC_037114.1', 'NC_037115.1', 'NC_037116.1', 'NC_037117.1', 'NC_037118.1', 'NC_037119.1', 'NC_037120.1', 'NC_037121.1', 'NC_037122.1', 'NC_037123.1', 'NC_037124.1', 'NC_037125.1', 'NC_037126.1', 'NC_037127.1', 'NC_037128.1', 'NC_037129.1', 'NC_037130.1')
    #     genomeSize <- 2425697331
    #     lineages <- list(
    #         "coho" = "coho",
    #         "chinook" = "chinook",
    #         "ACC" = c("coho", "chinook"),
    #         "rainbow_trout" = "rainbow_trout",
    #         "ACCR" = c("ACC", "rainbow_trout"),
    #         "arctic_char" = "arctic_char"
    #     )
    #     ancestral_lineage = list(
    #         "ACCRA" = c("ACCR", "arctic_char")
    #     )
    #     outgroups = "atlantic_salmon"
    # }
    # if (species == "lizards") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/lizards/anoCar2.simpleRepeat.txt.gz"
    #     rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/lizards/anoCar2.rmsk.txt.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/lizards/lizards.acccegilosv.GATKug.filtered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/lizards/anoCar2.fa.gz"
    #     dict_file <- "/data/smew1/rdavies/motifLossAnalysis/lizards/anoCar2.dict"
    #     chrlist <- convert_dict_to_chrlist(dict_file, min_size = 1e6) 
    #     genomeSize <- 1799143587
    #     lineages <- list(
    #         "A_grahami" = "A_grahami",
    #         "A_lineatopus" = "A_lineatopus",
    #         "AGL" = c("A_grahami", "A_lineatopus"),
    #         "A_valencienni" = "A_valencienni",
    #         "AGLV" = c("AGL", "A_valencienni"),
    #         "A_sagrei" = "A_sagrei",
    #         "AGLVS" = c("AGLV", "A_sagrei"),
    #         "A_evermanni" = "A_evermanni",
    #         "A_cristatellus" = "A_cristatellus",
    #         "AEC" = c("A_evermanni", "A_cristatellus"),            
    #         "A_angusticeps" = "A_angusticeps",
    #         "A_insolitus" = "A_insolitus",
    #         "A_cybotes" = "A_cybotes"
    #     )
    #     ancestral_lineage = list(
    #         "AGLVSECAIC" = c("AGLVS", "AEC", "A_angusticeps", "A_insolitus", "A_cybotes")
    #     )
    #     outgroups = c("A_occultus", "A_chlorocyanus")
    # }
    # if (species == "bats") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/bats/myoLuc2.simpleRepeat.txt.gz"
    #     rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/bats/myoLuc2.rmsk.txt.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/bats/bat.drbb.GATKug.filtered.vcf.gz"
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/bats/myoLuc2.fa.gz"
    #     dict_file <- "/data/smew1/rdavies/motifLossAnalysis/bats/myoLuc2.dict"
    #     chrlist <- convert_dict_to_chrlist(dict_file, min_size = 1e6)
    #     genomeSize <- 2034575300
    #     lineages <- list(
    #         "red_bat" = "red_bat",
    #         "davids_bat" = "davids_bat",
    #         "ARD" = c("red_bat", "davids_bat"),
    #         "brandts_bat" = "brandts_bat"
    #     )
    #     ancestral_lineage = list(
    #         "ARDB" = c("ARD", "brandts_bat")
    #     )
    #     outgroups = c("big_brown_bat")
    # }
    # if (species == "canidae") {
    #     nCores <- 16
    #     simpleRepeat_file <- "/data/smew1/rdavies/motifLossAnalysis/canidae/canFam3.simpleRepeat.txt.gz"
    #     rmask_file <- "/data/smew1/rdavies/motifLossAnalysis/canidae/canFam3.rmsk.txt.gz"
    #     vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/canidae/canidae.aaacddggir.GATKug.filtered.vcf.gz"
    #     ## vcf_file <- "/data/smew1/rdavies/motifLossAnalysis/canidae/canidae.acdggir.GATKug.filtered.vcf.gz"        
    #     reference <- "/data/smew1/rdavies/motifLossAnalysis/canidae/canFam3.fa.gz"
    #     chrlist <- paste0("chr", 1:37)
    #     genomeSize <- 2410976875
    #     lineages = list(
    #         "coyote"              = "coyote",
    #         "dog"                 = "dog",
    #         "ACD"                 = c("coyote", "dog"),
    #         "african_golden_wolf" = "african_golden_wolf",
    #         "ACDA"                = c("ACD", "african_golden_wolf"),
    #         "golden_jackal"       = "golden_jackal",
    #         "ACDAG"               = c("ACDA", "golden_jackal"),
    #         "dhole"               = "dhole",
    #         "ACDAGD"              = c("ACDAG", "dhole"),
    #         "african_hunting_dog" = "african_hunting_dog",
    #         "ACDAGDA"             = c("ACDAGD", "african_hunting_dog"),
    #         "andean_fox"          = "andean_fox",
    #         "ACDAGDAA"            = c("ACDAGDA", "andean_fox"),
    #         "red_fox"             = "red_fox"
    #     )
    #     ancestral_lineage = list(
    #         "ACDAGDAAR"           = c("ACDAGDAA", "red_fox")
    #     )
    #     outgroups = c("grey_fox", "island_fox")
    # }
    # return(
    #     list(
    #         # nCores = nCores,
    #         simpleRepeat_file = simpleRepeat_file,
    #         rmask_file  = rmask_file,
    #         # vcf_file = vcf_file,
    #         reference = reference,
    #         chrlist = chrlist,
    #         genomeSize = genomeSize,
    #         lineages  = lineages,
    #         lineages_to_build = lineages_to_build,
    #         ancestral_lineage = ancestral_lineage,
    #         outgroups = outgroups,
    #         vcf_load_split_num_files = vcf_load_split_num_files,
    #         Klist = Klist,
    #         cgte = cgte,
    #         rgte = rgte,
    #         n_extra_random_starts = n_extra_random_starts, ## since slow
    #         max_iters_atToGC = max_iters_atToGC, ## also since slow, who cares what parameters are! 
    #         use_gradient_and_hessian_for_ATGC_model_fitting = use_gradient_and_hessian_for_ATGC_model_fitting, ## disable as uses nlminb which is slow with these
    #         ancestral_map_window_size = ancestral_map_window_size,
    #         n_initial_atToGC_fitting_reps = n_initial_atToGC_fitting_reps,
    #         skip_at_to_gc_ci_fitting = skip_at_to_gc_ci_fitting
    #     )
    # )
}
