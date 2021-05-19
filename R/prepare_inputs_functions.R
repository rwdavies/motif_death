## sample_alias <- info[[i_species]][1]
## sample_name <- info[[i_species]][3]
## paper <- info[[i_species]][4]
## n_mapping_pieces <- as.integer(info[[i_species]][5])

## canids / foxes
get_canid_info <- function() {
## x <- read.table("PRJNA494815 (multiple canidae).txt", sep = "\t", header = TRUE)
## x[order(x[, "sra_bytes"]), c("sample_alias", "sample_title", "fastq_bytes", "sample_accession", "scientific_name")]
## x <- read.table("SRA307300 (dhole).txt", header = TRUE, sep = "\t")
## y <- read.table("PRJNA312115 (also urocyon).txt", header = TRUE, sep = "\t")
paper1 <- "PRJNA494815 (multiple canidae).txt"
## 
## gray fox,      GrayFox_GOGANRA_NPS_GF30           NA 18308882706;20098897733
## SAMN09516312 Urocyon cinereoargenteus
info <- list(
    c("Reef", "Vulpes vulpes", "red_fox", "PRJNA378561 (red fox).txt", 50),
    c("Lcu2_Pastora", "Lycalopex culpaeus", "andean_fox", "SRS523207 (andean fox).txt", 50),
    c("Sister 1", "Lycaon pictus", "african_hunting_dog", "PRJNA488046 (african wild dog).txt", 50),
    c("RUFZCHN00001", "Cuon alpinus", "dhole", "SRA307300 (dhole).txt", 50)    ,
    c("GoldenJackal_Syria", "Canis aureus", "golden_jackal", paper1, 50),
    c("AfricanGoldenWolf_Ethiopia", "Canis anthus", "african_golden_wolf", paper1, 50),
    c("cac", "Canis latrans", "coyote", "SRS661477 (coyote).txt", 50),
    c("GFO41F", "Urocyon cinereoargenteus", "grey_fox", "PRJNA312115 (also urocyon).txt", 50),
    c("SMI15F", "Urocyon littoralis", "island_fox", "PRJNA312115 (also urocyon).txt", 50),
    c("SAMD00009664", "Canis lupus familiaris", "dog", "PRJDB2266 (korean dog).txt", 50)
)
    return(info)
}


get_bird_info <- function() {
## birds
n_mapping_pieces <- 50
paper1 <- "PRJEB10586 (finch recomb paper).txt"
paper2 <- "PRJNA156703 (medium ground finch).txt"
paper3 <- "PRJNA212872 (golden collared manakin).txt"
paper4 <- "PRJNA212877 (rifleman).txt"
paper5 <- "PRJNA212900 (kea).txt"
paper6 <- "PRJNA212869 (american crow).txt"
info <- list(
    c("DBF", "Taeniopygia bichenovii", "double_barrelled_finch", paper1, 50),
    c("26881", "Taeniopygia guttata", "zebra_finch", paper1, 50),
    c("73958", "Poephila acuticauda hecki", "long_tailed_finch", paper1, 50),
    c("Geospiza_fortis reads", "Geospiza fortis", "medium_ground_finch", paper2, 50),
    c("BGI_N305", "Manacus vitellinus", "golden_collared_manakin", paper3, 50),
    c("BGI_N310", "Acanthisitta chloris", "rifleman", paper4, 50),
    c("BGI_N333", "Nestor notabilis", "kea", paper5, 50),
    c("BGI_N302", "Corvus brachyrhynchos", "american_crow", paper6, 50)
)
    return(info)
}


get_salmon_info <- function() {
## salmon
info <- list(
    c("OkisDH3", "Oncorhynchus kisutch", "coho", "PRJNA352719 (coho).txt", 100),
    c("Swanson DH", "Oncorhynchus mykiss", "rainbow_trout", "PRJNA335610 (rainbow trout).txt", 300),
    c("Chilliwack Gynogen DE9421", "Oncorhynchus tshawytscha", "chinook", "PRJNA416144 (chinook).txt", 100),
    c("Salp-Liver_IW2", "Salvelinus alpinus", "arctic_char", "PRJNA348349 (arctic char).txt", 100),
    c("1133201704016 Sally gDNA, muscle SAL280BP_PE_IL47-CWRU Sample", "Salmo salar", "atlantic_salmon", "PRJNA72713 (atlantic salmon).txt", 100)
)
## note - can get away with ~23 GB / 24 hours with 100 I think
## ideally would be more like ~12 hours
    return(info)
}

get_bat_info <- function() {
## bats
info <- list(
    c("PRJNA167910.BU_THK_EF1", NA, "big_brown_bat", "PRJNA72449 (big brown bat).txt", 100),
    c("Myotis brandtii WGS", NA, "brandts_bat", "PRJNA178678. (brandts).txt"),
    c("Myotis davidii", NA, "davids_bat", "PRJNA171994 (davids).txt"),
    c("Myotis rufoniger, Gosudonggul cave, Danyang, in South Korea", NA, "red_bat", "PRJNA381377 (red).txt")
)
##    c(NA, NA, "little_brown", "PRJNA277738 (little brown bat).txt"),
    return(info)
}


get_bear_info <- function() {
    paper1 <- "PRJNA169236 (3 bears).txt"
    paper2 <- "SRP064940 (giant panda).txt"
    info <- list(
        c("PB7 N7773", NA, "polar_bear", paper1),
        c("GRZ 100", NA, "grizzly_bear", paper1),
        c("ABC2 Lucky", NA, "black_bear", paper1),
        c("Generic sample from Ailuropoda melanoleuca", NA, "giant_panda", paper2)
    )
    return(info)
}

get_b <- function(paper, sample_alias) {
    b <- read.table(paper,  header = TRUE, sep = "\t")
    b <- b[(b[, "sample_alias"] == sample_alias), ]
    if (sample_alias == "Geospiza_fortis reads") {
        b <- b[
            (b[, "nominal_length"] >= 200) &
            (b[, "nominal_length"] <= 500) , ]
    } else if (sample_alias == "OkisDH3") {
        b <- b[grep("PE400", b[, "experiment_alias"]), ]
        b[, "nominal_length"] <- 400
    } else if (sample_alias == "Swanson DH") {
        b <- b[grep("Swanson_450bp", b[, "experiment_alias"]), ]
        b[, "nominal_length"] <- 450
    } else if (sample_alias == "Chilliwack Gynogen DE9421") {
        b <- b[grep("HI.4075", b[, "experiment_alias"]), ]
        b[, "nominal_length"] <- 250 ## unclear
    } else if (sample_alias == "Salp-Liver_IW2") {
        b <- b[grep("Illumina HiSeq 2000 sequencing; WGS ", b[, "experiment_title"]), ]
        b[, "nominal_length"] <- 250 ## unclear
    } else if (sample_alias == "1133201704016 Sally gDNA, muscle SAL280BP_PE_IL47-CWRU Sample") {
        b <- b[b[, "library_strategy"] == "WGS", ]
        ## hmm, just take hiseq
        b <- b[grep("Illumina HiSeq", b[, "experiment_title"]), ]        
        b <- b[b[, "nominal_length"] < 1000 , ]
    } else if (sample_alias == "PRJNA167910.BU_THK_EF1") {
        b <- b[b[, "nominal_length"] == 180, ]
        b <- b[b[, "sra_bytes"] > 33000000000, ]
        ## honestly, 2 should be enough
        b <- b[b[, "run_accession"] %in% c("SRR363679", "SRR363680"), ]
    } else if (sample_alias == "Myotis rufoniger, Gosudonggul cave, Danyang, in South Korea") {
        b <- b
    } else if (sample_alias == "Reef") {
        read_length <- b[, "base_count"] / b[, "read_count"] / 2
        b <- b[read_length == 150, ]
    } else if (sample_alias == "SAMD00009664") {
        b[, "library_name"] <- "lib1"
    } else if (sample_alias == "PB7 N7773") {
        b <- b[b[, "library_name"] == "232C", ]
        b <- b[b[, "run_accession"] == "SRR518683", ]
    } else if (sample_alias == "GRZ 100") {
        b <- b[b[, "run_accession"] == "SRR518713", ]
    } else if (sample_alias == "Generic sample from Ailuropoda melanoleuca") {
        ## b <- b[b[, "run_accession"] == "SRR518713", ]
    } else {
        b <- b[is.na(b[, "nominal_length"]) | b[, "nominal_length"] < 1000 , ]
    }
    b$read_length <- b$base_count / b$read_count / 2
    return(b)
}


