#!/usr/bin/env Rscript

library("HATBAG")

R_DIR <- Sys.getenv("R_DIR")
source(file.path(R_DIR, "run_all_functions.R"))
args <- commandArgs(trailingOnly = TRUE)

if (1 == 0) {

    setwd("/well/davies/users/dcc832/primates/")
    args <- c(
        "~/proj/motif_death/R/",
        "artiodactyla",
        "A",
        "2021_05_26",
        1,
        "vcf/bovidae.cbgogwd.GATKug.filtered.vcf.gz"
    )

}

print(args)
analysisDir <- getwd()

print(analysisDir)
R_DIR <- args[1]
species_order <- args[2]
to_run <- args[3]
outputDate <- args[4]
nCores <- as.integer(args[5])
vcf_file <- file.path(getwd(), args[6])
HATBAG_DIR <- args[7]
outputDir <- args[8]
callable_bed <- args[9]

HATBAG_params_1 <- get_params(outputDate)
HATBAG_params_2 <- get_per_species_params(species_order)
ndge <- 3
gcW2 <- 5000
# Klist <- 10 # 6 for test

message(paste0("HATBAG analyzing ", species_order))
message(paste0("HATBAG output in ", outputDir))

masterDirHDD <- file.path(outputDir, outputDate)
dir.create(masterDirHDD, recursive = TRUE)
system(paste0("cd ", HATBAG_DIR," && git log | head -n10 > ", outputDir, "/", outputDate, "/HATBAG_head_git_log.txt"))

# argh - hack for now - using too much ram
if (to_run == "E") {
    nCores <- max(c(1, floor(nCores / 2)))
    print(paste0("Lowering nCores to ", nCores))
}

# TODO: use fewer cores 
out <- do.call(
    HATBAG, 
    c(list(species=species_order, run=to_run, outputDate=outputDate, outputDir=outputDir, vcf_file=vcf_file, nCores=nCores, callable_bed=callable_bed, ndge=ndge, gcW2=gcW2),
    HATBAG_params_1,
    HATBAG_params_2
    )
)

for(i in 1:10) {
    gc(reset = TRUE)
}

quit()
    









### SWAP IN
##

    


    #############################TEMP
    if (i_species == 1) {
        similar_kmer_criterion <- c("cpg", "at")
        use_one_sided_pvalue <- FALSE
        outputDir <- paste0(analysisDir, species, "-without-mrle-two-sided-p/")
    } else if (i_species == 2) {
        similar_kmer_criterion <- c("cpg", "at")
        use_one_sided_pvalue <- TRUE
        outputDir <- paste0(analysisDir, species, "-without-mrle-one-sided-p/")
    } else if (i_species == 3) {
        similar_kmer_criterion <- c("cpg", "at", "mrle")
        use_one_sided_pvalue <- FALSE
        outputDir <- paste0(analysisDir, species, "-with-mrle-two-sided-p/")
    } else if (i_species == 4) {
        similar_kmer_criterion <- c("cpg", "at", "mrle")
        use_one_sided_pvalue <- TRUE
        outputDir <- paste0(analysisDir, species, "-with-mrle-one-sided-p/")
    }
    #############################TEMP    


## THIS CODE CAN BE USED TO REGENERATE QQ PLOTS

## -- temporary to QQ plotting
setwd(file.path(HATBAG_DIR, "HATBAG", "R"))
sapply(dir()[-grep("~", dir())], source)
masterDirHDD <- file.path(outputDir, outputDate)    
testingNames = c("losslin", "lossat", "gainlin", "gainat")
plottingNames = c("Loss Lineage","Loss AT to GC","Gain Lineage","Gain AT to GC")
load(snp_counts_filename(masterDirHDD))
load(mask_info_filename(masterDirHDD))

analyze_QQ_data(
    masterDirHDD,
    Klist = 10,
    testingNames,
    lineageNames,
    nLin,
    plottingNames,
    repeatNamesAll,
    species
)
## end of temporary for QQ plotting



## debug ancestral map plots


setwd(file.path(HATBAG_DIR, "HATBAG", "R"))
sapply(dir()[-grep("~", dir())], source)
library("testthat"); library("HATBAG")

make_and_plot_one_map_per_order(
    masterDirHDD,
    K,
    testingNames
)
system("rsync -av /data/smew1/rdavies/motifLossAnalysis//primates_nean/2018_03_21/E_cluster/ancestralMap.*.averaged.png florence:~/")

analysisDir <- "/data/smew1/rdavies/motifLossAnalysis/"
outputDate <- "2018_03_21"
species <- "primates_nean"
masterDirHDD <- file.path(analysisDir, species, outputDate)

test <- "losslin"
iLin <- 3
repeatName <- "nonRepeat"
motifNumber <- 1
chrlist <- paste0("chr", c(1:22))
        plot_an_ancestral_map(
            masterDirHDD = masterDirHDD,
            test = test,
            iLinDiscovery = iLin,
            iLinCluster = iLin,
            repeatName = repeatName,
            motifNumber = motifNumber,
            chrlist = chrlist
        )


o <- lapply(c(3, 6), function(iLin) {
    load(file = ancestral_map_filename_no_chr(masterDirHDD, test, iLin, iLin, repeatName, motifNumber))
    return(cbind(
    valL_new = unlist(sapply(ancestral_map, function(x) x[, "valL_new"])),
    valG_new = unlist(sapply(ancestral_map, function(x) x[, "valG_new"])),
    valAL_new = unlist(sapply(ancestral_map, function(x) x[, "valAL_new"])),
    valAG_new = unlist(sapply(ancestral_map, function(x) x[, "valAG_new"]))
    ))
})


o4 <- o[[1]]
colnames(o4) <- paste0(colnames(o4), ".AHN")
o6 <- o[[2]]
colnames(o6) <- paste0(colnames(o6), ".gorilla")

m <- cbind(o4, o6)
output <- array(NA, c(8, 8))
rownames(output) <- colnames(m)
colnames(output) <- colnames(m)
for(i in 1:8) {
    for(j in 1:8) {
        output[i, j] <- round(cor(m[, i], m[, j], use = "pairwise.complete.obs") ** 2, 3)
    }
}

## for AHN, semi-OK correlation between loss, AT to GC signal
## what is comparison like between maps
## what is comparison 
## very little correlation between gorilla, chimp loss map?
output
cor(valL_new, valAL_new, use = "pairwise.complete.obs") ** 2 ## better

colnames(ancestral_map[[1]])



## one per order
## ?average?

K <- 10
load(motif_results_dataframe_filename(masterDirHDD, K) )
mrd <- motif_results_dataframe
mrd_test <- mrd[mrd[, 1] == "losslin" & mrd[, "repeatName"] == "nonRepeat", ]
iLinDiscovery <- 3
iLinCluster <- 3
repeatName <- "nonRepeat"
motifNumber <- 1

## OK I think I can normalize across ALL species
## OK, start by getting all ORs
## then, for each, load up
## then make 1 plot per order, for discovery using "losslin", etc



for(motifNumber in 1:7) {
    ## more noise? 
    load(file = ancestral_map_filename_no_chr(masterDirHDD, test, iLinDiscovery, iLinCluster, repeatName, motifNumber))
    print(motifNumber)
    print(range(ancestral_map[[1]][, "valL_new"], na.rm = TRUE))
    print(range(ancestral_map[[1]][, "valG_new"], na.rm = TRUE))    
}




    chr,
    motif_results_dataframe,
    iLin,
    test,
    masterDirHDD,
    K,
    gcW2,
    tmpdir,
    motifSuperResults_all,
    nLin,
    ancNames_for_ref_building,
    lineageNames_single
