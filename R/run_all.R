library("HATBAG")

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
print(getwd())
R_DIR <- args[1]
species_to_run <- args[2]
to_run <- args[3]
outputDate <- args[4]
nCores <- as.integer(args[5])
vcf_file <- file.path(getwd(), args[6])
HATBAG_DIR <- args[7]

source(file.path(R_DIR, "run_all_functions.R"))

analysisDir <- getwd()


## outputDate <- "2018_04_27" ## no missing allowed, mrle >= 6, etc
## outputDate <- "2018_05_22" ## allow missing in 
## outputDate <- "2018_05_29" ## no missing allowed, set hets to NA, otherwise same as 2018_04_27
## outputDate <- "2018_05_30" ## full run. 1 missing lineage allowed + only 1 outgroup required, set hets to NA, mrle >= 5, mncdn <= 3, use mrle as similar k-mer
## outputDate <- "2018_07_18" ## no update from previous, just to get up to speed after long absence
## outputDate <- "2018_11_21" ## test to get running again



## source("~/proj/motif_death/R/run_all_functions.R")
## "ruminantia"
##run <- "ABCDEF"
## setwd("/well/davies/users/dcc832/primates/")
##analysisDir <- "/data/smew1/rdavies/motifLossAnalysis/"
##
##species_to_run <- c("mice", "primates_nean", "felidae", "primates", "ruminantia", "avian", "salmon")
## species_to_run <- c("mice", "hominoidea", "cercopithecidae", "primates_nean", "felidae", "ruminantia", "avian", "salmon", "lizards", "bats", "canidae")
## ## lizards did not work
## ## try to debug, then run bats
## i_species_to_run <- 11
## to_run <- c(run, run, run, run, run, run, run, run, run, run, "CDEF")

i_species <- 1
i_species_to_run <- 1


for(i_species in i_species_to_run) {
    
    run <- to_run[i_species]
    species <- species_to_run[i_species]
    ## specify some params that get modified often
    out <- get_params(outputDate)
    num_non_missing_outgroups_required <- out$num_non_missing_outgroups_required
    num_missing_lineages_allowed <- out$num_missing_lineages_allowed
    mrle <- out$mrle
    mncdnle <- out$mncdnle
    similar_kmer_criterion <- out$similar_kmer_criterion
    use_one_sided_pvalue <- out$use_one_sided_pvalue
    ## haven't seen fit to modify in a while
    callable_bed <- NULL    
    ndge <- 3
    gcW2 <- 5000
    Klist <- 10
    ##
    message(paste0("Analyzing ", species))
    ##
    out <- get_per_species_params(species) 
    ## nCores <- out$nCores ## here for nCores keep the value set above
    simpleRepeat_file <- out$simpleRepeat_file
    rmask_file <- out$rmask_file
    ## vcf_file <- out$vcf_file
    reference <- out$reference
    chrlist <- out$chrlist
    genomeSize <- out$genomeSize
    lineages <- out$lineages
    ancestral_lineage <- out$ancestral_lineage
    outgroups <- out$outgroups
    lineages_to_build <- out$lineages_to_build
    ## TODO: make sure this aligns with HATBAG_OUTPUT_DIR specified in Snakefile_reference_order
    outputDir <- file.path(analysisDir, "hatbag", species, "/")
    message(outputDir)
    
    ##
    if (1 == 1) {
        
        masterDirHDD <- file.path(outputDir, outputDate)
        dir.create(masterDirHDD, recursive = TRUE)
        system(paste0("cd ", HATBAG_DIR," && git log | head -n10 > ", outputDir, "/", outputDate, "/HATBAG_head_git_log.txt"))
        ##

        ## argh - hack for now - using too much ram
        if (run == "E") {
            nCores <- max(c(1, floor(nCores / 2)))
            print(paste0("Lowering nCores to ", nCores))
        }
        ## use fewer cores 
        ## 
        
        out <- HATBAG(
            species = species,
            run = run,
            outputDate = outputDate,
            outputDir = outputDir,
            vcf_file = vcf_file,
            reference = reference,
            chrlist = chrlist,
            genomeSize = genomeSize,
            rmask_file = rmask_file,
            # Klist = Klist,
            nCores = nCores,
            lineages = lineages,
            lineages_to_build = lineages_to_build,
            ancestral_lineage = ancestral_lineage,
            outgroups = outgroups,
            similar_kmer_criterion = similar_kmer_criterion,
            callable_bed = callable_bed,
            use_one_sided_pvalue = use_one_sided_pvalue,
            simpleRepeat_file = simpleRepeat_file,
            ndge = ndge,
            mrle = mrle,
            mncdnle = mncdnle,
            gcW2 = gcW2,
            num_non_missing_outgroups_required = num_non_missing_outgroups_required,
            num_missing_lineages_allowed = num_missing_lineages_allowed,
            vcf_load_split_num_files = out$vcf_load_split_num_files,
            Klist = out$Klist,
            cgte = out$cgte,
            rgte = out$rgte,
            n_extra_random_starts = out$n_extra_random_starts,
            max_iters_atToGC = out$max_iters_atToGC,
            use_gradient_and_hessian_for_ATGC_model_fitting = out$use_gradient_and_hessian_for_ATGC_model_fitting,
            ancestral_map_window_size = out$ancestral_map_window_size,
            n_initial_atToGC_fitting_reps = out$n_initial_atToGC_fitting_reps
        )

    }
    for(i in 1:10) {
        gc(reset = TRUE)
    }
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
