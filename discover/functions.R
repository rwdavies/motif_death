library("ape")
library("rotl")
library("rentrez")
library("taxize")
library("jsonlite")
library("parallel")
library("data.table")
library("rotl")

get_match_against_subfamilies <- function(subfamilies, genuses = NULL) {
    fields <- "&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created"
    if (is.null(genuses)) {
        genuses <- unique(unlist(lapply(subfamilies, function(subfamily) {
            lots_of_snakes <- taxonomy_subtree(tnrs_match_names(subfamily)[,"ott_id"])
            possibilities <- unique(sapply(strsplit(lots_of_snakes$edge_label, "_"), function(x) x[1]))
            genuses <- possibilities[-grep("(", possibilities, fixed = TRUE)]
            genuses <- genuses[-grep("'", genuses, fixed = TRUE)]
            return(genuses)
        })))
    }
    if (length(genuses) == 0) {
        message("Nothing found")
        return(NULL)
    }
    super_results <- mclapply(1:length(genuses), mc.cores = 48, function(iGenus) {
        ##
        genus <- genuses[iGenus]
        url <- paste0("https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-search/", genus, "?limit=1000")
        file <- tempfile()
        x <- system(paste0("curl -s ", url, " > ", file), intern = TRUE)
        x <- readLines(file, n = 1)
        if(x == "No results.") {
            unlink(file)
            return(NULL)
        }
        available <- fromJSON(file)
        unlink(file)
        message(paste0(genus, " has ", nrow(available)))
        ##
        ##
        ##
        iRow <- 1
        url <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree(", available[iRow, "taxId"], ")", fields)
        file <- tempfile()
        x <- system(paste0("curl -s ", shQuote(url), " > ", file), intern = TRUE)
        if (file.info(file)["size"] == 0) {
            return(NULL)
        }
        results <- read.table(file, sep = "\t", header = TRUE, comment.char="", quote = "")
        unlink(file)
        remove1 <- results[, "library_source"] %in% c("TRANSCRIPTOMIC", "METAGENOMIC")
        remove2 <- results[, "library_selection"] %in% c("Reduced Representation", "Restriction Digest", "Hybrid Selection", "ChIP", "DNase")
        remove3 <- results[, "library_strategy"] %in% "Hi-C"
        results <- results[!remove1 & !remove2 & !remove3, ]
        ## keep most things?
        ## remove3 <- results[, "instrument_model"]
        ## results <- results[grep("HiSeq", results[ ,"instrument_model"]), ]
        results <- results[results[ ,"read_count"] > 10e6, ]  ## honestly, prob too small?
        results
    })
    results <- data.frame(rbindlist(lapply(super_results[!sapply(super_results, is.null)], function(x) x)))
    ## remove duplicates?
    results <- unique(results)
    ## remove weird NA
    results <- results[!is.na(results[, "study_accession"]), ]
    ## ALSO - run accession should be unique, so just take one
    ## table(tapply(results[, "base_count"], INDEX = factor(results[, "run_accession"]), function(x) length(unique(x))))
    results <- results[match(unique(results[, "run_accession"]), results[, "run_accession"]), ]
    print(table(results[, "scientific_name"]))
    results
}



get_ucsc_assemblies <- function(motif_death_dir = "~/proj/motif_death/") {
    ##
    ## all assemlies
    ##
    ##system("cd ~/Downloads && curl -s https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt > UCSC_GI.assemblyHubList.txt")
    assemblies <- read.table(file.path(motif_death_dir, "discover/resources/UCSC_GI.assemblyHubList.txt"), skip = 11, sep = "\t", comment.char="@", header = TRUE, quote = "")
    assemblies[, "scientific.name"] <- gsub(" ", "_", assemblies[, "scientific.name"])
    assemblies$class2 <- NA
    ##
    ## annotate them some more
    ##
    for(what in c("primates", "mammals", "vgp")) {
        ## manually downloaded from e.g. https://hgdownload.soe.ucsc.edu/hubs/mammals/index.html
        ## copy and paste into txt file, fix header by removing line breaks
        primates <- read.table(file.path(motif_death_dir, paste0("discover/resources/UCSC", what, ".txt")), sep = "\t", quote = "", header = TRUE)
        primates$accession2 <- sapply(strsplit(primates[, "NCBI.assembly"], "_"), function(x) paste0(x[1], "_", x[2]))
        m <- match(primates$accession2, assemblies[, "X..accession"])
        assemblies[m, "class"] <- what
        if (what == "vgp") {
            assemblies[m[!is.na(m)], "class2"] <- primates[!is.na(m), "class.VGP.link"]
        }
    }
    return(assemblies)
}


investigate <- function(keyword) {
    ##
    ##
    a <- rotl::tnrs_match_names(keyword) ##
    b <- taxonomy_subtree(ott_id = a[, "ott_id"])
    result_phylo <- taxonomy_subtree(ott_id = a[, "ott_id"], output_format="phylo")
    ##
    ## get data
    ##
    subfamilies <- sapply(strsplit(b$edge_label, "_"), function(x) x[[1]])
    f <- function(wer, subfamilies) {
        x <- grep(wer, subfamilies, fixed = TRUE)
        if (length(x) > 0) {
            subfamilies <- subfamilies[-x]
        }
        subfamilies
    }
    subfamilies <- f("ott", subfamilies)
    subfamilies <- f("(", subfamilies)
    subfamilies <- unique(subfamilies)
    message(paste0("There are ", length(subfamilies), " subfamilies to investigate"))
    results <- get_match_against_subfamilies(genuses = subfamilies)
    ##
    ## get some info for the plots
    ##
    ##
    ## do plot
    ##
    make_informative_plot(results, result_phylo, keyword)
    return(
        list(
            results = results,
            result_phylo = result_phylo
        )
    )
}

make_informative_plot <- function(results, result_phylo, keyword, plotdir = "~/Downloads/") {
    ##
    ## get some information to plot about them
    ##
    mapped_per_species <- tapply(X = results[, "base_count"], INDEX = as.factor(results[, "scientific_name"]), FUN = sum, na.rm =TRUE) / 1e9
    N_per_species <- tapply(X = results[, "base_count"], INDEX = as.factor(results[, "scientific_name"]), FUN = length)
    stopifnot(sum(!(names(N_per_species) == names(mapped_per_species))) == 0) ## check they are done in the same way
    ##
    ## now, ideally, want to take the above newick tree
    ##   - colour reference genomes red
    ##   - colour available genomes blue
    ##   - possibly remove everything else?
    ##
    tree <- result_phylo
    tip.color <- rep("black", length(tree$tip.label))
    tip.label <- tree$tip.label
    keep_tip <- rep(FALSE, length(tip.label))
    ##
    ## do references here?
    ##
    for(i_what in 1:2) {
        ## so with the two versions, can do with, without subspecies
        to_compare <- sapply(strsplit(tip.label, "_"), function(x) {
            x <- x[!(substr(x, 1, 3) == "ott")]
            paste0(x[1:(i_what + 1)], collapse = "_")
        })
        species_have_refs <- which(!is.na(match(to_compare, assemblies[, "scientific.name"])))
        for(i_species in species_have_refs) {
            keep_tip[i_species] <- TRUE
            x <- assemblies[match(to_compare[i_species], assemblies[, "scientific.name"]), ]
            ## want the class
            tip.label[i_species] <- paste0(
                tip.label[i_species], " (REF, class=",x[, "class"], ")"
            )
            tip.color[i_species] <- "red"
        }
    }
    ##
    ## do found species
    ##
    x <- gsub(" ", "_", names(mapped_per_species))
    for(i_species in 1:length(mapped_per_species)) {
        species_name<- x[i_species]
        i_match <- grep(species_name, tip.label)[1] ## just keep first?
        if (is.na(i_match)) {
            warning(paste0("Could not match ", species_name))
        } else {
            keep_tip[i_match] <- TRUE
            if (tip.color[i_match] == "red") {
                tip.color[i_match] <- "green"
            } else {
                tip.color[i_match] <- "blue"
            }
            ## match uniquely?
            tip.label[i_match] <- paste0(
                tip.label[i_match],
                ": (N=", N_per_species[i_species], ", base_count(G)=",
                mapped_per_species[i_species], ")"
            )
        }
    }
    ##
    ## drop unnecessary things
    ##
    tree$tip.label <- tip.label
    tree2 <- drop.tip(tree, tip.label[!keep_tip])
    tip.color <- tip.color[match(tree2$tip.label, tree$tip.label)]
    height <- max(10, 3 * length(tree$tip.label) / 100)
    pdf(file.path(plotdir, paste0(keyword, ".pdf")), height = height, width = 20)
    plot(tree2, tip.color = "white", show.node.label = TRUE, show.tip.label = TRUE, srt = 90)
    par(new = TRUE)
    plot(tree2, tip.color = tip.color, show.node.label = FALSE)
    dev.off()
}



look_at_one_species_or_study <- function(results, study = NA, scientific_name = NA) {
    cols <- c("scientific_name", "instrument_platform", "instrument_model", "library_source", "library_selection", "fastq_bytes", "nominal_length", "first_created", "read_count", "study_accession", "run_accession", "library_strategy", "run_accession", "base_count", "experiment_title", "sample_alias")
    w1 <- rep(TRUE, nrow(results))
    w2 <- rep(TRUE, nrow(results))
    if (!is.na(study)) {
        w1 <- results[, "study_accession"] == study
    }
    if (!is.na(scientific_name)) {
        w2[] <- FALSE
        x <- grep(scientific_name, results[, "scientific_name"])
        if (length(x) > 0) {
            w2[x] <- TRUE
        } else {
            stop("Could not find that name")
        }
    }
    m <- results[w1 & w2, ]
    m <- m[order(-m[, "base_count"]), ]
    m <- cbind(m, run_total = cumsum(m[, "base_count"]) / 6e10)
    m[, c(cols, "run_total")]
}
