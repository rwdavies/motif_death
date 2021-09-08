library("ape")
library("rotl")
library("rentrez")
library("taxize")
## devtools::install_github("GuangchuangYu/treeio")


system("cd ~/Downloads && curl -O http://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.phast-242.nh")
original_tree <- scan("~/Downloads/241-mammalian-2020v2.phast-242.nh", what = "character")
tree <- original_tree

## try to get latin names
x <- unlist(strsplit(tree, "_"))

start <- gsub("(", "", sapply(strsplit(x, ":"), head, n = 1), fixed = TRUE)
end <- gsub("(", "", sapply(strsplit(x, ",", fixed = TRUE), tail, n = 1), fixed = TRUE)
species <- cbind(end[-length(end)], start[-1])
## manually fix dog
## Warning messages:
## 1: Some names were duplicated: ‘canis lupus’, ‘ceratotherium simum’.
## 2: lupus familiaris, simum cottoni are not matched
## ugh - whatever - fix manually
species[200, ] <- c("Canus", "lupus_familiaris")
species <- species[-201, ]
species[225, ] <- c("Ceratotherium", "simum_cottoni")
species <- species[-226, ]


##
## get rotl stuff
##
matched <- rotl::tnrs_match_names(paste0(species[, 1], " ", species[, 2]))

lineages <- taxonomy_taxon_info(matched[, "ott_id"], include_lineage=TRUE)

## get common name, etc?
## synonyms <- unlist(sapply(lineages, function(x) {
##     if ("synonyms" %in% names(x)) {
##         if (length(x[["synonyms"]]) > 0) {
##             x[["synonyms"]][[1]]
##         }
##     } else {
##         NA
##     }
## }))


##
##
##
more_details <- t(sapply(1:length(lineages), function(iSpecies) {
    x <- t(sapply(lineages[[iSpecies]][["lineage"]], I))
    w <- match(c("genus", "family", "suborder", "order") , unlist(x[, "rank"]))
    to_out <- array("NA", 4)
    to_out[!is.na(w)] <- unlist(x[w, c("unique_name")])
    sapply(strsplit(to_out, " ", fixed = TRUE), head, n = 1)
}))
## remove () stuff


##
## get taxize stuff - SLOOOW - API call one at a time!
##
common_names <- sci2comm(paste0(species[, 1], " ", species[, 2]))



##
## can I query for ref genomes?
##

entrez_db_summary("assembly")

entrez_db_searchable("assembly")

##
## get "best" match for each
##
assembly_info <- parallel::mclapply(paste0(species[, 1], " ", species[, 2]), mc.cores = 1, function(species_name) {
    print(species_name)
    x <- entrez_search(db = "assembly", term = species_name)
    m <- t(sapply(x$ids, function(id) {
        y <- entrez_summary(db = "assembly", id = id)
        y[c("speciesname", "assemblyname", "assemblyclass", "coverage", "contign50", "scaffoldn50")]
    }))
    if (length(m) == 0) {
        return(rep(NA, 6))
    }
    m[which.max(m[, "scaffoldn50"]), ]
})


assembly_info <- t(sapply(assembly_info, I)) ## I think

assembly_qual <- paste0(round(as.numeric(assembly_info[, "contign50"]) / 1e3), "_", round(as.numeric(assembly_info[, "scaffoldn50"]) / 1e3, 1))

## re-label names!
tree <- original_tree
suffixes <- apply(more_details, 1, paste0, collapse = "_")
for(i in 1:nrow(species)) {
    species_name <- paste0(species[i, 1], "_", species[i, 2])
    common_name <- gsub(" ", "_", common_names[[i]])
    tree <- gsub(species_name, paste0(assembly_qual[i], "_", common_name, "_", species_name, "_", suffixes[i]), tree)
}

pdf("~/Downloads/temp.pdf", height = 100, width = 25)
plot(ape::read.tree(text = tree))
dev.off()







##
## all assemlies
##
system("cd ~/Downloads && curl -s https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt > UCSC_GI.assemblyHubList.txt")
assemblies <- read.table("~/Downloads/UCSC_GI.assemblyHubList.txt", skip = 11, sep = "\t", comment.char="@", header = TRUE, quote = "")
assemblies[, "scientific.name"] <- gsub(" ", "_", assemblies[, "scientific.name"])



## alt version
for(what in c("primates", "mammals", "vgp")) {
    ## manually downloaded from e.g. https://hgdownload.soe.ucsc.edu/hubs/mammals/index.html
    ## copy and paste into txt file, fix header by removing line breaks
    primates <- read.table(paste0("~/Downloads/UCSC", what, ".txt"), sep = "\t", quote = "", header = TRUE)
    primates$accession2 <- sapply(strsplit(primates[, "NCBI.assembly"], "_"), function(x) paste0(x[1], "_", x[2]))
    assemblies[match(primates$accession2, assemblies[, "X..accession"]), "class"] <- what
}

table(assemblies[, "class"])






##
## Feliformia e.g. cats, meerkats
##
a <- rotl::tnrs_match_names("Herpestidae") ##
b <- taxonomy_subtree(ott_id = a[, "ott_id"])
b2 <- taxonomy_subtree(ott_id = a[, "ott_id"], output_format="phylo")
##
## get data
##
subfamilies <- sapply(strsplit(b$edge_label, "_"), function(x) x[[1]])
subfamilies <- subfamilies[-grep("ott", subfamilies)]
subfamilies <- subfamilies[-grep("(", subfamilies, fixed = TRUE)]
results <- get_match_against_subfamilies(genuses = subfamilies)

mapped_per_species <- tapply(X = results[, "base_count"], INDEX = as.factor(results[, "scientific_name"]), FUN = sum, na.rm =TRUE) / 1e9
N_per_species <- tapply(X = results[, "base_count"], INDEX = as.factor(results[, "scientific_name"]), FUN = length)
stopifnot(sum(!(names(N_per_species) == names(mapped_per_species))) == 0) ## check they are done in the same way

## now, as a start
## get how many Gbp exist, for that species

##
## now, ideally, want to take the above newick tree
##   - colour reference genomes red
##   - colour available genomes blue
##   - possibly remove everything else?
##
tree <- b2
tip.color <- rep("black", length(tree$tip.label))
tip.label <- tree$tip.label
keep_tip <- rep(FALSE, length(tip.label))
##
## do references here?
##
to_compare <- sapply(strsplit(tip.label, "_"), function(x) {
    paste0(x[1:2], collapse = "_")
})
species_have_refs <- which(!is.na(match(to_compare, assemblies[, "scientific.name"])))
for(i_species in species_have_refs) {
    keep_tip[i_species] <- TRUE
    ## AM HERE
    ## FINISH ADDING THESE TO THE PLOT
    ## CAN I ALSO LOOK AT QUALITY SOMEHOW?
    ## THEN MAKE THEM RED, THE REST BLUE
    x <- assemblies[match(to_compare[i_species], assemblies[, "scientific.name"]), ]
    ## want the class
    tip.label[i_species] <- paste0(
        tip.label[i_species], " (REF, class=",x[, "class"], ")"
    )
    tip.color[i_species] <- "red"
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
pdf("~/Downloads/feliformia.pdf", height = 20, width = 15)
plot(tree2, tip.color = tip.color)
dev.off()


