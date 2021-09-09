source("~/proj/motif_death/discover/functions.R")

assemblies <- get_ucsc_assemblies()

out <- investigate("Afrotheria")
##
results <- out$results

out <- investigate("Pholidota")

## manis looks promising?

## I think so? use pentadactyla as the outgroup, build the tree with the other two

chrom <- read.table("https://hgdownload.soe.ucsc.edu/goldenPath/manPen1//bigZips/manPen1.chrom.sizes")
## UGH, I dunno, too too many chromosomes to run, doesn't seem reasonable for now

out <- investigate("Testudines") ## turtles


out <- investigate("Emydidae") ## more specific turtles



## AM HERE
## FINISH UP TURTLES

## turtles
m <- look_at_one_species_or_study(results, scientific_name="Chrysemys picta")

Painted turtle, Chrysemys picta (outgroup) PRJNA589899, SRR11059563
Diamondback terrapin, Malaclemys terrapin terrapin PRJNA339452, SRR4048682, SRR4048684, SRR4048681, SRR4048683
red eared slider, Trachemys scripta elegans, PRJNA552319, SRR13043486

## so just three species?
## if mutations don't go so far back, could probably push it further

## not really?
yellow bellied slider, Trachemys scripta scripta



## reference
CAS_Tse_1.0
## can manually download simpleRepeat and repeatMasker from the usual place



## not sure about this one? Actinemys marmorata


1








##
## scratch
##

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








table(assemblies[, "class"])






