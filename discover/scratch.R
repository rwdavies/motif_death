motif_death_dir <- "~/proj/motif_death"
source(file.path(motif_death_dir, "discover/functions.R"))
assemblies <- get_ucsc_assemblies(motif_death_dir)


investigate_results <- investigate("Chondrichthyes", nCores = 48, plotdir = "~/proj/motif_death/discover")

results <- out$results


## get first set 
assemblies_to_consider <- get_refs_to_consider(result_phylo)







## https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/515/GCF_017639515.1_sCarCar2.pri/GCF_017639515.1_sCarCar2.pri_genomic.fna.gz

    
1

test <- fread("https://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt", data.table = FALSE)

## pseudo-algorithm




out <- investigate(keyword = "Canidae", nCores = 48)
results <- out$results

## optionally nuke afterwards
unlink(tempdir(), recursive = TRUE)
dir.create(tempdir(), showWarnings = FALSE)

## if this keeps crashing, just redo african wild dog SRR7874817 with SRR7874816
look_at_one_species_or_study(results, scientific_name = "ictus")


## scallops - possible? genome a bit small
xPecMax1.1
bay scallop, Argopecten irradians, PRJNA304363, SRR2976401, SRR2976403
great scallop, Pecten maximus, PRJEB22206, ERR3316025, ERR3316027, ERR3316026, ERR3316028
YessoEso, Mizuhopecten yessoensis, PRJNA259405, SRR1712896, SRR1712894
Azumapecten farreri aka Chlamys farreri, PRJNA185465, SRR2736013, I think OK
no plactopen, it is RADseq






## could try, using same one, against same reference?
head(look_at_one_species_or_study(results, scientific_name = "Rhinatrema bivittatum"), 1) ## two-lined caecilian, likely outgroup
head(look_at_one_species_or_study(results, scientific_name = "Geotrypetes seraphini"), 1)## Gaboon caecilian
head(look_at_one_species_or_study(results, scientific_name = "Microcaecilia unicolor"), 1) ## Microcaecilia unicolor

aGeoSer1.1
3,779,430,019
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/459/505/GCF_902459505.1_aGeoSer1.1/GCF_902459505.1_aGeoSer1.1_genomic.fna.gz





## should be OK?

## OK worth trying, genome looks fine
## Petromyzon marinus
## kPetMar1.pri
## https://www.ncbi.nlm.nih.gov/assembly/GCA_010993605.1

look_at_one_species_or_study(results, scientific_name = "Petromyzon")[1:5, ]

Sea lamprey, Petromyzon marinus, PRJNA385973, SRR5535434 (which is blood not sperm) (ref genome)
Brook lamprey,Lampetra planeri, PRJNA420358, SRR6329407, looks good
Far Eastern Brook Lamprey, Lethenteron reissneri, PRJNA558325, SRR9964061
Lethenteron camtschaticum - NO - not enough + is testis looks 


## not entirely sure it will work, divergence is high, potentially 40m years ago

look_at_one_species_or_study(results, scientific_name = "Lethenteron camtschaticum")
look_at_one_species_or_study(results, scientific_name = "Lethenteron reissneri")








## try the skates
## Amblyraja looks pretty good
## https://www.ncbi.nlm.nih.gov/assembly/GCF_010909765.2#/def
## check out the PDF
## assembly sAmbRad1.1.pri

## looks OK - sScyCan1.1
## 


Scyliorhinus

look_at_one_species_or_study(results, scientific_name = "Bathyraja")
")

Carcharodon carcharias - great white shark
Rhincodon typus - whale shark
diverged several hundred million years ago


## anolis looks fine, BUT, do we have enough?
m <- look_at_one_species_or_study(results, scientific_name="Anolis")
m <- m[m[, "base_count"] > 3e9, ]

## so these ones are well captured. can possibly do some more?



## not sure this will work...
## so these 4, potential ref fence lizard
## SceUnd_v1.1, looks OKx
Sceloporus grammicus ## NO
Sceloporus occidentalis ## Maybe
Sceloporus tristichus ## yes
Sceloporus undulatus ## yes
##


## here for frogs
https://www.sciencedirect.com/science/article/abs/pii/S105579031300300X


## Really not enough data
## reference genome aRanTem1.1, VGP primary
Rana dabieshanensis
Rana dalmatina
Rana muscosa
Rana temporaria

look_at_one_species_or_study(results, scientific_name="Rana temporaria")

##

## xenopus, doesn't look great. too diverged
## bufo bufo common toad


## frogs,


## possible amphibians?
## check out reference genomes?
1379 GCF_901001135.1 aRhiBiv1.1  Rhinatrema_bivittatum    two-lined caecilian
1382 GCF_901765095.1 aMicUni1.1 Microcaecilia_unicolor tiny Cayenne caecilian
1393 GCF_902459505.1 aGeoSer1.1  Geotrypetes_seraphini        Gabon caecilian
1416 GCF_905171765.1 aBufBuf1.1              Bufo_bufo            common toad
1417 GCF_905171775.1 aRanTem1.1        Rana_temporaria   European common frog



## manis looks promising?

## I think so? use pentadactyla as the outgroup, build the tree with the other two

chrom <- read.table("https://hgdownload.soe.ucsc.edu/goldenPath/manPen1//bigZips/manPen1.chrom.sizes")
## UGH, I dunno, too too many chromosomes to run, doesn't seem reasonable for now

out <- investigate("Testudines") ## turtles


## pangolins
chinese pangolin, Manis pentadactyla, PRJNA529540, SRR9018595 (outgroup)
Indian pangolin, Manis crassicaudata, PRJNA490788, SRR7874732
sunda pangolin, Manis javanica, PRJNA529540, SRR9018632


look_at_one_species_or_study(results, scientific_name="Manis pentadactyla")[1:4, ]






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






