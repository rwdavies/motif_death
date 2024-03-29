##
## snakes
##

## sub-families manually copied in
subfamilies <- c("Dipsadidae", "Sibynophiinae", "Natricinae", "Colubrinae", "Colubridae", "Ahaetuliinae", "Calamariinae", "Grayiinae", "Dipsadinae", "Pseudoxenodontinae")
results <- get_match_against_subfamilies(subfamilies)
results <- unique(results)

f <- function(study, sample = NA) {
    cols <- c("scientific_name", "instrument_platform", "instrument_model", "library_source", "library_selection", "fastq_bytes", "nominal_length", "first_created", "read_count", "study_accession", "run_accession", "library_strategy", "run_accession", "base_count")
    if (!is.na(sample)) {
        m <- results[results[, "study_accession"] == study & results[, "scientific_name"] == sample, ]
    } else {
        m <- results[results[, "study_accession"] == study & results[, "scientific_name"] == sample, ]
    }
    m <- m[order(-m[, "read_count"]), ]
    m <- cbind(m, run_total = cumsum(m[, "base_count"]) / 6e10)
    m[, c(cols, "run_total")]
}


## all in Family = Colubridae

## REFERENCE ONLY - Thamnophis elegans (what weve been using!)
f("PRJNA189551") ## Thamnophis sirtalis (what weve been using!) ["SRR786678", "SRR770194"] (or what was used previously)
f("PRJNA588151") ## Ptyas mucosa (asssembly - looks fine) ["SRR10412127", "SRR10412131" "SRR10412132" "SRR10412124" "SRR10412126" "SRR10412129" "SRR10412130" "SRR10412123"]
f("PRJNA268069") ## Pantherophis guttatus (asssembly - looks fine) ["SRR9596764", "SRR9596760", "SRR9596763", "SRR9596761"]
f("PRJNA587592") ## Pantherophis obsoletus, Western Rat Snake (assembly - looks fine) ["SRR10405238", "SRR10405237", "SRR10405232", "SRR10405231", "SRR10405221", "SRR10405230", "SRR10405229", "SRR10405222", "SRR10405225", "SRR10405234", "SRR10405233", "SRR10405223", "SRR10405236", "SRR10405235", "SRR10405226", "SRR10405224"]
## gold mine below
f("PRJNA473624", "Thermophis baileyi") ## Baileys snake ["SRR7286062"]
f("PRJNA473624", "Thermophis zhaoermii") ## Sichuan hot-spring keelback ["SRR7293238", "SRR7293239"]
f("PRJNA473624", "Pseudoxenodon macrops") ## Chinese False Cobra ["SRR7293240", "SRR7293241"]
f("PRJNA473624", "Thermophis shangrila") ## Shangrila hot-spring keelback ["SRR7293244", "SRR7293245"]
f("PRJNA473624", "Pseudoxenodon bambusicola") ## Bamboo False Cobra ["SRR7293247", "SRR7293246"]

## genome is ~2Gbp

## NO TO THESE
PRJNA209850, Lampropeltis getula, common kingsnake (not enough reads)
Masticophis flagellum (not enough reads)


## PREVIOUSLY USED
ijima_sea_snake Emydocephalus_ijimae Family = Elapidae
blue_ringed_sea_krait Family = 	Elapidae
slender_necked_sea_snake = Elapidae
yellow_lipped_sea_krait = Elapidae
common_viper = Viperidae
garter_snake = 	Colubridae (family Natricinae)


pdf("~/Downloads/snaketree.pdf", height = 300, width = 20)
tree <- ape::read.tree("~/Downloads/snakes.nh")
tip.color <- rep("black", length(tree$tip.label))
tip.color[grep("Thermophis", tree$tip.label)] <- "blue"
tip.color[grep("Thamnophis", tree$tip.label)] <- "red"
plot(tree, tip.color = tip.color)
dev.off()



m <- results[results[, "scientific_name"] %in% "Masticophis flagellum" , cols]

## PROBABLY NOT
PRJNA209850, Lampropeltis getula, Family = Colubridae (I think this is OK? see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3842982/)


head(results[, )])

## 1M reads is ~100Mbp, or maybe 200 Mbp
## clearly need more like 100M reads
## so only two studies with > 10M reads
results[results[, "read_count"] > 1e7, ]

## this does feel very very inefficient though hmm

## Lots of promise here..
results

## just seems like the one garter snake
## can I go up a level?








##
##
##























    results <- data.frame(rbindlist(lapply(w, function(iRow) {
        cbind(
            iGenus = iGenus,
            genus = genus,
            iRow = iRow,
            taxId = available[iRow, 1],
            scientificName = available[iRow, 2],
            displayName = available[iRow, 3],
            results_list[[iRow]][, c("scientific_name", "instrument_platform", "instrument_model", "library_source", "library_selection", "fastq_bytes", "nominal_sdev", "first_created", "read_count", "study_accession")]
        )
    })))



pdf("~/Downloads/snaketree.pdf", height = 300, width = 8)
plot(ape::read.tree("~/Downloads/snakes.nh"))
dev.off()



x <- data.table::fread("~/Downloads/3811890/S1 File.PHY", data.table = FALSE)
x <- x[-1, ]
x[, 2] <- as.character(x[, 2])
x2 <- t(sapply(x[, 2], function(a) unlist(strsplit(a, ""))))


## Look at some favourites?
f <- function(a, b) {
m <- table(
    x2[grep(a, x[, 1]), ],
    x2[grep(b, x[, 1]), ]
)
m <- m[-1, ]
m <- m[, -1]
good <- sum(diag(m))
bad <- sum(m) - sum(diag(m))
1 - bad / (good + bad)
}


snakes <- c("Thamnophis_sirtalis","Thamnophis_elegans","Pantherophis_guttatus","Ptyas_mucosa","Pantherophis_obsoletus","Thermophis_baileyi","Thermophis_zhaoermii","Pseudoxenodon_macrops", "Laticauda_laticaudata", "Laticauda_colubrina")#,"Thermophis_shangrila") #,"Pseudoxenodon_bambusicola")
d <- array("", c(length(snakes), length(snakes)))
colnames(d) <- snakes
rownames(d) <- snakes
for(i in 1:length(snakes)) {
    for(j in 1:length(snakes)) {
        a <- snakes[i]
        b <- snakes[j]
        d[i, j] <- round(f(a, b), 3)
    }
}













## data <- read.table("~/Downloads/results_sequence_tsv.txt", header = TRUE, sep = "\t")

## ## is this - everything?

## require(XML)
## data <- xmlParse("~/Downloads/ena_read_run_20210804-2032.xml")
## xml_data <- xmlToList(data)
## length(xml_data)


## ## convert to a data frame?
## ## then, um, save?



## devtools::install_github("https://github.com/cstubben/ENAbrowseR")
## library("ENAbrowseR")
## yp <- ena_search("tax_tree(632)", result= "assembly", drop = FALSE)




