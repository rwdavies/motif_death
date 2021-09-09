source("~/proj/motif_death/discover/functions.R")

assemblies <- get_ucsc_assemblies()

out <- investigate("Lagomorpha")
##
results <- out$results





##
## rabbits and friends
##
out <- investigate("Leporidae")
results <- out$results

## brush rabbit, Sylvilagus bachmani, PRJNA512907,  SRR12437581 <- this is the outgroup, probably
## european rabbit, Oryctolagus cuniculus, PRJEB28783, ERR2811813
## Granada hare, Lepus granatensis, PRJNA399194, SRR5949623
## Broom hare, Lepus castroviejoi, PRJNA562432, SRR10023741
## snowshoe hare, Lepus americanus, PRJNA420081, SRR6485265
## mountain hare, Lepus timidus, PRJNA399194, SRR5949626
## european hare, Lepus europaeus, PRJNA561428, SRR10011655, SRR10012548
## Black-tailed jackrabbit, Lepus californicus, PRJNA420081, SRR11020276





##
## pika
##

## see here https://doi.org/10.1093/molbev/msaa026
## Out of Tibet: Genomic Perspectives on the Evolutionary History of Extant Pikas
## Fig 2 for instance
out <- investigate("Lagomorpha")
results <- out$results


## so all the Ochotona are from two studies:
## PRJNA716776, which appears unpublished
## and PRJNA74593, which is the reference
## actually also PRJEB28783



## Illumina - each trio serves as the outgroup to the other?
## pika - Ochotona, (yes just Ochotona), SRR14101820, SRR14101857 ## conothoa
## american pika - Ochotona princeps, SRR413468, SRR413470 ## pika
## chinese red pika - Ochotona erythrotis, SRR14101824, SRR14101825 ## conothoa
## large eared pika - Ochotona macrotis, SRR14101814 ## conothoa
## mongolian pika - Ochotona pallasi, SRR14101830, SRR14101832 ## pika
## steppe pika - Ochotona pusilla, SRR14101831 ## pika


##

## BGIseq, but single end
## Ochotona huangensis
## Ochotona mantchurica
## Ochotona cansus, SRR14101836
## Ochotona roylei
## Ochotona thibetana

## not obviously enough?
## Ochotona coreana
## Ochotona curzoniae
## Ochotona dauurica
## Ochotona gaoligongensis
## Ochotona gloveri
## Ochotona himalayana
## Ochotona huanglongensis
## Ochotona muliensis
## Ochotona nubrica
## Ochotona sikimaria
## Ochotona thomasi

##
look_at_one_species_or_study(results, scientific_name = "Ochotona princeps")[1:3, ]

look_at_one_species_or_study(results, scientific_name = "Ochotona dauurica")

##
look_at_one_species_or_study(results, scientific_name = "Ochotona macrotis")


## seems like a lot?
m <- results[grep("Ochotona", results[, "scientific_name"]), ]
m <- m[m[, "instrument_platform"] != "BGISEQ", c("scientific_name", "base_count", "instrument_platform", "fastq_bytes")]
m <- m[m[, "scientific_name"] != "Ochotona princeps", ]
m <- cbind(m, n = sapply(strsplit(m[, "fastq_bytes"], ";"), length))
m[, 2] <- m[, 2] / 6e10
m$fastq_bytes <- NULL

m[order(-m[, 4], -m[, 2]), ]

1



## rabbits
Oryctolagus cuniculus
Oryctolagus cuniculus algirus
Oryctolagus cuniculus cuniculus

