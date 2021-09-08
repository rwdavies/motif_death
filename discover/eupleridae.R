source("~/proj/motif_death/discover/functions.R")

assemblies <- get_ucsc_assemblies()

##
## eupleridae (within feliformia but not feloidea)
## useful: https://www.biorxiv.org/content/10.1101/2020.10.05.326090v1.full.pdf
##
feliformia <- investigate("Feliformia")


## so here, basically everything that's not felidae, is worth trying
results <- feliformia$results

## OK, so use Meerkat reference
## GCF_006229205.1 aka meerkat_22Aug2017_6uvM2_HiC
##
## OK looks like there might be sadly many outgroups
## basically everything in hyaenidae will be outgroup
## then there is fossa
## then there are things in herpestidae like meercats and mongooses
##
##

## striped hyaena
look_at_one_species_or_study(results, scientific_name = "Hyaena hyaena")
## use: SRR5904112 (can also use SRR5904111)

## brown hyaena
look_at_one_species_or_study(results, scientific_name = "Parahyaena brunnea")
## use: SRR5886633

## aardwolf
look_at_one_species_or_study(results, scientific_name = "Proteles cristata")
## use: SRR13528972

## spotted hyaena
look_at_one_species_or_study(results, scientific_name = "Crocuta crocuta")[1:5, ]
## use: SRR9914662

## fossa
look_at_one_species_or_study(results, scientific_name = "Cryptoprocta ferox")
## use: SRR11097184 (can also use SRR11428701)

## banded mongoose
look_at_one_species_or_study(results, scientific_name = "Mungos mungo")
## just this one
## use: SRR7704821

## Common dwarf mongoose
look_at_one_species_or_study(results, scientific_name = "Helogale parvula")
## use: SRR7637809

## Meerkat
look_at_one_species_or_study(results, scientific_name = "Suricata suricatta")
##use: SRR11434616






