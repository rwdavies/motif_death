source("~/proj/motif_death/discover/functions.R")

assemblies <- get_ucsc_assemblies()

##
## eupleridae (within feliformia but not feloidea)
## useful: https://www.biorxiv.org/content/10.1101/2020.10.05.326090v1.full.pdf
##
canidae <- investigate("Canidae")
results <- canidae$results

## obviously use normal dog reference CanFam3.1
## then there are two groups that you can run, using the other as the outgroups
## group1 = the one with vulpes
## group2 = the one with canis

f <- function(x) {
    results[match(x, results[, "run_accession"]), "study_accession"]
}

## foxes and friends
## use Arctic fox reference ASM1834538v1 aka GCF_018345385.1 (I think?)
## else use regular dog reference? not that much saved for island fox?
## else can use otocyon as reference possibly? a bit better?
## outgroups should be two Urocyon's (I think)
## and or can also include Otocyon megalotis
Urocyon littoralis - Island fox, PRJNA478450, SRR7458264
Urocyon cinereoargenteus - Gray fox, PRJNA478450, SRR7458270
Otocyon megalotis - Bat eared fox, PRJNA681015, SRR13177425
(ignoring subspecies difference megalotis vs virgatus)
## Nyctereutes procyonoides - Common racoon dog ## NOT ELIGIBLE - not enough sequence
Vulpes lagopus - Arctic fox, PRJEB43377, ERR5417970
Vulpes zerda - Fennec fox, PRJNA648123, SRR14750511
Vulpes vulpes - Red fox, PRJNA681015, SRR13177416


## dogs and friends (Caninae)
## use normal dog reference CanFam3.1 (actually, newer one)
## use Chrysocyon brachyurus and Lycalopex culpaeus as outgroup to start
Chrysocyon brachyurus - Maned wolf, PRJNA512907, SRR13167987
Lycalopex culpaeus - Culpeo or Andean Fox, PRJNA232497, SRR1066702, SRR7107709, SRR1066703 (probably different animals? but I think OK?)
## unclear?
## going down
Canis mesomelas elongae - Black backed jackal, PRJEB31639, ERR3210523
Lycaon pictus - African wild dog, PRJNA488046, SRR7874817
Cuon alpinus lepturus - Dhole, PRJNA448733, SRR7107853
Canis lupaster - African wolf, PRJEB31639, ERR3245533, ERR3245532
Canis lupus familiaris - Dog, PRJEB39198, ERR4318109 ## chose one semi at random

## not used
## Canis adustus, Side-striped jackal NO, it is BGI, don't trust tech for this


