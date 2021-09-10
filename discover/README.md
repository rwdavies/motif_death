About
=====

Files in this folder are related to discovery new sets of species that could be suitable for motif_death/HATBAG

How to run
==========

Use code like the following to make a plot, and to investigate individual species
```
source("~/proj/motif_death/discover/functions.R")
assemblies <- get_ucsc_assemblies()
out <- investigate("Afrotheria")
results <- out$results
look_at_one_species_or_study(results, scientific_name="Chrysemys picta")
```

Note, to download tables, to get more options in the UCSC table browser, start in the assembly you want from the Genomes -> Other setting, then it will take you over to that species.

Samples to run / have been run
==============================

Todo, revisit birds, fish, possibly go further!

|Status|run name|Class, order, family|Common name|Notes|
|---|---|---|---|---|
|ready to go| anolis | Reptilia, Squamata | anolis | looks good, should be enough I hope|
|ready to go (need manual download)| emydidae | Reptilia, Testudines | turtles | nothing else looked good. if these turtles aren't doo distantly related, there are more than can be added (Emys, Terrapene, Actinemys)|
|ready to go|leporidae | Mammals, lagomorpha | rabits/hares|looks great|
|ready to go (need manual download)|ochotona | Mammals, lagomorpha | pikas|only using Illumina data. samples are from two phyla, use as outgroups to each other|
|unsure|eupleridae|Mammalia,Carnivora, Feliformia| meerkat |config made by Jhamat?|
|unsure|canis|Mammalia, Carnivora, Canidae| dogs | name should somehow be Cerdocyonina and Canina, but cannot get unique any better vs foxes|
|unsure|vulpes|Mammalia, Carnivora, Canidae | foxes | name should somehow be Vulpini and Urocyon, but cannot get unique any better vs dogs|




Amphibians
==========
Can't easily get Xenopus or Rana to work (not enough / too diverged), same with toad (bufo)
Couldn't obviously get caecilians to work either!



Reptiles - Lizards
==================

Should work fine with standard Anolis reference genome
bridled anole, Anolis frenatus, PRJNA400786, SRR8100042, SRR8100035 this is the outgroup
green anole, Anolis carolinensis, PRJNA381064, SRR5813770
grass anole, Anolis auratus, PRJNA400787, SRR8148313, SRR8148312
no other name, Anolis apletophallus, PRJNA400788, SRR8111633, SRR8111630, SRR8111632

Tried to get SceUnd_v1.1 to work, but not sure if occidentalis is enough
Sceloporus grammicus, NO
Sceloporus occidentalis, Maybe
Sceloporus tristichus, yes
Sceloporus undulatus, yes



Reptiles - Turtles
==================

Used

Painted turtle, Chrysemys picta (outgroup) PRJNA589899, SRR11059563

Diamondback terrapin, Malaclemys terrapin terrapin PRJNA339452, SRR4048682, SRR4048684, SRR4048681, SRR4048683

red eared slider, Trachemys scripta elegans, PRJNA552319, SRR13043486

Not used (could use if not too distantly related?)
Terrapene carolina triunguis, Emys orbicularis, Actinemys marmorata

Couldn't obviously get other turtles to work (I think?)


Pangolins
=========

So these below work, but the reference genome looks too discontinuous for this approach to work at this time


chinese pangolin, Manis pentadactyla, PRJNA529540, SRR9018595 (outgroup)
Indian pangolin, Manis crassicaudata, PRJNA490788, SRR7874732
sunda pangolin, Manis javanica, PRJNA529540, SRR9018632

There is a good sample Phataginus tricuspis but this looks too far diverged
https://www.pangolinsg.org/wp-content/uploads/sites/4/2018/04/Screenshot-7.png




Future / notes
==============


https://commons.wikimedia.org/wiki/File:An_evolutionary_tree_of_mammals.jpeg

Afrotheria, looks like elephants (3) and sea cows (3) are a possibility
Reference genomes for either don't look so great, the scaffolds aren't assembled into chromosomes, but could be fine? They are UCSC assemblies


Turtles, so far looks good for Emydidae
If they aren't too diverged / not too many SNPs, there are some other nearby turtles that could also potentially be used
