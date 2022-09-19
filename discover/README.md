About
=====

Files in this folder are related to discovery of new sets of species that could be suitable for motif_death/HATBAG. functions.R is used for investigating a specific taxon (eg. Testudines, Hominidae) for a suitable reference genome and short read data. progress_visualise.R is used for marking groups of species that have been investigated for this pipeline on a giant tree plot (eg. root: Chordata, leaf rank: order). 

How to use functions.R
======================

Use code like the following to make a plot, and to investigate individual species
```
motif_death_dir <- "~/proj/motif_death"
source(file.path(motif_death_dir, "discover/functions.R"))
assemblies <- get_ucsc_assemblies(motif_death_dir)
out <- investigate("Afrotheria")
results <- out$results
look_at_one_species_or_study(results, scientific_name="Chrysemys picta")
```

How to use progress_visualise.R
================================

Data on currently investigated groups of species is stored in `motif_death/progress_track.json`. 

Use code like this to make a plot:

```
motif_death_dir <- "~/proj/motif_death"
ncbi_names <- read.csv(file.path(motif_death_dir, "discover/resources/names.csv.gz"))
source(file.path(motif_death_dir, "discover/progress_visualise.R"))

make_plot(root="amphibia", leaf_rank="family", plot_dir="~/Downloads", plot_scale = 0.2, max_iter = 100)
```

Getting simpleRepeat / repeatMasker
===================================

Note, to download tables, to get more options in the UCSC table browser, start in the assembly you want from the Genomes -> Other setting, then it will take you over to that species.

Also, if you have to download multiple, see file merge_different_repeatMasker_files.sh as a way to do it.

Samples to run / have been run
==============================

<!--- marker -->
| Common name | Config name | Status | Notes |
|-|-|-|-|
| cod icefish | notothenioidei_20220803.json | running | Richard running, have an error 0 significant k-mers |
| mice | mus_20220825.json | run success |  |
| salmon | salmonidae_20220823.json | running |  |
| apes | hominidae_20220830.json | running |  |
| Eagles | aquila_20220902.json | revisit | Four species: Heliaca, [chrysaetos canadensis, chrysaetos chrysaetos], spilogaster, heliaca outgroup. Divergence maybe too high at 6% |
| hummingbirds | trochilidae_20220906.json | ready |  |
| cats | felidae_20220913.json | running |  |
| anoles | anolis_20210909.json | not feasible | Ref genome is 1Gbp data in chromosomes. Currently (Sep 2022) 4 species with sufficient data, but divergence at 50M years is too high. Previous notes: looks good, should be enough I hope. Tried to get SceUnd_v1.1 to work, but not sure if occidentalis is enough Sceloporus grammicus, NO Sceloporus occidentalis, Maybe Sceloporus tristichus, yes Sceloporus undulatus, yes |
| cows and similar | artiodactyla_20210708.json | revisit |  |
| finches and similar | avian_20210804.json | revisit |  |
| turtles | emydidae_20210909.json | run success | worked |
| caecilians | gymnophiona_20210929.json | revisit | 3 available, all VGP refs, not sure if will work, not sure about divergence |
| rabbits/hares | leporidae_20210909.json | revisit | looks great |
| pikas | ochotona_20210909.json | revisit | only using Illumina data. samples are from two phyla, use as outgroups to each other |
| lampreys | petromyzontidae_20210929.json | run failed | could be worth revisiting. Keyword: Petromyzontiformes. Four possibilities, can make a tree from 3 of them (the other not enough data). Quite possibly the sea lamprey is too far diverged against the others, but we'll see. Sea lamprey, Petromyzon marinus, PRJNA385973, SRR5535434 (which is blood not sperm) (ref genome) Brook lamprey,Lampetra planeri, PRJNA420358, SRR6329407, looks good Far Eastern Brook Lamprey, Lethenteron reissneri, PRJNA558325, SRR9964061 Lethenteron camtschaticum - NO - not enough + is testis looks |
| Pangolins | pholidota_20210909.json | not feasible | No good ref (Sep 2022). Previous notes: So these below work, but the reference genome looks too discontinuous for this approach to work at this time: chinese pangolin, Manis pentadactyla, PRJNA529540, SRR9018595 (outgroup) Indian pangolin, Manis crassicaudata, PRJNA490788, SRR7874732 sunda pangolin, Manis javanica, PRJNA529540, SRR9018632. There is a good sample Phataginus tricuspis but this looks too far diverged https://www.pangolinsg.org/wp-content/uploads/sites/4/2018/04/Screenshot-7.png |
| snakes | snakes_20210810.json | revisit | Snakes are very diverged. Our only good ref is in Colubridae which is ~40Mya. So I think stay within Colubridae. Unsure of divergences within Colubridae. There are not too many species with data, so could be that only available options are too diverged. Revisit sometime but low priority for now. |
| whale (minke ref) | whippomorpha_20210713.json | run success |  |
| whale (blue ref) | whippomorpha_bluewhale_ref_20210802.json | run success |  |
| Scallops | pectinidae_20210929.json | run failed | 4X coverage for bay scallop. Too diverged (~6-8%), and there are no other feasible scallops available. |
| dogs | canis_20210910.json | revisit | name should somehow be Cerdocyonina and Canina, but cannot get unique any better vs foxes. coverage OK though andean fox and dhole low at about 10 to 15 X coverage, could be worth adding more for a re-run |
| foxes | vulpes_20210910.json | run success | not much loss signal, wide AT to GC, likely like dogs which do not have PRDM9 and have recombination localize at promoters |
| meerkats | eupleridae_20210908.json | revisit | config made by Jhamat, I think? Also should probably be called Viverroidea. Note if this doesn't look good at the start, remove the hyaenidae species, and re-run with just fossa as the outgroup |
| Sharks |  | not feasible | Can't find anything obvious that works. Many sets of 3-5 species identified. However divergence very high between them (e.g. great white shark vs whale shark), and within sets of 3, can't find compatible high quality reference genome, AND 3 or more sets of high quality Illumina non-GAII sequences. Maybe in a few more years. |
| Starfish |  | not feasible | Genome probably too small at less than 400 Mbp |
| Molluscs |  |  | Worth considering octopuses, as well as oysters, they look viable. For octoposes, in particular https://www.ncbi.nlm.nih.gov/assembly/GCF_006345805.1 ASM634580v1 Muusoctopus leioderma Octopus rubescens Octopus bimaculoides Octopus vulgaris. |
| Amphibians |  |  | Can't easily get Xenopus or Rana to work (not enough / too diverged), same with toad (bufo) Couldn't obviously get caecilians to work either! |
| Sea squirts |  | not feasible | Genome too small at 100MB |
| Kiwi |  | not feasible | Don't have a good reference. The phylogenetic tree looks good though, we can take 4 different kiwis. (Sep 2022) |
| Duck |  | not feasible | Have good tufted duck ref. Divergence at 20M years feels too high |
| Penguins |  | not feasible | Don't have good ref (Sep 2022) |
| Sturgeons |  | not feasible | Not enough data (Sep 2022). |
| Chicken |  | not feasible | Chickens are too undiverged (<1%). But maybe could work if look outside just chickens? |
| Tinamou |  | not feasible | Think too diverged + no good ref. (Sep 2022) |
| Moa |  | not feasible | Not enough data (Sep 2022) |
| emu |  | not feasible | Not enough data (Sep 2022) |
| Rhea |  | not feasible | Not enough data (Sep 2022) |
<!--- marker -->

Future / notes
==============


https://commons.wikimedia.org/wiki/File:An_evolutionary_tree_of_mammals.jpeg

Afrotheria, looks like elephants (3) and sea cows (3) are a possibility
Reference genomes for either don't look so great, the scaffolds aren't assembled into chromosomes, but could be fine? They are UCSC assemblies


Turtles, so far looks good for Emydidae
If they aren't too diverged / not too many SNPs, there are some other nearby turtles that could also potentially be used
















