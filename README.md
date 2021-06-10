motif_death
===========

## About

This repository contains code to download, pre-process, map, and variant call related species for use in motif death research.

## Workthrough / explanation

### Set up notes
User specific paths are defined in `activate`. Please change as needed before running. You should then load them into the environment by running:

```
source ./activate
```

You need to install R. If running on BMRC (cluster), you can use an existing R installation, for example:

```
module load R/3.6.2-foss-2019b
```

Note, dependent R packages are not currently captured, and will cause program failure at some point for example R package `data.table`.

### Manual curation of input file and downloading reference material

Briefly, currently, identify a set of related species, appropriate for the model, and papers with their raw data, with data deposited in the EBI ENA. 

For example, for artiodactyla (includes ruminants and related species), this might include cows and goats, and more distantly related animals like okapi and giraffe.

For a given deposition project ID like PRJNA313910, we can then go to e.g. https://www.ebi.ac.uk/ena/browser/view/PRJNA313910 and manually download a text file with information like the download links to the fastq (raw data) files.

In `prepare_inputs.R` and `prepare_inputs_functions.R` we extract information about the raw data in a semi-automated fashion.

```
R -f R/prepare_inputs.R --args artiodactyla
```

### Downloading per-sample raw data

Analysis is done per-sample for all samples in an order using the `run.sh` file, which requires manual intervention to note how species are linked to an order. For example `goat` is in `artiodactyla`, specified in `run.sh`.

To use `run.sh` further requires a file in `reference_info`, here `Snakefile_reference_artiodactyla`, with manual specification of things like the download path for the reference genome, the species, etc. Most of these can probably be set more automatically.

Note that run.sh pulls sources `activate`, which contains some manual paths, specific to the current machine and user. Similarly, `Snakefile_programs` contains a mixture of manual paths, which should be moved to the activate script (a script containing machine and user specific paths), and some programs set depending on paths that depend on the activate script, which is OK. 

Eventually, `run.sh` is actually run, on a head or data node. It requires snakemake to be installed, currently done using `install_snakemake_2.sh` semi-interactively through conda
```
run.sh preprocess all local goat --dryrun ## to check what will be run
run.sh preprocess all local goat ## to run it
```
Note this might take from 30 minutes to a few hours of wall clock time to run per job, depending on size of individual fastqs.

### Mapping per-sample

Run for each species, jobs submitted to the cluster, this would look like the following

```
run.sh mapping all cluster goat --dryrun
run.sh mapping all cluster goat
```

This might take a few hours, or a few days to run through, depending on the sample being processed, and cluster availability. It first maps reads, then marks PCR duplicates, then runs indel realignment, then merges the bams together.

TODO: Once mapping complete and BAM OK, add in deletion of original fastq files (now done manually)

## Downstream variant calling per-order

Once all the samples from an order have been downloaded and mapped, it is time to call variants, and calculate depth of coverage, from which callable regions are derived, where "callable region" means regions where we can reasonably expect to call variants. Here we can perform this step using any of the names of the species in that order.

```
./run.sh downstream all cluster goat --dryrun
./run.sh downstream all cluster goat
```

Note that the number of cores being used per chromosome is a parameter set somewhere, one option is to make this parameter more visible, another is to break the region into chunks and call the genome in chunks them combine back together. This would run faster and more easily on the cluster but requires some coding.

TODO: `get_per_sample_average_cor.R` requires manual intervention to set list of chromosomes for new species, fix this so it runs automatically

## Treemix

The above should also run treemix (can be manually run using `./run.sh downstream treemix cluster goat`), which creates some plots in the analysis directory `treemix/` which shows the relationship between samples. This is pre-supposed in most cases, i.e. based on prior evidence, we assume a relationship between samples, though ideally we would get this from this plot. This is currently maually interpreted / assumed, though it would be good it this could be automated (at least, the within tree relationships - probably OK to assume the outgroup). The results from this (the interpreted results) are currently used in `R/run_all_functions.R` which should probably be automated.

## Analysis using HATBAG

HATBAG itself is a relatively OK piece of code. It can be run from the command line.

The current integration of this into the script is very hacky. It should work - it does work - but it is inelegant, and was put together now in a brute force fashion (previously this ran on a single compute server with 300+ GB of RAM, rather than on a cluster). Old code towards this goal (in `Snakefile_HATBAG`) still exists but isn't quite working.

This step requires manual downloading of files (explanation of how this was done, and how this might be automated, to be discussed in time) from UCSC as specified in `R/run_all_functions.R` for now (simpleRepeat file, rmask file). This also requires manual specification of paths and values in `R/run_all_functions.R`

This does currently run on the cluster using

```
./run.sh HATBAG HATBAG_HACK cluster goat
```

## Other TODO for Robbie

Clear out old files in `/well/davies/users/dcc832/primates/hatbag_OLD_TO_DELETE`
Make HATBAG use far far fewer temporary files in steps B -> C (and later?)
Check out Step C, gainat, nonRepeat, artiodactyla, took ~6 hours?
Also can Step B run better when using a cluster? Use far more jobs, 1 core each, split out types of test, run all in parallel?
Should cluster / snakemake use be incorporated into HATBAG itself? Do I use snakemake in tests for some repo? 

