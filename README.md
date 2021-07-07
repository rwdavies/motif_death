motif_death
===========

## About

This repository contains code to download, pre-process, map, and variant call related species for use in motif death research.

## Warning!

**Repo will only run through for test data in current state, as additional specifications for HATBAG in `R/run_all_functions.R` have not been done for other orders / species.**

## Testing

```
. activate
./scripts/test.sh
```

## Workthrough / explanation

### Set-up
All packages and user-specific paths are specified in `activate`. **A new user should define their own packages and paths there.**  Some notes on installations below:

1. For the version of R you are using (i.e. specified in `activate`), you need to install of [**HATBAG**](https://github.com/rwdavies/HATBAG). On the cluster, you can [follow instructions here](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/r-and-rstudio-on-the-bmrc-cluster) on installing R packages. This requires creating a `~/.Rprofile` file, downloading the HATBAG repo, then running the HATBAG installation as described in its' repo. 
    * Note: You can install once on rescomp 1 or 2, and ignore instructions on rescomp 3.
    * Note: you should be able to ignore pdflatex errors that arise during installation
    * Note: Make sure you add `options(bitmapType='cairo')` to `~/.Rprofile` as specified in cluster instructions ([more details](https://stackoverflow.com/questions/24999983/r-unable-to-start-device-png-capabilities-has-true-for-png))
1. You can install your own version of **Snakemake** and specify it in `activate`. See [Snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more instructions, which requires first installing conda or mamba, and then installing snakemake within an environment.

**Before running any other steps, ensure you run the following first:**
```
. activate
```

### A. Manual curation of input file and downloading reference material

Briefly, currently, identify a set of related species, appropriate for the model, and papers with their raw data, with data deposited in the EBI ENA. 

For example, for artiodactyla (includes ruminants and related species), this might include cows and goats, and more distantly related animals like okapi and giraffe.

For a given deposition project ID like PRJNA313910, we can then go to e.g. https://www.ebi.ac.uk/ena/browser/view/PRJNA313910 and manually download a text file with information like the download links to the fastq (raw data) files.

In `prepare_inputs.R` and `prepare_inputs_functions.R` we extract information about the raw data in a semi-automated fashion, and output `species_mapping_info/{species}.json` for each species.

```
R -f R/prepare_inputs.R --args artiodactyla
```

### B Motif Death & HATBAG run

For a given run, titled 'artiodactyla', consisting of several species e.g. goat, cow.
Once you have `species_mapping_info/{species}.json` for all species you are interested in (from step A), create `config/artiodactyla.json` with all arguments for motif death and HATBAG (you can reference .jsons for other orders in the same folder).

### B.1. Downloading files

You need internet access to download fastqs and reference files, so this needs to be run on a head node on the cluster ('local'). 
```
./run.sh config/artiodactyla.json download_all local --dryrun ## to check what will be run
./run.sh config/artiodactyla.json download_all local ## to run it
```
Note this might take from 30 minutes to a few hours of wall clock time to run per job, depending on size of individual fastqs.

### B.2 Running the pipeline

```
./run.sh config/artiodactyla.json download_all cluster --dryrun ## to check what will be run
./run.sh config/artiodactyla.json download_all cluster ## to run it
```

**Notes:**
Mapping: This might take a few hours, or a few days to run through, depending on the sample being processed, and cluster availability. It first maps reads, then marks PCR duplicates, then runs indel realignment, then merges the bams together.

Downstream: Once all the samples from an order have been downloaded and mapped, it is time to call variants, and calculate depth of coverage, from which callable regions are derived, where "callable region" means regions where we can reasonably expect to call variants.

Downstream: The above should also run treemix (can be manually run using `./run.sh config/artiodactyla.json treemix cluster`), which creates some plots in the analysis directory `treemix/` which shows the relationship between samples. This is pre-supposed in most cases, i.e. based on prior evidence, we assume a relationship between samples, though ideally we would get this from this plot. This is currently maually interpreted / assumed, though it would be good it this could be automated (at least, the within tree relationships - probably OK to assume the outgroup). The results from this (the interpreted results) are currently used in `R/run_all_functions.R` which should probably be automated.

HATBAG: HATBAG itself is a relatively OK piece of code. It can be run from the command line.  The current integration of this into the script is very hacky. It should work - it does work - but it is inelegant, and was put together now in a brute force fashion (previously this ran on a single compute server with 300+ GB of RAM, rather than on a cluster). Old code towards this goal (in `HATBAG.smk`) still exists but isn't quite working.  

## Appendix
Fixed directory and filenames specified in `config/filenames.json`.

## TODOs

* Mapping: Note that the number of cores being used per chromosome is a parameter set somewhere, one option is to make this parameter more visible, another is to break the region into chunks and call the genome in chunks them combine back together. This would run faster and more easily on the cluster but requires some coding.
* Once mapping complete and BAM OK, add in deletion of original fastq files (now done manually)
* `get_per_sample_average_cor.R` requires manual intervention to set list of chromosomes for new species, fix this so it runs automatically
* HATBAG requires manual downloading of files (explanation of how this was done, and how this might be automated, to be discussed in time) from UCSC as specified in `R/run_all_functions.R` for now (simpleRepeat file, rmask file). This also requires manual specification of paths and values in `R/run_all_functions.R`
* Clear out old files in `/well/davies/users/dcc832/primates/hatbag_OLD_TO_DELETE`
* Make HATBAG use far far fewer temporary files in steps B -> C (and later?)
* Check out Step C, gainat, nonRepeat, artiodactyla, took ~6 hours?
* Also can Step B run better when using a cluster? Use far more jobs, 1 core each, split out types of test, run all in parallel?
* Should cluster / snakemake use be incorporated into HATBAG itself? Do I use snakemake in tests for some repo? 

