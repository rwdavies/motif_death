motif_death
===========

## About

This repository contains code to download, pre-process, map, and variant call related species for use in motif death research.

## Workthrough / explanation

### Set-up
1. **Pre-requisites:** Most required packages are installed for shared use on the cluster (see `Snakefile_programs` for details). You **also need R and other packages**. If running on cluster, you can use existing installations, for example running this before anything else:
    ```
    module load R/3.6.0-foss-2018b
    module load SAMtools/1.9-foss-2018b
    module load HTSlib/1.9-foss-2018b
    export PATH=/apps/well/git/2.3.4/bin/:${PATH}
    ```
    You can put these lines in your `~/.bash_profile` and they will be run as part of the shell script to submit jobs (`run.sh`).

1. You need your own installation of [**HATBAG**](https://github.com/rwdavies/HATBAG). On the cluster, you can [follow instructions here](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/r-and-rstudio-on-the-bmrc-cluster) on installing R packages. This requires creating a `~/.Rprofile` file, downloading the HATBAG repo, then running the HATBAG installation as described in its' repo. 
    * Note: You can install once on rescomp 1 or 2, and ignore instructions on rescomp 3.
    * Note: you should be able to ignore pdflatex errors that arise during installation
    * Note: Make sure you add `options(bitmapType='cairo')` to `~/.Rprofile` as specified in cluster instructions ([more details](https://stackoverflow.com/questions/24999983/r-unable-to-start-device-png-capabilities-has-true-for-png))
1. You also need **Snakemake**.  On the BRMC (cluster), a fixed Snakemake installation is specified in the code and you don't need to install anything. Otherwise on local you can use `install_snakemake_2.sh` to semi-interactively install with Conda (see Snakemake docs for more info).
1. Once done, specify **user specific paths in `activate`. Ensure your HATBAG installation location is correctly captured.**

### A. Manual curation of input file and downloading reference material

Briefly, currently, identify a set of related species, appropriate for the model, and papers with their raw data, with data deposited in the EBI ENA. 

For example, for artiodactyla (includes ruminants and related species), this might include cows and goats, and more distantly related animals like okapi and giraffe.

For a given deposition project ID like PRJNA313910, we can then go to e.g. https://www.ebi.ac.uk/ena/browser/view/PRJNA313910 and manually download a text file with information like the download links to the fastq (raw data) files.

In `prepare_inputs.R` and `prepare_inputs_functions.R` we extract information about the raw data in a semi-automated fashion.

```
R -f R/prepare_inputs.R --args artiodactyla
```

Note: you have to run this from your project directory for environment variables (from `.activate`) to be read.

### B.1 Downloading per-sample raw data

Analysis is done per-sample for all samples in an order using the `run.sh` file, which requires manual intervention to note how species are linked to an order. For example `goat` is in `artiodactyla`, specified in `run.sh`.

To use `run.sh` further requires a file in `reference_info`, here `Snakefile_reference_artiodactyla`, with manual specification of things like the download path for the reference genome, the species, etc. Most of these can probably be set more automatically.

Note that run.sh pulls sources `activate`, which contains some manual paths, specific to the current machine and user. Similarly, `Snakefile_programs` contains a mixture of manual paths, which should be moved to the activate script (a script containing machine and user specific paths), and some programs set depending on paths that depend on the activate script, which is OK. 

Eventually, `run.sh` is actually run, on a head or data node. 
```
./run.sh preprocess all local goat --dryrun ## to check what will be run
./run.sh preprocess all local goat ## to run it
```
Note this might take from 30 minutes to a few hours of wall clock time to run per job, depending on size of individual fastqs.

### B.2 Download and prepare reference

On top of downloading files for each species above, you also need to download the reference genome and index it.
Run the below for any **one** of the species in your group (as reference is shared across all species in group).

```
./run.sh prep_reference all local goat --dryrun ## to check what will be run
./run.sh prep_reference all local goat ## to run it
```

### C. Mapping per-sample

Run for each species, jobs submitted to the cluster, this would look like the following

```
./run.sh mapping all cluster goat --dryrun
./run.sh mapping all cluster goat
```

This might take a few hours, or a few days to run through, depending on the sample being processed, and cluster availability. It first maps reads, then marks PCR duplicates, then runs indel realignment, then merges the bams together.

TODO: Once mapping complete and BAM OK, add in deletion of original fastq files (now done manually)

## D. Downstream variant calling per-order

Once all the samples from an order have been downloaded and mapped, it is time to call variants, and calculate depth of coverage, from which callable regions are derived, where "callable region" means regions where we can reasonably expect to call variants. Here we can perform this step using any of the names of the species in that order.

```
./run.sh downstream all cluster goat --dryrun
./run.sh downstream all cluster goat
```

Note that the number of cores being used per chromosome is a parameter set somewhere, one option is to make this parameter more visible, another is to break the region into chunks and call the genome in chunks them combine back together. This would run faster and more easily on the cluster but requires some coding.

TODO: `get_per_sample_average_cor.R` requires manual intervention to set list of chromosomes for new species, fix this so it runs automatically

## E. Treemix

The above should also run treemix (can be manually run using `./run.sh downstream treemix cluster goat`), which creates some plots in the analysis directory `treemix/` which shows the relationship between samples. This is pre-supposed in most cases, i.e. based on prior evidence, we assume a relationship between samples, though ideally we would get this from this plot. This is currently maually interpreted / assumed, though it would be good it this could be automated (at least, the within tree relationships - probably OK to assume the outgroup). The results from this (the interpreted results) are currently used in `R/run_all_functions.R` which should probably be automated.

## F. Analysis using HATBAG

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

