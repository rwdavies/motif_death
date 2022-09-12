motif_death
===========

## About

This repository contains code to download, pre-process, map, and variant call related species for use in motif death research.

# Table of contents
1. [Introduction](#paragraph-testing)
2. [Setup](#paragraph-setup)
3. [Running Motif Death and HATBAG](#paragraph-motifandhatbag)
    1. [Downloading files](#paragraph-download)
    2. [Running the pipeline](#paragraph-running)
4. [Appendix](#paragraph-appendix)
4. [TODOs](#paragraph-todos)

## Testing <a name="paragraph-testing"></a>

```
. activate
./scripts/test.sh --dryrun # check what will be run (will fail test)
./scripts/test.sh 4 # specify number of cores, default 8 if left blank
```


## Setup <a name="paragraph-setup"></a>
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

## Running Motif Death and HATBAG  <a name="paragraph-motifandhatbag"></a>

For a group of species (titled `SPECIES_ORDER` e.g., artiodactyla) you want to analyze, you should write a `config/{...}.json` (see other files in folder for reference):
* **Reference**: Start with finding a suitable reference for your group of species. You need the reference genome, as well as simpleRepeat and repeat mask tables.
    * Note: By assumption, one reference is used for all mapping within a given `SPECIES_ORDER`. This means once a species is mapped, it will not be remapped if used in a different run. To re-map species with a different reference, you should create a new config with a different `SPECIES_ORDER`. This will redownload fastqs in a different folder, so watch storage space.
* `SPECIES_LIST`: Specify your desired species, which shouldn't be too diverged (find a reference phylogenetic tree). For example, for artiodactyla (focusing on ruminants and related species), this might include cows and goats, and more distantly related animals like okapi and giraffe. Also include a more distant species as as outgroup, e.g. white tailed deer (as `TREEMIX_OUTGROUP`).
    * For each species, data must exist in the EBI ENA, with a project ID e.g. PRJNA313910.  Under this Project ID, choose the 'Run Accessions' (i.e. units) to use for each species.
    * TODO: add notes about what is good choice (high coverage, short read ...)
* `RUN_ID`: (e.g., 20200706) Coverage, VCF, Treemix, and HATBAG are all run specific, as you can have multiple runs within the same order using different species.
* `HATBAG_OUTPUT_DIR`: HATBAG has an additional ID parameter, so you can redo the same order and run with different HATBAG settings.
* `REF_NAME` and `REF_URL`: The reference URL (self-explanatory) and reference name (a short form description, usually the part before `.fa`).
* `CHR_LIST_ONLY_AUTOS` and `CHR_LIST`: a list of the chromosomes that will be analyzed. Despite the names, the current assumption in the code is that only non-sex chromosomes will be used.
* `TREEMIX_OUTGROUP` species to be used as the outgroup for Treemix and throughout.
* `HATBAG_PARAMS` specifically `genomeSize`, is the reference fasta genome length.


The resulting folder structure is below:

```
{ANALYSIS_DIR}
├── external # rmasks, simpleRepeats downloaded here
├── ref # references downloaded here
├── mapping
│   └── {SPECIES_ORDER}
│       └── {species} # fastq downloaded here and mapped
├── coverage
│   └── {SPECIES_ORDER}
│       └── {RUN_ID}
├── treemix
│   └── {SPECIES_ORDER}
│       └── {RUN_ID}
├── vcf
│   └── {SPECIES_ORDER}
│       └── {RUN_ID}
├── hatbag
│   └── {SPECIES_ORDER}
│       └── {RUN_ID}
│           └── {HATBAG_OUTPUT_DIR}
└── logs # snakemake logs
```

### Downloading files <a name="paragraph-download"></a>

You need internet access to download fastqs and reference files, so this needs to be run on a head node on the cluster ('local'). 
```
./run.sh config/artiodactyla.json download_all local --dryrun ## to check what will be run
./run.sh config/artiodactyla.json download_all local ## to run it
```
Note this might take from 30 minutes to a few hours of wall clock time to run per job, depending on size of individual fastqs.

### Running the pipeline <a name="paragraph-running"></a>

```
./run.sh config/artiodactyla.json all cluster --dryrun ## to check what will be run
./run.sh config/artiodactyla.json all cluster ## to run it
```

**Notes:**

Mapping: This might take a few hours, or a few days to run through, depending on the sample being processed, and cluster availability. It first maps reads, then marks PCR duplicates, then runs indel realignment, then merges the bams together.

Downstream: Once all the samples from an order have been downloaded and mapped, it is time to call variants, and calculate depth of coverage, from which callable regions are derived, where "callable region" means regions where we can reasonably expect to call variants.

Downstream: Treemix will run automatically (can also be manually run using `./run.sh config/artiodactyla.json treemix cluster`), which creates plots in `treemix/` showing the relationship between samples. If in `HATBAG_PARAMS` you do not specify `lineages`, it is automatically deduced from Treemix output.  However, outgroup is then assumed to be the same as `TREEMIX_OUTGROUP`.

HATBAG: HATBAG itself is a relatively OK piece of code. It can be run from the command line.  The current integration of this into the script is very hacky. It should work - it does work - but it is inelegant, and was put together now in a brute force fashion (previously this ran on a single compute server with 300+ GB of RAM, rather than on a cluster). Old code towards this goal (in `HATBAG.smk`) still exists but isn't quite working.  

## Appendix <a name="paragraph-appendix"></a>
Directory paths and filenames specified in `config/filenames.json`.

## Robbie rsync process <a name="paragraph-appendix"></a>
Now that BMRC has two factor authentication, it is cumbersome to rsync from individual directories. So we stage a copy over in Stats.

```
rsync -av --exclude=*RData --exclude=*gz /well/davies/shared/motif_death_analysis/hatbag smew:/data/smew1/rdavies/
rsync -av --exclude=*RData --exclude=*gz /well/davies/shared/motif_death_analysis/treemix smew:/data/smew1/rdavies/

```


## TODOs <a name="paragraph-todos"></a>

* Test that reference chromosome names match what are input in run config json.
* Split chunk_fastq into two jobs to run simultaneous
* delete fastqs
* Write up notes on wildcard constraints for species, chr, units
* Mapping: Note that the number of cores being used per chromosome is a parameter set somewhere, one option is to make this parameter more visible, another is to break the region into chunks and call the genome in chunks them combine back together. This would run faster and more easily on the cluster but requires some coding.
* Once mapping complete and BAM OK, add in deletion of original fastq files (now done manually)
* Clear out old files in `/well/davies/users/dcc832/primates/hatbag_OLD_TO_DELETE`
* Make HATBAG use far far fewer temporary files in steps B -> C (and later?)
* Check out Step C, gainat, nonRepeat, artiodactyla, took ~6 hours?
* Also can Step B run better when using a cluster? Use far more jobs, 1 core each, split out types of test, run all in parallel?
* Should cluster / snakemake use be incorporated into HATBAG itself? Do I use snakemake in tests for some repo? 

