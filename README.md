motif_death
===========

## About

This repository contains code to download, pre-process, map, and variant call related species for use in motif death research.

## Workthrough / explanation

### Manual curation of input file and downloading reference material

Briefly, currently, identify a set of related species, appropriate for the model, and papers with their raw data, with data deposited in the EBI ENA. 

For example, for artiodactyla (includes ruminants and related species), this might include cows and goats, and more distantly related animals like okapi and giraffe.

For a given deposition project ID like PRJNA313910, we can then go to e.g. https://www.ebi.ac.uk/ena/browser/view/PRJNA313910 and manually download a text file with information like the download links to the fastq (raw data) files.

In `prepare_inputs.R` and `prepare_inputs_functions.R` we extract information about the raw data in a semi-automated fashion.

```
R -f prepare_inputs.R --args artiodactyla
```

Note, dependent R packages are not currently captured, and will cause program failure at some point for example R package `data.table`.

### Downloading per-sample raw data

Analysis is done per-sample for all samples in an order using the `run.sh` file, which requires manual intervention to note how species are linked to an order. For example `goat` is in `artiodactyla`, specified in `run.sh`.

To use `run.sh` further requires `Snakefile_reference_artiodactyla`, with manual specification of things like the download path for the reference genome, the species, etc.

Note that run.sh pulls sources `activate`, which contains some manual paths, specific to the current machine and user. Similarly, `Snakefile_programs` contains a mixture of manual paths, which should be moved to the activate script, and some programs set depending on paths that depend on the activate script, which is OK. 

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

TODO: `get_per_sample_average_cor.R` requires manual intervention to set list of chromosomes for new species

## Analysis using HATBAG

TODO: get this to run on rescomp?

See on smew ~/proj/HOTSPOT_DEATH_ANALYSIS/

TODO (for Robbie): Clear out old files in `/well/davies/users/dcc832/primates/hatbag_OLD_TO_DELETE`b

