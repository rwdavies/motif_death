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

### Downloading data

Analysis is done using the `run.sh` file, which requires manual intervention to note how species are linked to an order. For example `goat` is in `artiodactyla`, specified in `run.sh`. To use `run.sh` further requires `Snakefile_reference_artiodactyla`, with manual specification of things like the path to the reference genome, the species, etc.

Eventually, `run.sh` is actually run, on a head or data node. It requires snakemake to be installed, currently done using `install_snakemake_2.sh` semi-interactively through conda
```
run.sh preprocess all local goat --dryrun ## to check what will be run
run.sh preprocess all local goat ## to run it
```
Note this might take from 30 minutes to a few hours of wall clock time to run per job, depending on size of individual fastqs.

### Mapping

Run for each species, jobs submitted to the cluster, this would look like the following

```
run.sh mapping all cluster goat
```

This might take one or a few days to run through to completion.

## Downstream variant calling

Here somehow

### Analysis using HATBAG

See on smew ~/proj/HOTSPOT_DEATH_ANALYSIS/