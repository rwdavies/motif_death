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

Run on a head or data node. 
```
run.sh preprocess all local goat
```

### Mapping

Run for each species, jobs submitted to the cluster, this would look like the following

```
run.sh mapping all local goat
```

## Downstream variant calling

Here somehow

### Analysis using HATBAG

See on smew ~/proj/HOTSPOT_DEATH_ANALYSIS/