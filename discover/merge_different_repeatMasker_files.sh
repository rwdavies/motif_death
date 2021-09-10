#!/bin/bash

cd /well/davies/shared/motif_death_analysis/external/

## vulpes 
ref_name="ASM1834538v1"
key="rmask"

ref_name="Dog10K_Boxer_Tasha"
key="rmask"

gunzip -c ${ref_name}.${key}.SINE.gz | head -n1 > ${ref_name}.rmsk

files=$(ls ${ref_name}.${key}*gz)
for f in $files
do
    echo $f
    gunzip -c ${f} | tail -n+2 >> ${ref_name}.rmsk
done

wc -l ${ref_name}.rmsk

sed -i '' -e '1s/#chrom/#genoName/' -e '1s/chromStart/genoStart/' -e '1s/chromEnd/genoEnd/' -e '1s/name/repName/' $ref_name.rmsk

gzip $ref_name.rmsk


## some R code to change chromosome names
a <- data.table::fread(cmd = "gunzip -c /well/davies/shared/motif_death_analysis/external/Dog10K_Boxer_Tasha.rmsk.gz", data.table = FALSE)w

