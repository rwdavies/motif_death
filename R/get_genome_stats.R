#!/usr/bin/env Rscript

library(argparse)

## simple script to look at how callable the reference is

parser <- ArgumentParser()
parser$add_argument("--ref_dir")
parser$add_argument("--ref")
parser$add_argument("--chr_prefix", type = "character", default = "")
parser$add_argument("chrlist", type= "character", nargs='*')

args <- parser$parse_args()
ref_dir <- args$ref_dir
ref <- args$ref
chr_prefix <- args$chr_prefix
chr_nums <- args$chrlist
chrlist <- paste0(chr_prefix, chr_nums)

## use /data/wildmice/ref$ cat NCBIM37_um.fa.amb to figure out the N's
amb <- read.table(file.path(ref_dir, paste0(ref, ".amb")))
amb2=amb[-1,]

## tricky 
ann <- read.table(file.path(ref_dir, paste0(ref, ".ann")),sep="\t")
to_start <- paste0("^0 ", chrlist, " ")
w <- sapply(to_start, function(x) grep(x, ann[, 1]))
m <- t(sapply(strsplit(as.character(ann[w + 1, 1]), " "), I))
annL3 <- data.frame(
    chr = chrlist,
    start = as.numeric(m[, 1]),
    length = as.numeric(m[, 2]),
    N = 0,
    stringsAsFactors = FALSE
)

## loop through each, partition appropriately
prev <- 0
for(ichr in 1:dim(annL3)[1]) {
    chr <- annL3[ichr,1]
    start <- annL3[ichr, "start"] + 1
    end <- start -1 + annL3[ichr, "length"]
    which <- amb2[,1]>=start & amb2[,1]<end
    annL3[ichr, "N"]=sum(amb2[which,2])
}

annL4 <- rbind(
    annL3,
    c("all", 0, sum(annL3[, "length"]), sum(annL3[, "N"]))
)

if (sum(is.na(annL4)) > 0)
    stop("Entries of output are NA")

write.table(
    annL4,
    file = file.path(ref_dir, paste0(ref, ".summary.txt")),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)
