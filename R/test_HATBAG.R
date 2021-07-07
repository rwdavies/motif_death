#!/usr/bin/env Rscript

library("testthat")

args <- commandArgs(trailingOnly = TRUE)
HATBAG_OUTPUT_DIR = args[1]

load(file.path(HATBAG_OUTPUT_DIR, "/D_summary/allp.losslin.K6.RData")) # path hardcoded in Snakefile_reference_test
# Assumes test k-mer length 6 and lost k-mers are "AGACAT", "GCGTCC"

m2.dropNA <- m2[(!is.na(m2[, "test_pop1"])) & (!is.na(m2[, "test_pop2"])), ]
SigKmers <- m2.dropNA[(m2.dropNA[, "test_pop1"] < (0.05 / 4^6)) | (m2.dropNA[, "test_pop2"] < (0.05 / 4^6)), ]
print(SigKmers)

test_bool <- setequal(rownames(SigKmers), c("AGACAT", "GCGTCC"))

if (as.numeric(test_bool) == 1) {
    cat("test passed; AGACAT and GCGTCC are the only significant k-mers.\n")
    quit(status = 0)
} else {
    cat("test failed; AGACAT and GCGTCC are not the only significant k-mers. See above.\n")
    quit(status = 1)
}
