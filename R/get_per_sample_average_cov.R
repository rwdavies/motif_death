##library("data.table")



## setwd("/well/davies/users/dcc832/primates/")
species <- commandArgs(trailingOnly = TRUE)[1]
ref_dir <- commandArgs(trailingOnly = TRUE)[2]
ref <- commandArgs(trailingOnly = TRUE)[3]
## species <- "chimp"; ref <- "hg38.fa";
## species <- "caroli"; ref <- "NCBIM37_um.fa";
## species <- "whitetaileddeer"; ref <- "bosTau8.fa"
## species <- "test_pop1"; ref_dir <- "ref/"; ref <- "ref.fa"
ref_summary_file <- paste0(ref_dir, ref, ".summary.txt")


get_chrlist <- function(ref) {
    if (ref == "bosTau8.fa") {
        chrlist <- 1:29
        chr_prefix <- "chr"
    } else if (ref == "hg38.fa") {
        chrlist <- 1:22
        chr_prefix <- "chr"
    } else if (ref == "NCBIM37_um.fa") {
        chrlist <- 1:19
        chr_prefix <- ""
    } else if (ref == "canFam3.fa") {
        chrlist <- 1:38
        chr_prefix <- "chr"
    } else if (ref == "ref.fa") { # Test
        chrlist <- 1:2
        chr_prefix <- ""
    } else {
        print(paste0("ref = ", ref))
        stop("ref not defined")
    }
    return(list(chrlist = chrlist, chr_prefix = chr_prefix))
}
    
chrlist <- get_chrlist(ref)$chrlist

RData_file_function <- function(species, chr)
    paste0("coverage/coverage.", species, ".chr", chr, ".RData")


rebuild <- FALSE
message("First pass")

chr <- chrlist[1]

depth_sum <- 0
for(chr in chrlist) {
    message(paste0(chr, ", ", date()))
    RData_file <- RData_file_function(species, chr)
    if (rebuild | (file.exists(RData_file) == FALSE)) {
        input_file <- paste0("coverage/coverage.", species, ".chr", chr, ".txt.gz")
        if (file.exists(input_file) == FALSE)
            stop(paste0("Cannot find file:", input_file))
        message("load")
        f <- tempfile()
        system(paste0("gunzip -c ", input_file, " | cut -f1,4 > ", f))
        nrows <- system(paste0("wc -l ", f), intern = TRUE)
        nrows <- as.numeric(strsplit(nrows, " ")[[1]][1])
        message(paste0("nrows = ", nrows))
        message(paste0("extract depth"))
        depth <- as.integer(system(paste0("cut -f2 ", f, " | sed 1d"), intern = TRUE))
        message(paste0("extract L"))        
        ## L <- as.integer(system(paste0("cut -f1 ", f, " | sed 1d | cut -c ", 5 + nchar(chr), "-"), intern = TRUE))
        L <- as.integer(system(paste0("cut -f1 ", f, " | sed 1d | cut -f2 --delimiter=':' "), intern = TRUE))
        ## data <- fread(
        ##     f, data.table = FALSE, nrows = nrows,
        ##     colClasses = c("character", "integer"),
        ##     skip = 1
        ## )
        ## unlink(f)
        ## message("extract")
        ## for(i in 1:10) {
        ##     gc(reset = TRUE)
        ## }
        ## L <- as.numeric(substr(data[, 1], 5 + nchar(chr), 100))
        ## gc(reset = TRUE); gc(reset = TRUE)
        ## depth <- data[, 2]
        ## rm(data)
        message("save")
        save(L, depth, file = RData_file, compress = FALSE) ## speed up since a temp fi
        for(i in 1:10) {
            gc(reset = TRUE)
        }
        ##unlink(input_file)        
    } else {
        message("load")
        load(RData_file)
    }
    depth_sum <- depth_sum + sum(as.numeric(depth), na.rm = TRUE)
    rm(L, depth)
    gc(reset = TRUE); gc(reset = TRUE)
}



ref_summary <- read.table(ref_summary_file, header = TRUE)
x <- ref_summary[ref_summary[, "chr"] == "all", ]
mappable_genome_size <- x["length"] - x["N"]
av_cov <- as.numeric(depth_sum / mappable_genome_size)
if (is.na(av_cov)) {
    stop("Average coverage is NA for an unknown reason")
}

cat(av_cov, file = paste0("coverage/average.", species, ".txt"))



message("Second pass")
for(chr in chrlist) {
    message(chr)
    message("load")
    load(file = paste0("coverage/coverage.", species, ".chr", chr, ".RData"))
    message("determine callable")
    callable <- (av_cov / 3) <= depth & depth <= (av_cov * 2)
    rm(depth); gc(reset = TRUE); gc(reset = TRUE)
    ## determine callable, uncallable regions
    message("generate regions")
    d <- diff(callable)
    starts <- c(1, which(d != 0) + 1)
    ends <- c(which(d != 0), length(callable))
    message("make bed and save")
    ## turn into bed file
    bed <- cbind(chr, L[starts] - 1, L[ends], callable[starts])
    callable_bed <- bed[bed[, 4] == 1, 1:3, drop = FALSE]
    write.table(
        callable_bed,
        file = paste0("coverage/coverage.", species, ".chr", chr, ".callableOnly.bed"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
    ##message("remove RData files")
    RData_file <- RData_file_function(species, chr)
    unlink(RData_file)
}

quit(status = 0)
