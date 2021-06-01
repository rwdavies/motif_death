## setwd("/Users/robert davies/Dropbox/Hotspot Death/Other ENA excel")

##
## workflow to add a species
## 1) identify one or more papers / ENA projects with data of interest
## 2) download relevant .tsv or .txt files from ENA with information about the raw data (e.g. for project ID PRJNA313910, go to https://www.ebi.ac.uk/ena/browser/view/PRJNA313910, select all columns, and manually download tsv file with information
## 3) manually make "get_<species>_info" with manually curated information from the above
##   a) make a list with one character vector entry per species in some arbitrary order
##   b) for each species the character vector has length 5 with the following elements in the following order
##      - name of the sample in the ENA file (sample_alias column)
##      - scientific name
##      - working name for the analysis (lower case, simple form)
##      - file to read in
##      - number of mapping pieces
##   c) manually edit "get_b" function that loads ENA file and sub-selects specific entries from ENA file, e.g. possibly taking a subset for space reasons, and only taking Illumina short read ideally high quality HiSeq PE 2x150 or similar (rather than mate pair, etc)
## 3 - appendix - a URL like this might work to download PRJ file https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA317745&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true

## 


R_DIR <- Sys.getenv("R_DIR")
ENA_DIR <- Sys.getenv("ENA_DIR")

if (1 == 0) {
    
    ANALYSIS_DIR='/well/myers/rwdavies/primates/'
    R_DIR="~/proj/motif_death/R/"
    ENA_DIR="/well/davies/shared/motif_death_resources/Other ENA excel/"
    group <- "artiodactyla"
    
}

source(file.path(R_DIR, "prepare_inputs_functions.R"))

## choose
group <- commandArgs(trailingOnly = TRUE)
if (group == "artiodactyla") {
    info <- get_artiodactyla_info()
} else {
    stop("fix this file!")
}
## info <- get_canid_info()
## info <- get_bird_info()
## info <- get_salmon_info()
## info <- get_bat_info()
## info <- get_bear_info()


i_species <- 1


lapply(1:length(info), function(i_species) {
    ##
    print(i_species)
    sample_alias <- info[[i_species]][1]
    sample_name <- info[[i_species]][3]
    paper <- info[[i_species]][4]
    n_mapping_pieces <- as.integer(info[[i_species]][5])
    ##
    b <- get_b(paper, sample_alias)
    b[, c(grep("run_accession", colnames(b)), grep("fastq", colnames(b)))] 
    ##
    print(sample_name)
    print(b[, c("instrument_model", "sample_accession", "nominal_length", "nominal_sdev", "read_length", "sample_alias", "scientific_name", "fastq_bytes")])
    ##    
    units <- lapply(1:nrow(b), function(i) {
        if (
        ("fastq_md5" %in% colnames(b)) &
        ("fastq_ftp" %in% colnames(b))        
        ) {
            a1 <- basename(strsplit(as.character(b[i, "fastq_ftp"]), ";")[[1]])
            a2 <- strsplit(as.character(b[i, "fastq_md5"]), ";")[[1]]
            md5_1 <- a2[grep("_1.fastq", a1)]
            md5_2 <- a2[grep("_2.fastq", a1)]            
        } else {
            md5_1 <- NA
            md5_2 <- NA            
        }
        x <- strsplit(as.character(b[i, "fastq_ftp"]), ";")[[1]]
        x2 <- sapply(x, basename)
        nominal_length <- b[i, "nominal_length"]
        if (is.na(nominal_length)) {
            nominal_length <- 200
        }
        lb <- b[i, "library_name"]
        if (is.na(lb) | lb == "") {
            lb <- paste0("lb", i)
        }
        o <- list(
            "1" = x[grep("_1.fastq.gz", x2)],
            "2" = x[grep("_2.fastq.gz", x2)],
            "md5_1" = md5_1,
            "md5_2" = md5_2,            
            "lb" = lb,
            "lb_insert_size" = as.character(nominal_length),
            "flowcell_barcode" = paste0("X", i),
            "flowcell_lane" = as.character(i)
        )
    })
    names(units) <- as.character(b[, "run_accession"])
    ##
    if (("fastq_bytes" %in% colnames(b)) & (FALSE == (sample_name %in% c("coho", "chinook")))) {
        max_in_gb <- max(as.numeric(unlist(strsplit(as.character(b[, "fastq_bytes"]), ";")))) / 1024 / 1024 / 1024
        n_mapping_pieces <- round((max_in_gb / 25) * 200)
        if ((n_mapping_pieces < 10) | (1000 < n_mapping_pieces)) {
            stop("bad assumption!")
        }
    }
    ##
    library("jsonlite")
    o <- list(
        "species" = sample_name,
        "n_mapping_pieces" = n_mapping_pieces,
        "platform" = "Illumina",
        "mapping_queue" = "short.qc",
        "units" = units
    )
    ##
    json_to_out <- toJSON(o, pretty = TRUE, auto_unbox = TRUE)
    print(    json_to_out)
    cat(json_to_out, file = file.path("species_mapping_info", paste0(sample_name, ".json")))
    ##
    gb <- sum(as.numeric(unlist(strsplit(as.character(b[, "fastq_bytes"]), ";")))) / 1024 / 1024 / 1024
    return(c(gb, sample_name))
})




quit()



## s <- data.table::fread("~/Downloads/speclist.txt", skip = 60, sep = "&", data.table = FALSE)
## ## re-format
## c <- c(grep(":", s[, 1]), nrow(s) + 1)
## out <- array("", c(length(c), 3))
## colnames(out) <- c("N", "S", "C")

## for(i_c in 1:length(c)) {
##     w <- (c[i_c]:(c[i_c + 1] - 1))
##     local <- s[w, ]
##     local[1] <- substr(local[1], 18, 1000)
##     ##a <- sapply(strsplit(s[w, ], ":"), function(x) x[length(x)])
##     b <- sapply(strsplit(local, "="), I)
##     b[1, ] <- gsub(" ", "", b[1, ])
##     for(col in c("N", "S", "C")) {
##         if (sum(b[1, ] == col) > 0) {
##             out[i_c, col] <- b[2, b[1, ] == col]
##         }
##     }
## }
## p


## g <- read.csv("~/Downloads/genomes.csv")
## ## add on matches if possible!
## g$commonS <- out[match(g[, 1], out[, "N"]), "S"]
## g$commonC <- out[match(g[, 1], out[, "N"]), "C"]

## g2 <- g[grep("Animals", g[, 2]), ]

##   Eukaryota;Animals;Amphibians                                                        2
##   Eukaryota;Animals;Birds                                                             9
##   Eukaryota;Animals;Fishes                                                           19
##   Eukaryota;Animals;Flatworms                                                         1
##   Eukaryota;Animals;Insects                                                          16
##   Eukaryota;Animals;Mammals                                                          32
##   Eukaryota;Animals;Other Animals                                                     3
##   Eukaryota;Animals;Reptiles                                                          2
##   Eukaryota;Animals;Roundworms                                                        6



## who <- "Eukaryota;Animals;Reptiles"
## who <- "Eukaryota;Animals;Insects"
## who <- "Eukaryota;Animals;Fishes"
## a <- g[g[, 2] == who & g[, "Chromosomes"] > 0, ]

## a[sort(unlist(apply(a, 2, function(x) grep("Salm", x)))), ]






## table(as.character(g2[, 2]))
## ## 
## g[grep("Animals", g[, 2]), ][1:10, ]


## system("curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR304/SRR304976/SRR304976.sra")

## system("curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR126/SRR1264540/SRR1264540.sra")
## ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1264540/SRR1264540_1.fastq.gz

