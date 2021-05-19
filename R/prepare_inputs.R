setwd("/Users/robert davies/Dropbox/Hotspot Death/Other ENA excel")

source("~/Dropbox/Hotspot Death/Other ENA excel/general_functions.R")

## choose
info <- get_canid_info()
info <- get_bird_info()
info <- get_salmon_info()
info <- get_bat_info()
info <- get_bear_info()

i_species <- 1

lapply(1:length(info), function(i_species) {
    ##
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
        if (is.na(lb)) {
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
    cat(json_to_out, file = paste0(sample_name, ".json"))
    ##
    gb <- sum(as.numeric(unlist(strsplit(as.character(b[, "fastq_bytes"]), ";")))) / 1024 / 1024 / 1024
    return(c(gb, sample_name))
})


system(paste0("rsync -av *.json rescomp:/users/flint/rwdavies/personal/proj/primates/"))



quit()



s <- data.table::fread("~/Downloads/speclist.txt", skip = 60, sep = "&", data.table = FALSE)
## re-format
c <- c(grep(":", s[, 1]), nrow(s) + 1)
out <- array("", c(length(c), 3))
colnames(out) <- c("N", "S", "C")

for(i_c in 1:length(c)) {
    w <- (c[i_c]:(c[i_c + 1] - 1))
    local <- s[w, ]
    local[1] <- substr(local[1], 18, 1000)
    ##a <- sapply(strsplit(s[w, ], ":"), function(x) x[length(x)])
    b <- sapply(strsplit(local, "="), I)
    b[1, ] <- gsub(" ", "", b[1, ])
    for(col in c("N", "S", "C")) {
        if (sum(b[1, ] == col) > 0) {
            out[i_c, col] <- b[2, b[1, ] == col]
        }
    }
}
p


g <- read.csv("~/Downloads/genomes.csv")
## add on matches if possible!
g$commonS <- out[match(g[, 1], out[, "N"]), "S"]
g$commonC <- out[match(g[, 1], out[, "N"]), "C"]

g2 <- g[grep("Animals", g[, 2]), ]

  Eukaryota;Animals;Amphibians                                                        2
  Eukaryota;Animals;Birds                                                             9
  Eukaryota;Animals;Fishes                                                           19
  Eukaryota;Animals;Flatworms                                                         1
  Eukaryota;Animals;Insects                                                          16
  Eukaryota;Animals;Mammals                                                          32
  Eukaryota;Animals;Other Animals                                                     3
  Eukaryota;Animals;Reptiles                                                          2
  Eukaryota;Animals;Roundworms                                                        6



who <- "Eukaryota;Animals;Reptiles"
who <- "Eukaryota;Animals;Insects"
who <- "Eukaryota;Animals;Fishes"
a <- g[g[, 2] == who & g[, "Chromosomes"] > 0, ]

a[sort(unlist(apply(a, 2, function(x) grep("Salm", x)))), ]






table(as.character(g2[, 2]))
## 
g[grep("Animals", g[, 2]), ][1:10, ]


system("curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR304/SRR304976/SRR304976.sra")

system("curl -O ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR126/SRR1264540/SRR1264540.sra")
ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/000/SRR1264540/SRR1264540_1.fastq.gz

