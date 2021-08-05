library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
config_json_path <- args[1]
order_csv_path <- args[2]

config <- fromJSON(config_json_path)

order_df <- data.frame()
for (species in names(config$SPECIES_LIST)) {
    species_prj <- config$SPECIES_LIST[[species]]$ENA_PRJ
    species_units <- config$SPECIES_LIST[[species]]$RUN_ACCESSION
    species_url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",species_prj,"&result=read_run&fields=run_accession,library_name,nominal_length,fastq_ftp,study_accession,scientific_name,instrument_platform,fastq_bytes,fastq_md5,sample_accession&format=tsv&download=true")
    species_table <- read.table(file = species_url, sep = '\t', header = TRUE)
    species_table$species <- species
    species_table <- species_table[species_table$run_accession %in% species_units, ]
    
    max_in_gb <- max(as.numeric(unlist(strsplit(as.character(species_table[, "fastq_bytes"]), ";")))) / 1024 / 1024 / 1024
    n_mapping_pieces <- round((max_in_gb / 25) * 200)
    if ((n_mapping_pieces < 10) | (1000 < n_mapping_pieces)) {
        stop("bad assumption!")
    }
    species_table$n_mapping_pieces <- n_mapping_pieces

    species_table[is.na(species_table$nominal_length), 'nominal_length'] <- 200

    order_df <- rbind(order_df, species_table)
}

order_df$flowcell_barcode <- "X1"
order_df$flowcell_lane <- "1"
order_df$mapping_queue <- "short.qc@@short.hge"

# Rename columns
names(order_df)[names(order_df) == "library_name"] <- "lb"
names(order_df)[names(order_df) == "nominal_length"] <- "lb_insert_size"
names(order_df)[names(order_df) == "run_accession"] <- "units"

write.csv(order_df, order_csv_path, row.names = FALSE)
message("Order dataframe saved in ", order_csv_path)
