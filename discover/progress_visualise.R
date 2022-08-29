# For the program to work you need to download a file from ncbi. 
# This file is too large to store on github (~200MB).
# This file contains name information for species, and it is used for getting
# a representative species of a taxon.

# How to download:
# Go to https://ftp.ncbi.nih.gov/pub/taxonomy/
# Download taxdmp.zip and extract
# Run the following in the taxdmp directory:
# cat names.dmp | sed 's/\t//g' | tr "|" "," > names.csv
# Copy names.csv to motif_death/discover/progress_visualise_data




# You only need to do the following step if you plan to do a large tree that is not
# cached. Check out motif_death/discover/progress_visualise_data/ncbi_downstream
# for currently saved api calls.

# How to set ncbi api key, for taxize::ncbi_downstream() to run faster:
# Go to: https://www.ncbi.nlm.nih.gov/account/
# Create an account. Click account name -> account settings to get a key
# Edit .Renviron (eg. usethis::edit_r_environ() )
# Add the following line: ENTREZ_KEY=api_key_here
# Restart R

root <- "Testudines"
leaf_rank <- "family" # Draw tree until this level. Use singular not plural
motif_death_dir <- "~/proj/motif_death/"

root_ott_id <- rotl::tnrs_match_names(root)$ott_id
root_rank <- attributes(rotl::tnrs_match_names(root))$original_response$results[[1]]$matches[[1]]$taxon$rank

#### Step 1: Read in progress_track.json ####

progress <- jsonlite::fromJSON(paste0(motif_death_dir, "progress_track.json"))$main

# Won't work if same rep used twice
rep_taxons <- rotl::tnrs_match_names(progress[,'representative_species'])
progress <- cbind(progress, rep_taxons)

root_parent <- datelife::get_ott_clade(ott_ids=progress[,'ott_id'],ott_rank=root_rank)[[root_rank]]
progress <- cbind(progress, root_parent)

#### Step 2: Get list of all taxons of leaf_rank from root ####
get_ncbi <- function(root, leaf_rank){
  filename <- paste0(tolower( paste0(root, '_', leaf_rank) ), ".RData")
  files <- list.files(path=paste0(motif_death_dir, "discover/progress_visualise_data/ncbi_downstream"))

  if (filename %in% files){
    print("get_ncbi: Saved api call exists")
    api_out <- readRDS( paste0(motif_death_dir, "discover/progress_visualise_data/ncbi_downstream/", filename))
  } else {
    want_save <- readline(prompt = "get_ncbi: No saved api call. Would you like to save? [y/n] ")
    root_ncbi_id <- taxize::get_ids(root, db='ncbi')$ncbi[[1]]
    api_out <- taxize::ncbi_downstream(root_ncbi_id, downto=leaf_rank)
    if (want_save == 'y'){
      saveRDS(api_out, file=paste0(motif_death_dir, "discover/progress_visualise_data/ncbi_downstream/", filename))
    }
  }

  return(rotl::tnrs_match_names(api_out[,'childtaxa_name'])$ott_id)
}

ncbi_ott_ids <- get_ncbi(root, leaf_rank)

#### Step 3: Create the tree (ape phylo object) ####

# Add ott_id of representative if it has the correct root_parent
valid <- progress[progress$root_parent == root_ott_id,]

all_ids <- c(ncbi_ott_ids, valid$ott_id)

# There are sometimes broken ott_ids so parse the error message and remove them
subtree <- tryCatch({
    rotl::tol_induced_subtree(all_ids)
}, error = function(e) {
  x <- e$message
  problems <- lapply(stringr::str_extract_all(x, "[0-9]+"), strtoi)[[1]]
  all_ids <- all_ids[ !(all_ids %in% problems) ]
  return(rotl::tol_induced_subtree(all_ids))
  }
)

#### Step 4: A function to convert scientific name to common name ####
# We use the ncbi taxonomy dump found here: https://www.ncbi.nlm.nih.gov/guide/taxonomy/
# The file is \t|\t delimited. Change this: cat names.dmp | sed 's/\t//g' | tr "|" "," > names.csv
# New file is motif_death/discover/progress_visualise_data/names.csv
# For the future: make into a csv with headers Scientific, Common and load that directly

# names2 is a cleaned up version of names
ncbi_names <- read.csv("~/proj/motif_death/discover/progress_visualise_data/names.csv")
ncbi_names <- ncbi_names[,c(1,2,4)]
colnames(ncbi_names) <- c('ncbi_id', 'name', 'name_type')
ncbi_names$name <- tolower(ncbi_names$name)

sci2comm_fast <- function(scientific_name){
  scientific_name <- tolower(scientific_name)
  if (scientific_name %in% ncbi_names$name){
    # Take first match
    id <- ncbi_names[ncbi_names$name == scientific_name,][1,"ncbi_id"]
    possible <- ncbi_names[ncbi_names$ncbi_id==id,]
    
    comm <- possible[grepl('common', possible$name_type, fixed=TRUE),]
    return(comm[1,"name"])
  } else {
    return(NA)
  }
}

#### Step 5: Function to convert ott_id into representative (or NA) ####

# cache = c("polar bear", "emerald rockcod")
# names(cache) = c("123456", "654321")
# Cache file: motif_death/discover/progress_visualise_data/representatives.RData
ott_id_to_rep <- function(ott_id, max_iter=100){
  ott_id <- as.character(ott_id)
  cache <- readRDS(paste0(motif_death_dir, "discover/progress_visualise_data/representatives.RData"))
  
  if (ott_id %in% names(cache)){
    print(paste("ott_id_to_rep: In cache, converted", ott_id, "to", cache[[ott_id]]))
    return(cache[[ott_id]])
  }else{
    rep <- tryCatch({
      tree <- rotl::taxonomy_subtree(strtoi(ott_id), output_format = 'phylo')
      labels <- tree$tip.label
      species <- c()
      for (i in labels){
        split <- strsplit(i, "_")[[1]]
        split <- split[-length(split)]
        name <- paste(split, collapse=" ")
        species <- c(species, name)
      }
      rep_candidate <- NA
      # Work with species = ['mus musculus domesticus', 'vulpes vulpes']
      for (i in species[1:min(length(species), max_iter)]){
        x <- sci2comm_fast(i)
        if (!is.na(x)){
          rep_candidate <- x
          break
        }
      }
      rep_candidate
    }, error = function(e){
      print(e)
      NA
    })
    
    cache <- c(cache, rep)
    names(cache)[length(names(cache))] <- ott_id
    saveRDS(cache, file=paste0(motif_death_dir, "discover/progress_visualise_data/representatives.RData"))
    print(paste("ott_id_to_rep: Not in cache, converted", ott_id, "to", rep))
    return(rep)
  }
}

#### Step 6: Function to convert ott id into new label and colour  ####
color_dict <- c('blue', 'orange', 'red', 'green', 'purple')
names(color_dict) <- c('ready to go', 'running', 'run failed', 'run success', 'not feasible')

convert_id <- function(ott_id){
  if (ott_id %in% progress$ott_id){
    row <- progress[progress$ott_id == ott_id,]
    return(c(row$common_name, color_dict[[row$status]]))
  }else{
    rep <- ott_id_to_rep(ott_id)
    return(c(rep, "black"))
  }
}

#### Step 7: Make plot ####

default_labels <- subtree$tip.label
# Use only default_ids
default_ids <- strtoi(stringr::str_extract(default_labels, "[0-9]+"))

new_labels <- rep(NA,length(default_labels))
new_colors <- rep(NA,length(default_labels))

for (i in 1:length(default_ids)){
  print(paste0("Converting id ", i, "/", length(default_ids), ": ", default_ids[i]))
  convert_id_out <- convert_id(default_ids[i])
  new_labels[i] <- paste0(default_labels[i], " (", convert_id_out[1], ")")
  new_colors[i] <- convert_id_out[2]
}

subtree$tip.label <- new_labels

scale <- 0.5*length(new_labels)
pdf(file = paste0('~/Downloads/', root, '_', leaf_rank, '.pdf'), width=1*scale, height=scale)
ape::plot.phylo(subtree, tip.color = new_colors)
dev.off()

