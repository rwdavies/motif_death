# library(stringr)
# library(rotl)
# library(datelife)
# library(rvest)
library(ggplot2)
library(ggtree)
library(ggimage)
# library(taxize)

# Note: Running on a large tree eg. Chordata, orders will be slow.
# Will implement writing api call result to Rdata.

# How to set ncbi api key, for taxize::ncbi_downstream() to run faster:
# Run in console: taxize::use_entrez()
# Create an account. Click account name -> account settings to get a key.
# Run in console: usethis::edit_r_environ()
# Add the following line: ENTREZ_KEY=api_key_here
# Restart R

root <- "Carnivora"

leaf_rank <- "family" # Draw tree until this level
motif_death_dir <- "~/proj/motif_death/"

root_ott_id <- rotl::tnrs_match_names(root)$ott_id
root_rank <- attributes(rotl::tnrs_match_names(root))$original_response$results[[1]]$matches[[1]]$taxon$rank

#### Step 1: Read in progress_track.csv ####

progress <- read.table(paste0(motif_death_dir,'progress_track.csv'), header = TRUE, sep="\t")
# For each config we will take the first run_accession of the first species as the representative
progress[,'run_accession'] <- NA

# Much faster to query all at once then loop over one by one
url <- "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query="
for (i in 1:length(progress[,'config_name'])){
  config <- rjson::fromJSON(file = paste0(motif_death_dir, "config/", progress[i, 'config_name']))
  run_accession <- config$SPECIES_LIST[1][[1]]$RUN_ACCESSION[1]
  progress[i,'run_accession'] <- run_accession
  if (i != length(progress[,'config_name'])){
    url <- paste0(url, "run_accession=%22", run_accession, "%22%20OR%20")
  }else{
    url <- paste0(url, "run_accession=%22", run_accession, "%22&fields=scientific_name")
  }
}

scientific_names <- read.table(url, header=TRUE, sep='\t')

df <- rotl::tnrs_match_names(scientific_names[,'scientific_name'])
rownames(df) <- scientific_names[,'run_accession']

leaf_parent <- datelife::get_ott_clade(ott_ids=df[,'ott_id'],ott_rank=leaf_rank)[[leaf_rank]]
root_parent <- datelife::get_ott_clade(ott_ids=df[,'ott_id'],ott_rank=root_rank)[[root_rank]]
df <- cbind(df, leaf_parent, root_parent)

df_sorted <- df[progress[,'run_accession'],]
progress <- cbind(progress, df_sorted)

#### Step 2: Get list of all taxons of leaf_rank from root ####
# Sources: taxize::ncbi_downstream, datelife::get_ott_children, progress[,'leaf_rank']

get_ncbi <- function(root, leaf_rank){
  root_ncbi_id <- taxize::get_ids(root, db='ncbi')$ncbi[[1]]
  ncbi_children <- taxize::ncbi_downstream(root_ncbi_id, downto=leaf_rank)
  return(rotl::tnrs_match_names(ncbi_children[,'childtaxa_name'])$ott_id)
}

get_otl <- function(root, leaf_rank){
  print(root)
  otl_children <- datelife::get_ott_children(ott_ids=root_ott_id, ott_rank=leaf_rank)[[root]]
  return(otl_children[otl_children[,'rank'] == leaf_rank,]$ott_id)
}

otl_ott_ids <- get_otl(root, leaf_rank)
ncbi_ott_ids <- get_ncbi(root, leaf_rank)


get_all_ids <- function(progress, ...){
  valid_configs <- progress[progress[,'root_parent'] == root_ott_id,]
  return(unique(c(valid_configs$ott_id, ...)))
}

all_ids <- get_all_ids(progress, ncbi_ott_ids, otl_ott_ids)


subtree <- tryCatch({
    rotl::tol_induced_subtree(all_ids)
}, error = function(e) {
  x <- e$message
  problems <- lapply(stringr::str_extract_all(x, "[0-9]+"), strtoi)[[1]]
  all_ids <- all_ids[ !(all_ids %in% problems) ]
  return(rotl::tol_induced_subtree(all_ids))
  }
)


tip.color <- rep("black", length(subtree$tip.label))

id_to_label_info <- function(ott_id){
  # Issue is multiple configs with same first run_accession - blue whale
  config <- progress[progress[,'ott_id'] == ott_id,][1,]
  colours <- c('blue', 'orange', 'green', 'red')
  names(colours) <- c('not_run', 'running', 'success', 'failed')
  return(c(config$config_name, colours[[config$status]]))
}


for (i in 1:length(tip.color)){
  id <- strtoi(stringr::str_extract(subtree$tip.label[i], "[0-9]+"))
  if ( id %in% progress$ott_id ){
    label_info <- id_to_label_info(id)
    print(subtree$tip.label[i])
    subtree$tip.label[i] <- paste(subtree$tip.label[i], label_info[1])
    
    tip.color[i] <- label_info[2]
  }
}

scale <- 0.2*length(tip.color)
pdf(file = paste0('~/Downloads/', root, '_', leaf_rank, '.pdf'), width=1*scale, height=scale)
ape::plot.phylo(subtree, tip.color = tip.color)
dev.off()

# # Download images for subtree$tip.label
# current_files <- list.files(path=img_dir)
# for (i in subtree$tip.label){
#   
#   tryCatch({
#     if( !(paste0(i, '.jpg') %in% current_files)){
#       name <- strsplit(i,'_')[[1]][1]
#       download_img(name, i)
#     }
#   },
#   error=function(e) {print(e)})
# }
# 
# # Input single string
# image_from_label <- function(lab){
#   images <- list.files(path=img_dir)
#   filename <- paste0(lab, '.jpg')
#   if (filename %in% images){
#     return(filename)
#   }else{
#     return('blank.jpg')
#   }
# }

# get_image_link <- function(keyword){
#   page <- rvest::read_html(paste0("https://en.wikipedia.org/wiki/", keyword))
#   img_container <- rvest::html_nodes(page, "tr:nth-child(2) img")[1]
#   img_link <- rvest::html_attr(img_container, "src")
#   return(substring(img_link, 3))  # lol engineering pi
# }
# 
# # Output: img_dir/filename.jpg
# download_img <- function(keyword, filename){
#   download.file(get_image_link(keyword), paste0(img_dir, filename, '.jpg'), mode= 'wb')
# }




