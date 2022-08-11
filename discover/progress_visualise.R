library(stringr)
library(rotl)
library(datelife)
library(rvest)
library(ggplot2)
library(ggtree)
library(ggimage)
library(taxize)
root <- "Chordata"

img_dir <- "~/r_projects/species_img/"
motif_config_dir <- "~/r_projects/motif_death/config/"



get_image_link <- function(keyword){
  page <- read_html(paste0("https://en.wikipedia.org/wiki/", keyword))
  img_container <- html_nodes(page, "tr:nth-child(2) img")[1]
  img_link <- html_attr(img_container, "src")
  return(substring(img_link, 3))  # lol engineering pi
}

# Output: img_dir/filename.jpg
download_img <- function(keyword, filename){
  download.file(get_image_link(keyword), paste0(img_dir, filename, '.jpg'), mode= 'wb')
}

#### get: all_ids, investigated_ids ####
root_id <- rotl::tnrs_match_names(root)[1, "ott_id"]
  
sub_taxons <- datelife::get_ott_children(ott_ids = root_id, ott_rank = "order")[[root]]
sub_taxons <- subset(sub_taxons, rank=="order")

# All orders from root
all_ids <- sub_taxons[,"ott_id"]

filenames <- list.files(path=motif_config_dir)
investigated_taxons <- lapply(strsplit(filenames, "_"), function(x) {x[1]})
investigated_taxons <- rotl::tnrs_match_names(investigated_taxons)
investigated_taxons <- investigated_taxons[!is.na(investigated_taxons[,"ott_id"]),]

# All orders we have investigated
investigated_ids <- datelife::get_ott_clade(ott_ids=investigated_taxons[,"ott_id"], ott_rank = "order")$order


# Our 'all_ids' is not actually everything - make sure at least investigated_ids are included
for (id in investigated_ids){
  if (!(id %in% all_ids)){
    all_ids <- c(all_ids, id)
  }
}

#### end ####

# Make a phylo object
subtree <- rotl::tol_induced_subtree(all_ids)

tip.color <- rep("black", length(subtree$tip.label))

for (i in 1:length(tip.color)){
  id <- strtoi(str_extract(subtree$tip.label[i], "[0-9]+"))
  if ( is.element(id, investigated_ids) ){
    tip.color[i] <- "green"
  }
}


# Download images for subtree$tip.label
current_files <- list.files(path=img_dir)
for (i in subtree$tip.label){
  
  tryCatch({
    if( !(paste0(i, '.jpg') %in% current_files)){
      name <- strsplit(i,'_')[[1]][1]
      download_img(name, i)
    }
  },
  error=function(e) {print(e)})
}

# Input single string
image_from_label <- function(lab){
  images <- list.files(path=img_dir)
  filename <- paste0(lab, '.jpg')
  if (filename %in% images){
    return(filename)
  }else{
    return('blank.jpg')
  }
}

height = 1*length(subtree$tip.label)
pdf(file = paste0("~/Downloads/", root, ".pdf"), height = height, width=height)
ggtree(subtree) + 
  xlim(0, 0.24*height) +
  geom_tiplab(geom='label', color=tip.color) +
  geom_tiplab(aes(image=paste0(img_dir, sapply(label, image_from_label))), geom='image', offset=1, size=.01)
  geom_nodelab()
dev.off()

# 
# # Make plot
# height = 1*length(subtree$tip.label)
# pdf(file = paste0("~/Downloads/", root, ".pdf"), height = height, width = height)
# 
# ape::plot.phylo(subtree, tip.color = "white", show.node.label = TRUE, show.tip.label = TRUE, srt = 90)
# par(new = TRUE)
# ape::plot.phylo(subtree, tip.color=tip.color)
# 
# dev.off()

