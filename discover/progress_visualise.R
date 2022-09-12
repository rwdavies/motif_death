validate_json <- function(){
  progress_track <- jsonlite::fromJSON(file.path(motif_death_dir, "progress_track.json"))
  progress <- progress_track$main
  
  # Check for duplicated rep species
  rep_species <- progress$representative_species
  rep_species_ne <- rep_species[rep_species != ""]
  if (TRUE %in% duplicated(rep_species_ne)){
    print("validate_json: Duplicated representative species present")
    print(rep_species_ne[duplicated(rep_species_ne)])
  }
  
  # Check for invalid statuses
  status <- progress$status
  
  invalid_status <- status != "" & !(status %in% progress_track$color_dict$names) 
  
  if (TRUE %in% invalid_status){
    print("validate_json: Invalid status present")
    print(progress[invalid_status,])
  }
  
  # Check for reps that don't match. Could just rely on rotl warning, but that feels hacky
  rep_taxons <- rotl::tnrs_match_names(rep_species_ne)
  unmatched <- is.na(rep_taxons$ott_id)
  if (TRUE %in% unmatched){
    print("validate_json: Unmatched representative species")
    print(rep_taxons[unmatched,])
  }
  
  # Check for reps that don't have status
  no_status <- !(rownames(progress)[rep_species != ""] %in% rownames(progress)[status != ""])
  if (TRUE %in% no_status){
    print("validate_json: There are representative species with no status")
    progress[rownames(progress)[rep_species != ""][no_status],]
  }
  
}

# Global variables available: motif_death_dir, ncbi_names
make_plot <- function(root, leaf_rank, plot_dir="~/Downloads", plot_scale=0.2, max_iter=100){
  #### Step 1: Read in and prune progress_track.json ####
  progress_track <- jsonlite::fromJSON(file.path(motif_death_dir, "progress_track.json"))

  color_dict <- progress_track$color_dict$colors
  names(color_dict) <- progress_track$color_dict$names
  
  progress <- progress_track$main
  progress <- progress[,colnames(progress) != 'notes']
  progress <- progress[progress$representative_species != "",]
  
  rep_taxons <- rotl::tnrs_match_names(progress[,'representative_species'])
  # Won't work if same rep used twice
  progress <- cbind(progress, ott_id = rep_taxons$ott_id)
  
  tax_info <- rotl::taxonomy_taxon_info(progress$ott_id, include_lineage = TRUE)
  lineage_info <- rotl::tax_lineage(tax_info)
  
  root_rank <- attributes(rotl::tnrs_match_names(root))$original_response$results[[1]]$matches[[1]]$taxon$rank
  
  for (i in 1:length(progress[,1])){
    lin <- lineage_info[[as.character(progress[i,'ott_id'])]]
    root_parent <- lin[tolower(lin$rank)==tolower(root_rank),][1,]$ott_id
    leaf_parent <- lin[tolower(lin$rank)==tolower(leaf_rank),][1,]$name
    progress[i,'root_parent'] <- root_parent
    progress[i,'leaf_parent'] <- leaf_parent
  }
  
  message("make plot: dropping entries with NA as root_parent")
  print(progress[is.na(progress$root_parent),])
  progress <- progress[!is.na(progress$root_parent),]
  
  root_ott_id <- rotl::tnrs_match_names(root)$ott_id
  progress <- progress[progress$root_parent == root_ott_id,]
  
  #### Step 2: Get list of all taxons of leaf_rank from root ####
  downstream_ott_ids <- get_downstream(root, leaf_rank)
  
  #### Step 3: Create the tree (ape phylo object) ####
  all_ids <- c(downstream_ott_ids, progress$ott_id)
  
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
  
  #### Step 4: Make plot ####
  default_labels <- subtree$tip.label
  # Conversion input is id only
  default_ids <- strtoi(stringr::str_extract(default_labels, "[0-9]+"))
  
  new_labels <- rep(NA, length(default_labels))
  new_colors <- rep(NA, length(default_labels))
  
  for (i in 1:length(default_ids)){
    id <- default_ids[i]
    message("make_plot: Converting id ", i, "/", length(default_ids), ": ", id)
    if (id %in% progress$ott_id){
      row <- progress[progress$ott_id == id, ]
      new_labels[i] <- paste0(row$leaf_parent, " (", row$common_name, ")")
      new_colors[i] <- color_dict[[row$status]]
    } else {
      lab <- default_labels[i]
      lab <- strsplit(lab, '_')[[1]]
      if (length(lab) >= 2){
        lab <- lab[-length(lab)]
      }
      lab <- paste(lab, collapse=" ")
      new_labels[i] <- paste0(lab, " (", ott_id_to_rep(id), ")")
      new_colors[i] <- 'black'
    }
  }
  
  subtree$tip.label <- new_labels
  
  scale <- plot_scale*length(new_labels)
  pdf(file = file.path(plot_dir, paste0(root, '_', leaf_rank, '.pdf')), width=scale, height=scale)
  ape::plot.phylo(subtree, tip.color = new_colors)
  dev.off()

  
}

get_downstream <- function(root, leaf_rank){
  cache <- readRDS(file.path(motif_death_dir, "discover/resources/get_downstream.RData"))
  entry_name <- tolower(paste0(root, '_', leaf_rank))
  
  if (entry_name %in% names(cache)){
    message("get_downstream: Saved api call exists")
    output <- cache[[entry_name]]
  } else {
    message("get_downstream: No saved api call exists, will save. This step may take some time")

    root_itis_id <- taxize::get_ids(root, db='itis')$itis
    api_out <- taxize::itis_downstream(root_itis_id, downto=leaf_rank)
    output <- rotl::tnrs_match_names(api_out$taxonname)$ott_id  

    output <- output[!is.na(output)]
    cache[[entry_name]] <- output
    saveRDS(cache, file.path(motif_death_dir, "discover/resources/get_downstream.RData"))
  }
  message("get_downstream: Output is")
  print(output)

  return(output)
}


#### Step 4: A function to convert scientific name to common name ####
sci2comm_fast <- function(scientific_name){
  scientific_name <- tolower(scientific_name)
  if (scientific_name %in% ncbi_names$name){
    # Take first match
    id <- ncbi_names[ncbi_names$name == scientific_name,][1,"ncbi_id"]
    possible <- ncbi_names[ncbi_names$ncbi_id==id,]
    
    comm <- possible[possible$name_type == 'comm',]
    return(comm[1,"name"])
  } else {
    return(NA)
  }
}


#### Step 5: Function to convert ott_id into representative (or NA) ####
# cache = c("polar bear", "emerald rockcod")
# names(cache) = c("123456", "654321")
# Cache file: motif_death/discover/resources/ott_id_to_rep.RData
ott_id_to_rep <- function(ott_id, max_iter=100, force_redo_cache=FALSE){
  ott_id <- as.character(ott_id)
  cache <- readRDS(file.path(motif_death_dir, "discover/resources/ott_id_to_rep.RData"))
  
  if (ott_id %in% names(cache) & !force_redo_cache){
    message("ott_id_to_rep: In cache, converted ", ott_id, " to ", cache[[ott_id]])
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
      for (i in sample(species)[1:min(length(species), max_iter)]){
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
    
    cache[[ott_id]] <- rep
    saveRDS(cache, file=file.path(motif_death_dir, "discover/resources/ott_id_to_rep.RData"))
    message("ott_id_to_rep: Not in cache, converted ", ott_id, " to ", rep)
    return(rep)
  }
}
