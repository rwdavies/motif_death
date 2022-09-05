motif_death_dir <- '~/proj/motif_death'

json <- jsonlite::fromJSON(file.path(motif_death_dir, "progress_track.json"))
json <- json$main
json <- json[,c(2,1,4,5)]

readme <- file.path(motif_death_dir, 'discover/README.md')
readme_tmp <- file.path(motif_death_dir, 'discover/README.md.tmp')

raw_text <- readChar(readme, file.info(readme)$size)

split <- strsplit(raw_text, '<!--- marker -->', fixed=TRUE)[[1]]

top <- split[1]
top <- substr(top, 1, nchar(top)-1) # Get rid of a \n
middle <- c("| Common name | Config name | Status | Notes |", "|-|-|-|-|")
bottom <- split[3]
bottom <- substr(bottom, 2, nchar(bottom))

for (i in 1:length(json[,1])){
  line <- paste(json[i,], collapse=' | ')
  line <- paste("|", line, "|")
  middle <- c(middle, line)
}

file <- file(readme_tmp)
writeLines(c(top, '<!--- marker -->', middle ,'<!--- marker -->', bottom), file)
close(file)

system(paste("mv", readme_tmp, readme ))
