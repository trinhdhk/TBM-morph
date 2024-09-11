library(tidytable)

# Helper functions ----
## Function to parse subcortical volume data ----
read_volume <- function(file){
  cat(file, '\n')
  xml <- XML::xmlParse(file) |> XML::xmlToList()
  info <- strsplit(basename(file), split='_')
  tissues <- xml$labels |> 
    unname() |> 
    purrr::transpose() |> 
    do.call(data.table, args=_) |>
    mutate(across(.fns=unlist))
  cbind(subject = sapply(info, "[[", 1),
        date = sapply(info, "[[", 2) |> anytime::anydate(),
        file=basename(file)|>gsub('\\.xml', '', x=_), tissues)
}

read_volumes <- function(xmldir){
  timepoint=vector('list',2)
  
  xmls = dir(xmldir,  full.names=TRUE)
  do.call(rbind,  lapply(xmls, read_volume))
}

## Function to parse intra-crainial data ----
read_btiv <- function(file){
  tab = read.table(file, sep = " ")[,1:3]
  tab = tab|> as_tidytable()
  info <- strsplit(tab[[1]], split='_')
  colnames(tab) = c('file', 'bTIV_voxel', 'bTIV_volume')
  tab$file <- gsub(".nii.gz", "", tab$file)
  cbind(subject=sapply(info, "[[", 1),
        date = sapply(info, "[[", 2) |> anytime::anydate(),
        tab) |> as_tidytable()
}


vols <- read_volumes('data/mri/volumes')
bTIV <- read_btiv('data/metadata/bTIV.txt')
vols <- left_join(vols, bTIV)
vols <- vols[!subject %in% c('215', '2061')] |> arrange(subject, number, name)
vols <- vols[, .(subject, date, 
                 number, name, 
                 vol = as.numeric(volumeCat)/bTIV_volume)
             ][, name := gsub('(Left\\s)|(Right\\s)', '', x=name, perl=T)] |> 
  unique()
# Remove the non-brain
vols <- vols[!grepl('Non-?Brain', name)
             ][!name %in% c('vessel', 
                            'Ventricular Lining', 
                            'Lesion',
                            "Cerebellum Exterior", 
                            "Cerebral Exterior")]
vols <- vols |>
  mutate(across(c(subject, number), as.integer)) |>
  arrange(subject, date, number) |>
  dt(, 
     data.table::let(
       number = number[[1]],
       vol = sum(vol)
  ), by=list(subject, date, name)) |>
  unique()

write.csv(vols, 'data/imported/volumetry.csv',
          row.names = FALSE)
saveRDS(vols, 'data/imported/volumetry.RDS')
