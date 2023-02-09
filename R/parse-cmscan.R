#functions for parsing and manipulating cmscan data

#' read cmscan data
#' function to load in cmscan data into tibble format
#' @param path
#'
#' @return
#' @export
#'
#' @examples
#' read_cmscan('~/Downloads/rna_infernal/TARA_ANE_2018_MAG_00101.tbl')
#'
read_cmscan <- function(inPath, fileExt = 'tbl') {
  if(file.exists(inPath)){
  name <- gsub(paste0('.',fileExt),'',basename(inPath))
  cols <- c(
    "idx","target_name","target_accession",
    "query_name","query_accession","clan_name",
    "mdl","mdl_from","mdl_to",
    "seq_from","seq_to","strand",
    "trunc","pass","gc",
    "bias","score","E_value",
    "inc","olp","anyidx",
    "afrct1","afrct2","winidx",
    "wfrct1","wfrct2","description_of_target"
  )

  suppressWarnings(
    suppressMessages(
      {cm = readr::read_table(inPath, comment = "#", col_names = cols)}
    )
  )
  cm <- tibble::tibble(cbind(name=name, cm))
  message('Loading ', nrow(cm) ,' cmscan hits',name )
  return(cm)
  }else{
    message(paste0('!!! failed to read cm file - file does not exist, please check the path: ',inPath))
  }
}




#' Multiread cmscan
#'
#' @param dirPath
#' @param fileExt
#'
#' @return
#' @export
#'
#' @examples
#' cm.multi <- multiread_cmscan(dirPath = '~/Downloads/rna_infernal/', fileExt = 'tbl')
multiread_cmscan <- function(dirPath, fileExt = 'tbl'){
  files <- list.files(full.names = TRUE, path = dirPath , pattern = paste0('*',fileExt))
  cms.multi <- c()
  for (f in files){
    cm <- read_cmscan(f)
    cms.multi <- rbind(cms.multi,cm)
  }
  return(tibble::tibble(cms.multi))
}


#Function to filter on targetname
cm_filter <- function(cm, query='rRNA'){
  tmp <- cm %>% dplyr::filter(grepl(query, target_name))
  return(tmp)
}



cm %>% dplyr::rename(start = seq_from,
                     end = seq_to)



# RUN EXAMPLES
#read cmscan data by pointing to file

cm <- read_cmscan(inPath = '/Users/fdb/DATA/Genomes/mrum-genome.tblout')
cm <- read_cmscan(inPath = '~/Downloads/rna_infernal/TARA_ANE_2018_MAG_00101.tbl')

#load multiple cmscan data files in single tibble
cm.multi <- multiread_cmscan(dirPath = '~/Downloads/rna_infernal/', fileExt = 'tbl')

#filter for rRNA or SSU/LSU
cm_filter(cm.multi, query='rRNA') %>% dplyr::select(name, target_name)
cm_filter(cm.multi, query='SSU') %>% dplyr::select(name, target_name)
cm_filter(cm.multi, query='LSU') %>% dplyr::select(name, target_name)
cm_filter(cm.multi, query='LSU')  %>% filter(name == 'Metagenome_centric_SAG_TOSAG00_12')


cm_filter(cm, query='LSU_rRNA_bacteria|SSU_rRNA_bacteria') %>% dplyr::select(name, target_name,seq_from, seq_to)
