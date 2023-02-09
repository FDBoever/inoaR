#parse-barrnap
read_barrnap <- function(inPath, fileExt = 'tbl') {
  if(file.exists(inPath)){
    name <- gsub(paste0('.',fileExt),'',basename(inPath))
    bn <- read.delim(inPath,header=FALSE,skip=1)
    colnames(bn) = c('contig','method','type','start','end','eval','strand','n','annotation')
    bn <- bn %>% tibble::tibble() %>%
      dplyr::mutate(contig=as.character(contig)) %>%
      dplyr::mutate(length = end - start) %>%
      tidyr::separate(annotation,sep=';',into=c('Name','product','note')) %>%
      dplyr::mutate(Name=gsub('Name=','',Name),
                    product=gsub('product=','',product),
                    note=gsub('note=','',note))%>%
      mutate(rna_type=ifelse(grepl("16S|18S",Name),'SSU',
                             ifelse(grepl("23S|28S",Name),'LSU','other'))) %>%
      mutate(name=name)
    message('Loading ', nrow(bn) ,' barrnap hits, ',name )
    return(bn)
  }else{
    message(paste0('!!! failed to read barrnap file - file does not exist, please check the path: ',inPath))
  }
}


multiread_barrnap <- function(dirPath, fileExt = 'txt'){
  files <- list.files(full.names = TRUE, path = dirPath , pattern = paste0('*',fileExt))
  brnp.multi <- c()
  for (f in files){
    if(length(readLines(file(f)))>1){
      brnp <- read_barrnap(f)
      brnp.multi <- rbind(brnp.multi,brnp)
      }
    }
  return(tibble(brnp.multi))
}


read_barrnap_bedtools <- function(path){
  bn <- read.delim(path,header=FALSE)
  colnames(bn) = c('source','contig','method','type','start','end','eval','strand','n','annotation','sequence')
  bn <- bn %>% tibble::tibble() %>%
    dplyr::mutate(sequence=as.character(sequence),
                  source=as.character(source),
                  contig=as.character(contig)) %>%
    dplyr::mutate(length = end - start) %>%
    tidyr::separate(annotation,sep=';',into=c('Name','product','note')) %>%
    dplyr::mutate(Name=gsub('Name=','',Name),
                  product=gsub('product=','',product),
                  note=gsub('note=','',note))%>%
    mutate(rna_type=ifelse(grepl("16S|18S",Name),'SSU',
                           ifelse(grepl("23S|28S",Name),'LSU','other')))
  return(bn)
}



# RUN EXAMPLES
#read barrnap data by pointing to file
brnp  <- read_barrnap(inPath='/Users/sa01fd/DATA/rna-miner/metagenme_assembly/out_barnap/final_assembly.fasta.barnap.txt')


#load multiple barrnap data files in single tibble
brnp.multi <- multiread_barrnap(dirPath = '~/Downloads/rna_barnap/',
                                fileExt = 'tbl')


