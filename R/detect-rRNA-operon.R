#detect-rRNA-operon
#load("./R/parse-cmscan.R")

#Compare Barrnap to cmscan


# Speficic code for Hackl2020 example dataset
#cmscan -Z 5.874406
#        --cut_ga
#        --rfam
#        --nohmmonly
#        --tblout GCA_008330645.1_CrBVI_genomic.cmscan.tblout
#        --fmt 2
#        --clanin ~/DATA/annotation_databases/Rfam/Rfam.clanin ~/DATA/annotation_databases/Rfam/Rfam.cm
#        /Users/fdb/DATA/EukGenome_database/Hackl2020/ncbi_dataset/data/GCA_008330645.1/GCA_008330645.1_CrBVI_genomic.fna

# barrnap --kingdom euk
#         --lencutoff 0.6
#         --reject 0.5 /Users/fdb/DATA/EukGenome_database/Hackl2020/ncbi_dataset/data/GCA_008330645.1/GCA_008330645.1_CrBVI_genomic.fna
#         > GCA_008330645.1_CrBVI_genomic.barrnap.txt

brnp  <- read_barrnap(inPath='/Users/fdb/DATA/EukGenome_database/Hackl2020/GCA_008330645.1_CrBVI_genomic.barrnap.txt')



brnp %>%
  dplyr::filter(contig == 'VLTN01000011.1') %>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin =  start , xmax =  end , y = contig, fill = Name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3")

cm <- read_cmscan(inPath = '/Users/fdb/DATA/EukGenome_database/Hackl2020/GCA_008330645.1_CrBVI_genomic.cmscan.tblout')

cm %>%
  dplyr::filter(query_name == 'VLTN01000011.1') %>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin =  seq_from , xmax =  seq_to , y = query_name, fill = target_name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3") + ggplot2::xlim(c(526647,536831))



cm %>%
  dplyr::filter(query_name == 'CM017891.1') %>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin =  seq_from , xmax =  seq_to , y = query_name, fill = target_name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3") + ggplot2::xlim(c(0,10000))








cm.rrna <- cm_filter(cm, query='rRNA')
prox.clust <- cm_proximity_clusters(feature.tbl=cm.rrna, distance = 10000)


out<- c()
for(p.cl in unique(prox.clust$prox.cluster)){

  for(gne in c('SSU', 'LSU')){
    tmp <- prox.clust %>%
      dplyr::filter(prox.cluster == p.cl) %>%
      dplyr::filter(grepl(gne,target_name)) %>%
      dplyr::arrange(E_value,desc(score)) %>%
      dplyr::slice(1)

    out <- rbind(out, tmp)
  }

}





cm_filter_kingdom <- function(cm){


}


brnp  <- read_barrnap(inPath='/Users/fdb/DATA/Genomes/fdb33_barnap.txt')
cm <- read_cmscan(inPath = '/Users/fdb/DATA/Genomes/mrum-genome.tblout')
cm_filter(cm, query='LSU_rRNA_bacteria|SSU_rRNA_bacteria') %>%
  dplyr::select(name, target_name,seq_from, seq_to)%>%
  mutate(length=seq_to-seq_from)


brnp.multi <- multiread_barrnap(dirPath = '~/DATA/EukGenome_database/analysis_out/rna_barrnap', fileExt = 'txt')
cm.multi <- multiread_cmscan(dirPath = '~/DATA/EukGenome_database/analysis_out/rna_infernal', fileExt = 'tbl')


brnp.multi %>% group_by(name) %>% tally() %>% arrange(desc(n))

brnp.multi %>% filter(name=='TARA_ARC_108_MAG_00310.fa.euk.barnap.txt') %>%
  group_by(contig) %>%
  tally() %>%
  arrange(desc(n))

brnp.multi %>% filter(contig=='TARA_SOC_28_MAG_00049_extended_000000299199')



brnp.multi %>%
  dplyr::filter(contig == 'TARA_SOC_28_MAG_00049_extended_000000007781') %>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin =  start , xmax =  end , y = contig, fill = Name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3")







#as in current implementation of run_brnp
brnp.multi <- brnp.multi %>% dplyr::filter(grepl('euk',name))


rna.merged <-  merge_barrnap_cm(brnp.multi, cm.multi)

rna.merged %>%
  ggplot2::ggplot(ggplot2::aes(product, length,color=method)) +
  ggplot2::geom_point() + ggplot2::coord_flip()




merge_barrnap_cm = function(brnp, cm){
  tmp.cm <- cm_filter(cm, query='LSU_rRNA|SSU_rRNA') %>%
    dplyr::rename(start = seq_from,
                  end = seq_to,
                  contig=query_name) %>%
    dplyr::mutate(length=abs(end-start),
                  method='cmscan_RFAM') %>%
    dplyr::select(name, contig,method, target_name, strand, start, end, length) %>%
    mutate(product = ifelse(target_name=='LSU_rRNA_archaea','',
                            ifelse(target_name=='SSU_rRNA_bacteria','16S ribosomal RNA',
                                   ifelse(target_name=='LSU_rRNA_bacteria','23S ribosomal RNA',
                                          ifelse(target_name=='SSU_rRNA_eukarya','18S ribosomal RNA',
                                                 ifelse(target_name=='LSU_rRNA_eukarya','28S ribosomal RNA',target_name
                                                 )))))) %>%
  select(-target_name)

  tmp.brnp <- brnp %>% select(name, contig, method, strand, start, end, length,product)
  tmp.out <- rbind(tmp.brnp, tmp.cm)
  return(tmp.out)
}

merge_barrnap_cm(brnp, cm) %>%
  ggplot2::ggplot(ggplot2::aes(product, length,color=method)) +
  ggplot2::geom_point() + ggplot2::coord_flip()



# RUN EXAMPLES
#read cmscan data by pointing to file
cm <- read_cmscan(inPath = '~/Downloads/rna_infernal/TARA_ANE_2018_MAG_00101.tbl', fileExt = 'tbl')
cm <- read_cmscan(inPath = '~/Downloads/primary_assembly_spades.tblout', fileExt = 'tblout')



#load multiple cmscan data files in single tibble
cm <- multiread_cmscan(dirPath = '~/Downloads/rna_infernal/', fileExt = 'tbl')

#filter for rRNA or SSU/LSU
cm.rrna <- cm_filter(cm, query='rRNA')
cm.rrna %>% dplyr::select(name, target_name)


#sort cm based on query_name and start position of genes
cm <- cm_sort(cm)

#Orient CMs strands
cm_oriented <- cm_orient(cm, seqs.tbl)
cm_oriented <- sortregions(cm_oriented)

#Orient sequences strands
seqs.oriented.tbl <- seq_cm_revcomp(cm, seqs.tbl)


#filter cm for eukaryotic rDNA operon detection, including default lenght filter
#f_euk_cm <- filter_euk_rDNA(cm = cm_oriented, option = 'length_filter')


f_euk_cm <- hardFilter(cm)


#keep contigs that containg at least lsu or ssu
f_euk_cm <- strict_lsu_ssu(f_euk_cm)


#run this only when filtered, or change it to include filtering step
conkingom.tbl <- cm_consensus_kingdom(f_euk_cm)

euk.contig.names <- conkingom.tbl %>%
  dplyr::filter(kingdom=='euk') %>%
  dplyr::select(query_name) %>%
  dplyr::pull()

prok.contig.names <- conkingom.tbl %>%
  dplyr::filter(kingdom=='prok') %>%
  dplyr::select(query_name) %>%
  dplyr::pull()

#part of polishing
f_euk_cm <- remove_5_8s(f_euk_cm, prok.contig.names)



selected_contigs <- f_euk_cm %>%
  dplyr::select(query_name) %>%
  unique() %>%
  pull()


#f_prok_cm <- filter_prok_rDNA(cm = cm_oriented, option = 'length_filter')
sortregions(f_euk_cm) %>% data.frame %>%
  dplyr::filter(query_name %in% euk.contig.names) %>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin = start, xmax = end, y = query_name, fill = target_name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3")


cm %>%  data.frame %>%
  dplyr::filter(query_name %in% c('s1412.ctg001879l')) %>%
  filter(seq_from<200000)%>%
  dplyr::mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin =  seq_from , xmax =  seq_to , y = query_name, fill = target_name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3")


#f_prok_cm <- filter_prok_rDNA(cm = cm_oriented, option = 'length_filter')
sortregions(f_euk_cm) %>% data.frame %>% mutate(direction=ifelse(strand=='+',1,0)) %>%
  ggplot2::ggplot(ggplot2::aes(xmin = start, xmax = end, y = query_name, fill = target_name)) +
  gggenes::geom_gene_arrow() +
  ggplot2::scale_fill_brewer(palette = "Set3")


prox.clust <- cm_proximity_clusters(feature.tbl=cm.rrna, distance = 1000)

prox.clust  %>%
  group_by(prox.cluster, query_name) %>%
  tally() %>%
  arrange(desc(n))



prox.clust %>% filter(prox.cluster == 'cluster286')
prox.clust %>%  filter(prox.cluster == 'cluster286') %>% select(target_name, name,seq_from,seq_to,score ,E_value)

cm_filter(prox.clust, query='eukarya')%>%  filter(prox.cluster == 'cluster69') %>% select(target_name, name,seq_from,seq_to,score ,E_value)
cm_filter(prox.clust, query='bacteria')%>%  filter(prox.cluster == 'cluster69') %>% select(target_name, name,seq_from,seq_to,score ,E_value)




cm_proximity_clusters
cm_proximity_clusters <- function(feature.tbl, distance = 1000){
  #use start and end
  f.tbl <- feature.tbl %>%
    dplyr::mutate(start=seq_from, end=seq_to)

  contigs <- f.tbl %>%
    dplyr::select(query_name) %>%
    dplyr::pull() %>%
    as.character() %>% unique()

  tracker <- 0 #stores the cluster number
  output <- NULL
  for(contig in contigs){
    #subset table
    sub.tbl <-  f.tbl %>% dplyr::filter(query_name == contig )

    glued.tbl <- sub.tbl %>%
      dplyr::mutate(diff = start - dplyr::lag(end,default=dplyr::first(end))) %>%
      dplyr::mutate(glued = ifelse(diff < distance, TRUE, FALSE))# %>%
    #dplyr::select(seqnames, locus_tag,start, end , diff, glued)

    glued.tbl[1,"glued"] <- FALSE

    #logic, for each false, make new group
    clusters <- c()
    for(o in glued.tbl$glued){
      if(!o){
        tracker <- tracker+1
        clusters <- c(clusters,paste0('cluster',tracker))
      }else{
        clusters <- c(clusters,paste0('cluster',tracker))
      }
    }

    glued.tbl <- glued.tbl %>%
      dplyr::mutate(prox.cluster = clusters)

    output <- rbind(output, glued.tbl)

  }
  return(output)
}



#
# bn_proximity_clusters <- function(feature.tbl, distance = 1000){
#   f.tbl <- feature.tbl
#
#   contigs <- f.tbl %>%
#     dplyr::select(contig) %>%
#     dplyr::pull() %>%
#     as.character() %>% unique()
#
#   tracker <- 0 #stores the cluster number
#   output <- NULL
#   for(cg in contigs){
#     sub.tbl <-  f.tbl %>% dplyr::filter(contig == cg )
#
#     glued.tbl <- sub.tbl %>%
#       dplyr::mutate(diff = start - lag(end,default=first(end))) %>%
#       dplyr::mutate(glued = ifelse(diff < distance, TRUE, FALSE))
#
#     glued.tbl[1,"glued"] <- FALSE
#
#     clusters <- c()
#     for(o in glued.tbl$glued){
#       if(!o){
#         tracker <- tracker+1
#         clusters <- c(clusters,paste0('cluster',tracker))
#       }else{
#         clusters <- c(clusters,paste0('cluster',tracker))
#       }
#     }
#
#     glued.tbl <- glued.tbl %>%
#       dplyr::mutate(prox.cluster = clusters)
#     output <- rbind(output, glued.tbl)
#   }
#   return(output)
# }


#------------------------------------------------------#
#   Functions
#-------------------------------------------------------#
cm2granges <- function(cm){
  gr_df <- cm %>%
    #odd, need to look into this more carefully
    dplyr::mutate(start = ifelse(strand=='+', seq_from, seq_to),
                  end = ifelse(strand=='+',seq_to, seq_from)) %>%
    dplyr::mutate(seqnames=query_name) %>%
    dplyr::select(seqnames,target_name, start,end,strand) %>%
    data.frame()

  gr <- GenomicRanges::makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  return(gr)
}


out_fasta_lsussu <- function(cm, seqs, type='SSU',analysisName){
  granges <- cm2granges(cm)
  selected_features <-  granges[grepl(type,granges$target_name),]
  feat_seq <- BSgenome::getSeq(seqs, selected_features)
  names(feat_seq) = paste0(
    as.character(data.frame(selected_features[grepl(type,selected_features$target_name),])$seqnames),
    data.frame(selected_features[grepl(type,selected_features$target_name),])$start,
    data.frame(selected_features[grepl(type,selected_features$target_name),])$target_name#,
  )
  Biostrings::writeXStringSet(feat_seq, paste0("/Users/sa01fd/DATA/rna-miner/GRanges_",type,"_",analysisName,".fasta"),append=FALSE)
}


remove_5_8s <- function(cm, prok.contig.names){
  out <- cm %>%
    mutate(prok = query_name %in% prok.contig.names,
           remove = ifelse(target_name=="5_8S_rRNA" & prok==TRUE, TRUE, FALSE)) %>%
    filter(remove==FALSE)
  return(out)

}


#sortregions
#=================================================================================#
#' Function that sets start to 0, and use majority rule direction (most on + strand, reverse visual)
#'
#' @param ann.tbl
#' @param majority TRUE when you want to use majority rule
#' @param OGs specify OG strand to define direction
#'
#' @return
#' @export
#'
#' @examples
#' otpt <- sortregions(ann.tbl=otpt)
sortregions <- function(ann.tbl,majority=TRUE,OGs=NULL){
  otpt <- NULL
  ann.tbl <- ann.tbl %>% dplyr::mutate(start=seq_from, end=seq_to)

  for(contig in unique(as.character(ann.tbl$query_name))){
    sb <- ann.tbl %>%
      dplyr::filter(query_name == contig)

    if(!is.null(OGs)){
      dirctn <- sb %>%
        dplyr::filter(target_name == OGs) %>%
        dplyr::select(strand) %>%
        dplyr::pull() %>%
        dplyr::head(n=1) %>%
        dplyr::as.character()

      dirctn <- ifelse(dirctn=='+',1,-1)
    }
    if(majority==TRUE){
      dirctn <- sb %>%
        dplyr::mutate(direction=ifelse(strand=="+",1,-1)) %>%
        dplyr::select(direction) %>%
        sum()
    }

    if(dirctn <0){
      maxend <- sb %>%
        dplyr::select(end) %>%
        max()

      sb <- sb %>%
        dplyr::mutate(start=maxend-start) %>%
        dplyr::mutate(end=maxend-end)
    }

    #General, all trimm to 0 start
    minstart <- sb %>%
      dplyr::select(start) %>% min()

    sb <- sb %>%
      dplyr::mutate(start=start-minstart) %>%
      dplyr::mutate(end=end-minstart)

    otpt <- rbind(otpt,sb)
  }
  return(otpt)
}


path<-sina.path
read.sina <- function(path){
  out <- read.delim(path,sep=',')
  out$name <- as.character(out$name)
  return(out)
}



cm_sort <- function(cm){
  out <- cm %>% group_by(query_name) %>% dplyr::arrange(seq_from) %>% ungroup()
  return(out)
}



hardFilter <- function(cm){
  # NEED TO TAKE INTO ACCOUNT OPERON, MULTIPLE OPERONS IN ONE CONTIG WONT YET WORK
  # take highest scoring SSU and LSU per contig
  cm <- cm_proximity_clusters(cm,distance=5000)

  filtered_lsu_ssu <- cm %>%
    dplyr::filter(grepl("LSU|SSU",target_name)) %>%
    dplyr::mutate(type=ifelse(grepl("LSU",target_name),'LSU','SSU')) %>%
    dplyr::group_by(type, query_name,prox.cluster) %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-type)

  non_lsu_ssu <- cm %>% filter(!grepl("LSU|SSU",target_name))

  out <- rbind(filtered_lsu_ssu, non_lsu_ssu)

  out <- out %>% dplyr::mutate(length=abs(end-start)) %>%
    dplyr::filter(target_name == "SSU_rRNA_bacteria" & length > 1200 |
                    target_name == "LSU_rRNA_bacteria" & length > 2100 |
                    target_name == "5S_rRNA" & length > 100 |
                    target_name == "tRNA" |
                    target_name == "SSU_rRNA_eukarya" & length > 1400 |
                    target_name == "LSU_rRNA_eukarya" & length > 2800 |
                    target_name == "5_8S_rRNA" & length > 150 |
                    target_name == "mtPerm_5S")

  return(out)
}


strict_lsu_ssu <- function(cm){
  contigs <- cm %>%
    dplyr::filter(grepl("LSU|SSU",target_name)) %>%
    dplyr::select(query_name) %>%
    dplyr::pull()
  out <- cm %>%
    dplyr::filter(query_name %in% contigs)
  return(out)
}



cm_orient <- function(cm, seqs.tbl){
  #determine majority strand orientation
  strand_major.count <- cm %>%
    dplyr::group_by(query_name, strand) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::group_by(query_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  n_r <- strand_major.count %>% filter(strand=='-') %>% tally() %>% pull()
  message(n_r,' contigs detected, and reversed')
  #change orientation in CM table
  out <- cm %>%
    dplyr::left_join(seqs.tbl %>%
                       dplyr::select(name, length), by=c('query_name'='name')) %>%
    dplyr::left_join(strand_major.count %>%
                       dplyr::mutate(majority_strand = strand) %>%
                       dplyr::select(query_name,majority_strand),by=c('query_name'='query_name'))#%>%
  #dplyr::mutate(start = ifelse(majority_strand=="-",length-seq_from,seq_from),
  #              end = ifelse(majority_strand=="-",length-seq_to,seq_to)) %>% select(-length)

  return(out)
}



seq_cm_revcomp <- function(cm, seqs.tbl){
  #determine majority strand orientation
  strand_major.count <- cm %>%
    dplyr::group_by(query_name, strand) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::group_by(query_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  sq.tb <- seqs.tbl %>% left_join(strand_major.count, by=c('name'='query_name')) %>% data.frame()

  cnt <- 0
  #sl <- c()
  for(i in 1:nrow(sq.tb)){
    if(!is.na(sq.tb[i,'strand']=='-')){
      if(sq.tb[i,'strand']=='-'){
        message('reverse complement of ',sq.tb[i,'name'])
        sq.tb[i,'sequence'] = as.character(Biostrings::reverseComplement(Biostrings::DNAString(sq.tb[1,'sequence'])))
        cnt <- cnt+1
      }
    }
  }
  if(cnt>1){
    message('Changed orientation of ', cnt, 'sequences')
  }else{
    message('All sequences were kept in original orientation')
  }
  sq.tb <- sq.tb %>% tibble::tibble()
  return(sq.tb)
}



gff_parse_seq <- function(euk_cm, seqs.tbl,start='seq_from',end='seq_to'){
  out.t <- c()
  for(qr_name in unique(euk_cm$query_name)){
    selected_sequence <- seqs.tbl %>%
      dplyr::filter(name==qr_name) %>%
      dplyr::select(sequence) %>%
      dplyr::pull()

    s <- euk_cm %>%
      dplyr::filter(query_name == qr_name) %>%
      data.frame()

    seq_tr <- c()
    for(i in 1:nrow(s)){
      seq_tr <- c(seq_tr, substr(selected_sequence, s[i,start], s[i,end]))
    }

    out <- cbind(s,sequence=seq_tr) %>%
      dplyr::mutate(lenght_c=nchar(seq_tr))

    out.t <- rbind(out.t, out)

  }
  out.t <- out.t %>%
    tibble::tibble() %>%
    dplyr::mutate(sequence=as.character(sequence))
  return(out.t)
}


#example conkingom.tbl <- cm_consensus_kingdom(f_euk_cm)
cm_consensus_kingdom <- function(cm){
  #based on SSU and LSU sequences
  out <- cm %>%
    dplyr::filter(grepl("SSU|LSU",target_name)) %>%
    dplyr::mutate(kingdom = ifelse(grepl('eukar',target_name),'euk',
                                   ifelse(grepl('bact',target_name),'prok','unkown'))) %>%
    dplyr::select(target_name, query_name, kingdom) %>%
    dplyr::group_by(query_name,kingdom) %>%
    dplyr::tally() %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n))
  print(out %>% group_by(kingdom) %>% tally())
  return(out)
}



#---------------------------------#
detect_ITS <- function(cm, t.name='ITS1', markers=c('SSU_rRNA_eukarya','5_8S_rRNA')){
  out.t <- c()
  for(qr_name in unique(cm$query_name)){
    message('::',qr_name)

    s <- cm %>%
      dplyr::filter(query_name == qr_name) %>%
      data.frame()

    if(markers[1] %in% s$target_name & markers[2] %in% s$target_name ){
      if(all(table(s$target_name) < 2 )){
        message('single targets detected')
      }else{
        message('multiple targets detected')
        print(table(s$target_name))
      }

      if (markers[1] %in% s$target_name & markers[2] %in% s$target_name) {
        message(markers[1], 'and', markers[2] ,' detected')
      }else{
        message("zero seqs passed rDNA opern search. check parameters")
      }

      r5_8S_end <- s[s$target_name==markers[1],'seq_to']
      lsu_start <- s[s$target_name==markers[2],'seq_from']

      if(lsu_start > r5_8S_end){
        out.its <-
          rbind(c(target_name = t.name,
                  query_name=qr_name,
                  seq_from=r5_8S_end+1,
                  seq_to=lsu_start-1,
                  strand='+',
                  length=lsu_start-r5_8S_end-2)) %>%
          data.frame() %>%
          dplyr::mutate(seq_from=as.numeric(as.character(seq_from)),
                        seq_to=as.numeric(as.character(seq_to)),
                        length=as.numeric(as.character(length)))
      }


      out.t <- rbind(out.t, out.its)
    }else{
      message('markers not detected in ',qr_name)
    }
  }

  return(out.t)
}



#-------------------------------------------------
detect_stem <- function(cm, stem_length , target_name = 'ITS2', markers = c('5_8S_rRNA','LSU_rRNA_eukarya')){
  out.t <- c()
  for(qr_name in unique(cm$query_name)){
    message('::',qr_name)

    s <- cm %>%
      dplyr::filter(query_name == qr_name) %>%
      data.frame()

    if(all(table(s$target_name) < 2 )){
      message('single targets detected')
    }else{
      message('multiple targets detected')
      print(table(s$target_name))
    }

    if (markers[1] %in% s$target_name & markers[2] %in% s$target_name) {
      message(markers[1], 'and',markers[] ,' detected')
    }else{
      message("zero seqs passed rDNA opern search. check parameters")
    }

    r5_8S_end <- s[s$target_name==markers[1],'seq_to']
    lsu_start <- s[s$target_name==markers[2],'seq_from']

    if(lsu_start > r5_8S_end){
      out.its <-
        rbind(c(target_name = paste0('pre_',target_name),
                query_name=qr_name,
                seq_from=r5_8S_end+1-25,
                seq_to=r5_8S_end,
                strand='+',
                length=25),
              c(target_name = paste0('post_',target_name),
                query_name=qr_name,
                seq_from=lsu_start+1,
                seq_to=lsu_start+1+25,
                strand='+',
                length=25)) %>%
        data.frame() %>%
        dplyr::mutate(seq_from=as.numeric(as.character(seq_from)),
                      seq_to=as.numeric(as.character(seq_to)),
                      length=as.numeric(as.character(length)))
    }

    out.t <- rbind(out.t, out.its)
  }

  return(out.t)
}
