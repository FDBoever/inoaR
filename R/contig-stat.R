# loading packages
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))


#' read_assembly
#' Function to load assembly and convert it to a contigs.tbl, including length, GC, percent A, T, G, C as well as the sequence
#' @param filePath path to the assembly fasta file
#'
#' @return tibble with contig stats
#' @export
#'
#' @examples
#' read_assembly(filePath = '~/DATA/EukGenome_database/Labarre2021/12430790/MAST-11-sp1.fasta')
read_assembly <- function(filePath){
  contigs <- Biostrings::readDNAStringSet(filePath)

  contigs.tbl <- tibble::tibble(name=names(contigs),
                                length = nchar(contigs),
                                GC = c(Biostrings::letterFrequency(contigs, letters = "GC", as.prob = TRUE)),
                                pA = c(Biostrings::letterFrequency(contigs, letters = "A", as.prob = TRUE)),
                                pT = c(Biostrings::letterFrequency(contigs, letters = "T", as.prob = TRUE)),
                                pG = c(Biostrings::letterFrequency(contigs, letters = "G", as.prob = TRUE)),
                                pC = c(Biostrings::letterFrequency(contigs, letters = "C", as.prob = TRUE)),
                                sequence=paste(contigs))
  return(contigs.tbl)
}

#' Summarise Assembly
#' Function to summarise assemblies (multi-fasta)
#' contining stats like, total lenght, n contigs, n contigs longer than given lenght, N50, L50 etc
#' @param contigs.tbl contigs.tbl generated via read_assembly function
#' @param name name of the assembly
#'
#' @return tibble with summary statistics
#' @export
#'
#' @examples
#' summarise_assembly(contigs.tbl = cntgs, name= 'MAST-11-sp1')
#'
summarise_assembly <- function(contigs.tbl, name){
  c_ls <- contigs.tbl$length
  c_ls <- rev(sort(c_ls))
  sum_ls <- sum(c_ls)
  cumsum_ls <- cumsum(c_ls)
  contig.summary <- tibble::tibble(
    name = name,
    total_length = sum(c_ls),
    n_contigs = length(c_ls),
    n_2p5kb = length(c_ls[c_ls>2500]),
    n_5kb = length(c_ls[c_ls>5000]),
    n_10kb = length(c_ls[c_ls>10000]),
    n_20kb =length(c_ls[c_ls>20000]),
    n_50kb = length(c_ls[c_ls>50000]),
    n_100kb = length(c_ls[c_ls>100000]),
    N25 = c_ls[which(cumsum_ls - sum_ls* .25 >=0)[1]],
    N50 = c_ls[which(cumsum_ls - sum_ls* .5 >=0)[1]],
    N75 = c_ls[which(cumsum_ls - sum_ls* .75 >=0)[1]],
    N90 = c_ls[which(cumsum_ls - sum_ls* .9 >=0)[1]],
    L25 = which(cumsum_ls - sum_ls* .25 >=0)[1],
    L50 = which(cumsum_ls - sum_ls* .5 >=0)[1],
    L75 = which(cumsum_ls - sum_ls* .75 >=0)[1],
    L90 = which(cumsum_ls - sum_ls* .9 >=0)[1],
    median_lenght = median(c_ls),
    mean_length = mean(c_ls),
    shortest = min(c_ls),
    longest = max(c_ls)
  )
  return(contig.summary)
}


#' Summarise assemblies
#' Function that summarises multiple assemblies
#' @param dirPath path to the directory containing fasta files
#' @param fileExt file extention of the fasta files (fasta, faa, fn)
#'
#' @return tibble with summarised assemblies (rows)
#' @export
#'
#' @examples
#' assemblies.summary <- summarise_assemlies(dirPath = '~/DATA/EukGenome_database/Labarre2021/12430790/', fileExt = 'fasta')

summarise_assemlies <- function(dirPath, fileExt = 'fasta'){
  files <- list.files(full.names = TRUE, path = dirPath , pattern = paste0('*',fileExt))
  assemblies.stat <- c()
  for (f in files){
    name <- gsub(paste0('.',fileExt),'',basename(f))
    cntgs <- read_assembly(f)
    assemblies.stat <- rbind(assemblies.stat,
                           summarise_assembly(cntgs, name))
  }
  return(assemblies.stat)
}




assemblies.summary <- summarise_assemlies(dirPath = '~/DATA/EukGenome_database/Labarre2021/12430790/',
                    fileExt = 'fasta')

