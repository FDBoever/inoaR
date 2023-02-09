# ccap-annotate-fasta.R
#----------------------------------------------------------#
# Script that takes ccap-db-tools fasta files with accession numbers as headers and enriches the headers with CCAP strain numbers and species names

#dependencies
# library(dplyr)
# library(seqinr)




# Functions
#----------------------------------------------------------#
#' ccap_annotate_fasta
#' save annotated output file as .named. that includes ccap strain number and species names
#' @param inPath path to the fasta file
#' @param metadataPath path to the metadata file (could be linked internally via sql)
#' @param fastaExtension extension of fasta file, could be '.fasta' (default), '.fa', '.fna', '.faa' etc..
#'
#' @return doesnt return anything, it saves an update of the file
#' @export
#'
#' @examples
#' ccap_annotate_fasta(inpath='~/DATA/phyloCCAP/alignments/rRNA_euk_all_recent_update.5_8S.MAFFT.aln.fasta',
#'                     metadataPath='~/DATA/phyloCCAP/tblSequences.txt',
#'                     fastaExtension = 'fasta')
#'
ccap_annotate_fasta <- function(inPath, metadataPath, fastaExtension='fasta'){
  metadata = read.delim(metadataPath,sep='\t',header = TRUE) %>% tibble::as_tibble()
  inFasta <- seqinr::read.fasta(inPath)

  outPath <- paste0(dirname(inPath),'/',gsub(pattern = paste0("\\.",fastaExtension,"$"), "",  basename(inPath)),'.named.',fastaExtension)

  headers <- names(inFasta)
  out <- metadata %>%
    dplyr::filter(accession_version %in% headers) %>%
    dplyr::select(consensus, NCBIaccession , accession_version, organism) %>%
    dplyr::mutate(combined = paste(NCBIaccession, consensus, organism, sep='_'))

  annotated_names <- gsub(' ', '_',out$combined)
  names(annotated_names) <- out$accession_version

  seqinr::write.fasta(sequences = inFasta,
                      names = unname(annotated_names[names(inFasta)]),
                      file.out = outPath)
}



# Example Usage
#----------------------------------------------------------#

ccap_annotate_fasta(inPath='~/DATA/phyloCCAP/reference_trees/eukaryotes/alignments/rRNA_euk_all_recent_update.18S.MAFFT.aln.fasta',
                     metadataPath='~/DATA/phyloCCAP/tblSequences.txt',
                     fastaExtension = 'fasta')

ccap_annotate_fasta(inPath='~/DATA/phyloCCAP/alignments/rRNA_euk_all_recent_update.5_8S.MAFFT.aln.fasta',
                    metadataPath='~/DATA/phyloCCAP/tblSequences.txt',
                    fastaExtension = 'fasta')

ccap_annotate_fasta(inPath='~/DATA/phyloCCAP/rRNA_euk_all_recent_update.28S.MAFFT.aln.fasta',
                    metadataPath='~/DATA/phyloCCAP/tblSequences.txt',
                    fastaExtension = 'fasta')

ccap_annotate_fasta(inPath='~/DATA/phyloCCAP/alignments/rRNA_euk_all_recent_update.5_8S.MAFFT.structure.aln.fasta',
                    metadataPath='~/DATA/phyloCCAP/tblSequences.txt',
                    fastaExtension = 'fasta')

