#PR2

library(dplyr)

install.packages(devtools)
remotes::install_github("pr2database/pr2database")

#--------------------------------------------------#

library("pr2database")

pr2 <- pr2_database()

pr2_ostreo <- pr2 %>% dplyr::filter(genus == "Pycnococcus")
pr2_ostreo <- pr2_ostreo %>% dplyr::select( genbank_accession, species,
                                            pr2_sample_type, gb_strain, gb_clone,
                                            pr2_latitude, pr2_longitude,
                                            sequence_length, sequence, reference_sequence  )

pr2 %>% dplyr::select(kingdom,supergroup,division,class, order) %>% unique()%>% arrange(division) %>% data.frame()

pr2 %>% dplyr::filter(division == "Chlorophyta") %>%
  dplyr::select(kingdom,supergroup,division,class, order, family,genus,species) %>%
  unique()

pr2_chlor <- pr2 %>% dplyr::filter(division == "Chlorophyta")

pr2_chlor <- pr2_chlor %>% dplyr::filter(sequence_length >1399) %>%
  dplyr::filter(pr2_sample_type %in% c('culture','isolate'))


#pr2_chlor <- pr2_chlor %>% dplyr::select( genbank_accession, species,
#                                            pr2_sample_type, gb_strain, gb_clone,
#                                            pr2_latitude, pr2_longitude,
#                                            sequence_length, sequence, reference_sequence  )


seq_chlor <- Biostrings::DNAStringSet(pr2_chlor$sequence)

names(seq_chlor) <- paste(pr2_chlor$genbank_accession, pr2_chlor$species,
                           "strain",pr2_chlor$gb_strain,
                           "clone",pr2_chlor$gb_clone,
                           sep="|")

Biostrings::writeXStringSet(seq_chlor, "~/DATA/pr2_chlor.fasta", width = 80)

#--------------------#

pr2_chlor_genus <- pr2_chlor %>% dplyr::group_by(genus) %>% dplyr::slice(5)
seq_chlor <- Biostrings::DNAStringSet(pr2_chlor_genus$sequence)

names(seq_chlor) <- paste(pr2_chlor_genus$genbank_accession, pr2_chlor_genus$species,
                          "strain",pr2_chlor_genus$gb_strain,
                          "clone",pr2_chlor_genus$gb_clone,
                          sep="|")

names(seq_chlor) <- paste(pr2_chlor_genus$genbank_accession, pr2_chlor_genus$species,
                          pr2_chlor_genus$gb_strain,pr2_chlor_genus$order,pr2_chlor_genus$family,
                          sep="|")

Biostrings::writeXStringSet(seq_chlor, "~/DATA/pr2_chlor_genus.fasta", width = 80)


