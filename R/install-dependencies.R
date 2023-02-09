#install-dependencies

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

install.packages('readr',dependencies=TRUE)
