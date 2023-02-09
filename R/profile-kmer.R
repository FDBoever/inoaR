suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(kmer))
suppressMessages(library(rgr))



k <-args[1]
k=5



#=======================================
# Calculate k-mer frequencies
#======================================
message(paste0(':: calculating ',k,'-mer frequencies'))

# K-mer frequency using kmer package
# Prepare sequences as input for kmer::kcount function (requires a bit of an odd input format)
# anyway, this is how you can do it
seq <- contigs.tbl$sequence # nothing else then a vector of sequences
contigs.list <- stringr::str_extract_all(seq, stringr::boundary("character"))

#k- is set as input varameter now
kmer.matrix <- kmer::kcount(contigs.list,k = k)
rownames(kmer.matrix) <- contigs.tbl$name


#=======================================#
# Centered log-ratio (CLR) transformation
#=======================================#

# given compositional nature of k-mer frequencies, people often use CLR transformation

#---- calculate clr manually?
# interesting souce for compositional data transformations --> https://github.com/ggloor/book/blob/master/05-transforms.Rmd
#kmer.dat.prop <- t(apply(kmer.matrix+1, 1, function(x) x/sum(x)))
#kmer.dat.clr <- t(apply(kmer.dat.prop, 1, function(x) log(x) - mean(log(x))))

#----- calculate clr using rgr package
# I add +1 to the matrix to avoid log(0) for zero counts... this is needed to not get Inf values
message(paste0(':: Centered log-ratio (CLR) transformation'))
kmer.clr <- rgr::clr(t(kmer.matrix+1),ifwarn = FALSE) %>% t()

# filter out kmers with zero-variance
kmer.clr.filtered <- kmer.clr[ , which(apply(kmer.clr, 2, var) != 0)]


#==================================================
# Jaccard Distance
#kmder.dist <- dist(kmer.matrix, method = "binary")

# NMDS/PCoA
#df.cmdscale <- cmdscale(kmder.dist) %>% data.frame
#colnames(df.cmdscale) <- c('PCOA1', 'PCOA2')
#df.cmdscale %>% ggplot2::ggplot(aes(PCOA1,PCOA2))+geom_point()



