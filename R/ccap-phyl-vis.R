#Libraries
if(!require(stringr)) install.packages("stringr", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ape)) install.packages("ape", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(phytools)) install.packages("phytools", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ggrepel)) install.packages("ggrepel", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org",dependencies=TRUE)
#if(!require(devtools)) install.packages("devtools",dependencies=TRUE)
#install.packages("devtools") # if you have not installed "devtools" package
if(!require(remotes)) install.packages("remotes", repos = "http://cran.us.r-project.org",dependencies=TRUE)
if(!require(ggtree)) remotes::install_github("YuLab-SMU/ggtree")

#library(stringr)
#library(ape)
#library(phytools)
#library(ggrepel)
#library(ggtree)
#library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
analysisName = args[1]
#print(args[2])
#print(args[3])

# In bash use this
#	find . -type f -name "*.tree" | xargs -I{} cat "{}" > trees.tree


#
#Acanthamoeba


########################################################################################
#~/DATA/phyloCCAP/CCAP_treboux_Darienko/CCAP_treboux_Darienko.18S.quality.txt

#analysisName = 'Darienko_2019'
#analysisName = 'CCAP_treboux_Darienko'
analysisName = 'Acanthamoeba'
analysisName = 'Chaetoceros'
analysisName = 'Isochrysis'
analysisName = 'Botryococcus'

analysisName = 'Porphyridium'




dir_path <- '~/DATA/phyloCCAP/allGenera/'
dir_path <- '~/DATA/phyloCCAP/'

tree_path = paste(c(dir_path , analysisName, '/trees/', analysisName,'.18S.MAFFT.trimal.aln_FastTree.tree'),collapse='')



# load the tree
tree = ape::read.tree(file= tree_path)

# PARSE THE OUTGROUP IF KNOWN AND USE HERE
# OR IF NOT KNOWN, USE MIDPOINT ROOT
mdTree = ape::ladderize(phytools::midpoint.root(tree),right=TRUE)

#--Load metadata------------------
metadata_path = paste(c(dir_path, analysisName,'/', analysisName,'.insd.metadata.txt'),collapse='')

NCBImetadata = read.delim(metadata_path,sep='\t',header = FALSE)
NCBImetadata  = NCBImetadata %>% unique()
rownames(NCBImetadata) = NCBImetadata$V10
#rownames(NCBImetadata) = gsub("\\.1","",rownames(NCBImetadata))

mdTree$tip.label = gsub("\\.1","",mdTree$tip.label)
keepRecord = mdTree$tip.label
mdTree$tip.label = paste(NCBImetadata[mdTree$tip.label,'V12'], NCBImetadata[mdTree$tip.label,'V11'],sep=' - ')


#-----------visualise tree---------------------
ccap_or_not = grepl(" CCAP", mdTree $tip.label)
ccap_color = c('TRUE'='red','FALSE'='black')
ccap_color[as.character(ccap_or_not)]


tipsize = 2
cols <- c("100" = "black", "95-99" = "grey20", "85-95" = "grey40", "75-85" = "grey", "<75" = "lightgrey")
p <- ggtree::ggtree(mdTree)
p$data$bootstrap <- as.numeric(p$data$label)
p$data$Source = grepl(" CCAP",p$data$label)
p <- p + geom_point2(aes(subset=!is.na(bootstrap) & as.numeric(bootstrap) == 1,label = "100", fill="100") ,pch=21) +
  geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) >0.95 & bootstrap < 1, label="95-99", fill="95-99"), pch=21) +
  geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) >0.85 & bootstrap < 0.95, label="85-95", fill="85-95"), pch=21) +
  geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) > 0.75 & bootstrap < 0.85, label="75-85", fill="75-85") ,pch=21) +
  geom_point2(aes(subset=!is.na(bootstrap) & as.numeric( bootstrap) <0.75 , label="<75", fill="<75") ,pch=21) +
  geom_treescale(x=1, y=10) +
  geom_tiplab(aes(color=Source),size=2 ) +
  xlim(0,0.5) +
  scale_fill_manual(name = "bootstrap", breaks=c("100", "95-99", "85-95", "75-85", "<75"), values = cols, labels = c("100", "95-99", "85-95", "75-85", "<75"))+
  theme(legend.position="right")+
  scale_color_manual(values=c('black','red'))

p


ggsave(filename = paste0('~/DATA/phyloCCAP/tree_',analysisName,'_testingCCAP-phyl-vis.pdf') , plot = p , device='pdf' , height = 25, width = 10)

