if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sangerseqR")
library(sangerseqR, quietly=TRUE)
library(Biostrings, quietly=TRUE)




hetab1 <- read.abif("~/DATA/Sanger_sequencing/FDBoever-17-10-2022-1/9E-G07-191022-05-18.ab1")
str(hetab1, list.len=20)


homosangerseq <- sangerseq(hetab1)

## --------------------------------------------------------------------------
#from a sequence file object
homosangerseq <- sangerseq(homoscf)

#directly from the file
hetsangerseq <- readsangerseq(system.file("extdata",
                                          "heterozygous.ab1",
                                          package="sangerseqR"))
str(hetsangerseq)

## --------------------------------------------------------------------------
#default is to return a DNAString object
Seq1 <- primarySeq(homosangerseq)
reverseComplement(Seq1)

#can return as string
primarySeq(homosangerseq, string=TRUE)

## --------------------------------------------------------------------------
chromatogram(hetsangerseq, width=200, height=2, trim5=50, trim3=100,
             showcalls='both', filename="chromatogram.pdf")

## --------------------------------------------------------------------------
hetcalls <- makeBaseCalls(hetsangerseq, ratio=0.33)
hetcalls

## --------------------------------------------------------------------------
chromatogram(hetcalls, width=100, height=2, trim5=50, trim3=100,
             showcalls='both', filename="chromatogram2.pdf")

## --------------------------------------------------------------------------
ref <- subseq(primarySeq(homosangerseq, string=TRUE), start=30, width=500)
hetseqalleles <- setAllelePhase(hetcalls, ref, trim5=50, trim3=300)
hetseqalleles

## --------------------------------------------------------------------------
pa <- pairwiseAlignment(primarySeq(hetseqalleles)[1:400],
                        secondarySeq(hetseqalleles)[1:400],
                        type="global-local")
writePairwiseAlignments(pa)





library(sangeranalyseR)
A_chloroticaFFN <- file.path("~/DATA/Sanger_sequencing/untitled folder/9E-191022-05-F.ab1",
                             "Allolobophora_chlorotica",
                             "ACHLO",
                             "9E-191022-05-F.ab1")
sangerReadF <- sangeranalyseR::SangerRead(readFeature = "Forward Read",
                          readFileName          = A_chloroticaFFN,
                          geneticCode           = GENETIC_CODE,
                          TrimmingMethod        = "M1",
                          M1TrimmingCutoff      = 0.0001,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE)

# using `new` method to create SangerRead instance
sangerReadF <- new("SangerRead",
                   readFeature           = "Forward Read",
                   readFileName          = A_chloroticaFFN,
                   geneticCode           = GENETIC_CODE,
                   TrimmingMethod        = "M1",
                   M1TrimmingCutoff      = 0.0001,
                   M2CutoffQualityScore  = NULL,
                   M2SlidingWindowSize   = NULL,
                   baseNumPerRow         = 100,
                   heightPerRow          = 200,
                   signalRatioCutoff     = 0.33,
                   showTrimmed           = TRUE)
