
library(Rsamtools)

#pileupsD <- list.files("~/genomedk/DELFI_data/RawData", recursive = T, full.names = T) # only BAM
pileupsD <- list.files("~/genomedk/DELFI_data/RawData", recursive = T, full.names = T, pattern = "bam$")

# extract the list of 27 CRC and 43 CTL from 
#sL <- read.csv("~/genomedk/matovanalysis/DELFI_analysis/R/sample_reference.csv")
#iC <- sL[sL[,6]=="Colorectal Cancer",3] # list of CRC patient numbers/names,
#iH <- sL[sL[,6]=="Healthy",3] # list of healthy/control samples, 215

# match the list from the pileup file and load the 70 BAM files

# test, load a BAM file and extract FRL364bp, also 205 and 198bp.
which <- GRanges(seqnames = c("seq1"),ranges = IRanges(c(364),c(364)))
what <- c("rname", "strand", "pos", "qwidth","seq")
param <- ScanBamParam(which = which, what = what)
bam <- scanBam(pileupsD[1], param=param)

# test load BAM
#bf <- BamFile(pileupsD[1])
#seqnames(seqinfo(bf))
