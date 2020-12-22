
library(Rsamtools)

#pileupsD <- list.files("~/genomedk/DELFI_data/RawData", recursive = T, full.names = T) # only BAM
#[204] "/home/alexandre/genomedk/DELFI_data/RawData/EGAF00002727334/PGDX5881P_WGS_processed_downsamp.bam"       
pileupsD <- list.files("~/genomedk/DELFI_data/RawData", recursive = T, full.names = T, pattern = "bam$")

# extract the list of 27 CRC and 43 CTL from 
#sL <- read.csv("~/genomedk/matovanalysis/DELFI_analysis/R/sample_reference.csv")
#iC <- sL[sL[,6]=="Colorectal Cancer",3] # list of CRC patient numbers/names,
#iH <- sL[sL[,6]=="Healthy",3] # list of healthy/control samples, 215

# match the list from the pileup file and load the 70 BAM files

# test, load a BAM file and extract FRL364bp, also 205 and 198bp.
which <- GRanges(seqnames = c("chr17"),ranges = IRanges(c(65000000),c(70000000)))
what <- c("rname", "strand", "pos", "qwidth","seq")
param <- ScanBamParam(which = which, what = what)
param <- ScanBamParam(which=which, what=scanBamWhat())
bam <- scanBam(pileupsD[204], param=param)

frl <- bam$`chr17:65000000-70000000`$isize # insert size, how is it related to fragment length?
which(frl==364) # 16
#     66   1364   9295  11447  21434  24023  40216  58545  62288  62766  64194  66174  95151  96881 108029 110268

bam$`chr17:65000000-70000000`$mpos # starting positions; how to find the ending positions?

# test load BAM
#bf <- BamFile(pileupsD[1])
#seqnames(seqinfo(bf))
