library(stringr)
library(devtools)
library (LaplacesDemon)
library(FNN)
library("RColorBrewer")
library(ROCR)
library("ggplot2")
library(gplots)
library(spatstat)
library("transport")
library(readxl)
setwd ('G:\\DELFI_data/Derived/fragment_length_in_bins')
setwd ('~/genomedk/DELFI_data/Derived/fragment_length_in_bins')
require("reticulate")
py_install("pandas")
py_install("https://github.com/Hogfeldt/ctDNAtool")
source_python("~/genomedk/matovanalysis/DELFI_analysis/python/pickle_reader.py")
#pickle_data <- read_pickle_file("~/genomedk/DELFI2/Workspaces/per_and_elias/delfi2_length_5Mbp/DL000978HLQ0_100AM.pickle")
##excel_sheets(path = "~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx")

#del2 <- read_excel("~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx", sheet = 1) 
#end2 <- read_excel("~/genomedk/DELFI2/RawData/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo.xlsx", sheet = 2)# 1069 samples: 169 colon, 100 rectum, 800 control.
del2 <- read_excel("~/genomedk/DELFI2/Workspaces/matov/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo2.xlsx", sheet = 1) 
end2 <- read_excel("~/genomedk/DELFI2/Workspaces/matov/201217_Delfi2_fastq_and_sample_manifest_updated_batchinfo2.xlsx", sheet = 2)# 1069 samples: 169 colon, 100 rectum, 800 control.


col_list <- which(end2$diagnostic_group=="Colon cancer") # 169 samples
col1 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="I") # 19
col2 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="II") # 65
col3 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="III") # 42
col4 <- which(end2$diagnostic_group=="Colon cancer" & end2$crc_stage=="IV") # 43
col0 <- which(end2$diagnostic_group=="Adenoma colon" & end2$adenoma_risk=="HIGH") # 46 col Adenoma HIGH
colA <- which(end2$diagnostic_group=="Adenoma colon" & end2$adenoma_risk=="LOW") # 53 col Adenoma LOW
ctl1 <- which(end2$diagnostic_group=="No comorbidity-no finding") # 284
ctl2 <- which(end2$diagnostic_group=="Comorbidity-no finding") # 191
ctl3 <- which(end2$diagnostic_group=="Other finding") # 325
rec_list <- which(end2$diagnostic_group=="Rectal cancer") # 100 samples
rec1 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="I") # 29
rec2 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="II") # 25
rec3 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="III") # 31
rec4 <- which(end2$diagnostic_group=="Rectal cancer" & end2$crc_stage=="IV") # 15
rec0 <- which(end2$diagnostic_group=="Adenoma rectum" & end2$adenoma_risk=="HIGH") # 21 rec Adenoma HIGH
recA <- which(end2$diagnostic_group=="Adenoma rectum" & end2$adenoma_risk=="LOW") # 14 rec Adenoma LOW

e2 <- data.frame(end2) # 1803   27 - stage and disease info comes from this file (sheet 2)
d2 <- data.frame(del2) # 765  10  - DELFI names (sheet 1)
m2 <- merge(e2,d2,by="SampleID") # 681  36 - only 681 are members of both groups


#bFr <- read.delim("filtered_window_data_5MB_1_700_m5000.txt")
#bf <- read.table("filtered_window_data_5MB_1_700_m5000.txt")
bFr <- read.table("filtered_window_data_5MB_1_700_m5000.txt", header = TRUE)
typeof(bFr) # list
lookUp <- read.table("5MB_map_file.txt")

# originial paper has 504 bins per patient sample, which covers 2.52G of genome
s1 = unlist(bFr[1,3:702]) # 700 (fragments) length histogram values for one 5mps bin
# up to 49 bins per chromosome, some get very few
s2 = unlist(bFr[2,3:702])
hist(s1[s1>100])
hist(s2[s2>100])

#sL <- read.csv("U:\\Documents/R/delfi_scripts-master/sample_reference.csv")
sL <- read.csv("~/genomedk/matovanalysis/DELFI_analysis/R/sample_reference.csv")

iL <- sL[sL[,6]=="Lung Cancer",3] # list  
iD <- sL[sL[,6]=="Cholangiocarcinoma",3] # list  
iB <- sL[sL[,6]=="Breast Cancer",3] # list  
iG <- sL[sL[,6]=="Gastric cancer",3] # list  
iC <- sL[sL[,6]=="Colorectal Cancer",3] # list of CRC patient numbers/names,
iC1 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="I",3] # list of CRC patient numbers/names,
iC2 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="II",3] # list of CRC patient numbers/names,
iC3 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="III",3] # list of CRC patient numbers/names,
iC4 <- sL[sL[,6]=="Colorectal Cancer" & sL[,7]=="IV",3] # list of CRC patient numbers/names,
iO <- sL[sL[,6]=="Ovarian Cancer",3] # list of ovarian cancer patient numbers/names,
iP <- sL[sL[,6]=="Pancreatic Cancer",3] # list  
iH <- sL[sL[,6]=="Healthy",3] # list of healthy/control samples, 215

length(bFr[bFr[,1]==iC[1],1])# 555 bins, i.e number of rows with CRC FRs for sample 1

# number of fragments per patient
sum(crc1[,3:702])#54026559

nbLCC <- length(iL)
lcc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbLCC) {
  #i= 1
  lcc1 <- bFr[bFr[,1]==iL[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    lcc <- lcc1 
  } else {
    lcc <- rbind(lcc , lcc1)# agggregate merged results over all CRC samples
  }
}

nbDCC <- length(iD)
dcc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbDCC) {
  #i= 1
  dcc1 <- bFr[bFr[,1]==iD[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    dcc <- dcc1 
  } else {
    dcc <- rbind(dcc , dcc1)# agggregate merged results over all CRC samples
  }
}

nbBCC <- length(iB)
bcc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbBCC) {
  #i= 1
  bcc1 <- bFr[bFr[,1]==iB[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    bcc <- bcc1 
  } else {
    bcc <- rbind(bcc , bcc1)# agggregate merged results over all CRC samples
  }
}

nbGCC <- length(iG)
gcc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbGCC) {
  #i= 1
  gcc1 <- bFr[bFr[,1]==iG[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    gcc <- gcc1 
  } else {
    gcc <- rbind(gcc , gcc1)# agggregate merged results over all CRC samples
  }
}

nbPCC <- length(iP)
pcc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbPCC) {
  #i= 1
  pcc1 <- bFr[bFr[,1]==iP[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    pcc <- pcc1 
  } else {
    pcc <- rbind(pcc , pcc1)# agggregate merged results over all CRC samples
  }
}

nbOVC <- length(iO)
ovc <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbOVC) {
  #i= 1
  ovc1 <- bFr[bFr[,1]==iO[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    ovc <- ovc1 
  } else {
    ovc <- rbind(ovc , ovc1)# agggregate merged results over all CRC samples
  }
}

nbCRC <- length(iC)
crc <-matrix(,nrow=555,ncol=702)
c <-matrix(,nrow=nbCRC ,ncol=700)
nFc  <- vector()
# for all patients
for(i in 1:nbCRC) {
  #i= 1
    crc1 <- bFr[as.character(bFr[,1])==as.character(iC[i]),]
    print(dim(crc1))
    
    if (i != 3){
    nFc[i] <- sum(crc1[,3:702])
    }
    
    if (i == 10){
      crcMax <- crc1[,3:702]
    }
    
    if (i == 16){
      crcMin <- crc1[,3:702]
    }
    
    c1 <- colSums(crc1[,3:702])
    c [i,] <- c1
    # distribution for a bin crc1[1,3:702]
    if (i==1){
      crc <- crc1 
      crcB <- t(crc1[,3:702])
    } else {
      crc <- rbind(crc , crc1)# agggregate merged results over all CRC samples
      if (i != 3){
        crcB <- rbind(crcB , t(crc1[,3:702]))# agggregate merged results over all CRC samples
        #ctl2 <- cbind(ctl , ctl1)# agggregate merged results over all CRC samples
      }
    }
}

nbC1 <- length(iC1)
cs1 <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbC1) {
  #i= 1
  cs11 <- bFr[bFr[,1]==iC1[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs1 <- cs11 
  } else {
    cs1 <- rbind(cs1 , cs11)# agggregate merged results over all CRC samples
  }
}

nbC2 <- length(iC2)
cs2 <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbC2) {
  #i= 1
  cs21 <- bFr[bFr[,1]==iC2[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs2 <- cs21 
  } else {
    cs2 <- rbind(cs2 , cs21)# agggregate merged results over all CRC samples
  }
}

nbC3 <- length(iC3)
cs3 <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbC3) {
  #i= 1
  cs31 <- bFr[bFr[,1]==iC3[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs3 <- cs31 
  } else {
    cs3 <- rbind(cs3 , cs31)# agggregate merged results over all CRC samples
  }
}

nbC4 <- length(iC4)
cs4 <-matrix(,nrow=555,ncol=702)
# for all patients
for(i in 1:nbC4) {
  #i= 1
  cs41 <- bFr[bFr[,1]==iC4[i],]
  # distribution for a bin crc1[1,3:702]
  if (i==1){
    cs4 <- cs41 
  } else {
    cs4 <- rbind(cs4 , cs41)# agggregate merged results over all CRC samples
  }
}

cS <-c("PGDX5881P"  ,"PGDX5882P",  "PGDX5883P1", "PGDX5884P"  ,"PGDX5886P"  ,"PGDX5889P"  ,"PGDX5890P1" ,"PGDX5891P"  ,"PGDX5892P",  "PGDX5894P1",
       "PGDX5969P1" ,"PGDX5895P1",
"PGDX5899P1", "PGDX5900P1" ,"PGDX5903P1", "PGDX5963P" , "PGDX5964P" , "PGDX5965P" , "PGDX5966P",  "PGDX5967P"  ,"PGDX5968P",
"PGDX6249P1", "PGDX8820P1",
"PGDX8823P" , "PGDX8825P1", "PGDX8828P",  "PGDX8829P1")
dimnames(c)[[1]]<-cS
cF <- c(191, 192 ,194 ,195, 196, 197, 198, 200, 201, 203 ,204, 205 ,206 ,357, 358 ,362 ,364 ,367, 368 ,370)
c20<- c[,cL]
dimnames(c20)[[2]]<-cF

nbCTL <- length(iH) / 5 # devide by 5 until laptop comes
ctl <-matrix(,nrow=555,ncol=702)
ctlB <-matrix(,nrow=700,ncol=555)
h <-matrix(,nrow=nbCTL ,ncol=700)
nFh  <- vector()
#nFh <- c(nFh, 1:nbCTL)
# for all patients
for(i in 1:nbCTL) {
  #i= 1
  ctl1 <- bFr[as.character(bFr[,1])==as.character(iH[i]),]
  
  print(dim(ctl1))
  
  if (i == 10){
    ctlMax <- ctl1[,3:702]
  }
  
  if (i == 11){
    ctlMin <- ctl1[,3:702]
  }
  
  
  if (i != 9){
    nFh[i] <- sum(ctl1[,3:702])
  }
  
  h1 <- colSums(ctl1[,3:702])
  h[i,] <- h1
  # distribution for a bin ctl1[1,3:702]
  if (i==1){
    ctl <- ctl1 
    ctlB <- t(ctl1[,3:702])
  } else {
    ctl <- rbind(ctl , ctl1)# agggregate merged results over all CRC samples
    if (i != 9){
    ctlB <- rbind(ctlB , t(ctl1[,3:702]))# agggregate merged results over all CRC samples
    #ctl2 <- cbind(ctl , ctl1)# agggregate merged results over all CRC samples
    }
  }
}
for(i in (nbCTL+1):(nbCTL+nbCTL+1)) {
  #i= 1
  ctl1 <- bFr[bFr[,1]==iH[i],]
  
  print(dim(ctl1))
  
  # distribution for a bin ctl1[1,3:702]
  if (i==(nbCTL+1)){
    ctl <- ctl1 
  } else {
    ctl <- rbind(ctl , ctl1)# agggregate merged results over all CRC samples
  }
}
#compare all CTL vs all CRC samples in terms of overall number of fragments, i.e. coverage
nFh[9]<-nFh[43]
nFh<-nFh[1:42]
nFc[3]<-nFc[27]
nFc<-nFc[1:26]
plot(nFh, col="green",ylim=range(c(23000000,86000000)))
par(new = TRUE)
plot(nFc, col="red",ylim=range(c(23000000,86000000)))

#hS <-c("PGDX5881P"  ,"PGDX5882P",  "PGDX5883P1", "PGDX5884P"  ,"PGDX5886P"  ,"PGDX5889P"  ,"PGDX5890P1" ,"PGDX5891P"  ,"PGDX5892P",  "PGDX5894P1",
       #"PGDX5969P1" ,"PGDX5895P1",
       #"PGDX5899P1", "PGDX5900P1" ,"PGDX5903P1", "PGDX5963P" , "PGDX5964P" , "PGDX5965P" , "PGDX5966P",  "PGDX5967P"  ,"PGDX5968P",
       #"PGDX6249P1", "PGDX8820P1",
       #"PGDX8823P" , "PGDX8825P1", "PGDX8828P",  "PGDX8829P1")
dimnames(h)[[1]]<-hS
hF <- c(191, 192 ,194 ,195, 196, 197, 198, 200, 201, 203 ,204, 205 ,206 ,357, 358 ,362 ,364 ,367, 368 ,370)
h20<- h[,cL]
dimnames(h20)[[2]]<-hF

# for all bins
# extract FR histograms and assign to a small matrix
# repeat with Healthy
# plug two matrices into KLD


length(unique(sName))# 533 patients, 295815 bins in total, which is 555 bin per patient, which covers 2.775G of genome

ln <- as.integer(unlist(lcc[,3:702]))# convert to numbers
lnm <- matrix(ln, ncol = dim(lcc)[1], byrow = (dim(lcc)[2]-2)) # convert back to matrix form

dn <- as.integer(unlist(dcc[,3:702]))# convert to numbers
dnm <- matrix(dn, ncol = dim(dcc)[1], byrow = (dim(dcc)[2]-2)) # convert back to matrix form

bn <- as.integer(unlist(bcc[,3:702]))# convert to numbers
bnm <- matrix(bn, ncol = dim(bcc)[1], byrow = (dim(bcc)[2]-2)) # convert back to matrix form

gn <- as.integer(unlist(gcc[,3:702]))# convert to numbers
gnm <- matrix(gn, ncol = dim(gcc)[1], byrow = (dim(gcc)[2]-2)) # convert back to matrix form

pn <- as.integer(unlist(pcc[,3:702]))# convert to numbers
pnm <- matrix(pn, ncol = dim(pcc)[1], byrow = (dim(pcc)[2]-2)) # convert back to matrix form

on <- as.integer(unlist(ovc[,3:702]))# convert to numbers
onm <- matrix(on, ncol = dim(ovc)[1], byrow = (dim(ovc)[2]-2)) # convert back to matrix form

cn <- as.integer(unlist(crc[,3:702]))# convert to numbers
cnm <- matrix(cn, ncol = dim(crc)[1], byrow = (dim(crc)[2]-2)) # convert back to matrix form
cnmD<-cnm*2 # mimic double coverage
cnmH<-cnm/2 # mimic half coverage

cn1 <- as.integer(unlist(cs1[,3:702]))# convert to numbers
cnm1 <- matrix(cn1, ncol = dim(cs1)[1], byrow = (dim(cs1)[2]-2)) # convert back to matrix form

cn2 <- as.integer(unlist(cs2[,3:702]))# convert to numbers
cnm2 <- matrix(cn2, ncol = dim(cs2)[1], byrow = (dim(cs2)[2]-2)) # convert back to matrix form

cn3 <- as.integer(unlist(cs3[,3:702]))# convert to numbers
cnm3 <- matrix(cn3, ncol = dim(cs3)[1], byrow = (dim(cs3)[2]-2)) # convert back to matrix form

cn4 <- as.integer(unlist(cs4[,3:702]))# convert to numbers
cnm4 <- matrix(cn4, ncol = dim(cs4)[1], byrow = (dim(cs4)[2]-2)) # convert back to matrix form

hn <- as.integer(unlist(ctl[,3:702]))# convert to numbers
hnm <- matrix(hn, ncol = dim(ctl)[1], byrow = (dim(ctl)[2]-2)) # convert back to matrix form

#compare all bins CTL vs all bins CRC per FRL
plot(rowSums(hnm), col="green",ylim=range(c(0,55000000)))
par(new = TRUE)
plot(rowSums(cnm), col="red",ylim=range(c(0,55000000)))

#KL.divergence(t(cnm), t(hnm), k = 10, algorithm=c("kd_tree", "cover_tree", "brute")) # from FNN; input must be matrices
# 504.4316 507.1018 508.6900 509.2213 509.0071 508.4851 507.8492 507.7506 507.1370 506.2944

#KLD(t(cnm), t(hnm)) # from LaplacesDemon
write.csv(t(v364[sample(1:14430,14430)]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_c364_14430.csv')
write.csv(t(hnm[364,]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_h364_23310.csv')

write.csv(t(lnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_lcc.csv')
write.csv(t(dnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_dcc.csv')
write.csv(t(bnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_bcc.csv')
write.csv(t(gnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_gcc.csv')
write.csv(t(pnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_pcc.csv')
write.csv(t(onm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ovc.csv')
write.csv(t(cnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc.csv')
write.csv(t(cnmD),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcD.csv')# double
write.csv(t(cnmH),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcH.csv')# half- doesnt work well
write.csv(crcMax,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcMax.csv')# max coverage crc pt
write.csv(crcMin,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcMin.csv')# min coverage crc pt
write.csv(ctlMax,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlMax.csv')# max coverage crc pt
write.csv(ctlMin,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlMin.csv')# min coverage crc pt
write.csv(t(hnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl.csv')
write.csv(t(hnm),'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl2.csv')
ctlBB<-unname(ctlB)
write.csv(ctlBB,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctlB.csv')
crcBB<-unname(crcB)
write.csv(crcBB,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crcB.csv')
write.csv(hM364,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl364.csv')
write.csv(cM364,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc364.csv')
write.csv(hM205,'tyG:\\matovanalysis/DELFI_analysis/python/delfi1_ctl205.csv')
write.csv(cM205,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc205.csv')
write.csv(hM198,'G:\\matovanalysis/DELFI_analysis/python/delfi1_ctl198.csv')
write.csv(cM198,'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc198.csv')
write.csv(t(cnm1),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc1.csv')
write.csv(t(cnm2),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc2.csv')
write.csv(t(cnm3),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc3.csv')
write.csv(t(cnm4),'G:\\matovanalysis/DELFI_analysis/python/delfi1_crc4.csv')

write.csv(t(v364[sample(1:14430,100)]),'~/genomedk/matovanalysis/DELFI_analysis/python/delfi1_v364_100.csv')

k_ch2 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCh2.csv')
k_cmaxmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmaxmin.csv')

kcmax_hmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmaxCRCmax.csv')

k_hmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmax.csv')
k_hmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCTLmin.csv')

k_cbin2 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_bins364.csv')

k_cbin <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRCbins.csv')
k_cmax <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmax.csv')
k_cmin <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCmin.csv')
k_c12 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC2.csv')
k_c13 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC3.csv')
k_c14 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1CRC4.csv')
k_c23 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2CRC3.csv')
k_c24 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2CRC4.csv')
k_c34 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC3CRC4.csv')
k_cD <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCd.csv')
k_cH <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRCh.csv')
k_c1 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC1.csv')
k_c2 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC2.csv')
k_c3 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC3.csv')
k_c4 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC4.csv')
k_c364 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC364.csv')
k_c205 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC205.csv')
k_c198 <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC198.csv')
k_l <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceLCC.csv')
k_cl <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_LCC.csv')
k_d <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceDCC.csv')
k_cd <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_DCC.csv')
k_b <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceBCC.csv')
k_cb <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_BCC.csv')
k_cg <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_GCC.csv')
k_g <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceGCC.csv')
k_p <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergencePCC.csv')
k_c <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC.csv') # CRC  
k_cR <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceCRC_Reverse.csv') # CRC  REVERSE

k_o <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceOVC.csv')
k_co <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_OVC.csv')
k_cp <- read.csv('G:\\matovanalysis/DELFI_analysis/python/KLdivergenceCRC_PCC.csv')

kcbin2<-k_cbin2[2:701,2]
plot(kcbin2)
kcmaxmin<-k_cmaxmin[2:701,2]
plot(kcmaxmin)
kch2<-k_ch2[2:701,2]
plot(kch2)
khmaxcmax<-kcmax_hmax[2:701,2]
plot(khmaxcmax)
khmin<-k_hmin[2:701,2]
plot(khmin)
khmax<-k_hmax[2:701,2]
plot(khmax)
#kcbin<-k_cbin[2:701,2]
plot(k_cbin)
kcmax<-k_cmax[2:701,2]
plot(kcmax)
kcmin<-k_cmin[2:701,2]
plot(kcmin)
kcD<-k_cD[2:701,2]
plot(kcD)
kc23<-k_c23[2:701,2]
plot(kc23)
kc24<-k_c24[2:701,2]
plot(kc24)
kc34<-k_c34[2:701,2]
plot(kc34)
kc14<-k_c14[2:701,2]
plot(kc14)
kc13<-k_c13[2:701,2]
plot(kc13)
kc12<-k_c12[2:701,2]
plot(kc12)
kc1<-k_c1[2:701,2]
plot(kc1)
kc2<-k_c2[2:701,2]
plot(kc2)
kc3<-k_c3[2:701,2]
plot(kc3)
kc4<-k_c4[2:701,2]
plot(kc4)
kc364<-k_c364[2:701,2]
plot(kc364)
kc205<-k_c205[2:701,2]
plot(kc205)
kc198<-k_c198[2:701,2]
plot(kc198)
kl<-k_l[2:701,2]
plot(kl)
kcl<-k_cl[2:701,2]
plot(kcl)
kd<-k_d[2:701,2]
plot(kd)
kcd<-k_cd[2:701,2]
plot(kcd)
kb<-k_b[2:701,2]
plot(kb)
kcb<-k_cb[2:701,2]
plot(kcb)
kcg<-k_cg[2:701,2]
plot(kcg)
kg<-k_g[2:701,2]
plot(kg)
kcp<-k_cp[2:701,2]
plot(kcp)
kp<-k_p[2:701,2]
plot(kp)
kc<-k_c[2:701,2]
plot(kc)
kcr<-k_cR[2:701,2]
plot(kcr)
colorKLD <- rowSums(hnm)/rowSums(cnm)
div <- data.frame(kc) 
div$clr = colorKLD
ve <- seq(1, 700, by=1)
cK <- ggplot(div, aes(x= ve, y = kc, color = clr) )  + geom_point()
cK+scale_color_gradientn(colours = rainbow(5))
ko<-k_o[2:701,2]
plot(ko)
kco<-k_co[2:701,2]
plot(kco)

plot(khmin, ylim=range(c(0,4)), col="green", main = "Divergence [in bits] from Healthy of healthy with highest and lowest coverage", type="l")
par(new = TRUE)
plot(khmax, ylim=range(c(0,4)), col="red", type="l")
legend(500,3,legend=c("CTL pt w 64mln fragments", "CTL w 28mln fragments"),col=c("green","red"),lty=1:1, cex=1.0)

plot(kcmin, ylim=range(c(0,4)), col="orange", main = "Divergence [in bits] from Healthy of individual samples", type="l")
par(new = TRUE)
plot(kcmax, ylim=range(c(0,4)), col="red", type="l")
par(new = TRUE)
plot(khmin, ylim=range(c(0,4)), col="green", type="l")
par(new = TRUE)
plot(khmax, ylim=range(c(0,4)), col="blue", type="l")
legend(420,3,legend=c("CRC pt w 23mln fragments", "CRC pt w 85mln fragments","CTL pt w 28mln fragments", "CTL w 64mln fragments"),col=c("orange","red","green","blue"),lty=1:1, cex=1.0)

plot(kcmin, ylim=range(c(0,4)), col="green", main = "Divergence [in bits] from Healthy of individual samples", type="l")
par(new = TRUE)
plot(kcmax, ylim=range(c(0,4)), col="red", type="l")
legend(500,3,legend=c("CRC pt w 23mln fragments", "CRC pt w 85mln fragments"),col=c("green","red"),lty=1:1, cex=1.0)

plot(kc23, ylim=range(c(0,2.4)), col="green", main = "Divergence [in bits] from CRC II of CRC III,IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc24, ylim=range(c(0,2.4)), col="blue", type="l")
par(new = TRUE)
plot(kc34, ylim=range(c(0,2.4)), col="red", type="l")
legend(500,1.6,legend=c("CRC Stage II->III", "CRC Stage II->IV", "CRC Stage III->IV"),col=c("green","blue","red"),lty=1:1, cex=1.0)

plot(kc12, ylim=range(c(0,2.2)), col="green", main = "Divergence [in bits] from CRC I of CRC II-IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc13, ylim=range(c(0,2.2)), col="blue", type="l")
par(new = TRUE)
plot(kc14, ylim=range(c(0,2.2)), col="red", type="l")
legend(500,1.6,legend=c("CRC Stage I->II", "CRC Stage I->III", "CRC Stage I->IV"),col=c("green","blue","red"),lty=1:1, cex=1.0)

plot(kc1, ylim=range(c(0,5)), col="green", main = "Divergence [in bits] from healthy of CRC Stages I-IV per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kc2, ylim=range(c(0,5)), col="blue", type="l")
par(new = TRUE)
plot(kc3, ylim=range(c(0,5)), col="orange", type="l")
par(new = TRUE)
plot(kc4, ylim=range(c(0,5)), col="red", type="l")
legend(500,2.6,legend=c("CRC Stage I", "CRC Stage II", "CRC Stage III", "CRC Stage IV"),col=c("green","blue","orange","red"),lty=1:1, cex=1.0)

plot(kc, ylim=range(c(0,2.7)), col="green", main = "Divergence [in bits] from healthy for cancers per fragment length [in bp]", type="l")
par(new = TRUE)
plot(ko, ylim=range(c(0,2.7)), col="blue", type="l")
par(new = TRUE)
plot(kb, ylim=range(c(0,2.7)), col="orange", type="p")
par(new = TRUE)
plot(kp, ylim=range(c(0,2.7)), col="red", type="l")
par(new = TRUE)
plot(kg, ylim=range(c(0,2.7)), col="brown", type="l")
par(new = TRUE)
plot(kl, ylim=range(c(0,2.7)), col="pink", type="l")
par(new = TRUE)
plot(kd, ylim=range(c(0,2.7)), col="purple", type="l")
legend(500,2.6,legend=c("CRC", "Ovarian", "Breast", "Pancreatic", "Gastric", "Lung", "Bile duct"), col=c("green","blue","orange","red","brown","pink","purple"),lty=1:1, cex=1.0)

plot(kcd, ylim=range(c(0,1.7)), col="green", main = "Divergence [in bits] from CRC of other cancers per fragment length [in bp]", type="l")
par(new = TRUE)
plot(kco, ylim=range(c(0,1.7)), col="blue", type="l")
par(new = TRUE)
plot(kcb, ylim=range(c(0,1.7)), col="orange", type="p")
par(new = TRUE)
plot(kcp, ylim=range(c(0,1.7)), col="red", type="l")
par(new = TRUE)
plot(kcg, ylim=range(c(0,1.7)), col="brown", type="l")
par(new = TRUE)
plot(kcl, ylim=range(c(0,1.7)), col="purple", type="l")
legend(500,1.7,legend=c("Bile duct", "Ovarian", "Breast", "Pancreatic", "Gastric", "Lung"), col=c("green","blue","orange","red","brown","purple"),lty=1:1, cex=1.0)

sum(s1[which(kc>1)])/sum(s1) # % fragments belonging to KLD>1 based on one bin analysis only

sum(cnm[which(kc>2.8),])/sum(cnm) # 0.0001387829
sum(cnm[which(kc>2),])/sum(cnm) # 0.06610582
sum(cnm[which(kc>1),])/sum(cnm) # 0.1513125

# compare Healthy and CRC number of fragments per bin for FR=364 and 222
hist(cnm[364,], xlim = c(0,70), breaks = 20)
hist(hnm[364,], xlim = c(0,70), breaks = 20)

# 1st most discriminate bin for 205 bp
hgA <- hist(cM205[,369],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hM205[,369], xlim = c(0,70),breaks = 20, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,80), ylim = c(0,4)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,80), ylim = c(0,4)) # Add 2nd histogram using different color

# 1st most discriminate bin for 198 bp
hgA <- hist(cM198[,447], breaks = 10 , plot = FALSE) # Save first histogram data
hgB <- hist(hM198[,447], breaks = 25, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,600), ylim = c(0,9)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,600), ylim = c(0,9)) # Add 2nd histogram using different color

# 1st most discriminate bin for 364 bp
hgA <- hist(cM364[,308],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hM364[,308], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,40), ylim = c(0,6)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,40), ylim = c(0,6)) # Add 2nd histogram using different color

# 2nd most discriminate bin for 364 bp
hgA <- hist(cM364[,519],xlim = c(0,70), breaks = 10 , plot = FALSE) # Save first histogram data
hgB <- hist(hM364[,519], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,50), ylim = c(0,6)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,50), ylim = c(0,6)) # Add 2nd histogram using different color

hgA <- hist(ctlMin[,364],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(crcMin[,364],xlim = c(0,70), breaks = 20 , plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm[364,],xlim = c(0,70), breaks = 50, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[364,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm3[210,],xlim = c(0,70), breaks = 50, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[210,], xlim = c(0,70),breaks = 50, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,70), ylim = c(0,2200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,70), ylim = c(0,2200)) # Add 2nd histogram using different color

hgA <- hist(cnm[205,],breaks = 30, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[205,], breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,560),  ylim = c(0,3800)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,560),  ylim = c(0,3800)) # Add 2nd histogram using different color

hgA <- hist(cnm[168,],xlim = c(0,15000), breaks = 90, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[168,], xlim = c(0,15000),breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,4000), ylim = c(0,3200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,4000), ylim = c(0,3200)) # Add 2nd histogram using different color

hgA <- hist(cnm[164,],xlim = c(0,15000), breaks = 90, plot = FALSE) # Save first histogram data
hgB <- hist(hnm[164,], xlim = c(0,15000),breaks = 60, plot = FALSE) # Save 2nd histogram data
plot(hgA, col = rgb(1,0,0,1/10),xlim = c(0,4200), ylim = c(0,3200)) # Plot 1st histogram using a transparent color
plot(hgB, col = rgb(0,1,0,1/10), add = TRUE,xlim = c(0,4200), ylim = c(0,3200)) # Add 2nd histogram using different color

hist(hnm[168,])
hist(cnm[168,])

#SHOW A HEATMAP OF HOW NBS FOR 364 vary per bin (555bins) per CRC patient or normal
# 555x27 and 555x43
cM364 <- matrix(cnm[364,], ncol = 555, byrow = 26) #  
hM364 <- matrix(hnm[364,], ncol = 555, byrow = 42) # 
cM205 <- matrix(cnm[205,], ncol = 555, byrow = 26) #  
hM205 <- matrix(hnm[205,], ncol = 555, byrow = 42) # 
cM198 <- matrix(cnm[198,], ncol = 555, byrow = 26) #  
hM198 <- matrix(hnm[198,], ncol = 555, byrow = 42) # 
# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC364<-apply(cM364,2,var)
vH364<-apply(hM364,2,var)
mC364<-apply(cM364,2,mean)
mH364<-apply(hM364,2,mean)

plot(vC364, ylim=range(c(0,200)), col="red")
par(new = TRUE)
plot(vH364, ylim=range(c(0,200)), col="green")
plot(mC364, ylim=range(c(0,30)), col="red")
par(new = TRUE)
plot(mH364, ylim=range(c(0,30)), col="green")

dV364 <- vH364-vC364
dM364 <- mH364-mC364
plot(dV364)
which(dV364==max(dV364))# for FRL364, CRL is more dispersed in bin #475
which(dV364==min(dV364)) # for FRL364, CRC is more dispersed in bin #146
rV364 <- vH364/vC364
which(rV364>8) # bins # 74 220 282 373 444 453 455 495 498 525 553
which(rV364>9) # bins # 74 220 453 495
plot(dM364)
which(dM364==max(dM364))# for FRL364, CRL is more dispersed in bin #475
which(dM364==min(dM364)) # for FRL364, CRC is more dispersed in bin #146

#SHOW A HEATMAP OF HOW NBS FOR 202 vary per bin (555bins) per CRC patient or normal
# 555x27 and 555x43
cM205 <- matrix(cnm[205,], ncol = 555, byrow = 26) #  
hM205 <- matrix(hnm[205,], ncol = 555, byrow = 42) # 
# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC205<-apply(cM205,2,var)
vH205<-apply(hM205,2,var)

plot(vC205, ylim=range(c(0,5000)), col="red")
par(new = TRUE)
plot(vH205, ylim=range(c(0,5000)), col="green")

dV205 <- vH205-vC205
plot(dV205)
which(dV205==max(dV205))# for FRL364, CRL is more dispersed in bin #475
which(dV205==min(dV205)) # for FRL364, CRC is more dispersed in bin #146
rV205 <- vH205/vC205

cM198 <- matrix(cnm[198,], ncol = 555, byrow = 26) #  
hM198 <- matrix(hnm[198,], ncol = 555, byrow = 42) # 
# check variance across bins/genomic positions to find highest variability/difference CRC to CTL
vC198<-apply(cM198,2,var)
vH198<-apply(hM198,2,var)

plot(vC198, ylim=range(c(0,12000)), col="red")
par(new = TRUE)
plot(vH198, ylim=range(c(0,12000)), col="green")

dV198 <- vH198-vC198
plot(dV198)
which(dV198==max(dV198))# for FRL364, CRL is more dispersed in bin #475
which(dV198==min(dV198)) # for FRL364, CRC is more dispersed in bin #146

# mix for clustering test
c364 <- cnm[364,]
h364 <- hnm[364,]

ze36 = c (1,  16, 14,  9,   12,  4,  4,  6,  7,   14, 15 , 9, 2, 2, 3, 5, 6, 7, 8, 9, 9, 11, 14, 19, 19 ,2, 4, 11, 14, 9, 14, 1, 1,
        1, 1 ,1 ,2 ,2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10,
        10, 10, 11, 11, 12, 12, 14, 14, 13, 14, 15, 16, 1, 16, 17, 17, 18, 18, 1, 2, 2, 21, 4, 9, 11, 1, 1, 1, 1, 1, 1, 1, 2, 2,2,  
         3, 3, 3, 3, 3, 3, 4, 4,  4, 4, 4, 4, 4, 4 ,4, 5, 5, 6, 6, 6, 6,6, 7,  7, 7, 7, 8, 8, 8,9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
        10, 11, 12, 12, 12, 12, 12,13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 17, 17, 17, 18, 19, 19, 20, 20, 21, 21, 9)
 
cL <- which(kc>2.38) # 191 192 194 195 196 197 198 200 201 203 204 205 206 357 358 362 364 367 368 370

#cT <- crc[91:100,]# select 50 random bins of 14430
cT <- crc[,3:702]
#cT[,1]="CRC"

#hT <- ctl[91:100,]# select 50 random bins of 23310
hT <- ctl[,3:702]
#hT[,1]="HLT"

cT20 <- cT[,cL]# pick all bins corresponding to the top 20
hT20 <- hT[,cL]# pick all bins corresponding to the top 20

cT1 <- cT[,364]
hT1<-hT[,364]
cT2 <- cT[,205]
hT2<-hT[,205]
cT3 <- cT[,198]
hT3<-hT[,198]

cv364 <- matrix(cT1, ncol = 555, byrow = 26) #
cv205 <- matrix(cT2, ncol = 555, byrow = 26) #  
cv198 <- matrix(cT3, ncol = 555, byrow = 26) #  
indx<-which(kc364>1.5) #  genomic bins: 49  64  90 130 132 209 216 267 308 324 508 519
cv2 <-matrix(,nrow=26,ncol=length(indx))
cv2<-cv364[,indx]
#cv2[,1]<-cv364[,308]
#cv2[,2]<-cv364[,519]
#cv2[,3]<-cv205[,282]
#cv2[,4]<-cv205[,369]
#cv2[,5]<-cv198[,36]
#cv2[,6]<-cv198[,447]
hv364 <- matrix(hT1, ncol = 555, byrow = 42) #  
hv205 <- matrix(hT2, ncol = 555, byrow = 42) # 
hv198 <- matrix(hT3, ncol = 555, byrow = 42) # 
hv2 <-matrix(,nrow=42,ncol=length(indx))
hv2<-hv364[,indx]
#hv2[,1]<-hv364[,308]
#hv2[,2]<-hv364[,519]
#hv2[,3]<-hv205[,282]
#hv2[,4]<-hv205[,369]
#hv2[,5]<-hv198[,36]
#hv2[,6]<-hv198[,447]

dimnames(cv2)[[1]] <- unique(crc[,1]) 
#dimnames(cv2)[[2]] <-  c("308","519","282","369","36","447")

hcTop20 <- rbind(cv2 , hv2)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) # 3 CTL mixed with CRC: #36, 55, 44, i.e: CTL #10, 29, 18. 
# hm$rowInd CRC 1-26, CTL 27-68 
#(bottom) 25 24 10 22  4 23 44 61 20  6 15 21 11 12 55 13 14 68 36 52  7 19 16  5  8  3  2 17  9 18  1 54 48 51 63 
# 67 28 56 49 62 50 53 33 58 47 45 66 29 30 26 64 46 59 57 42 60 65 27 31 34 32 40 41 38 43 37 35 39 (top end)
# hm$colInd order of 12 bins: 4  6 10  8  7  9  3  2 12  5  1 11

distr <- dist(df)
hcr <- hclust(distr)
hcr$order # 38 41 37 43 35 39 27 31 34 32 40 50 53 47 45 66 33 58 62 49 56 28 67 63 51 48 54 29 30 26 64 57 42 60 65 46 59  
#8  5 16 19 17  2  3  1  9 18 14 13 55 11 12 10 24 25 23 4 22 15 21  6 20 44 61  7 52 36 68

d2 <- dist(hcTop20,method = "euclidean", diag = FALSE, upper = TRUE)
c2 <- hclust(d2, method = "ward.D2", members = NULL)

#heatmap(x, Rowv = NULL, Colv = if(symm)"Rowv" else NULL,
#        distfun = dist, hclustfun = hclust,
#        reorderfun = function(d, w) reorder(d, w),
#        add.expr, symm = FALSE, revC = identical(Colv, "Rowv"),
#        scale = c("row", "column", "none"), na.rm = TRUE,
#        margins = c(5, 5), ColSideColors, RowSideColors,
#        cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
#        labRow = NULL, labCol = NULL, main = NULL,
#        xlab = NULL, ylab = NULL,
#        keep.dendro = FALSE, verbose = getOption("verbose"), …)

heatmap.2(df)

y <- df
## Row- and column-wise clustering
hr <-hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <-hclust(as.dist(1-cor(y, method="spearman")), method="complete")
## Tree cutting
mycl <-cutree(hr, h=max(hr$height)/1.5); mycolhc <-rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)]
## Plot heatmap
mycol <-colorpanel(40, "darkblue", "yellow", "white") 
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none")




dNc= c("CRC1","CRC2","CRC3","CRC4","CRC5","CRC6","CRC7","CRC8","CRC9","CRC10")#,"CRC11","CRC12","CRC13","CRC14","CRC15","CRC16",
       #"CRC17","CRC18","CRC19",
      #"CRC20","CRC21","CRC22","CRC23","CRC24","CRC25","CRC26","CRC27","CRC28","CRC29","CRC30","CRC31","CRC32","CRC33","CRC34",
      #"CRC35","CRC36","CRC37","CRC38","CRC39",
      #"CRC40","CRC41","CRC42","CRC43","CRC44","CRC45","CRC46","CRC47","CRC48","CRC49","CRC50","CRC51")

cpts <-unique(crc[,1])
cpts <- unique(crc$sample)[26]

cPtL <- crc [,1] == unique(crc$sample)[1]

crcP1 <- crc[,crc[,1]==cpts[1]]




dimnames(cT20)[[1]] <- unique(crc[,1]) # 27 CRC patients id numbers

c20[3,]<-c20[27,]# 3rd CRC is all zeros
c20<-c20[1:26,]

h20[9,]<-h20[43,]# 9th Healthy is all zeros
h20<-h20[1:42,]

hcTop20 <- rbind(c20 , h20)

#c21<-transpose(c20)# dimension names get confused

df<-scale(hcTop20)
#df<-scale(cM364)
#df<-scale(hM364)


col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)

heatmap(df, scale = "none", col =  col)

cT1<-c20[,17]
hT1<-h20[,17]
#rocCH <- append(c20,h20)
TP <- matrix(rep(1, 26*12), ncol = 12, byrow = 26)
TN <- matrix(rep(0, 42*12), ncol = 12, byrow = 42)
rocTF <- rbind(TP , TN)

#rocTF <- append(rep(1, 26),integer(42))
#pred <- prediction(rocCH, rocTF)

pred <- prediction(hcTop20, rocTF)
perf<-performance(pred,"tpr", "fpr")
plot(perf)
auc<- performance(pred,"auc")
auc

sName = bFr[,1] # list of sample names
typeof(sName) # [1] "character"
#list
#typeof(aux)
#[1] "character
#as.numeric(str_extract_all(string, aux)[[1]])
#as.numeric(strsplit(string, "aux")[[1]][-1])

cM364[,308]
#  5  7  8 12  5 12 17  2  7  8 20 17 13 17  7  1  9 12  1 12  8 10  6 10  5 25
hM364[,308]
# 31 20 28 20 33 23 22 19 41 15 33 41 41 28 34 22 35 13 18 22 21 18 14 19 19  5 23 14 15 17 22 22 24 16 11 19 22 23 18 21 13 15
mean(cM364[,308])
# 9.846154
mean(hM364[,308])
# 22.14286
mean(cM364[,519])
# 8.538462
mean(hM364[,519])
# 20.85714

cM364[,519]
# 15 15  8  7  4 12  8  6  8  6  7 11 11  5 13  2  9 11  4  9  6  3 11  8  8 15
hM364[,519]
# 23 16 16 17 31 23 17 42 31 18 43 34 49 34 35 17 29 16 20 14 21 16 22 26 16 10 20 10 10 18  6 17 18 17 15 15 16 13 18 23 10 14

cM198[,447]
# 132 215 186 227 121 294 201 114 226 175 189 193 182 194 163  72 107 192 110 185 197 129 185 210 162 255
hM198[,447]
# 422 365 323 305 426 365 277 431 508 208 412 431 371 379 440 361 535 203 226 251 379 249 274 385 282 214 359 284 200 293 255 226 254 282 145 232 305 216 284 325 246 184
wasserstein1d(cM364[,519], hM364[,519])
wasserstein(cM364,hM364,p=1)