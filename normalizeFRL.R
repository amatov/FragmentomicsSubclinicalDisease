ctl1D2N = array(0, dim=c(74,574,499))
j=1
for (i in 1:74  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  
  ctl1D2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  
  # save 74 individual profiles as xls files
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1N_individual', i,'.csv')
  auxSAVE <-ctl1D2N[j,,] #matrix( unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))),nrow=574,ncol=499)
  
  write.csv(unlist(auxSAVE),sa_name)
  j=j+1
}

colD2N = array(0, dim=c(79,574,499))
j=1
for (i in 1:79  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  
  colD2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  
  # save 79 individual profiles as xls files
  sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_colN_individual', i,'.csv')
  auxSAVE <-colD2N[j,,] #matrix( unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))),nrow=574,ncol=499)
  
  write.csv(unlist(auxSAVE),sa_name)
  j=j+1
}
recD2N = array(0, dim=c(50,574,499))
j=1
for (i in 1:50  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listREC[i]], header = TRUE) # sample per sample, file per file. 
  
  recD2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  
  # save 79 individual profiles as xls files
  #sa_name <- paste0('~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_colN_individual', i,'.csv')
  #auxSAVE <-colD2N[j,,] #matrix( unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))),nrow=574,ncol=499)
  
  #write.csv(unlist(auxSAVE),sa_name)
  j=j+1
}
crcD2N <- abind(colD2N,recD2N,along=1)
# convert to 2D format
ctl1D22N = array(0, dim=c(74*574,499))
colD22N = array(0, dim=c(79*574,499))
for (i in 1:499) {
  #i=1
  auxCTL <- ctl1D2N[,,i]
  ctl1D22N[,i] <- auxCTL 
  auxCOL <- colD2N[,,i]
  colD22N[,i] <- auxCOL
}
# save cohorts to disk
write.csv(colD22N*100000000000000,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_colN_all79.csv')
write.csv(ctl1D22N*100000000000000,'~/genomedk/matovanalysis/DELFI_analysis/python/delfi2_ctl1N_74.csv')

# read and plot KLD vector for FRL1-499bp
k_d2_colN <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/KLdivergenceD2_COL79_CTL1_norm.csv')
kd2_colN<- k_d2_colN[2:575,2]
plot(kd2_colN)

# UMIIMPROVE
pileupsUI <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/57PRE5Mb", recursive = T, full.names = T, pattern = "tsv")
nbUMII <- length(pileupsUI)
umiiN = array(0, dim=c(nbUMII,595,499))
j=1
for(i in 1:nbUMII) {
  #i= 2
  print(i)
  auxFR <- read.table( pileupsUI[i], header = TRUE) # sample per sample, file per file. 
  
  umiiN[j,,] <- as.integer(unlist(auxFR[,2:500]))/ sum(unlist(auxFR[,2:500])) 
  j=j+1
}
dim(umiiN) # 17850x499 for 30PreOps, 33320x499 for 56 PreOps
# UMICRUK
pileupsUC <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/CRUK5Mb", recursive = T, full.names = T, pattern = "tsv")
nbUMIC <- length(pileupsUC)
umicN = array(0, dim=c(nbUMIC,595,499))
j=1
for(i in 1:nbUMIC) {
  #i= 1
  print(i)
  auxFR <- read.table( pileupsUC[i], header = TRUE) # sample per sample, file per file. 
  
  umicN[j,,] <- as.integer(unlist(auxFR[,2:500]))/ sum(unlist(auxFR[,2:500]))
  j=j+1
}
dim(umicN) # 40460 x 499 for  68 CRUK PreOps
#UMISEQ
pileupsU <- list.files("~/genomedk/PolyA/faststorage/BACKUP/N140_Targeting/specs/umiseq_paper/divergence/data/45PON5Mb", recursive = T, full.names = T, pattern = "tsv")
nbUMI <- length(pileupsU)
umiN = array(0, dim=c(nbUMI,595,499))
j=1
for(i in 1:nbUMI) {
  #i= 2
  print(i)
  auxFR <- read.table( pileupsU[i], header = TRUE) # sample per sample, file per file. 
  
  umiN[j,,] <- as.integer(unlist(auxFR[,2:500])) / sum(unlist(auxFR[,2:500]))
  j=j+1
}
dim(umiN) # 17850   499 for 30PONs, 26775   499 for 45 PONs
