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
