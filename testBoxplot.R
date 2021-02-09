FrCTL1D2 <- vector(); Rctl1D2<-vector()
for (i in 1:74){
  Rctl1D2[i]<-sum(ctl1D2[i,,]<150)/sum(ctl1D2[i,,])
  FrCTL1D2[i] <- sum(ctl1D2[i,,])
}
FrColD2 <- vector(); RcolD2<-vector()
for (i in 1:79){
  RcolD2[i]<-sum(colD2[i,,]<150)/sum(colD2[i,,])
  FrColD2[i] <- sum(colD2[i,,])
}

FrCTL1D1 <- vector(); Rctl1D1<-vector()
for (i in 1:43){
  auxFR <- read.table(pileupsD1[unlist(listD1CTL215)[i]], header = TRUE) # sample per sample, file per file. 
  Rctl1D1[i]<-sum(auxFR[,2:500]<150)/sum(auxFR[,2:500])
  FrCTL1D1[i] <- sum(auxFR[,2:500])
}
FrCrcD1 <- vector(); RcrcD1<-vector()
for (i in 1:27){
  auxFR <- read.table(pileupsD1[listD1CRC27[i]], header = TRUE) # sample per sample, file per file.
  RcrcD1[i]<-sum(auxFR[,2:500]<150)/sum(auxFR[,2:500])
  FrCrcD1[i] <- sum(auxFR[,2:500])
}
erPQ <- list(pShFrD1crc=RcrcD1,pShFrD1ctl= Rctl1D1,pShFrD2cc=RcolD2,pShFrD2ctl=Rctl1D2)
boxplot(erPQ,notch = TRUE,horizontal = TRUE,border = "brown",col = c("blue","green","brown","orange"))

erPQ <- list(nFrD1crc=FrCrcD1,nFrD1ctl=FrCTL1D1,nFrD2cc=FrColD2,nFrD2ctl=FrCTL1D2)
boxplot(erPQ,notch = TRUE,horizontal = TRUE,border = "brown",col = c("blue","green","brown","orange"))