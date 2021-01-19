ctl1D23 = array(0, dim=c(74,499))
for (i in 1:74  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  ctl1D23[j,] <- unlist(colSums(auxFR[,2:500])) / sum(unlist(auxFR[,2:500]))
  
}
colD23 = array(0, dim=c(79,499))
for (i in 1:79  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  colD23[j,] <- unlist(colSums(auxFR[,2:500]))  / sum(unlist(auxFR[,2:500]))
}
x <- rbind(ctl1D23, colD23)
y = rep(0, 153)
y[75:153]=1
# x 499 by number of samples, both crc and ctl, 
# y vector 1s for crc and 0s for ctl

cvm = cv.glmnet(x, y, family = "binomial", alpha=1, nfolds=10) 

plot(cvm)

plot(cvm$glmnet.fit) 