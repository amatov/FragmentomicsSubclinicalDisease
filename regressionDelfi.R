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



# Splitting the data into test and train
samplesTr = array(0, dim=c((37+40),499))
#samplesTr <- rbind(ctl1D23[1:37,], colD23[1:40,])
samplesTr <- rbind(ctl1D2[,,195][1:37,], colD2[,,195][1:40,])
selectionTr = rep(0, 77)
selectionTr[41:77]=1

samplesTe = array(0, dim=c((37+39),499))
#sampplesTe<- rbind(ctl1D23[38:74,], colD23[41:79,])
samplesTe <- rbind(ctl1D2[,,195][38:74,], colD2[,,195][41:79,])
selectionTe= rep(0, 76)
selectionTe[40:76]=1 

cvm = cv.glmnet(samplesTr, selectionTr, family = "binomial", alpha=1, nfolds=10) 

# identifying best lamda
best_lam <- cvm$lambda.min
best_lam # 0.02529132
# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(samplesTr, selectionTr, alpha = 1, lambda = best_lam)
pred <- predict(lasso_best, s = best_lam, newx = samplesTe)

plot(cvm)
plot(cvm$glmnet.fit)

final <- cbind(selectionTe, pred)
final

