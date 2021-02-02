ctl1D23 = array(0, dim=c(74,499))
for (i in 1:74  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCTL1[i]], header = TRUE) # sample per sample, file per file. 
  ctl1D23[i,] <- unlist(colSums(auxFR[,2:500])) / sum(unlist(auxFR[,2:500]))
  
}
colD23 = array(0, dim=c(79,499))
for (i in 1:79  ) {
  print(i)
  #i=2
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  colD23[i,] <- unlist(colSums(auxFR[,2:500]))  / sum(unlist(auxFR[,2:500]))
}
#x <- rbind(ctl1D23[1:37,], colD23[1:37,])
#y = rep(0, 153)
#y[75:153]=1

samplesTr = array(0, dim=c((37+40),499))
samplesTr <- rbind(ctl1D23[1:37,], colD23[1:40,])
selectionTr = rep(0, 77)
selectionTr[38:77]=1

samplesTe = array(0, dim=c((37+39),499))
samplesTe <- rbind(ctl1D23[38:74,], colD23[41:79,])
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

cvm = cv.glmnet(samplesTr, selectionTr, family = "binomial", alpha=1, nfolds=10) 

# identifying best lamda
best_lam <- cvm$lambda.min
best_lam 
# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(samplesTr, selectionTr, alpha = 1, lambda = best_lam)
pred <- predict(lasso_best, s = best_lam, newx = samplesTe)

plot(cvm)
plot(cvm$glmnet.fit)

final <- cbind(selectionTe, pred)
final
################################################################
colD2N = array(0, dim=c(129,574,499))
j=1
for (i in 1:129  ) {
  print(i)
  #i=2
  if (i<80){
  auxFR <- read.table(pileupsD2[listCOL[i]], header = TRUE) # sample per sample, file per file. 
  } else {
    k <- i-79
  auxFR <- read.table(pileupsD2[listREC[k]], header = TRUE) # sample per sample, file per file. 
  }
  colD2N[j,,] <- unlist(auxFR[,2:500])/ sum(unlist(auxFR[,2:500]))
  j=j+1
}
# Splitting the data into test and train
samplesTr = array(0, dim=c((37+40),length(ind195)))
#samplesTr <- rbind(ctl1D23[1:37,], colD23[1:40,])
samplesTr <- rbind(ctl1D2N[,ind195,195][1:37,], colD2N[,ind195,195][1:40,])
selectionTr = rep(0, 77)
selectionTr[38:77]=1

samplesTe = array(0, dim=c((37+39),length(ind195)))
#sampplesTe<- rbind(ctl1D23[38:74,], colD23[41:79,])
samplesTe <- rbind(ctl1D2N[,ind195,195][38:74,], colD2N[,ind195,195][41:79,])
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

cvm = cv.glmnet(samplesTr, selectionTr, family = "binomial", alpha=1, nfolds=10) 

# identifying best lamda
best_lam <- cvm$lambda.min
best_lam 
# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(samplesTr, selectionTr, alpha = 1, lambda = best_lam)
pred <- predict(lasso_best, s = best_lam, newx = samplesTe)

plot(cvm)
plot(cvm$glmnet.fit)

final <- cbind(selectionTe, pred)
final
plot(final)
########################################################################################

# Splitting the NORMALIZED data into test and train
samplesTr = array(0, dim=c((37+40),574*2))
ctlTr <- rbind(t(ctl1D2N[,,195][1:37,]), t(ctl1D2N[,,137][1:37,]))
#ctlTr <- rbind(ctlTr, t(ctl1D2[,,365][1:37,]))
colTr <- rbind(t(colD2N[,,195][1:40,]), t(colD2N[,,137][1:40,]))
#colTr <- rbind(colTr, t(colD2[,,365][1:40,]))
samplesTr <- rbind(t(ctlTr),t(colTr))# dim(samplesTr) 77 574 w 37 ctl and 40 cc
selectionTr = rep(0, 77)
selectionTr[38:77]=1

samplesTe = array(0, dim=c((37+39),574*2))
ctlTe <- rbind(t(ctl1D2N[,,195][38:74,]), t(ctl1D2N[,,137][38:74,]))
#ctlTe <- rbind(ctlTe, t(ctl1D2[,,365][38:74,]))
colTe <- rbind(t(colD2N[,,195][41:79,]), t(colD2N[,,137][41:79,]))
#colTe <- rbind(colTe, t(colD2[,,365][41:79,]))
samplesTe <- rbind(t(ctlTe), t(colTe))
selectionTe= rep(0, 76)
selectionTe[38:76]=1 

cvm = cv.glmnet(samplesTr, selectionTr, family = "binomial", alpha=1, nfolds=10) 

# identifying best lamda
best_lam <- cvm$lambda.min
best_lam 
# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(samplesTr, selectionTr, alpha = 1, lambda = best_lam)
pred <- predict(lasso_best, s = best_lam, newx = samplesTe)

plot(cvm)
plot(cvm$glmnet.fit)

final <- cbind(selectionTe, pred)
final
