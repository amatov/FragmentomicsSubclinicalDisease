
grays = rgb(red = 0:255/255, blue = 0:255/255, green = 0:255/255)
#heatmap(ctl1D2[74,,],Rowv=NA,Colv=NA,col=grays, scale = "none")
#heatmap(col4D2[23,,],Rowv=NA,Colv=NA,col=grays, scale = "none")
#heatmap(col1D2[8,,],Rowv=NA,Colv=NA,col=grays, scale = "none")
#heatmap(umiiD2[1,,],Rowv=NA,Colv=NA,col=grays, scale = "none")
#heatmap(umicD2[1,,],Rowv=NA,Colv=NA,col=grays, scale = "none")
heatmap(umiD2[1,,],Rowv=NA,Colv=NA,col=grays, scale = "none")

# all stages 79 colon cancers vs 74 control1 no comorbidity
samplesTr = array(0, dim=c((69+64),574,499))
samplesTr[1:69,,] <- colD2[1:69,,]
samplesTr[70:133,,] <- ctl1D2[1:64,,]
selectionTr = rep(1, 133)
selectionTr[70:133]=0
write.csv(samplesTr,'~/genomedk/matovanalysis/DELFI_analysis/python/samplesTr.csv')
write.csv(selectionTr,'~/genomedk/matovanalysis/DELFI_analysis/python/selectionTr.csv')
samplesTe = array(0, dim=c((10+10),574,499))
samplesTe[1:10,,] <- colD2[70:79,,]
samplesTe[11:20,,] <- ctl1D2[65:74,,]
selectionTe = rep(1, 20)
selectionTe[11:20]=0
write.csv(samplesTe,'~/genomedk/matovanalysis/DELFI_analysis/python/samplesTe.csv')
write.csv(selectionTe,'~/genomedk/matovanalysis/DELFI_analysis/python/selectionTe.csv')
#########################################################
# all stages 79 colon cancers, 20 colon adeH, 28 colon adeL, 74 control1 no comorbidity
samplesTr4 = array(0, dim=c((69+10+18+64),574,499))
samplesTr4[1:69,,] <- colD2[1:69,,] # class 4
samplesTr4[70:79,,] <- col0D2[1:10,,] # class 3
samplesTr4[80:97,,] <- colAD2[1:18,,] # class 2
samplesTr4[98:161,,] <- ctl1D2[1:64,,] # class 1
selectionTr4 = rep(4, 161)
selectionTr4[70:79]=3
selectionTr4[80:97]=2
selectionTr4[98:161]=1
write.csv(samplesTr4,'~/genomedk/matovanalysis/DELFI_analysis/python/samplesTr4.csv')
write.csv(selectionTr4,'~/genomedk/matovanalysis/DELFI_analysis/python/selectionTr4.csv')
npySave("~/genomedk/matovanalysis/DELFI_analysis/python/samplesTr4.npy", samplesTr4)
sTr4 <- read.csv('~/genomedk/matovanalysis/DELFI_analysis/python/samplesTr4.csv')

samplesTe4 = array(0, dim=c((10+10+10+10),574,499))
samplesTe4[1:10,,] <- colD2[70:79,,]
samplesTe4[11:20,,] <- col0D2[11:20,,] # class 3
samplesTe4[21:30,,] <- colAD2[19:28,,] # class 2
samplesTe4[31:40,,] <- ctl1D2[65:74,,]
selectionTe4 = rep(4, 40)
selectionTe4[11:20]=3
selectionTe4[21:30]=2
selectionTe4[31:40]=1
write.csv(samplesTe4,'~/genomedk/matovanalysis/DELFI_analysis/python/samplesTe4.csv')
write.csv(selectionTe4,'~/genomedk/matovanalysis/DELFI_analysis/python/selectionTe4.csv')
