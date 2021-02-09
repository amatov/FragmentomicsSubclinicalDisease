
grays = rgb(red = 0:255/255, blue = 0:255/255, green = 0:255/255)
heatmap(ctl1D2[2,,],Rowv=NA,Colv=NA,col=grays, scale = "none")

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
