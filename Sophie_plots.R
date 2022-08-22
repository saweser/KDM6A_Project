###
rm(list=ls());
load('Sophie.RData');

##### Quality Control ######################################################################################
qcReport(RGset, sampNames = pd$Sample_Name,sampGroups = pd$Sample_Name, pdf = "qcReport.pdf") # macht Probleme wenn er kein Sample_Group hat 

# Density Plot #
densityPlot(RGset, sampGroups = pd$Sample_Group, main = "Beta", xlab = "Beta",legend = TRUE,pal = brewer.pal(3,"Dark2"))
# pal = heat.colors(29,alpha = 1)
# pal = rainbow(29)
# sampGroup weglassen oder = NULL zeigt alle werte
# sampGroup=pd$Sample_Name zeigt nur 8 samples an....wieso?--> weil er als default pal=brewer.pal(8,"Dark2") nimmt und diese palette hat nicht mehr als 8 farben
# andere farbpallette nehmen !Palettes (grDevices)
# legend kann man nicht auf 2 columns ändern! also lieber legend =FALSE
# uses getBeta to get the beta values.

# Density Bean Plot # 
par(oma=c(0,3,0,0))
densityBeanPlot(RGset, sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)

# MDS Plot
mdsPlot(Mset.raw, numPositions = 1000, sampGroups = pd$Sample_Group, sampNames = NULL,pch = 16)

# QC plot #
qc <- getQC(Mset.raw)
head(qc)
plotQC(qc,badSampleCutoff = 10.5)
# haben log2 x und y achsen scala, rechnen werte +1 
# default für badSampleCutoff = 10.5......aber gerade ist bei über 13??

# eigener QC plot #
qc.data.frame<-data.frame(qc)
qc.data.frame<-qc.data.frame+1
p<-ggplot(qc.data.frame,aes(y=uMed, x=mMed))
p.point<-p+geom_point(size=4,shape=1) +ylab("uMed (log2)")+xlab("mMed (log2)")+ xlim(8, 14) +ylim(8,14)
p.point+coord_trans(x="log2")

# quality strip plot #
controlStripPlot(RGset, controls = c("BISULFITE CONVERSION I","BISULFITE CONVERSION II"),sampNames = pd$Sample_Name)

# plot betas by type?
# plotCpg?

# shinyMethyl #
summary <- shinySummarize(RGset)
runShinyMethyl(summary)
# wieso keine zuordnung zu chips?
########################################################################################################


##### Cell Type Composition #####
require(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(RGset)
View(cellCounts)
#####################################################################################################

##### Predicted Sex #####
predictedSex <- getSex(gset.funnorm, cutoff = -2)$predictedSex
head(predictedSex)
plotSex(getSex(gset.funnorm, cutoff = -2),id=pd$Sex)
# used id=pd$Sex um predicted mit sex aus sample sheet zu vergleichen.

plotSex(RGset)
addSex(gset.funnorm)

pd2<-pData(gset.funnorm)
pd
pd2$predictedSex
pd2$Sex
