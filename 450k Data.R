
#################### 450k Data Sophie ########################################################################
rm(list=ls());  # empty workspace

require(minfi)
library(shinyMethyl)
library(RColorBrewer)
library(ggplot2)
library(reshape)
setwd("/Users/Bria/Documents/KUM/R/Sophie")
#############################################################################################################

##### Read Data ##########################################################################################
baseDir<-file.path(getwd(),"idat")
targets<-read.450k.sheet(baseDir)
targets$Sample_Group<-c(rep("SPRYD2",12),rep("PDX",3),rep("Cell lines",14))
list.files(baseDir)
list.files(file.path(baseDir, "200397860014_idat"))
RGset <- read.450k.exp(targets = targets)
pd <- pData(RGset)
############################################################################################################

##### Preprocess #####
Mset.raw <- preprocessRaw(RGset)
Mset.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "no", reference = 2)
gset.funnorm <- preprocessFunnorm(RGset)
# i would used funnorm--> paper form kasper hansen 2014
#############################################################################################################

##### What normalization is better? ########################################################################
Mset.raw.beta<-getBeta(Mset.raw)
Mset.norm.beta<-getBeta(Mset.norm)
gset.funnorm.beta<-getBeta(gset.funnorm)

#### convert data frame to long format ######
Mset.raw.beta.long<-melt(Mset.raw.beta)
colnames(Mset.raw.beta.long)<-c("CpG","sample","value")
Mset.norm.beta.long<-melt(Mset.norm.beta)
colnames(Mset.norm.beta.long)<-c("CpG","sample","value")
gset.funnorm.beta.long<-melt(gset.funnorm.beta)
colnames(gset.funnorm.beta.long)<-c("CpG","sample","value")
##### boxplots ######
p<-ggplot(gset.funnorm.beta.long,aes(x=sample, y=value,fill=sample))       
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")
p.box+labs(y="beta-value")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(legend.position="none")

##### Quality Control ######################################################################################
qcReport(RGset, sampNames = pd$Sample_Name,sampGroups = pd$Sample_Name, pdf = "qcReport.pdf") # macht Probleme wenn er kein Sample_Group hat 

# Density Plot #
densityPlot(RGset, sampGroups = pd$Sample_Group, main = "Beta", xlab = "Beta",legend = TRUE,pal = brewer.pal(3,"Dark2"))
# pal = heat.colors(29,alpha = 1)
# pal = rainbow(29)
# sampGroup weglassen oder = NULL zeigt alle werte
# sampGroup=pd$Sample_Name zeigt nur 8 samples an....wieso?--> weil er als default pal=brewer.pal(8,"Dark2") nimmt und diese palette hat nicht mehr als 8 farben
# andere farbpallette nehmen !Palettes (grDevices)
# legend kann man nicht auf 2 columns ??ndern! also lieber legend =FALSE
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
# default f??r badSampleCutoff = 10.5......aber gerade ist bei ??ber 13??

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

##### Remove Probes with SNP ###########################################################################
snps <- getSnpInfo(gset.funnorm)
head(snps,10)
# Probe = SNPs present inside probe body 
# CpG at CpG interrogation
# SBE at single nucleotide extension
# rs = names of SNPs 
# maf = minor allele frequency of the SNPs
# drop probes with SNP at CpG interrogation or at single nucleotide extension
gset.funnorm <- dropLociWithSnps(gset.funnorm, snps=c("SBE","CpG"), maf=0)
#  drop the probes for any minor allele frequency maf=0 (maf is minor allele frequency)
####################################################################################################

##### Differentially Methylated Positions: dmpFinder ################################################
mset<-Mset.norm[,c(19,20)]        # only for two cellines SJ-AML4 and SJ-AML5
# geht genau so wie mit columns (also mset$Sample_Names anzeigen lassen und dann z??hlen welche column)
# also pehno has to be adapted
pd_2<-pd[c(19,20),]
pd_2$Sample_Group=c("High","Low")
table(pd_2$Sample_Group)
M <- getBeta(mset, type="Illumina")
sum(is.na(beta))
dmp <- dmpFinder(M, pheno=pd_2$Sample_Group, type="categorical")
# groups to be compared are difned in pheno=pd$Sample_Group.....here i want to compare the Names
# kasper Hansen sagt beta values sind besser
# betaThreshold damit keine 0 werte auftauche
#####################################################################################################

##### Differentially Methylated Regions: bumphunter #################################################
pheno <- pData(gset.funnorm[,c(19,20)])$Sample_Name
designMatrix <- model.matrix(~ pheno)
# Run the algorithm with B = 0 permutation on the Beta-values
# with a medium difference cutoff, say 0.2
# (which corresponds to 20% difference on the Beta-values):
dmrs <- bumphunter(], design = designMatrix, cutoff = 0.3, B=0, type="Beta")
# If number of candidate bumps is large, say > 30000, increase the cutoff 
# becuase: the most of the additional candidate regions found by lowering the cutoff 
# will be found nonsignificant after permutation scheme
# therefore time can be saved by being more stringent on the cutoff
# run algorithm with large number of permutations, B = 1000:
dmrs <- bumphunter(gset.funnorm, design = designMatrix, cutoff = 0.3, B=1000, type="Beta")
# Since permutation scheme is expensive, parallel computation is implemented in bumphunter (foreach)
# computing time was too long for my PC evan for cutoff= 0.8
names(dmrs)
head(dmrs$table, n=3)
######################################################################################################

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

# ????? getSex macht bei position 20 M, aber addSexmacht bei 20 F ?????
# ???? How to exclude Sex chromosomes?


source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

detach("package:IlluminaHumanMethylation450kanno.ilmn12.hg19", unload=TRUE)

.default.450k.annotation <- "ilmn12.hg19"
readGEORawFile("GSE35069_Matrix_signal_intensities.txt", sep = "\t", Uname = "Unmethylated Signal", 
               Mname = "Methylated signal", row.names = 1, pData = NULL,
               array = "IlluminaHumanMethylation450k",
               annotation = .default.450k.annotation,
               mergeManifest = FALSE, showProgress = TRUE) 

# change sep to what is used in file
# chnage Uname and Mname to what is used in file
# annotation: the .defualt.450k.annotation will us the default annotation
# default annotation should be IlluminaHumanMethylation450kanno.ilmn12.hg19







