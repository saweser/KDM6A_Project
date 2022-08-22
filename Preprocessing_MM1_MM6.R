########## AML_triplets_Preprocessing #############################################################
# Preprocessing of raw idat files
# check for normalization methods
# quality control plots

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/MM1_vs_MM6")
###########################################################################################################
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(ggplot2)
library(shinyMethyl)
library(reshape)
library(gridExtra)
library(colorspace)

source("/data/htp/A07/AML_triplets/pcaFunction.R")
############################################################################################################

##### Read Data ###########################################################################################
baseDir<-file.path(getwd(),"/idat")
targets<-read.metharray.sheet(baseDir) 
# reorder columns
targets$Sample_Well<-NULL
targets$Sample_Plate<-NULL
targets$Sample_Group<-NULL
targets$Pool_ID<-NULL

# add column for type (MM1, MM6)
targets$Sample_Type<-""
for (i in 1:nrow(targets)){                            
  if (grepl("MM_1",targets[i, "Sample_Name"])){           
    targets[i, "Sample_Type"]<-"MM1"                
  } else if( grepl("MM_6",targets[i, "Sample_Name"])){    
    targets[i, "Sample_Type"]<-"MM6"            
  }
}
targets<-targets[,c(1,5,2,3,4)]

RGset <- read.metharray.exp(targets = targets)
pd <- pData(RGset)
###########################################################################################################

##### Quality Control #####################################################################################
### calculate p-values ###
# generate detection p-value for every CpG in every sample
# minfi: compare total signal (M + U) for each probe to background signal(from negative control probes)
# Very small p-values: reliable signal 
# large p-values (>0.01): generally poor quality signal
detP <- detectionP(RGset)
# plot mean detection p-value for each sample
detP_mean<-data.frame(colMeans(detP))
detP_mean$sample<-rownames(detP_mean)
detP_mean$Sample_Type<-targets$Sample_Type

p<-ggplot(detP_mean,aes(x=sample,y=colMeans.detP.,fill=Sample_Type))
p.bar<-p+geom_bar(stat="identity")
p.bar+ labs(y="Detection p-value")+
  theme(legend.position="right")+
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x=element_blank())+
  scale_fill_manual(values = c("skyblue3","plum3"))

ggsave("Detection-p_MM1_MM6.eps",plot = p.bar, width= 105, height= 75, units= "mm", device = "eps",path = "/data/htp/A07/MM1_vs_MM6/plots")

### minfi qc report 
qcReport(RGset, sampNames=targets$Sample_Name, sampGroups=targets$Sample_Type, pdf="/data/htp/A07/MM1_vs_MM6/plots/qcReport.pdf")

### Density Plot 
pdf("plots/Density_MM1_MM6.pdf", width=8, height=6)
densityPlot(RGset, sampGroups = targets$Sample_Type, xlab = "Beta",legend = FALSE,pal = c("skyblue3","plum3"))
legend("topright", legend = levels(factor(targets$Sample_Type)),text.col=c("skyblue3","plum3"),bty = "n")
dev.off()
# producing a density plot that looks like this is not possible with ggplot
# can calculate density (density(x$y) )and then use geom_line 
# for ggplot would have to calculate beta before (getBeta) that function rounds the values 
# everything below 0=0 and everything above 1 will be 1
# so range 0-1
# the densityPlot dunction gets the RGset! it will not round the values like that so it plots also values above and below that range (0-1)

### Density Bean Plot
pdf("plots/Density_Beanplot_MM1_MM6.pdf", width=8, height=6)
par(oma=c(0,3,0,0))
densityBeanPlot(RGset, sampGroups = targets$Sample_Name, sampNames = pd$Sample_Name, pal = colorspace::rainbow_hcl(12))
dev.off()

### shinyMethyl
summary <- shinySummarize(RGset)
runShinyMethyl(summary)
##########################################################################################################

##### remove samples with high detection p-value ########################################################
##########################################################################################################

##### Normalization ######################################################################################
Mset.raw <-preprocessRaw(RGset)
Mset.noob<- preprocessNoob(RGset)
Mset.swan <-preprocessSWAN(RGset)
Mset.illumina <-preprocessIllumina(RGset,bg.correct = TRUE, normalize="controls")
Gset.quantile <- preprocessQuantile(RGset)
Gset.funnorm <-preprocessFunnorm(RGset)
##########################################################################################################

##### Plot Normalizations ##############################################################################
# get beta values
Mset.raw.beta<-getBeta(Mset.raw)
Mset.noob.beta<-getBeta(Mset.noob)
Mset.swan.beta<-getBeta(Mset.swan)
Mset.illumina.beta<-getBeta(Mset.illumina)
Gset.quantile.beta<-getBeta(Gset.quantile)
Gset.funnorm.beta<-getBeta(Gset.funnorm)
# get M values
Mset.raw.m<-getM(Mset.raw)
Mset.noob.m<-getM(Mset.noob)
Mset.swan.m<-getM(Mset.swan)
Mset.illumina.m<-getM(Mset.illumina)
Gset.quantile.m<-getM(Gset.quantile)
Gset.funnorm.m<-getM(Gset.funnorm)

#### convert data frame to long format ######
Mset.raw.beta.long<-melt(Mset.raw.beta)
Mset.noob.beta.long<-melt(Mset.noob.beta)
Mset.swan.beta.long<-melt(Mset.swan.beta)
Mset.illumina.beta.long<-melt(Mset.illumina.beta)
Gset.quantile.beta.long<-melt(Gset.quantile.beta)
Gset.funnorm.beta.long<-melt(Gset.funnorm.beta)
Mset.raw.m.long<-melt(Mset.raw.m)
Mset.noob.m.long<-melt(Mset.noob.m)
Mset.swan.m.long<-melt(Mset.swan.m)
Mset.illumina.m.long<-melt(Mset.illumina.m)
Gset.quantile.m.long<-melt(Gset.quantile.m)
Gset.funnorm.m.long<-melt(Gset.funnorm.m)

colnames(Mset.raw.beta.long)<-c("CpG","sample","value")
colnames(Mset.noob.beta.long)<-c("CpG","sample","value")
colnames(Mset.swan.beta.long)<-c("CpG","sample","value")
colnames(Mset.illumina.beta.long)<-c("CpG","sample","value")
colnames(Gset.quantile.beta.long)<-c("CpG","sample","value")
colnames(Gset.funnorm.beta.long)<-c("CpG","sample","value")
colnames(Mset.raw.m.long)<-c("CpG","sample","value")
colnames(Mset.noob.m.long)<-c("CpG","sample","value")
colnames(Mset.swan.m.long)<-c("CpG","sample","value")
colnames(Mset.illumina.m.long)<-c("CpG","sample","value")
colnames(Gset.quantile.m.long)<-c("CpG","sample","value")
colnames(Gset.funnorm.m.long)<-c("CpG","sample","value")

##### boxplots beta-values ######
p1<-ggplot(Mset.raw.beta.long,aes(x=sample, y=value))       
p1<-p1+geom_boxplot(stat="boxplot",position="dodge",fill="plum")+
  labs(y="beta-value",title="preprocessRaw")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

p2<-ggplot(Mset.noob.beta.long,aes(x=sample, y=value))       
p2<-p2+geom_boxplot(stat="boxplot",position="dodge",fill="yellow")+
  labs(y="beta-value",title="preprocessNoob")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

p3<-ggplot(Mset.swan.beta.long,aes(x=sample, y=value))       
p3<-p3+geom_boxplot(stat="boxplot",position="dodge",fill="lightblue")+
  labs(y="beta-value",title="preprocessSwan")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

p4<-ggplot(Mset.illumina.beta.long,aes(x=sample, y=value))       
p4<-p4+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise")+
  labs(y="beta-value",title="preprocessIllumina")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

p5<-ggplot(Gset.quantile.beta.long,aes(x=sample, y=value))       
p5<-p5+geom_boxplot(stat="boxplot",position="dodge",fill="pink")+
  labs(y="beta-value",title="preprocessQuantile")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

p6<-ggplot(Gset.funnorm.beta.long,aes(x=sample, y=value))       
p6<-p6+geom_boxplot(stat="boxplot",position="dodge",fill="indianred")+
  labs(y="beta-value",title="preprocessFunnorm")+
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=10))+
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title.y = element_text(size=12,color="black"))+ 
  theme(panel.grid.major.x=element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position="none")

grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2)
ggsave("Normalization_beta_MM1_MM6.eps",plot = grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2),device = "eps",path = "/data/htp/A07/MM1_vs_MM6/plots")
#######################################################################################################

##### visualise what the data looks like before and after normalisation ###############################
pdf("plots/Raw_Normalized_MM1_MM6.pdf", width=8, height=4)
par(mfrow=c(1,2))
densityPlot(RGset, sampGroups=pd$Sample_Type, main="Raw", legend=FALSE, pal = brewer.pal(3,"Set1"))
legend("top", legend = levels(factor(targets$Sample_Type)),text.col=brewer.pal(3,"Set1"),bty = "n")
densityPlot(getBeta(Gset.quantile), sampGroups=targets$Sample_Type,main="Normalized", legend=FALSE, pal = brewer.pal(3,"Set1"))
legend("top", legend = levels(factor(targets$Sample_Type)),text.col=brewer.pal(3,"Set1"),bty = "n")
dev.off()
# legend made separately because can then set where to put it ("top": top middle)
# but be careful with the colors, need to be the same in plot and legend!
##########################################################################################################

##### Data exploration ###################################################################################
##### PCA: before filtering ############################################################################################

Gset.quantile.m<-getM(Gset.quantile)
# add a column with slide and array to compare the Mvalues data frame to, to make sure the order is the same
targets<-targets %>%
  unite(Slide_Array,Slide,Array,sep="_",remove=FALSE)
# check if samples are in same order
summary(colnames(Gset.quantile.m)==targets$Slide_Array)
colnames(Gset.quantile.m) <- targets$Sample_Name
rownames(targets)<-targets$Sample_Name

PCA <- pcaFunction(mat = Gset.quantile.m, inf = targets, ngenes = 500, col = "Sample_Type") + 
  ggtitle("Mvalues")+
  theme(legend.title = element_blank())

ggsave("PCA_before_filtering_MM1_MM6.pdf", plot = PCA, width= 150, height= 100, units= "mm",device = "pdf",path = "/data/htp/A07/MM1_vs_MM6/plots")
#########################################################################################################################

##### get the 450k annotation data #########################################################################
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
############################################################################################################

##### Filtering out poor performing probes ####################################################################################
# filter out probes that have failed in one or more samples based on detection p-value
# ensure probes are in the same order 
detP <- detP[match(featureNames(Gset.quantile),rownames(detP)),]
# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(Gset.quantile)
# detP<0.01 will return logical for each value
# rowSums will then count how many values per row are true
# and that number needs to be ncol(mSetSq), so 10 in this case
# means it will exclude probes that have a detP value <0.01 in at least one sample
table(keep) 
Gset.quantile.keep <- Gset.quantile[keep,]
##########################################################################################################

### filter out probes on sex chromosomes #################################################################
# keep2 <- !(featureNames(Gset.quantile.keep) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep2) # filtered out 11261 probes
# Gset.quantile.keep <- Gset.quantile.keep[keep2,]

# leave in because is cell lines?
###########################################################################################################

#### leave in SNPs !!!!!! ######################################################################################
# mSetSqFlt <- dropLociWithSnps(mSetSqFlt) # function to drop SNPs
#########################################################################################################

### filter out cross reactive probes ####################################################################
# list: Chen(2013) Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray
# can be obtained from website www.sickkids.cs --> research--> Weksberg Lab--> publications
xReactiveProbes <- read.csv("48639-non-specific-probes-Illumina450k.csv" , stringsAsFactors=FALSE)
keep3 <- !(featureNames(Gset.quantile.keep) %in% xReactiveProbes$TargetID)
table(keep3)
Gset.quantile.keep <- Gset.quantile.keep[keep3,]
##########################################################################################################

##### PCA after filtering ################################################################################
Mvalues.keep<-getM(Gset.quantile.keep)

# check if samples are in same order
summary(colnames(Mvalues.keep)==targets$Slide_Array)
colnames(Mvalues.keep) <- targets$Sample_Name

PCA <- pcaFunction(mat = Mvalues.keep, inf = targets, ngenes = 500, col = "Sample_Type") + 
  ggtitle("Mvalues")+
  theme(legend.title = element_blank())

ggsave("PCA_after_filtering_MM1_MM6.pdf", plot = PCA, width= 150, height= 100, units= "mm",device = "pdf",path = "/data/htp/A07/MM1_vs_MM6/plots")
#############################################################################################################

##### save objects for downstream analysis ###############################################################
save(list=c("detP","pd","targets","ann450k","Gset.quantile.keep"),
     file="Preprocessed_Data_MM1_MM6.RData");
##########################################################################################################








