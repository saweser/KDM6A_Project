##### Probe-wise differential methylation analysis #######################################################

rm(list=ls());  # empty workspace
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
library(limma)
setwd("/data/htp/A07/MM1_vs_MM6")

# load workspace
load('Preprocessed_Data_MM1_MM6.RData');
source("/data/htp/A07/AML_triplets/plotCpg_Bria.R")
#########################################################################################################

##### calculate M and beta values ########################################################################
Mvalues <- getM(Gset.quantile.keep)
bvalues <- getBeta(Gset.quantile.keep)
##########################################################################################################

##### Design and Contrast Matrix ########################################################################
design <- model.matrix(~Sample_Type, targets) # rownames=samples,prints "1" for every "type" that is true for that sample
# model.matrix braucht nur den vector hier Sample_Type, hier nimmt er in aus targets, man kann ihm auch nur den vector geben.
colnames(design) <- c("Intercept","MM6")
# rownames(design) <- paste($condition, inf$Genotype, sep = "_") this is from beate, but makes not much sense here because do not have second covariate

fit <- lmFit(Mvalues, design) # fit linear model
fit2 <- eBayes(fit) # apply empirical Bayes smoothing to the standard errors.

# One can check whether the variance is now constant over mean expression levels
# by plotting the models residual variances (sigma) against average expression values with plotSA()

pdf("plots/SA_plot_MM1_MM6.pdf", width=4, height=3)
plotSA(fit2, xlab="Average log-expression", ylab="log2(sigma)",
       zero.weights=FALSE, pch=16, cex=0.2)
dev.off()
####################################################################################################

###### DMP tables ################################################################################################
# ordered by B-statistic by default, using topTable
# order by p-value: sort.by="p"; (ordering to B-statistic should be identical)
ann450kSub <- ann450k[match(rownames(Mvalues),ann450k$Name),c(1:4,12:19,24:ncol(ann450k))]
# make sure rownames of Mvalues and annotation are the same!
# take out columns of annotation that are not needed c(1:4,12:19,24:ncol(ann450k)), these are the columns from @listData
DMPs <- topTable(fit2, number =Inf, coef=2, confint = TRUE, adjust.method = "BH", genelist=ann450kSub)
DMPs[DMPs==""]  <- NA 

# separate the column UCSC_RefGene_Name so that it just shows one gene name only
DMPs<-DMPs %>%
  separate(UCSC_RefGene_Name, c("UCSC_RefGene_Name","UCSC_RefGene_Name_2"), ";",extra = "merge", fill = "right")

##### Save DMPs full table .csv #####################################################################################
write.table(DMPs, file="/data/htp/A07/MM1_vs_MM6/tables/DMPs_MM1_MM6_full.csv", sep=",", row.names=FALSE)
###################################################################################################################

# only significant ones
DMPs<-DMPs[DMPs$adj.P.Val<0.05,]
###################################################################################################################

##### How many Positions (CpGs) are differentially methylated ? ##################################################
# how many positions
length(rownames(DMPs))
# how many up regulated?
sum(DMPs$logFC>0)
# how many down regulated?
sum(DMPs$logFC<0)
###################################################################################################################

##### Save DMP tables as .csv #####################################################################################
write.table(DMPs, file="/data/htp/A07/MM1_vs_MM6/tables/DMPs_MM1_MM6.csv", sep=",", row.names=FALSE)
###################################################################################################################

##### Islands, Shelfs, Shores #####################################################################################
# plot number DMPS on Ilands, Shelfs, Shores
p<-ggplot(DMPs,aes(x=Relation_to_Island, fill=Relation_to_Island))
p.hist<-p+geom_histogram(stat = "count")
p.hist+
  labs(title="Differentially methylated positions")+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())
  
# plot all the cpgs on the chip
AllProbes<-read.csv("/data/htp/A07/AML_triplets/HumanMethylation450_15017482_v1-2.csv", skip = 7, header=T, strip.white = T, as.is = T)
AllProbes<-AllProbes[1:485577,]

p<-ggplot(AllProbes,aes(x=Relation_to_UCSC_CpG_Island, fill=Relation_to_UCSC_CpG_Island))
p.hist<-p+geom_histogram(stat = "count")
p.hist+
  labs(title="All probes")+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())

DMPs_plot<-as.data.frame(dplyr::count(DMPs,Relation_to_Island))
All<- as.data.frame(dplyr::count(AllProbes,Relation_to_UCSC_CpG_Island))

All$Relation_to_UCSC_CpG_Island[1]<-"OpenSea"
All<-All[order(All$Relation_to_UCSC_CpG_Island),]
DMPs_plot<-DMPs_plot[order(DMPs_plot$Relation_to_Island),]
colnames(All)<-c("Relation_to_Island","n_probes")

DMPs_plot <- merge(DMPs_plot, All,by="Relation_to_Island")
DMPs_plot$percentage<-(DMPs_plot$n/DMPs_plot$n_probes)

# plot the percentage
p<-ggplot(DMPs_plot,aes(x=Relation_to_Island, y=percentage, fill=Relation_to_Island))
p.bar<-p+geom_bar(stat = "identity")+
  labs(title="All probes")+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())
p.bar

ggsave("Percentage_Cpgs_MM1_MM6.eps",plot = p.bar, width= 105, height= 75, units= "mm", device = "eps",path = "/data/htp/A07/MM1_vs_MM6/plots")
###################################################################################################################

##### plot top 4 DMPs ############################################################################################
### Comparison 1: Diagnosis-Relapse
pdf("plots/Top4_DMPs_MMM1_MM6.pdf", width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
plotCpg_Bria(bvalues, cpg=rownames(DMPs)[1:4], pheno=targets$Sample_Type, ylab = "Beta values",mainPrefix = DMPs$UCSC_RefGene_Name[1:4],mainSuffix = DMPs$chr[1:4])
title(main="MM1 vs MM6", outer=T)
dev.off()
#########################################################################################################

# make own plotCpg with ggplot so can add mean, or median?

# dataframe consisting of bvalues and targets
long <- melt(bvalues)
colnames(long)<-c("probe","Slide_Array","value")
df<-left_join(targets,long,by="Slide_Array")

# plot the first 4 DMPs
p1<-ggplot(df[df$probe==rownames(DMPs)[1],], aes(x=Sample_Type, y=value, colour=Sample_Type))
p1<-p1+geom_point(position="jitter")+
labs(title=c(DMPs$UCSC_RefGene_Name[1]))+
theme(legend.position="none")+
theme(axis.title.x = element_blank())+
theme(axis.ticks.x=element_blank())+
theme(panel.grid.minor = element_blank()) +                                 
stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")

p2<-ggplot(df[df$probe==rownames(DMPs)[2],], aes(x=Sample_Type, y=value, colour=Sample_Type))
p2<-p2+geom_point(position="jitter")+
  labs(title=c(DMPs$UCSC_RefGene_Name[2]))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(panel.grid.minor = element_blank()) +                                  
  stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")

p3<-ggplot(df[df$probe==rownames(DMPs)[3],], aes(x=Sample_Type, y=value, colour=Sample_Type))
p3<-p3+geom_point(position="jitter")+
  labs(title=c(DMPs$UCSC_RefGene_Name[3]))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(panel.grid.minor = element_blank()) +                                  
  stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")

p4<-ggplot(df[df$probe==rownames(DMPs)[4],], aes(x=Sample_Type, y=value, colour=Sample_Type))
p4<-p4+geom_point(position="jitter")+
  labs(title=c(DMPs$UCSC_RefGene_Name[4]))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(panel.grid.minor = element_blank()) +                                  
  stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")
  
grid.arrange(p1,p2,p3,p4, ncol=2)
ggsave("DMPs_ggplot_MM1_MM6.eps",plot = grid.arrange(p1,p2,p3,p4, ncol=2), width= 105, height= 75, units= "mm", device = "eps",path = "/data/htp/A07/MM1_vs_MM6/plots")
#####################################################################################################

# try adding the mean as a line not as a point
geom_line(stat = "hline", yintercept = "mean" )
stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "crossbar")
geom_errorbar(stat = "hline", yintercept = "mean",
              width=0.8,aes(ymax=..y..,ymin=..y..))
geom_errorbar(stat = "identity", yintercept = "mean", width=0.6)
stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point")
##############################################################################









but you can use dplyr syntax i.e. df %>% left_join(a,df,by="probeid)


you can use directly like that in ggplot data

dplyr is a pipe : pipe a into b a %>% b
if a is a data frame it can go in the dplyr function: 
like count(df,column1) or df %>% count(column1)


the first one actually makes more sense when you do long%>%left_join(targets,.,by="Slide_Array")
. means current input

blah<-left_join(targets,long,by="Slide_Array")
 blah2<-long%>%left_join(targets,long,by="Slide_Array")
                                                  





