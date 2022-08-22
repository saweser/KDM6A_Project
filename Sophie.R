### empty workspace
rm(list=ls());

#################### 450k Data Sophie ########################################################################
require(minfi)
library(shinyMethyl)
library(RColorBrewer)
library(ggplot2)
setwd("/home/weser/450k_Data/Sophie")
#############################################################################################################

##### Read Data ##########################################################################################
baseDir<-file.path(getwd(),"idat")
targets<-read.450k.sheet(baseDir)
targets$Sample_Group<-c(rep("D2",12),rep("PDX",3),rep("AML Cell Line",14))
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


##### Remove Probes with SNP ###########################################################################
snps <- getSnpInfo(gset.funnorm)
#head(snps,10)
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
mset<-Mset.norm[,c(19,20,21,22)]        # only for two cellines SJ-AML4 and SJ-AML5
# geht genau so wie mit columns (also mset$Sample_Names anzeigen lassen und dann zählen welche column)
# also pehno has to be adapted
pd_2<-pd[c(19,20,21,22),]
pd_2$Sample_Group=c("High","Low")
table(pd_2$Sample_Group)
M <- getBeta(mset, type="Illumina")
#sum(is.na(beta))
dmp <- dmpFinder(M, pheno=pd_2$Sample_Group, type="categorical")
# groups to be compared are difned in pheno=pd$Sample_Group.....here i want to compare the Names
# kasper Hansen sagt beta values sind besser
# betaThreshold damit keine 0 werte auftauche
#####################################################################################################
aa<-5
##### Differentially Methylated Regions: bumphunter #################################################
pheno <- pData(gset.funnorm[,c(19,20)])$Sample_Name
designMatrix <- model.matrix(~ pheno)
# Run the algorithm with B = 0 permutation on the Beta-values
# with a medium difference cutoff, say 0.2
# (which corresponds to 20% difference on the Beta-values):
dmrs0 <- bumphunter(gset.funnorm[,c(19,20)], design = designMatrix, cutoff = 0.3, B=0, type="Beta")
# If number of candidate bumps is large, say > 30000, increase the cutoff 
# becuase: the most of the additional candidate regions found by lowering the cutoff 
# will be found nonsignificant after permutation scheme
# therefore time can be saved by being more stringent on the cutoff
# run algorithm with large number of permutations, B = 1000:
dmrs <- bumphunter(gset.funnorm[,c(19,20)], design = designMatrix, cutoff = 0.3, B=1000, type="Beta")
# Since permutation scheme is expensive, parallel computation is implemented in bumphunter (foreach)
# computing time was too long for my PC evan for cutoff= 0.8

#names(dmrs)
#head(dmrs$table, n=3)

######################################################################################################
## save full workspace
save(list=ls(all.names=TRUE), file="Sophie.RData");


# ????? getSex macht bei position 20 M, aber addSexmacht bei 20 F ?????
# ???? How to exclude Sex chromosomes?