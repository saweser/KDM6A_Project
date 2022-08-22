
# Settings
rm(list=ls());  # empty workspace
library(powsimRDev)
library(tidyverse)

setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("counts.RData")
##################################################################################################################

# Input
# The count table and annotation table.

cnts <- counts_mus
annot <- read.table(file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_mus.txt", header=T, row.names = 1, stringsAsFactors = F)

# use the 2A3 cell line with WT as baseline.
annot.2A3 <- annot %>% tibble::rownames_to_column() %>% dplyr::filter(cell_line=="2A3", state=="WT")
cnts.2A3 <- cnts[,annot.2A3$rowname]

# Estimation

# We will estimate the negative binomial distribution parameters for bulk RNAseq data, using TMM as normalisation.

estparam.2A3 <- estimateParam(countData = cnts.2A3, 
                              batchData = NULL, 
                              spikeData = NULL, 
                              spikeInfo = NULL, 
                              Lengths = NULL, 
                              MeanFragLengths = NULL, 
                              Distribution = 'NB',
                              RNAseq = 'bulk', 
                              normalisation = 'TMM', 
                              sigma = 1.96, 
                              NCores = NULL, verbose = TRUE)

plotParam(estParam.out = estparam.2A3, annot=T)


# Simulation
# The differential expression setup:
  
lfc.foo <- function(x) rnorm(x, sd = 1)
desetup.2A3 <- DESetup(ngenes = 10000, 
                       nsims = 25, 
                       p.DE = 0.1, 
                       p.B = NULL, 
                       pLFC = lfc.foo, 
                       bLFC = NULL, 
                       sim.seed = 43872)

# Combine the DE settings with the estimated parameters:
  
simsetup.2A3 <- SimSetup(desetup = desetup.2A3, 
                         params = estparam.2A3, 
                         spike = NULL, 
                         size.factors = 'equal',
                         downsample = FALSE)

# Running the DE simulations with limma-voom and TMM normalisation.


simres.2A3 <- simulateDE(n1 = c(3,6,9,15,20), n2=c(3,6,9,15,20), 
                         sim.settings = simsetup.2A3, 
                         DEmethod = "limma-voom", 
                         normalisation = "TMM", 
                         Preclust = FALSE, 
                         Preprocess = NULL, 
                         GeneSelect = NULL,
                         DimReduce = NULL, 
                         ClustMethod = NULL, 
                         Pipeline = NULL, 
                         spikeIns = FALSE, 
                         NCores = 1, 
                         verbose = TRUE)

save.image(file="powsimR_2A3.RData")

# Evaluation

# Calculating the error rates stratified by log2 fold changes. There are also other options available, e.g. stratify by mean expression or filter out lowly expressed genes.

evalDE.2A3 <- evaluateDE(simRes = simres.2A3, alpha.type = "adjusted", MTC='BH', alpha.nominal = 0.05,
                         stratify.by = 'mean', filter.by = 'none', strata.filtered = 1, target.by = 'lfc', delta=0)

plotEvalRes(evalRes = evalDE.2A3, rate = 'marginal', quick = T, annot = T)

plotEvalRes(evalRes = evalDE.2A3, rate = 'stratified', quick = T, annot = T)



### check KDM6A
kdm6a<-cnts.2A3[rownames(cnts.2A3)=="ENSMUSG00000037369",]
kdm6a<-as.numeric((kdm6a))
mean(kdm6a)

# knocdown
a<-cnts[rownames(cnts)=="ENSMUSG00000037369",]
a<-a[,grepl("2A3A9",colnames(a))]
a<-as.numeric(a)
mean(a)
