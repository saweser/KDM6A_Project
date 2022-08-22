
# I will provide you with example code to conduct posteriori power analysis. Please change the parameters to your own needs, especially the log fold changes, percentage DE, sample sizes.

# Settings
rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
library(powsimRDev)
library(tidyverse)

# Input
# The count table and annotation table.
load("counts.RData")


annot <- anno %>% tibble::rownames_to_column() %>% dplyr::filter(AraC=="No")
cnts <- counts[,annot$rowname]

# Estimation
# We will estimate the negative binomial distribution parameters for bulk RNAseq data, using TMM as normalisation.
estparam <- estimateParam(countData = cnts, 
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

plotParam(estParam.out = estparam, annot=T)


# Simulation
# The differential expression setup:
lfc.foo <- function(x) rnorm(x, sd = 1)
desetup <- DESetup(ngenes = 10000, 
                       nsims = 25, 
                       p.DE = 0.1, 
                       p.B = NULL, 
                       pLFC = lfc.foo, 
                       bLFC = NULL, 
                       sim.seed = 43872)

# Combine the DE settings with the estimated parameters:
simsetup <- SimSetup(desetup = desetup, 
                         params = estparam, 
                         spike = NULL, 
                         size.factors = 'equal',
                         downsample = FALSE)

#Running the DE simulations with limma-voom and TMM normalisation.
simres <- simulateDE(n1 = c(3,6,9,15,20), n2=c(3,6,9,15,20), 
                         sim.settings = simsetup, 
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

# save.image(file="powsimR_HEK.RData")
```

# Evaluation
# Calculating the error rates stratified by log2 fold changes. There are also other options available, e.g. stratify by mean expression or filter out lowly expressed genes.
evalDE <- evaluateDE(simRes = simres, alpha.type = "adjusted", MTC='BH', alpha.nominal = 0.05,
                         stratify.by = 'mean', filter.by = 'none', strata.filtered = 1, target.by = 'lfc', delta=0)

plotEvalRes(evalRes = evalDE, rate = 'marginal', quick = T, annot = T)

plotEvalRes(evalRes = evalDE, rate = 'stratified', quick = T, annot = T)


