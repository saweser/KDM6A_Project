########## run bumphunter ######################################################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/MM1_vs_MM6")
library(doParallel)
library(minfi)

load('Preprocessed_Data_MM1_MM6.RData');
################################################################################################################

##### set up design matrix #################################################
# bumphunter can only hunt bumps for ONE coefficient
# set the coeffcienct in the design matrix with coef=2
# make design with remission in Intercept

design<- stats::model.matrix(~Sample_Type, data=targets)
#################################################################################################

##### run bumphunter ############################################################################
registerDoParallel(cores = 4)
dmrs <- bumphunter(Gset.quantile.keep, design = design, cutoff = 0.08, B=1000, 
                     type="Beta",coef=2, verbose=TRUE)

# save full workspace
save(list=ls(all.names=TRUE), file="bumphunter_MM1_MM6_test.RData");

