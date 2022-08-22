annotation<-getAnnotation(gset.funnorm)
Rs<-annotation$Probe_rs
head(Rs)

sum(is.na(Rs))
length(Rs)

same for Mset.norm

green<-getGreen(RGset)
red<-getRed(RGset)
beta_RG<-getBeta(RGset)
beta_Meth<-getBeta(Mset.norm)
View(beta_Meth)
View(beta_RG)

sind rs probes in RGset, Mset und gset??

grep_beta<-beta_RG[grep("^rg",row.names(beta_RG)),]
grep_Meth<-beta_Meth[grep("^rg",row.names(beta_Meth)),]          
View(grep_beta)
View(grep_Meth)

grep_green<-green[grep("^rg",row.names(green)),]
View(grep_green)

SNP_RG<-getSnpBeta(RGset)
SNP_Meth<-getSnpBeta(Mset.norm)

snps <- getSnpInfo(gset.funnorm)
head(snps,10)
