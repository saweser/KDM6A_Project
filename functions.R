
########## Functions #####################################################################################

### calculate UPM from UMI counts #######################################################################
calc_UPM <-
  function (mat) 
  {
    sums <- colSums(mat)
    upm<-t(t(mat)/sums* 1e+06) # transform, because R will calculate row wise
    return(upm)
  }
#########################################################################################################

#### PCA function ######################################################################################
pcaFunction<-function(mat, inf, genes, ngenes, col="treatment", label="sid",shape ) { 
  
  if(!missing(genes)){
    vsdMat <- mat[genes,rownames(inf)]
  }else{
    vsdMat<- mat[,rownames(inf)]
  }
  rowVar<-apply(mat,1,var)
  
  if(!missing(ngenes)){
    mv500<-order(rowVar,decreasing = T)[1:ngenes]
  }else{
    mv500<-(rowVar>0)
  }
  if(missing(shape)){
    shape=col
  }
  
  pc<-prcomp(t(vsdMat[mv500,]),scale=T)
  pc.sum<-summary(pc)$importance
  varExp<-round(pc.sum[2,]*100,2)
  pcs<-data.frame(pc$x, inf)
  pcs$sid<-rownames(pcs)
  pp<-ggplot(pcs,aes_q(x=quote(PC1),y=quote(PC2),col=as.name(col)))
  if(!missing(label)){
    pp<-pp+geom_text(aes(label=as.name(label)),size=3)
  }else{
    if(missing(shape)){
      pp<-pp+geom_point(aes(shape=as.name(shape)),alpha=0.5,size=3)
    }else{
      
      pp<-pp+geom_point(alpha=0.5,size=3)
    }
  }
  pp+xlab(paste("PC1 (",varExp[1],"%)"))+
    ylab(paste("PC2 (",varExp[2],"%)"))
}
################################################################################################################

#### Boxplot for normalization #################################################################################
normBoxplot<-function(mm,title){
  df<- data.frame(mm) %>% gather(key=sample_name,value="raw",1:dim(mm)[2])
  box<-ggplot(df,aes(x=sample_name,y=raw))+
    geom_boxplot()+ylab(expression(log[2](counts)))+ 
    coord_flip()+
    theme_grey()+
    theme(axis.title.y=element_blank())
  if(!missing(title)){
    box+ggtitle(title)
  }else{
    box
  }
}
###############################################################################################################

#### get geneIDs from the ENSEMBL IDs #############################################################################
# use biomaRT::getBM
# insert die gene IDs in the DE data frame
# rownames need to be the ensembl IDs

# human
getGeneID_hsa<-function(mat){
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  hgnc <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),filters = 'ensembl_gene_id', values = rownames(mat), mart = ensembl)
  mat$ensembl_gene_id<-rownames(mat)
  mat<-left_join(mat,hgnc,by="ensembl_gene_id")
  rownames(mat)<-mat$ensembl_gene_id
  mat$ensembl_gene_id<-NULL
  refcols<-"external_gene_name"
  mat<-mat[,c(refcols,dplyr::setdiff(names(mat),refcols))]
  
}
# mouse
getGeneID_mus<-function(mat){
  ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  hgnc <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),filters = 'ensembl_gene_id', values = rownames(mat), mart = ensembl)
  mat$ensembl_gene_id<-rownames(mat)
  mat<-left_join(mat,hgnc,by="ensembl_gene_id")
  rownames(mat)<-mat$ensembl_gene_id
  mat$ensembl_gene_id<-NULL
  refcols<-"external_gene_name"
  mat<-mat[,c(refcols,dplyr::setdiff(names(mat),refcols))]
  
}
# mouse does not have hgnc_symbol, use external_gene_name isnted , for human its the same
####################################################################################################################

