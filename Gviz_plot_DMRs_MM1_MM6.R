########## Plot the DMRs with the Gviz package ####################################################################################

# introduction from Kaspar Daniel HAnsen (https://github.com/hansenlab/tutorial.450k/blob/master/vignettes/methylation450k.Rmd)
# input is DMR table from bumphunter!
# tracks and tables from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
################################################################################################################################

rm(list=ls());  # empty workspace

library(Gviz)
setwd("/data/htp/A07/MM1_vs_MM6")
load("bumphunter_MM1_MM6.RData")
################################################################################################################################

#### subset only significant DMRs #########################################################################################
dmrs_table<-subset(dmrs$table,fwer<=0.05)
#############################################################################################################################

#### get gene names ################################################################################################################
dmrs_table$chr<-gsub(dmrs_table$chr,pattern = "chr",replacement = "")

mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh=37)
results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                 filters = c("chromosome_name", "start", "end"), values = list(dmrs_table[,1], dmrs_table[,2], 
                                                                               dmrs_table[,3]), mart = mart)


ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)


# useMarts does not take the version, useEnsembl does take GRCh=37!


##### Specifiy which DMR to plot ##############################################################################################
dmr <- dmrs_table[1,] # which DMR table and which DMR to use
################################################################################################################################

##### General Stuff #############################################################################################################
genome <- "hg19"
chrom <- dmr$chr
start <- dmr$start
end <- dmr$end
minbase <- start - 0.25 * (end - start)
maxbase <- end   + 0.25 * (end - start)
pal <- c("plum", "midnightblue","seagreen1")
################################################################################################################################

##### Building the tracks ########################################################################################################
### Chromosome Ideogram
# red lines is current location
# centromer is red shape
iTrack <- IdeogramTrack(genome = genome, 
                        chromosome = dmr$chr, name = "")
### Genome Axis Track: 
# is relative to the other tracks, adds a genomic axis
gTrack <- GenomeAxisTrack(col = "black", cex = 1, 
                          name = "", fontcolor = "black")

### Annotation Track
# showing DMR position
dmrTrack <- AnnotationTrack(start = start, end = end, genome = genome, chromosom = chrom, 
                            name = "DMR")

### Biomart GENE REGION TRACK
biomTrack <- BiomartGeneRegionTrack(genome = "hg19",chromosome = chrom, start = minbase, end = maxbase,
                                    name = "ENSEMBL", transcriptAnnotation = "name2",
                                    collapseTranscripts = FALSE)

### methylation data track
gr <- granges(Gset.quantile.keep, use.names = FALSE) # makes GRanges object
beta<-as.data.frame(getBeta(Gset.quantile.keep))     # takes beta values out, will be plotted in the methTrack 
values(gr) <- cbind(values(gr), beta)                # add beta values to the metadata of the GRanges object

methTrack <- DataTrack(range = gr,
                       groups = targets$Sample_Type,
                       genome = genome,
                       chromosome = chrom, 
                       ylim = c(-0.05, 1.05),             # limits of the y axis
                       col = pal,                         # colors to be used
                       type = "boxplot",                  # plot type
                       do.out = FALSE,                    # remove outliers
                       box.width=8,                       # size of the box
                       name = "DNA Meth.\n(beta value)",  # name of the track
                       background.panel = "white",        # color of background
                       legend = TRUE,                     # draw a legend for the beta values
                       cex.title = 0.8,                   # text of the track title
                       cex.axis = 0.8,                    # size of the axis text 
                       cex.legend = 0.8,                  # size of the legend text
                       grid = FALSE)                       # gridlines in methTrack


# boxplots does not work when cpgs are too far part! with and without groups its not working
# box.width makes boxes bigger! then it works!!!!!
# still there is that dotted line now
# frame around box plot looks like it is the size of the IRanges.....
# in our case when it is just one base...this will be very small
# why only one base per cpG

#### Plot the tracks #################################################################################
# plot with methylation track
plotTracks(list(iTrack, gTrack, methTrack, dmrTrack, biomTrack), 
           from = minbase, 
           to = maxbase, 
           add53 = TRUE,            # makes 3`and 5' marks on the genome Axis
           add35 = TRUE,            # makes 3`and 5`marks on the genome axis
           sizes = c(1,1,4,1,1.5))   # relative vertical sizes for the tracks

# plot without methylation track
plotTracks(list(iTrack,gTrack,dmrTrack,biomTrack),
           add53 = TRUE,
           add35 = TRUE,
           sizes = c(1,1,1.5,2))

# because methylation track needs a lot of space scale will be diffrent 
# so to see the gen names one has to leave methylation track out ?
######################################################################################################

## Annotation of the DMR using bumphunter::matchGenes

islands <-dmrs_table

library(bumphunter)

library("TxDb.Hsapiens.UCSC.hg19.knownGene")

genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tab<- matchGenes(islands,genes,type="fiveprime",promoterDist = 2500, skipExons = FALSE, verbose=TRUE)

write.table(tab, file="/data/htp/A07/MM1_vs_MM6/Annotated_DMRs.csv", sep=",", row.names=FALSE)







##### bumphunter ############################################################################################
L: The number of genomic locations for which to simulate data
clusterL: The number of locations in each cluster.

# make a cutoff on the fwer!!!!what cutoff and why?
#  what is area?
# what cutoff does bumphunter take by default?

# cutoff in bumphunter used to kick out dmrs before bumphunter, will be anyway not significant
# on what exactely is the cutoff?

# The FWER is the probability of making at least one type I error in the family
# should be<0.05 meaning that the probability of making an error is below 5 %
# error oh what?--> of the multiple testing?--> that one probe is wrong?

# area: explained in "package bumphunter"
############################################################################################################

###### The GRanges object ##################################################################################
# The GRanges class is a container for the genomic locations and their associated annotations.

### How to set it up ?
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
# or get ranges out of range based object
gr<-granges(GenomicRatioSet)

#### metadata of GRanges
values(gr)
mcols(gr)

### show colnames of metadata
colnames(values(gr))
colnames(mcols(gr))

### add columns to metadata
values(gr) <- cbind(values(gr), beta)

### show ranges and metadata columns
head(gr)
############################################################################################################









