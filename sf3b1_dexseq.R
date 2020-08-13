#Exon counts created on compute1 as described here:
#/storage1/fs1/mgriffit/Active/kommagani_SplicingAnalysis/get_DEXseq_exon_counts.txt

#Set up parallel compute
#library("BiocParallel")
#BPPARAM = MulticoreParam(6)

#Read in the exon count data
inDir = "/Users/mgriffit/Google Drive/Manuscripts/SF3B1_Kommagani/dex-seq-analysis/exon_counts"
countFiles = list.files(inDir, pattern="exon_counts.tsv$", full.names=TRUE)
basename(countFiles)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

#Set output dir
outDir = "/Users/mgriffit/Google Drive/Manuscripts/SF3B1_Kommagani/dex-seq-analysis/exon_counts/results/"
setwd(outDir)

#Prepare a sample table
sampleTable = data.frame(
  row.names = c( "control_rep1", "control_rep2", "control_rep3", 
                 "sf3b1_rep1", "sf3b1_rep2", "sf3b1_rep3" ),
  condition = c("control", "control", "control",  
                "sf3b1", "sf3b1", "sf3b1" ),
  libType = c( "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end" ) )

library( "DEXSeq" )

dxd = DEXSeqDataSetFromHTSeq( countFiles, sampleData=sampleTable, 
                              design= ~ sample + exon + condition:exon,
                              flattenedfile=flattenedFile )

#Inspect the example data
colData(dxd)
head( counts(dxd), 5 )

#Note that the number of columns is 12, the first six (we have six samples) 
# corresponding to the number of reads mapping to out exonic regions and the 
# last six correspond to the sum of the counts mapping to the rest of the 
# exons from the same gene on each sample.
split( seq_len(ncol(dxd)), colData(dxd)$exon )
head( featureCounts(dxd), 5 )
head( rowRanges(dxd), 3 )
sampleAnnotation( dxd )

#Perform normalization
dxd = estimateSizeFactors( dxd )

#Perform dispersion estimation
dxd = estimateDispersions( dxd )
#dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
 
#As a shrinkage diagnostic, the DEXSeqDataSet use the method plotDispEsts() 
# that plots the per-exon dispersion estimates versus the mean normalised count,
# the resulting fitted values and the a posteriori (shrinked) dispersion 
# estimates (Figure 1). 
plotDispEsts( dxd )

#Test for differential expression 
#dxd = testForDEU( dxd, BPPARAM=BPPARAM )
dxd = testForDEU( dxd )

#Get fold change values for differential exon usage
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
#dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

#Get the final results
dxr1 = DEXSeqResults( dxd )
dxr1

#Show column descriptions
mcols(dxr1)$description

#How many exonic regions are significant with a false discovery rate of 10%:
table ( dxr1$padj < 0.1 )  #61583

#We may also ask how many genes are affected:
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) ) #10766

#To see how the power to detect differential exon usage depends on the number 
# of reads that map to an exon, a so-called MA plot is useful, which plots the 
# logarithm of fold change versus average normalized count per exon and marks 
# by red color the exons which are considered significant; here, the exons with 
# an adjusted p values of less than 0.1 (There is of course nothing special 
# about the number 0.1, and you can specify other thresholds in the call to plotMA().
plotMA( dxr1, cex=0.8 )

#Pull some examples with strong pvalues and fold changes
x = dxr1[(which (dxr1$padj < 0.001 & abs(dxr1$log2fold_sf3b1_control) > 5)),]

#Try some visualizations
plotDEXSeq( dxr1, "ENSG00000073711", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr1, "ENSG00000132915", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Different visualization options ...

#Show transcripts
plotDEXSeq( dxr1, "ENSG00000073711", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Show counts from individual samples
plotDEXSeq( dxr1, "ENSG00000073711", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# DEXSeq is designed to find changes in relative exon usage, i.e., changes in 
# the expression of individual exons that are not simply the consequence of 
# overall up- or down-regulation of the gene. To visualize such changes, it is 
# sometimes advantageous to remove overall changes in expression from the plots. 
# Use the option splicing=TRUE for this purpose.
plotDEXSeq( dxr1, "ENSG00000073711", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Combine the options above
plotDEXSeq( dxr1, "ENSG00000073711", displayTranscripts=TRUE, expression=FALSE, norCounts=TRUE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#Create plots for the top X by p-value


#Filter down to examples that meet certain criteria to find good examples?


#Create a volcano plot


#Create a heatmap using exon values from a subset of genes (e.g. known targets for SF3B1 splicing)
outDir2 = "/Users/mgriffit/Google Drive/Manuscripts/SF3B1_Kommagani/dex-seq-analysis/exon_counts/results/auto/"
setwd(outDir2)
x = dxr1[(which (dxr1$padj < 0.001 & abs(dxr1$log2fold_sf3b1_control) > 2)),]
o = order(x$padj)
y = x[o[1:1000],]
z = unique(y$groupID)[1:100]
for (i in 1:100){
  gene = z[i]
  name = paste(gene, ".pdf", sep="")
  pdf(file=name)
  plotDEXSeq( dxr1, gene, expression=FALSE, norCounts=TRUE, splicing=TRUE,
              legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()
}
setwd(outDir)

#Write out the significant results
filtered_data = dxr1[(which (dxr1$padj < 0.01 & abs(dxr1$log2fold_sf3b1_control) > 2)),]
dim(filtered_data)
length(unique(filtered_data[,1]))

library("xlsx")
write.xlsx(as.data.frame(filtered_data[,1:12]), file="DEXseq_Significant_SF3B1vsControl.xlsx", sheetName = "Sheet1",col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(as.data.frame(filtered_data[,1:12]), file="DEXseq_Significant_SF3B1vsControl.tsv", quote=FALSE, sep="\t", row.names=FALSE)

#Create a web report
outDir3 = "/Users/mgriffit/Google Drive/Manuscripts/SF3B1_Kommagani/dex-seq-analysis/exon_counts/results/html/"
setwd(outDir3)
DEXSeqHTML( dxr1, FDR=0.01, color=c("#FF000080", "#0000FF80") )
setwd(outDir)
