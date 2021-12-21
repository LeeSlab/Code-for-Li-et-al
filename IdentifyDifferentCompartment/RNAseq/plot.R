
MyReadFiles <- function(inputfiles1,inputfiles2,COMPtype){
  tempFile1 <- MyReadDelim(inputfiles1,sep='\t',header=F)
  tempFile2 <- MyReadDelim(inputfiles2,sep='\t',header=F)
  tempFile <- rbind(tempFile1,tempFile2)
  tempFile <- tempFile[!duplicated(tempFile$V5),]
  rownames(tempFile) <- tempFile$V5
  tempFile$geneType <- geneInfo[rownames(tempFile),'V8']
  tempFile$COMPtype <- COMPtype
  return(tempFile)
}

geneinProR1A <- MyReadDelim('./geneIn_ProR1_compartmentA.txt',sep='\t',header=F)
rownames(geneinProR1A) <- geneinProR1A$V5
geneinProR1B <- MyReadDelim('./geneIn_ProR1_compartmentB.txt',sep='\t',header=F)
rownames(geneinProR1B) <- geneinProR1B$V5
geneinOrthoR1A <- MyReadDelim('./geneIn_OrthoR1_compartmentA.txt',sep='\t',header=F)
rownames(geneinOrthoR1A) <- geneinOrthoR1A$V5
geneinOrthoR1B <- MyReadDelim('./geneIn_OrthoR1_compartmentB.txt',sep='\t',header=F)
rownames(geneinOrthoR1B) <- geneinOrthoR1B$V5

geneinProR2A <- MyReadDelim('./geneIn_ProR2_compartmentA.txt',sep='\t',header=F)
rownames(geneinProR2A) <- geneinProR2A$V5
geneinProR2B <- MyReadDelim('./geneIn_ProR2_compartmentB.txt',sep='\t',header=F)
rownames(geneinProR2B) <- geneinProR2B$V5
geneinOrthoR2A <- MyReadDelim('./geneIn_OrthoR2_compartmentA.txt',sep='\t',header=F)
rownames(geneinOrthoR2A) <- geneinOrthoR2A$V5
geneinOrthoR2B <- MyReadDelim('./geneIn_OrthoR2_compartmentB.txt',sep='\t',header=F)
rownames(geneinOrthoR2B) <- geneinOrthoR2B$V5



myfilter <- function(filetemp,COM){
  filetemp <- filetemp[filetemp$V7 == 'protein_coding',]
  filetemp$V8 <- COM
  return(filetemp)
}

geneinProR1A <- myfilter(geneinProR1A,'A')
geneinProR1B <- myfilter(geneinProR1B,'B')
geneinOrthoR1A <- myfilter(geneinOrthoR1A,'A')
geneinOrthoR1B <- myfilter(geneinOrthoR1B,'B')

geneinProR2A <- myfilter(geneinProR2A,'A')
geneinProR2B <- myfilter(geneinProR2B,'B')
geneinOrthoR2A <- myfilter(geneinOrthoR2A,'A')
geneinOrthoR2B <- myfilter(geneinOrthoR2B,'B')

tpm_strandSpecific <- read.delim2('../../../../../../20181124_RNAseq_chain_specificity_LD/htseq/tpm.rev.csv',sep=',')
rownames(tpm_strandSpecific)<-tpm_strandSpecific$gene_id
head(tpm_strandSpecific)


plotData <- data.frame(TPM_value = c(tpm_strandSpecific[rownames(geneinProR1A),'ProR1'],
                                     tpm_strandSpecific[rownames(geneinProR1B),'ProR1'],
                                     tpm_strandSpecific[rownames(geneinProR2A),'ProR2'],
                                     tpm_strandSpecific[rownames(geneinProR2B),'ProR2'],
                                     tpm_strandSpecific[rownames(geneinOrthoR1A),'OrthoR1'],
                                     tpm_strandSpecific[rownames(geneinOrthoR1B),'OrthoR1'],
                                     tpm_strandSpecific[rownames(geneinOrthoR2A),'OrthoR2'],
                                     tpm_strandSpecific[rownames(geneinOrthoR2B),'OrthoR2']),
                       sample = c(rep('ProR1',length(c(rownames(geneinProR1A),rownames(geneinProR1B)))),
                                  rep('ProR2',length(c(rownames(geneinProR2A),rownames(geneinProR2B)))),
                                  rep('OrthoR1',length(c(rownames(geneinOrthoR1A),rownames(geneinOrthoR1B)))),
                                  rep('OrthoR2',length(c(rownames(geneinOrthoR2A),rownames(geneinOrthoR2B))))),
                       COMP = c(geneinProR1A$V8, geneinProR1B$V8, geneinProR2A$V8, geneinProR2B$V8,
                                geneinOrthoR1A$V8,geneinOrthoR1B$V8,geneinOrthoR2A$V8,geneinOrthoR2B$V8))
library(ggplot2)
plotData$TPM_value <- as.numeric(plotData$TPM_value)
plotData$sample <- factor(plotData$sample, levels = c('ProR1','ProR2','OrthoR1','OrthoR2'))
p <- ggplot(plotData, aes(sample,log10(TPM_value+1),fill=COMP)) +
  geom_boxplot(outlier.alpha=0) +
  ylab("log10(TPM+1)") + 
  ylim(0,3) +
  theme_wf

ggsave('./geneExpInCompAorB.pdf',plot = p,width = 5,height = 5)
