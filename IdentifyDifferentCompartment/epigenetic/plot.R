

temp_colNames <- c('chr','start','end',
                   'H3K27ac_ProR1','H3K27ac_ProR2','H3K27ac_OrthoR1','H3K27ac_OrthoR2',
                   'H3K9me3_ProR1','H3K9me3_ProR2','H3K9me3_OrthoR1','H3K9me3_OrthoR2')

histoneInProR1A <- MyReadDelim('./histoneIn_ProR1_compartmentA.txt',sep='\t',header=F)
colnames(histoneInProR1A) <- temp_colNames
histoneInProR1A$COMP <- 'A'
histoneInProR1A$sample <- 'ProR1'
histoneInProR1B <- MyReadDelim('./histoneIn_ProR1_compartmentB.txt',sep='\t',header=F)
colnames(histoneInProR1B) <- temp_colNames
histoneInProR1B$COMP <- 'B'
histoneInProR1B$sample <- 'ProR1'
histoneInOrthoR1A <- MyReadDelim('./histoneIn_OrthoR1_compartmentA.txt',sep='\t',header=F)
colnames(histoneInOrthoR1A) <- temp_colNames
histoneInOrthoR1A$COMP <- 'A'
histoneInOrthoR1A$sample <- 'OrthoR1'
histoneInOrthoR1B <- MyReadDelim('./histoneIn_OrthoR1_compartmentB.txt',sep='\t',header=F)
colnames(histoneInOrthoR1B) <- temp_colNames
histoneInOrthoR1B$COMP <- 'B'
histoneInOrthoR1B$sample <- 'OrthoR1'

histoneInProR2A <- MyReadDelim('./histoneIn_ProR2_compartmentA.txt',sep='\t',header=F)
colnames(histoneInProR2A) <- temp_colNames
histoneInProR2A$COMP <- 'A'
histoneInProR2A$sample <- 'ProR2'
histoneInProR2B <- MyReadDelim('./histoneIn_ProR2_compartmentB.txt',sep='\t',header=F)
colnames(histoneInProR2B) <- temp_colNames
histoneInProR2B$COMP <- 'B'
histoneInProR2B$sample <- 'ProR2'
histoneInOrthoR2A <- MyReadDelim('./histoneIn_OrthoR2_compartmentA.txt',sep='\t',header=F)
colnames(histoneInOrthoR2A) <- temp_colNames
histoneInOrthoR2A$COMP <- 'A'
histoneInOrthoR2A$sample <- 'OrthoR2'
histoneInOrthoR2B <- MyReadDelim('./histoneIn_OrthoR2_compartmentB.txt',sep='\t',header=F)
colnames(histoneInOrthoR2B) <- temp_colNames
histoneInOrthoR2B$COMP <- 'B'
histoneInOrthoR2B$sample <- 'OrthoR2'

mergeAll_histone <- rbind(histoneInProR1A,histoneInProR1B,histoneInOrthoR1A,histoneInOrthoR1B,
                          histoneInProR2A,histoneInProR2B,histoneInOrthoR2A,histoneInOrthoR2B)
mergeAll_histone$H3K27ac_Pro <- rowMeans(mergeAll_histone[,c("H3K27ac_ProR1","H3K27ac_ProR2")])
mergeAll_histone$H3K27ac_Ortho <- rowMeans(mergeAll_histone[,c("H3K27ac_OrthoR1","H3K27ac_OrthoR2")])
mergeAll_histone$H3K9me3_Pro <- rowMeans(mergeAll_histone[,c("H3K9me3_ProR1","H3K9me3_ProR2")])
mergeAll_histone$H3K9me3_Ortho <- rowMeans(mergeAll_histone[,c("H3K9me3_OrthoR1","H3K9me3_OrthoR2")])


mergeAll_histone <- mergeAll_histone[,c("chr","start","end","COMP","sample",
                                        "H3K27ac_Pro","H3K27ac_Ortho","H3K9me2_Pro","H3K9me2_Ortho","H3K9me3_Pro","H3K9me3_Ortho",
                                        "H4k20me1_Pro","H4k20me1_Ortho")]

head(mergeAll_histone)


plotH3K27ac <- data.frame(Signal=c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$H3K27ac_Pro,
                                   mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$H3K27ac_Ortho),
                          sample=rep(c("ProR1","ProR2","OrthoR1","OrthoR2"),each=dim(mergeAll_histone)[1]/4),
                          COMPtype = c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$COMP,
                                       mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$COMP),
                          Sampletype = c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$sample,
                                         mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$sample))
plotH3K27ac$Signal <- as.numeric(plotH3K27ac$Signal)
plotH3K27ac$sample <- factor(plotH3K27ac$sample,levels = c("ProR1","ProR2","OrthoR1","OrthoR2"))
(plotk27ac <- ggplot(plotH3K27ac, aes(sample,log2(Signal+1),fill=COMPtype)) +
  geom_boxplot(outlier.alpha=0) +
  ggtitle("h3k27ac") +
  ylim(0,0.5) +
  theme_wf)


plotH3K9me3 <- data.frame(Signal=c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$H3K9me3_Pro,
                                   mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$H3K9me3_Ortho),
                          sample=rep(c("ProR1","ProR2","OrthoR1","OrthoR2"),each=dim(mergeAll_histone)[1]/4),
                          COMPtype = c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$COMP,
                                       mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$COMP),
                          Sampletype = c(mergeAll_histone[mergeAll_histone$sample %in% c('ProR1','ProR2'),]$sample,
                                         mergeAll_histone[mergeAll_histone$sample %in% c('OrthoR1','OrthoR2'),]$sample))
plotH3K9me3$Signal <- as.numeric(plotH3K9me3$Signal)
plotH3K9me3$sample <- factor(plotH3K9me3$sample,levels = c("ProR1","ProR2","OrthoR1","OrthoR2"))
(plotk9me3 <- ggplot(plotH3K9me3, aes(sample,log2(Signal+1),fill=COMPtype)) +
  geom_boxplot(outlier.alpha=0) +
  ggtitle("h3k9me3") +
  ylim(0,1.5) +
  theme_wf)

