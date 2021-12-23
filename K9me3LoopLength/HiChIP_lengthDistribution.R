
#title#######################################################
## 20210226
## Fan Wu
## Test loop length
#title#######################################################


rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("~/Documents/work/20180429_RBC_development_LD/20201119_HiChIP_terminal_K9me3_LD/")
source('../../public/WFMyFunction.R')

## using loops length > 5kb

MyTreatH3K9HiChIP <- function(filePath, sampleName){
  temp_file <- read.delim2(filePath,header=F,sep=' ')
  temp_file$sample <- sampleName
  temp_file$length <- (temp_file$V6+temp_file$V5)/2-(temp_file$V2+temp_file$V3)/2
  temp_file$Peak1 <- as.integer((temp_file$V2+temp_file$V3)/2)
  temp_file$Peak2 <- as.integer((temp_file$V6+temp_file$V5)/2)
  temp_file$CPM <- temp_file$V8/sum(temp_file$V8)*1e6
  temp_file <- temp_file[temp_file$length >= 5000,]
  return(temp_file)
}

K9me3_ProR1 <- MyTreatH3K9HiChIP('./HiCH_Pro_K9me3_R1_callLoopByHichipper/Pro_K9me3_R1.intra.loop_counts.bedpe', "Pro_K9me3_R1")
K9me3_ProR2 <- MyTreatH3K9HiChIP('./HiCH_Pro_K9me3_R2_callLoopByHichipper/Pro_K9me3_R2.intra.loop_counts.bedpe', "Pro_K9me3_R2")
K9me3_OrthoR1 <- MyTreatH3K9HiChIP('./HiCH_Ortho_K9me3_R1_callLoopByHichipper/Ortho_K9me3_R1.intra.loop_counts.bedpe', "Ortho_K9me3_R1")
K9me3_OrthoR2 <- MyTreatH3K9HiChIP('./HiCH_Ortho_K9me3_R2_callLoopByHichipper/Ortho_K9me3_R2.intra.loop_counts.bedpe', "Ortho_K9me3_R2")


write.table(K9me3_ProR1, 'K9me3_ProR1_MT5kb.bed', sep='\t', col.names = F, row.names = F, quote = F)
write.table(K9me3_ProR2, 'K9me3_ProR2_MT5kb.bed', sep='\t', col.names = F, row.names = F, quote = F)
write.table(K9me3_OrthoR1, 'K9me3_OrthoR1_MT5kb.bed', sep='\t', col.names = F, row.names = F, quote = F)
write.table(K9me3_OrthoR2, 'K9me3_OrthoR2_MT5kb.bed', sep='\t', col.names = F, row.names = F, quote = F)


library(ggplot2)
plotData <- rbind(K9me3_ProR1,K9me3_ProR2,K9me3_OrthoR1,K9me3_OrthoR2)
theme_wf <- theme(panel.grid.major =element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 1, colour = "black")) +
  theme(plot.title = element_text(size = 25, colour = "black", face = 'bold',hjust=0.5)) +
  theme(axis.title.x = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black")) + 
  theme(legend.text = element_text(size = 25, colour = "black")) +
  theme(legend.title = element_text(size = 25, colour = "black")) + 
  theme(axis.ticks.x = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 18, colour = "black",vjust = 1,hjust = 1,angle = 30),
        axis.text.y = element_text(size = 18, colour = "black"))


plotData$sample <- factor(plotData$sample, levels = unique(plotData$sample))


pdf('Test_HiChIP_loop_length_distribution_5kb.pdf',6,5)
(p1 <- ggplot(plotData, aes(x=length, color = sample)) + 
    # geom_histogram(bins=100,fill="white", aes(y=..density..), position="identity", alpha=0.5) +
    geom_density(alpha=0.6, size = 1) +
    labs(title="Loop length distribution", x = "Loop length", y = "Density") +
    # geom_vline(xintercept=1e6, color="red",linesample="dashed") +
    # geom_vline(xintercept=2e6, color="blue",linesample="dashed") +
    # xlim(0,500000) +
    scale_x_log10() +
    theme_wf)

(p2 <- ggplot(plotData, aes(x=sample, y=length, color = sample)) + 
    geom_boxplot(outlier.shape=NA) +
    labs(title="Loop length", x = "H3K9me3", y = "Size of loops") +
    scale_y_log10() +
    geom_hline(yintercept=1e6, color="blue",linesample="dashed") +
    # ylim(0,1e+07) +
    # geom_vline(xintercept=2e6, color="blue",linesample="dashed") +
    theme_wf)
dev.off()

