


#####################################################################################
## 2019-11-17
## Fan Wu
## calculate the boundary strength
#####################################################################################

rm(list=ls())

setwd("~/boundary_compare")
source("WFFunction.R")

options(digits=6)
options(stringsAsFactors = F)

suppressPackageStartupMessages({
  library(ggplot2)
  library(gplots)
})


MyCallBoundary <- function(chr = "chr1"){
  inputFile1name <- paste0("HiCProR1_",chr,".boundaries")
  inputFile2name <- paste0("HiCProR2_",chr,".boundaries")
  inputFile3name <- paste0("HiCOrthoR1_",chr,".boundaries")
  inputFile4name <- paste0("HiCOrthoR2_",chr,".boundaries")
  
  inputFile1 <- read.delim2(inputFile1name, sep="\t", header = T)
  inputFile2 <- read.delim2(inputFile2name, sep="\t", header = T)
  inputFile3 <- read.delim2(inputFile3name, sep="\t", header = T)
  inputFile4 <- read.delim2(inputFile4name, sep="\t", header = T)
  
  inputFile1$type <- "ProHiCR1"
  inputFile2$type <- "ProHiCR2"
  inputFile3$type <- "OrthoHiCR1"
  inputFile4$type <- "OrthoHiCR2"
  
  plotfile <- rbind(inputFile1,inputFile2,inputFile3,inputFile4)
  plotfile$insulationScore <- as.numeric(plotfile$insulationScore)
  plotfile$type <- factor(plotfile$type, level= c("ProHiCR1","ProHiCR2","OrthoHiCR1","OrthoHiCR2"))
  
  p <- ggplot(plotfile, aes(x = type, y=insulationScore, colour = type)) +
    labs(y = "Boundary Strength") +
    stat_boxplot(geom = "errorbar", width=0.2) +
    geom_boxplot() + 
    ggtitle(chr) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


p1 <- MyCallBoundary("chr1")
p2 <- MyCallBoundary("chr2")
p3 <- MyCallBoundary("chr3")
p4 <- MyCallBoundary("chr4")
p5 <- MyCallBoundary("chr5")

p6 <- MyCallBoundary("chr6")
p7 <- MyCallBoundary("chr7")
p8 <- MyCallBoundary("chr8")
p9 <- MyCallBoundary("chr9")
p10 <- MyCallBoundary("chr10")

p11 <- MyCallBoundary("chr11")
p12 <- MyCallBoundary("chr12")
p13 <- MyCallBoundary("chr13")
p14 <- MyCallBoundary("chr14")
p15 <- MyCallBoundary("chr15")

p16 <- MyCallBoundary("chr16")
p17 <- MyCallBoundary("chr17")
p18 <- MyCallBoundary("chr18")
p19 <- MyCallBoundary("chr19")
p20 <- MyCallBoundary("chr20")

p21 <- MyCallBoundary("chr21")
p22 <- MyCallBoundary("chr22")
p23 <- MyCallBoundary("chrX")


pdf("boundaryStrength.pdf",14,17)
plot_grid(p1,p2,p3,p4,p5,p6,p7,
          p8,p9,p10,p11,p12,p13,p14,p15,p16,
          nrow = 4, ncol = 4)
plot_grid(p17,p18,p19,p20, p21,p22,
          nrow = 4, ncol = 4)
dev.off()



## calculate the power
MyPowOneVar(groups = 4, 
            means = c(0.452, 0.459, 0.368, 0.143),
            sd = mean(0.184,0.187,0.137,0.143), 
            sig.level = 0.01, 
            power = 0.90)



## visualize the summary of infomation
library(FSA)
Summarize(insulationScore ~ type, data=plotfile, digits=3)


MyCalSig <- function(plotfile){
  model = lm(insulationScore ~ type, data=plotfile)
  # summary(model)
  
  # library(car)
  # Anova(model, type="II")
  # anova(model)
  # summary(model)

  # library(agricolae)
  # (HSD.test(model, "type"))
  # (LSD.test(model, "type", alpha = 0.01, p.adj="fdr"))
  
  library(multcomp)
  mc = glht(model,mcp(type = "Tukey"))
  mcs = summary(mc, test=adjusted("single-step"))
  return(mcs)
}



mergeAllChr <- rbind(p1$data,p2$data,p3$data,p4$data,p5$data,p6$data,p7$data,p8$data,p9$data,p10$data,
                     p11$data,p12$data,p13$data,p14$data,p15$data,p16$data,p17$data,p18$data,p19$data,p20$data,
                     p21$data,p22$data)

pdf("boundaryStrengthAllChr.pdf",8,7.5)
ggplot(mergeAllChr, aes(x = type, y=insulationScore, colour = type)) +
  labs(y = "Boundary Strength") +
  stat_boxplot(geom = "errorbar", width=0.2) +
  geom_boxplot() + 
  ggtitle("All Chromosomes") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


MyCalSig(mergeAllChr)


## calculate the test value 
plotfile <- mergeAllChr
model = lm(insulationScore ~ type, data=plotfile)
library(multcomp)
mc = glht(model,mcp(type = "Tukey"))
mcs = summary(mc, test=adjusted("single-step"))
# return(mcs)
mcs$test$pvalues
mcs




