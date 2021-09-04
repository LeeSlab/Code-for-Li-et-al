
## title #####plot the DLR and ICF
## 2020-12-20
## Fan Wu
## DLR&ICF


# Part1 DLR -------------------------------------------------------------------

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/DLR&ICF")

## DLR
MyRead <- function(path,ss,tt,
                   ...){
  temp_file <- read.delim2(path,sep='\t',header=F, skip = 1)
  temp_file$sample <- ss
  temp_file$type <- tt
  temp_file$V4 <- as.numeric(temp_file$V4)
  return(temp_file)
}

ProR1 <- MyRead('./homerProR1.DLR.bedGraph','ProE R1','DLR')
ProR2 <- MyRead('./homerProR2.DLR.bedGraph','ProE R2','DLR')
OrthoR1 <- MyRead('./homerOrthoR1.DLR.bedGraph','OrthoE R1','DLR')
OrthoR2 <- MyRead('./homerOrthoR2.DLR.bedGraph','OrthoE R2','DLR')
head(ProR1)

(median(ProR1$V4)+median(ProR2$V4))/2
(median(OrthoR1$V4)+median(OrthoR2$V4))/2

mergeDLR <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeDLR$sample <- factor(mergeDLR$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))

library(Rmisc)
sum = summarySE(mergeDLR, measurevar="V4", groupvars=c("sample","type"))
sum

library(ggplot2)
theme_wf <- theme(panel.grid.major =element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 0.7, colour = "black")) +
  theme(plot.title = element_text(size = 15, colour = "black", face = 'bold',hjust=0.5)) +
  theme(axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black")) + 
  theme(axis.ticks.x = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"))


pd = position_dodge(.2)
(p1 <- ggplot(sum, aes(x=sample,
                       y=V4)) +
    # geom_errorbar(aes(ymin=V4-se,
    #                   ymax=V4+se),
    #               width=.2, size=0.7) +
    geom_bar(stat="identity", width=0.5) +
    # geom_point(shape=15, size=4) +
    theme_bw() +
    theme(axis.title.y = element_text(vjust= 1.8),
          axis.title.x = element_text(vjust= -0.5),
          axis.title = element_text(face = "bold")) +
    scale_color_manual(values = c("red","red", "blue","blue")) +
    # ylim(0.3,0.5) +
    ylab("Boundary strength") +
    xlab("Compartment") + 
    theme(axis.title.x = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(size = 20, colour = "black"),
          axis.text.x = element_text(size = 13, colour = "black"),
          axis.text.y = element_text(size = 13, colour = "black")))

(p1 <- ggplot(sum, aes(V4)) +
    geom_bar())

ggplot(sum,
       aes(x=sample,y=V4))+
  geom_bar(stat="identity")


(p <- ggplot(mergeDLR,aes(x=sample,y=V4,color=sample))+
  geom_boxplot(na.rm = T, outlier.shape = NA) + 
  ylab("DLR(distal to local)") + 
  # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
  ylim(-3,3.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme_wf)


pdf('DLR.pdf',5,6.5)
print(p)
print(p)
dev.off()

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16
wilcox.test(ProR1$V4, OrthoR1$V4)
wilcox.test(ProR1$V4, OrthoR1$V4)


# Part2 ICF -------------------------------------------------------------------

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/DLR&ICF")

## ICF
MyRead <- function(path,ss,tt,
                   ...){
  temp_file <- read.delim2(path,sep='\t',header=F, skip = 1)
  temp_file$sample <- ss
  temp_file$type <- tt
  temp_file$V4 <- as.numeric(temp_file$V4)
  return(temp_file)
}

ProR1 <- MyRead('./homerProR1.ICF.bedGraph','ProE R1','ICF')
ProR2 <- MyRead('./homerProR2.ICF.bedGraph','ProE R2','ICF')
OrthoR1 <- MyRead('./homerOrthoR1.ICF.bedGraph','OrthoE R1','ICF')
OrthoR2 <- MyRead('./homerOrthoR2.ICF.bedGraph','OrthoE R2','ICF')
head(ProR1)

mergeICF <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeICF$sample <- factor(mergeICF$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))

library(Rmisc)
sum = summarySE(mergeICF, measurevar="V4", groupvars=c("sample","type"))
sum

library(ggplot2)
theme_wf <- theme(panel.grid.major =element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 0.7, colour = "black")) +
  theme(plot.title = element_text(size = 15, colour = "black", face = 'bold',hjust=0.5)) +
  theme(axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black")) + 
  theme(axis.ticks.x = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"))


(p <- ggplot(mergeICF,aes(x=sample,y=V4,color=sample))+
    geom_boxplot(na.rm = T, outlier.shape = NA) + 
    ylab("ICF(inter-chromosome frequency)") + 
    # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
    ylim(0,0.75) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme_wf)


pdf('ICF.pdf',5,6.5)
print(p)
dev.off()

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16



# Part3 DLR in A and B -------------------------------------------------------------------

MyRead <- function(path,ss,tt,
                   ...){
  temp_file <- read.delim2(path,sep='\t',header=F, skip = 1)
  temp_file$sample <- ss
  temp_file$type <- tt
  temp_file$V4 <- as.numeric(temp_file$V4)
  return(temp_file)
}

ProR1 <- MyRead('./DLR_in_compartment/ProR1DLRInProR1.A.bed','ProE R1','DLR')
ProR2 <- MyRead('./DLR_in_compartment/ProR2DLRInProR2.A.bed','ProE R2','DLR')
OrthoR1 <- MyRead('./DLR_in_compartment/OrthoR1DLRInOrthoR1.A.bed','OrthoE R1','DLR')
OrthoR2 <- MyRead('./DLR_in_compartment/OrthoR2DLRInOrthoR2.A.bed','OrthoE R2','DLR')
head(ProR1)

median(OrthoR2$V4)

(median(ProR1$V4)+median(ProR2$V4))/2
(median(OrthoR1$V4)+median(OrthoR2$V4))/2

mergeDLR <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeDLR$sample <- factor(mergeDLR$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))

library(Rmisc)
sum = summarySE(mergeDLR, measurevar="V4", groupvars=c("sample","type"))
sum

library(ggplot2)
theme_wf <- theme(panel.grid.major =element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 0.7, colour = "black")) +
  theme(plot.title = element_text(size = 15, colour = "black", face = 'bold',hjust=0.5)) +
  theme(axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black")) + 
  theme(axis.ticks.x = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"))



(p <- ggplot(mergeDLR,aes(x=sample,y=V4,color=sample))+
    geom_boxplot(na.rm = T, outlier.shape = NA) + 
    ylab("DLR(distal to local)") + 
    ggtitle('DLR in comaprtment A') +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
    ylim(-3,3.5) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme_wf)

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16


ProR1 <- MyRead('./DLR_in_compartment/ProR1DLRInProR1.B.bed','ProE R1','DLR')
ProR2 <- MyRead('./DLR_in_compartment/ProR2DLRInProR2.B.bed','ProE R2','DLR')
OrthoR1 <- MyRead('./DLR_in_compartment/OrthoR1DLRInOrthoR1.B.bed','OrthoE R1','DLR')
OrthoR2 <- MyRead('./DLR_in_compartment/OrthoR2DLRInOrthoR2.B.bed','OrthoE R2','DLR')
head(ProR1)


(median(ProR1$V4)+median(ProR2$V4))/2
(median(OrthoR1$V4)+median(OrthoR2$V4))/2


mergeDLR <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeDLR$sample <- factor(mergeDLR$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))

(p1 <- ggplot(mergeDLR,aes(x=sample,y=V4,color=sample))+
    geom_boxplot(na.rm = T, outlier.shape = NA) + 
    ylab("DLR(distal to local)") + 
    ggtitle('DLR in comaprtment B') +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
    ylim(-3,3.5) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme_wf)

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16

pdf('boxplot_DLRInCompAB.pdf',5,6.5)
print(p)
print(p1)
dev.off()




# Part4 ICF in A and B -------------------------------------------------------------------

ProR1 <- MyRead('./ICF_in_compartment/ProR1ICFInProR1.A.bed','ProE R1','ICF')
ProR2 <- MyRead('./ICF_in_compartment/ProR2ICFInProR2.A.bed','ProE R2','ICF')
OrthoR1 <- MyRead('./ICF_in_compartment/OrthoR1ICFInOrthoR1.A.bed','OrthoE R1','ICF')
OrthoR2 <- MyRead('./ICF_in_compartment/OrthoR2ICFInOrthoR2.A.bed','OrthoE R2','ICF')
head(ProR1)

mergeDLR <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeDLR$sample <- factor(mergeDLR$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))



(p <- ggplot(mergeDLR,aes(x=sample,y=V4,color=sample))+
    geom_boxplot(na.rm = T, outlier.shape = NA) + 
    ylab("ICF(inter-chromosome frequency))") + 
    ggtitle('ICF in comaprtment A') +
    # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
    ylim(0,0.75) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme_wf)

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16


ProR1 <- MyRead('./ICF_in_compartment/ProR1ICFInProR1.B.bed','ProE R1','ICF')
ProR2 <- MyRead('./ICF_in_compartment/ProR2ICFInProR2.B.bed','ProE R2','ICF')
OrthoR1 <- MyRead('./ICF_in_compartment/OrthoR1ICFInOrthoR1.B.bed','OrthoE R1','ICF')
OrthoR2 <- MyRead('./ICF_in_compartment/OrthoR2ICFInOrthoR2.B.bed','OrthoE R2','ICF')
head(ProR1)

mergeDLR <- rbind(ProR1, ProR2, OrthoR1, OrthoR2)
mergeDLR$sample <- factor(mergeDLR$sample, levels = c('ProE R1', 'ProE R2', 'OrthoE R1', 'OrthoE R2'))

(p1 <- ggplot(mergeDLR,aes(x=sample,y=V4,color=sample))+
    geom_boxplot(na.rm = T, outlier.shape = NA) + 
    ylab("ICF(inter-chromosome frequency))") + 
    ggtitle('ICF in comaprtment B') +
    # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey56") +
    ylim(0,0.75) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme_wf)

t.test(ProR1$V4, OrthoR1$V4) ## p-value < 2.2e-16
t.test(ProR2$V4, OrthoR2$V4) ## p-value < 2.2e-16

pdf('boxplot_ICFInCompAB.pdf',5,6.5)
print(p)
print(p1)
dev.off()





