
### -------- 0. analysis TAD score in Ti4h & Nu ----------------------------
## Fan
## 20220403

rm(list=ls())
options(stringsAsFactors = F)
options(scipen=2)

setwd("~/Documents/work/20180429_RBC_development_LD/20220223_HiC_NuTi_terminal_LD/TADscore/")
source('../../../public/WFMyFunction.R')


### -------- 1. read files -------------------------


MyReadFolder <- function(folder){
  tempFiles <- MyReadDelim(folder, sep='\t', header =F)
  tempFiles <- tempFiles[,c(1:5,8)]
  colnames(tempFiles) <- c('chr','start','end','binS','binE','TADscore')
  return(tempFiles)
}

folders <- list.files('./', pattern = 'bed')
readFileList <- lapply(folders, MyReadFolder)
names(readFileList) <- c('Nuclei_R1', 'Nuclei_R2', 
                         'Ti4h_R1', "Ti4h_R2", 
                         "Ti4hCon_R1", "Ti4hCon_R2")


NuR1 <- readFileList[[1]]
NuR2 <- readFileList[[2]]
Ti4hR1 <- readFileList[[3]]
Ti4hR2 <- readFileList[[4]]
TiConR1 <- readFileList[[5]]
TiConR2 <- readFileList[[6]]

NuR1 <- NuR1[!is.na(NuR1$TADscore),]
NuR2 <- NuR2[!is.na(NuR2$TADscore),]
Ti4hR1 <- Ti4hR1[!is.na(Ti4hR1$TADscore),]
Ti4hR2 <- Ti4hR2[!is.na(Ti4hR2$TADscore),]
TiConR1 <- TiConR1[!is.na(TiConR1$TADscore),]
TiConR2 <- TiConR2[!is.na(TiConR2$TADscore),]

NuR1$sample <- 'Nuclei R1'
NuR2$sample <- 'Nuclei R2'
Ti4hR1$sample <- '4h Ti R1'
Ti4hR2$sample <- '4h Ti R2'
TiConR1$sample <- '4h control R1'
TiConR2$sample <- '4h control R2'

mergeAllTADscore <- rbind(NuR1, NuR2, Ti4hR1, Ti4hR2, TiConR1, TiConR2)
head(mergeAllTADscore)
mergeAllTADscore$sample <- factor(mergeAllTADscore$sample, 
                                  levels = c("4h control R1", "4h control R2",
                                             "4h Ti R1", "4h Ti R2",
                                             "Nuclei R1", "Nuclei R2"))


### -------- 2. plot the boxplot -------------------------
library(ggplot2)
mergeAllTADscore$binLength <- mergeAllTADscore$binE - mergeAllTADscore$binS + 1
mergeAllTADscore <- mergeAllTADscore[mergeAllTADscore$binLength > 0,]

(p1 <- ggplot(mergeAllTADscore,
             aes(sample,TADscore,fill=sample)) +
  geom_boxplot() +
  theme_wf)

ggsave('./TADscore_compare.pdf', plot= p1, width = 9.5, height = 6)

t.test(mergeAllTADscore[mergeAllTADscore$sample%in%c('4h Ti R1'),]$TADscore,
       mergeAllTADscore[mergeAllTADscore$sample%in%c('4h control R1'),]$TADscore)$p.value
t.test(mergeAllTADscore[mergeAllTADscore$sample%in%c('4h Ti R2'),]$TADscore,
       mergeAllTADscore[mergeAllTADscore$sample%in%c('4h control R2'),]$TADscore)$p.value


t.test(mergeAllTADscore[mergeAllTADscore$sample%in%c('Nuclei R1'),]$TADscore,
       mergeAllTADscore[mergeAllTADscore$sample%in%c('4h control R1'),]$TADscore)$p.value
t.test(mergeAllTADscore[mergeAllTADscore$sample%in%c('Nuclei R2'),]$TADscore,
       mergeAllTADscore[mergeAllTADscore$sample%in%c('4h control R2'),]$TADscore)$p.value


