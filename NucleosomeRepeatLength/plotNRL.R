##############################################################################

rm(list = ls())
setwd("~/Documents/work/20180429_RBC_development_LD/20190506_Mnaseq-seq_terminal_LD/NucTools/51PE_merge")

library("pacman")
pacman::p_load( ggplot2, grid, plyr, zoo, optparse )

##############################################################################
##############################################################################

source("./../150PE_merge/function.R")

############################
# read command line options
############################
peaks_tweek<-1

out=NULL
dir="~/Documents/work/20180429_RBC_development_LD/20190506_Mnaseq-seq_terminal_LD/NucTools/51PE_merge"
info=NULL;maxX=2000;maxY=NULL;minY=NULL;span=0.05;peakW=20

# check if required arguments specified
wd <- dir
xlim<-maxX
span<-span
w<-peakW
sample_info <- paste("span=",span, ", window=",w,sep="")

subtitle.NRL<-c("Nucleosome repeat length")




input="ProMnaseq.merge.last.extSE.NRL.bed"
sample="ProMnaseqMergeLastExtSE"
out<-paste(sample,"_span_",span,"_win_",w,".png",sep="")
df.NRL<-load_data(wd,input)
dim(df.NRL)

minY<-min(df.NRL$frequency)
maxY<-max(df.NRL$frequency)
ylim<-c(minY,maxY)
tt1<-plot_approx(df.NRL[1:xlim,], sampleID, subtitle.NRL, sample_info, peaks_tweek, wd,output.NRL, ylim, xlim, span, w)




input="OrthoMnaseq.merge.last.extSE.NRL.bed"
sample="OrthoMnaseqMergeLastExtSE"
out<-paste(sample,"_span_",span,"_win_",w,".png",sep="")
df.NRL<-load_data(wd,input)
dim(df.NRL)

minY<-min(df.NRL$frequency)
maxY<-max(df.NRL$frequency)
ylim<-c(minY,maxY)
tt2<-plot_approx(df.NRL[1:xlim,], sampleID, subtitle.NRL, sample_info, peaks_tweek, wd,output.NRL, ylim, xlim, span, w)

tt1[[1]]$type <- "ProE R1"
tt2[[1]]$type <- "OrthoE R1"

tt1[[2]]$type <- "ProE R1"
tt2[[2]]$type <- "OrthoE R1"


tt1[[1]]$loess_line2 <- tt1[[1]]$loess_line/sum(tt1[[1]]$loess_line)
tt2[[1]]$loess_line2 <- tt2[[1]]$loess_line/sum(tt2[[1]]$loess_line)

NRL_all.df <- rbind(tt1[[1]], tt2[[1]])
NRL_all.df$type <- factor(NRL_all.df$type, levels = c("ProE R1", "OrthoE R1"))

reg_all.df <- rbind(tt1[[2]], tt2[[2]])
reg_all.df$type <- factor(reg_all.df$type, levels = c("ProE R1", "OrthoE R1"))
# ps2 <- data.frame(x = c(0,0), y=c(0,0),type=c("ProMnase", "OrthoMnase"))
# reg_all.df <- rbind(reg_all.df, ps2)

library(ggplot2)

reg.equation1<-lm_eqn(tt1[[2]])
reg.equation2<-lm_eqn(tt2[[2]])

NRL_all.df <- NRL_all.df[NRL_all.df$distance > 100,]

source('../../../../public/WFMyFunction.R')
MyWriteTab(NRL_all.df, 'NRL_all.csv',sep=',')

pdf("00_Pro_Ortho_merge_again_extSE_NRL.pdf",8,7)
ggplot(NRL_all.df, aes(x = distance, y = loess_line/10^4, colour = type)) +
  geom_line(size=0.6) +
  xlim(0,2000) +
  labs(x = "Distance(bp)", y="Counts(x 10^4)", title="Nucleosome Repeat Length")

ggplot(NRL_all.df, aes(x = distance, y = loess_line2*10^4, colour = type)) +
  geom_line(size=0.6) +
  xlim(0,2000) +
  ylim(4.75,5.5) +
  labs(x = "Distance(bp)", y="Normalized Counts Density(10^-4)", title="Nucleosome Repeat Length")
  
ggplot(reg_all.df, aes(x = x, y = y, colour = type)) +
  geom_point() +
  # geom_line() +
  geom_smooth(se = FALSE, method = "lm", formula=y~x-1, size = 0.5) +
  scale_x_continuous(limits = c(1,8),breaks = c(1:8)) +
  labs(y = "Distances (bp)", x = "Peak Number") + 
  geom_text(x = 1.5, y = 1734, hjust=0, label = reg.equation1, parse = T, size=6, show.legend = F) +
  geom_text(x = 1.5, y = 1600, hjust=0, label = reg.equation2, parse = T, size=6, show.legend = F)
dev.off()



round(mean(reg_all.df[reg_all.df$type=='ProE R1','nrl']))
## 194
round(mean(reg_all.df[reg_all.df$type=='OrthoE R1','nrl']))
## 192


