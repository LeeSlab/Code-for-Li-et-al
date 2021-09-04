## calculate the TAD score

args=(commandArgs(TRUE))
file=as.character(args[1])
file2=as.character(args[2])
sample=as.character(args[3])

TADScore <- function(hic_matrix, TAD_boundaries, chr = "chr", resolution = 20000) {
  hic_matrix <- as.matrix(hic_matrix)
  total_bins <- dim(hic_matrix)[1]
  boundary_num <- dim(TAD_boundaries)[1]
  # define a TAD
  TAD_start <- TAD_boundaries$binEnd[-boundary_num] +1 # limf: +1把boundary的位置全部去掉
  TAD_end <- TAD_boundaries$binStart[-1]
  IntraTAD <- function(hic_matrix, start_bin, end_bin) {
    sub_TAD_Matrix <- hic_matrix[start_bin:end_bin, start_bin:end_bin]
    IntraTAD <- sum(sub_TAD_Matrix, na.rm = T) + sum(diag(sub_TAD_Matrix), na.rm = T)
    return(IntraTAD/2)
  }
  # limf: edit InterTAD calculation
  # InterTAD <- function(hic_matrix, start_bin, end_bin) {
  #   InterTAD <- sum(hic_matrix[start_bin:end_bin, ], na.rm = T)
  # }
  InterTAD <- function(hic_matrix, start_bin, end_bin) {
    InterTAD <- sum(hic_matrix[start_bin:end_bin, ], na.rm = T) - sum(hic_matrix[start_bin:end_bin, start_bin:end_bin], na.rm = T)
  }
  
  TAD_score <- numeric()
  for ( i in 1:length(TAD_start)) {
    intra_TAD_interactions <- IntraTAD(hic_matrix, start_bin = TAD_start[i], end_bin = TAD_end[i])
    inter_TAD_interactions <- InterTAD(hic_matrix, start_bin = TAD_start[i], end_bin = TAD_end[i])
    TAD_score <- c(TAD_score, intra_TAD_interactions/(intra_TAD_interactions + inter_TAD_interactions))
  }
  # limf:
  # output_bed <- cbind(chr = chr, chr_start = (TAD_start - 1) * resolution, chr_end = TAD_end * resolution, TAD_score = TAD_score)
  output_bed <- cbind(chr = chr, chr_start = TAD_start * resolution, 
                      chr_end = TAD_end * resolution, 
                      binStart = TAD_start,
                      binEnd = TAD_end,
                      TAD_score = TAD_score)
  output_bed <- as.data.frame(output_bed)
  return(output_bed)
}

filePath=paste0('terminal_cellCycle_LD_MergeAll/',sample,'/hic_results/matrix/',sample,'/iced/20000/')
# file='dense.matrix'
# file2='bin.bed'

inputMatrix <- read.delim(paste0(filePath,file),header=F,sep='\t')
inputTADBoundary <- read.delim(file2,header=F,sep='\t')
colnames(inputTADBoundary) <- c('chr','binStart','binEnd')
chromosome <- as.character(inputTADBoundary[1,1])

outPutBed <- TADScore(inputMatrix, inputTADBoundary, chromosome)
write.table(outPutBed,paste0(sample,'_20000_iced',chromosome,'_TADscore.bed'),
  sep='\t',quote=F,row.names=F,col.names=F)


