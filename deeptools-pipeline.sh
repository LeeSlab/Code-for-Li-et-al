# deeptools pipeline:

## you can process bwfile form rawdata from GSA

## bwfile = CTCF/SMC3 files
## bedFile = example-peak-bed

computeMatrix reference-point --referencePoint center -b 3000 -a 3000 -S ${bwFile}/sample1.bw s${bwFile}/ample2.bw -R ${example-peak-bed} --binSize 10 --skipZeros -o sample.matrix.gz
plotProfile -m sample.matrix.gz -out ${outName}_Profile.pdf --samplesLabel name1 name2 --dpi 300 --yMin ${YMIN} --yMax ${YMAX}



## CTCF/SMC3 colocalize with GATA1
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 --binSize 10 -R ${example-peak-bed}/GATA1.peak.bed -S ${bwFile}/ProCTCF_R1.bw ${bwFile}/ProCTCF_R2.bw --skipZeros -o matrix.gz
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 --binSize 10 -R ${example-peak-bed}/GATA1.peak.bed -S ${bwFile}/OrthoCTCF_R1.bw ${bwFile}/OrthoCTCF_R2.bw --skipZeros -o matrix.gz 

computeMatrix reference-point --referencePoint center -b 3000 -a 3000 --binSize 10 -R ${example-peak-bed}/GATA1.peak.bed -S ${bwFile}/Pro_SMC3_R1.bw ${bwFile}/Pro_SMC3_R2.bw --skipZeros -o matrix.gz 
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 --binSize 10 -R ${example-peak-bed}/GATA1.peak.bed -S ${bwFile}/Ortho_SMC3_R1.bw ${bwFile}/Ortho_SMC3_R2.bw --skipZeros -o matrix.gz

plotProfile -m ${matrixFile} -out ${outName}_Profile.pdf --samplesLabel ${Name}R1 ${Name}R2 --dpi 300



