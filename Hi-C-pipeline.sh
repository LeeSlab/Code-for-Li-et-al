# Hi-C pipeline:

## HiC-Pro: https://github.com/nservant/HiC-Pro
HiC-Pro -c HiCProConfForterminal.txt -i inputFileName -o outputName

sort -R sample_allValidPairs -T tmp > tmp/sample_allValidPairs.random
head -n 584592707 tmp/sample_allValidPairs.random > sample_allValidPairs.consistent
## 584592707 for four samples

HiC-Pro_2.10.0/bin/utils/hicpro2juicebox.sh -i ${sample} -g hg19 -j ./software/Juicer_CM_tools/juicer_tools.jar
bash cal_resolution.sh ${sample}

## Call TAD：
## cword: https//github.com/dekkerlab/cworld-dekker
i=1..22,X
perl matrix2insulation.pl -i input_matrix -is 1000000 -ids 200000 -im mean -bmoe 3 -nt 0.25 -v
perl insulation2tad.pl -i sample_chr${i}.is1000001.ids200001.insulation -b sample_chr${i}.is1000001.ids200001.insulation.boundaries -o sample_20000_TAD_chr${i} --mbs 0.25
perl matrix2compartment.pl -i sample_20000_iced_${chro}.input.matrix --et


## homer http://homer.ucsd.edu/homer/interactions2/index.html
makeTagDirectory ${TagFile} -format HiCsummary ${validPair} -tbp 1 -genome hg19
echo "makeTagDirectory finshed"
runHiCpca.pl ${pcaout} ${TagFile} -res 50000 -window 100000 -genome hg19 -cpu 4
echo "runHiCpca finshed"
findTADsAndLoops.pl find ${TagFile} -cpu 4 -res 3000 -window 15000 -genome hg19
echo "findTADsAndLoops finshed"
findHiCCompartments.pl ${pcaout}.PC1.txt > ${pcaout}.compartments.txt 
echo "findHiCCompartments.pl finshed"
findHiCCompartments.pl ${pcaout}.PC1.txt -opp > ${pcaout}.compartments.opp.txt
analyzeHiC ${TagFile} -res 5000 -window 15000 -nomatrix -compactionStats auto -cpu 10

## arrowhead:https://github.com/aidenlab/juicer/wiki/Arrowhead
java -jar juicer_tools_1.11.04_jcuda.0.8.jar arrowhead -c -k KR -m 2000 -r 10000 ${hic_file} arrowhead/${name}_contact_list


## pairQC: https://github.com/4dn-dcic/pairsqc
## boundary strength
awk '{FS="\t";OFS="\t"} {print $1,$2,$3,$5,$6,$4,$7}' ${sample}_validPairs > ${sample}.4DNpairs
sort -k2,2 -k4,4 -k3,3n -k5,5n ${sample}.4DNpairs -T ./tmp/ > ${sample}.4DNpairs.sort
rm ${sample}.4DNpairs
./bgzip ${sample}.4DNpairs.sort
./pairix -p pairs ${sample}.4DNpairs.sort.gz
python ./pairsqc/pairsqc.py -p ${sample}.4DNpairs.sort.gz -c ./hg19.chrom.sizes -t P -O ${name}
Rscript ./pairsqc/plot.r 4 ${name}_report/

## 
