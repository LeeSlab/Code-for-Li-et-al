# ChIPseq pipeline:

# Read mapping:
bowtie2 -p 20 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x bowtie2.hg19.index -1 sample.1.fq.gz -2 sample.2.fq.gz -S sample_cutrun.sam
samtools view -bhS sample_cutrun.sam > sample_cutrun.bam
samtools sort -O bam sample_cutrun.bam -o sample.sorted.bam
picard MarkDuplicates I=sample.sorted.bam O=sample.sort_marked.bam REMOVE_DUPLICATES=true
samtools view -h -q 30 sample.sort_marked.bam | samtools sort -O BAM -@ 20 -o - > sample.sort_rmDup_q30.bam
samtools index sample.sort_rmDup_q30.bam

# Peak calling:
samtools view -h sample.bam > sample.sam

# 1. For transcription factors: 120bp
grep "@" sample.sam > sample.1_120.header
grep -v "@" sample.sam | awk ' $9 <= 120 && $9 >= 1 || $9 >= -120 && $9 <= -1 ' - >> sample.1_120.header
samtools view -S -b -o sample.1_120.bam sample.1_120.header

# 2. For histone modifications: 150-500bp
grep "@" sample.sam > sample.150_500.header
grep -v "@" sample.sam | awk ' $9 <= 500 && $9 >= 150 || $9 >= -500 && $9 <= -150 ' - >> sample.150_500.header
samtools view -S -b -o sample.150_500.bam sample.150_500.header

macs2 callpeak -t sample.bam -c input.bam -q 0.05 -f BAMPE -g hs -B --outdir MACS2/name -n name 2>./MACS2/name.macs2.log

