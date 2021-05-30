#!/bin/sh

#QC & trimming


for fq in `ls *_1.fq`
do
  fastq="${fq%_1.*}";
  trim_galore --no_report_file --cores 8 -o -r1 ${fastq}_1.fq -r2 ${fastq}_2.fq;
done

#mapping

idx=${0%/*}
for fq in `ls *_1.fq`
do
  bamfile="${fq%_1.*}";
  echo $bamfile ;
  bowtie2 -p 8 -x ${idx}/bt2Build_mm10/mm10 -1 ${bamfile}_1.fq -2 ${bamfile}_2.fq \
  | samtools view -bSu -F 4 - \
  | pv \
  | samtools sort -@ 4 -o ${bamfile}.bam -
  samtools index -b ${bamfile}.bam
done

#convert wig generated by MEDIPS to bigwig

chrom=${0%/*}/mm10.chrom.sizes.txt
for wig in `ls *.wig`; do
  bw=${wig%.*}.bw
  wigToBigWig $wig $chrom $bw
done

#computeMatrix

#5mC
computeMatrix reference-point -p 10 \
-bs 100 -a 2500 -b 2500 --referencePoint center \
--skipZeros \
-o 5mC_narrowPeak.txt.gz \
-R ~/ngs/reference/GSE90893_ESCs_histone_ATAC-seq/*.narrowPeak \
~/ngs/reference/GSE90893_ESCs_histone_ATAC-seq/*.broadPeak \
-S ~/ngs/reference/GSE72481_WT_DIP/WT_5mC.bw \
~/ngs/mydata/bigwig/1KO_5mC_G1.bw \
~/ngs/mydata/bigwig/1KO_5mC_G2.bw \

#5hmC
computeMatrix reference-point -p 10 \
-bs 100 -a 2500 -b 2500 --referencePoint center \
--skipZeros \
-o 5hmC_narrowPeak.txt.gz \
-R ~/ngs/reference/GSE90893_ESCs_histone_ATAC-seq/*.narrowPeak \
~/ngs/reference/GSE90893_ESCs_histone_ATAC-seq/*.broadPeak \
-S ~/ngs/reference/GSE72481_WT_DIP/WT_5hmC.bw \
~/ngs/mydata/bigwig/1KO_5hmC_G1.bw \
~/ngs/mydata/bigwig/1KO_5hmC_G2.bw \