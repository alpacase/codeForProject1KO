#!/bin/sh

# QC and trimming
trim_galore --cores 10 --paired --dont_gzip -o ./ *.fq
mkdir report_trimGalore
mv *.txt ./report_trimGalore

# mapping and bam sort
mkdir report_bam
mkdir bam

for fq in $(ls *_1_val_1.fq); do
        filename="${fq%_1_val_1.fq}"
        echo $filename
        hisat2 -p 8 --dta -x ~/ngs_analysis/resource/ht2Build_mm10/mm10 -1 ${filename}_1_val_1.fq -2 ${filename}_2_val_2.fq |
                samtools view -bSu -F 4 - |
                pv |
                samtools sort -@ 4 -o ./bam/${filename}.bam -
        samtools index ./bam/${filename}.bam
        samtools idxstats ./bam/${filename}.bam >./report_bam/${filename}_stats.txt
done

mkdir fastq_trimmed
mv *_1_val_1.fq *_2_val_2.fq ./fastq_trimmed
mkdir fastq_raw
mv *.fq fastq_raw

# create bigwig files
for bam in $(ls *.bam); do
        name=${bam%%.*}
        bamCoverage -b $bam -o ../bigwig/${name}.bw --normalizeUsing RPKM -p 10 -v -bs 1
done

# homer
# makeTagDirectory
for f in $(ls *.bam); do
        filename="${f%.*}"
        makeTagDirectory ../homer/${filename} $f
done

# makeBedGraph
sample=$(
        cat <<__EOF__
d0_1KO
d0_CI6
d0_WT
d2_1KO
d2_CI6
d2_WT
d4_1KO
d4_CI6
d4_WT
__EOF__
)

for i in ${sample[@]}; do
        makeUCSCfile $i -fragLength given -o auto
done

# bedGraph to BigWig
for i in $(ls *.bedGraph); do
        name=${i%.*}
        bedGraphToBigWig $i mm10.chrom.sizes.txt ${name}.bw
done

# calculating gene expression
analyzeRepeats.pl rna mm10 -strand both -count exons -d */ -rpkm >expression_rpkm.txt

# calculating repeat expression
analyzeRepeats.pl reference/mm10/rmsk.gtf mm10 -d */ -rpkm >expression_repeat_rpkm.txt
