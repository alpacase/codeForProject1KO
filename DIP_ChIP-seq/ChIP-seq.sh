#!/bin/sh

# ChIP-seq
# download fastq file
SRR=$(
  cat <<srr
SRR5077640
SRR5077644
SRR5077648
SRR5077660
SRR5077672
srr
)
for srr in ${SRR[@]}; do
  fasterq-dump -p -v -e 6 -t ./ -O ./ $srr
done

# quality check & trimming
trim_galore --no_report_file --cores 8 *.fq

# mapping
for fq in $(ls *.fq); do
  bamfile="${fq%.*}"
  bowtie2 -p 10 -x resource/bt2Build_mm10/mm10 -U ${bamfile}.fq |
    samtools view -bSu -F 4 -@ 10 - |
    pv |
    samtools sort -@ 10 -o ${bamfile}_sorted.bam -
  samtools index -b ${bamfile}_sorted.bam
done

# peak calling
# narrow peak
macs2 callpeak -f BAM -g mm -q 0.01 -c GSM2417124_ESC_ChIPinput_sorted.bam -n GSM2417096_ESC_H3K27ac -t GSM2417096_ESC_H3K27ac_sorted.bam
macs2 callpeak -f BAM -g mm -q 0.01 -c GSM2417124_ESC_ChIPinput_sorted.bam -n GSM2417092_ESC_H3K9ac -t GSM2417092_ESC_H3K9ac_sorted.bam

# broad peak
macs2 callpeak -f BAM -g mm -q 0.01 --broad -c GSM2417124_ESC_ChIPinput_sorted.bam -n GSM2417100_ESC_H3K27me3 -t GSM2417100_ESC_H3K27me3_sorted.bam
macs2 callpeak -f BAM -g mm -q 0.01 --broad -c GSM2417124_ESC_ChIPinput_sorted.bam -n GSM2417112_ESC_H3K9me3 -t GSM2417112_ESC_H3K9me3_sorted.bam

# ATAC-seq
# download fastq file
prefetch SRR5077752
fasterq-dump -v SRR5077752/*.sra

# QC and trimming
trim_galore --no_report_file --cores 8 *.fastq

# mapping
bowtie2 -p 6 -x ~/ngs/ref_genome/bt2Build_mm10/mm10 -U *_trimmed.fq |
  samtools view -bSu -F 4 -@ 6trim - |
  pv |
  samtools sort -@ 6 -o atac-seq.bam -
samtools index -b atac-seq.bam

# peak calling
macs2 callpeak -f BAM -g mm -q 0.01 --nomodel --call-summits --extsize 150 --shift -75 -n atac-seq -t atac-seq.bam
