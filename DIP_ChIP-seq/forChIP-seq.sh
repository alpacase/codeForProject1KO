#!/bin/sh

# downroad fastq file
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
trim_galore --no_report_file --cores 8 -o *.fq

# mapping
for fq in $(ls *.fq); do
  bamfile="${fq%.*}"
  bowtie2 -p 10 -x resource/bt2Build_mm10/mm10 \
    -U ${bamfile}.fq |
    samtools view -bSu -F 4 -@ 10 - |
    pv |
    samtools sort -@ 10 -o ${bamfile}_sorted.bam -
  samtools index -b ${bamfile}_sorted.bam
done

# peak calling
# narrow peak
macs2 callpeak -f BAM -g mm -q 0.01 \
  -c GSM2417124_ESC_ChIPinput_sorted.bam \
  -n GSM2417096_ESC_H3K27ac \
  -t GSM2417096_ESC_H3K27ac_sorted.bam
macs2 callpeak -f BAM -g mm -q 0.01 \
  -c GSM2417124_ESC_ChIPinput_sorted.bam \
  -n GSM2417092_ESC_H3K9ac \
  -t GSM2417092_ESC_H3K9ac_sorted.bam

# broad peak
macs2 callpeak -f BAM -g mm -q 0.01 --broad \
  -c GSM2417124_ESC_ChIPinput_sorted.bam \
  -n GSM2417100_ESC_H3K27me3 \
  -t GSM2417100_ESC_H3K27me3_sorted.bam
macs2 callpeak -f BAM -g mm -q 0.01 --broad \
  -c GSM2417124_ESC_ChIPinput_sorted.bam \
  -n GSM2417112_ESC_H3K9me3 \
  -t GSM2417112_ESC_H3K9me3_sorted.bam

# ATAC-seq
wget 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2417nnn/GSM2417076/suppl/GSM2417076_ESC_ATAC.bed.gz'
tar -xzf GSM2417076_ESC_ATAC.bed.gz
