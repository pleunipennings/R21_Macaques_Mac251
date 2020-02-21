#!/bin/bash

#This script for the data from Mac251 R21 project

mkdir ProcessedData/SAMPLE

#0 adapter trimming
~/Dropbox/bbmap/bbduk.sh in1=FASTQFile1 in2=FASTQFile2  out=ProcessedData/SAMPLE/SAMPLE_adp.trimmed.fastq ref=~/Dropbox/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=ProcessedData/SAMPLE/stats30_0.txt

#1 Trim reads at both ends at Score<30
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/SAMPLE/SAMPLE_adp.trimmed.fastq out=ProcessedData/SAMPLE/SAMPLE_trimmed.q30.fastq qtrim=rl trimq=30 stats=ProcessedData/SAMPLE/stats30_1.txt

#2. Kmer filtering
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/SAMPLE/SAMPLE_trimmed.q30.fastq out=ProcessedData/SAMPLE/SAMPLE_unmatched.q30.fq ref=~/Dropbox/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=ProcessedData/SAMPLE/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/SAMPLE/SAMPLE_unmatched.q30.fq out=ProcessedData/SAMPLE/SAMPLE_clean.q30.fq maq=30 stats=ProcessedData/SAMPLE/stats30_3.txt

#4. Align the file using bwa to the reference 
~/Dropbox/bwa-0.7.17/bwa index -p SIV -a is OriginalData/SIV_Env.fasta

~/Dropbox/bwa-0.7.17/bwa mem -t 4 -k 15 -a SIV ProcessedData/SAMPLE/SAMPLE_clean.q30.fq  > ProcessedData/SAMPLE/SAMPLE_BWAmapped.sam

#5. convert sam to bam
~/src/samtools/samtools view -S -b ProcessedData/SAMPLE/SAMPLE_BWAmapped.sam > ProcessedData/SAMPLE/SAMPLE_BWAmapped.bam


#6. sort the bam file
~/src/samtools/samtools sort  ProcessedData/SAMPLE/SAMPLE_BWAmapped.bam -o  ProcessedData/SAMPLE/SAMPLE_BWA_sort.bam

#7. index the bam file
~/src/samtools/samtools index  ProcessedData/SAMPLE/SAMPLE_BWA_sort.bam  ProcessedData/SAMPLE/SAMPLE_BWA_sort.bam.bai

rm -r ProcessedData/SAMPLE/SAMPLE_BWAmapped.bam
rm -r ProcessedData/SAMPLE/SAMPLE_BWAmapped.sam
rm -r ProcessedData/SAMPLE/SAMPLE_clean.q30.fq
rm -r ProcessedData/SAMPLE/SAMPLE_unmatched.q30.fq
rm -r ProcessedData/SAMPLE/SAMPLE_trimmed.q30.fastq
rm -r ProcessedData/SAMPLE/SAMPLE_adp.trimmed.fastq

