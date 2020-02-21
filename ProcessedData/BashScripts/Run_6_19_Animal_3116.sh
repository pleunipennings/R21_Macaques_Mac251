#!/bin/bash

#This script for the data from Mac251 R21 project

mkdir ProcessedData/Run_6_19_Animal_3116

#0 adapter trimming
~/Dropbox/bbmap/bbduk.sh in1=/Volumes/PleuniBackup/MacaquesSIVmac251/mac251_Run_6-138794658/FASTQ_Generation_2019-07-22_03_39_46Z-190187015/19_L001-ds.c641f57035e84fb6a05e814d1dc87322/19_S18_L001_R1_001.fastq.gz in2=/Volumes/PleuniBackup/MacaquesSIVmac251/mac251_Run_6-138794658/FASTQ_Generation_2019-07-22_03_39_46Z-190187015/19_L001-ds.c641f57035e84fb6a05e814d1dc87322/19_S18_L001_R2_001.fastq.gz  out=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_adp.trimmed.fastq ref=~/Dropbox/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=ProcessedData/Run_6_19_Animal_3116/stats30_0.txt

#1 Trim reads at both ends at Score<30
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_adp.trimmed.fastq out=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_trimmed.q30.fastq qtrim=rl trimq=30 stats=ProcessedData/Run_6_19_Animal_3116/stats30_1.txt

#2. Kmer filtering
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_trimmed.q30.fastq out=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_unmatched.q30.fq ref=~/Dropbox/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=ProcessedData/Run_6_19_Animal_3116/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
~/Dropbox/bbmap/bbduk.sh in=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_unmatched.q30.fq out=ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_clean.q30.fq maq=30 stats=ProcessedData/Run_6_19_Animal_3116/stats30_3.txt

#4. Align the file using bwa to the reference 
~/Dropbox/bwa-0.7.17/bwa index -p SIV -a is OriginalData/SIV_Env.fasta

~/Dropbox/bwa-0.7.17/bwa mem -t 4 -k 15 -a SIV ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_clean.q30.fq  > ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.sam

#5. convert sam to bam
~/src/samtools/samtools view -S -b ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.sam > ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.bam


#6. sort the bam file
~/src/samtools/samtools sort  ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.bam -o  ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWA_sort.bam

#7. index the bam file
~/src/samtools/samtools index  ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWA_sort.bam  ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWA_sort.bam.bai

rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.bam
rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_BWAmapped.sam
rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_clean.q30.fq
rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_unmatched.q30.fq
rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_trimmed.q30.fastq
rm -r ProcessedData/Run_6_19_Animal_3116/Run_6_19_Animal_3116_adp.trimmed.fastq

