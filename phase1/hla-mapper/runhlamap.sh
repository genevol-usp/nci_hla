#!/usr/bin/env bash

SAMPLE=66K00003 
CPUS=8

BAM=/scratch/vitor/mapping/${SAMPLE}_Aligned.sortedByCoord.out.bam
HLAMAPPER=/home/ecastelli/hla-mapper/hla-mapper
DB=/home/ecastelli/hla-mapper/hla-mapper_db_004.1_HLA
BWA=/home/ecastelli/bwa/bwa-0.7.17/bwa
STAR=$HOME/Libraries/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
OUT=$HOME/nci_hla/phase1/hla-mapper/results/$SAMPLE

mkdir -p $OUT

$HLAMAPPER rna bam=$BAM sample=$SAMPLE db=$DB output=$OUT threads=$CPUS bwa=$BWA star=$STAR
