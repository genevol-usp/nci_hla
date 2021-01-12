#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=64Gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1
#PBS -N hlamapper
#PBS -j oe
#PBS -o $PBS_JOBNAME.log

cd $PBS_O_WORKDIR

SAMPLE=$(awk "FNR==$PBS_ARRAYID" samples.txt)
CPUS=$PBS_NUM_PPN

BAM=/scratch/vitor/mapping/66K00003_Aligned.out.bam
HLAMAPPER=/home/ecastelli/hla-mapper/hla-mapper
DB=/home/ecastelli/hla-mapper/hla-mapper_db_004.1_HLA
BWA=/home/ecastelli/bwa/bwa-0.7.17/bwa
STAR=/home/vitor/Libraries/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
OUT=/home/vitor/nci_hla/phase1/hla-mapper/results/$SAMPLE

mkdir -p $OUT

$HLAMAPPER rna bam=$BAM sample=$SAMPLE db=$DB output=$OUT threads=$CPUS bwa=$BWA star=$STAR
