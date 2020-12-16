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

FASTQ1=( $(ls -v /raid/genevol/nci_rnaseq/phase1/fastq/${SAMPLE}*R1_001.fastq.gz) )
FASTQ2=( $(ls -v /raid/genevol/nci_rnaseq/phase1/fastq/${SAMPLE}*R2_001.fastq.gz) )
OUT=/home/vitor/nci_hla/phase1/hla-mapper/results/$SAMPLE
DB=/home/ecastelli/hla-mapper/hla-mapper_db_004.1_HLA
HLAMAPPER=/home/ecastelli/hla-mapper/hla-mapper
BWA=/home/ecastelli/bwa/bwa-0.7.17/bwa
STAR=/home/vitor/Libraries/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

mkdir -p $OUT

if [ "${#FASTQ1[@]}" -eq 1 ]; then
    
    $HLAMAPPER rna r1=$FASTQ1 r2=$FASTQ2 sample=$SAMPLE db=$DB output=$OUT threads=$CPUS bwa=$BWA star=$STAR

elif [ "${#FASTQ1[@]}" -gt 1 ]; then

    FQTMP1=/scratch/vitor/${SAMPLE}_R1.fastq
    FQTMP2=/scratch/vitor/${SAMPLE}_R2.fastq

    zcat ${FASTQ1[@]} > $FQTMP1 &
    zcat ${FASTQ2[@]} > $FQTMP2
    wait
    
    $HLAMAPPER rna r1=$FQTMP1 r2=$FQTMP2 sample=$SAMPLE db=$DB output=$OUT threads=$CPUS bwa=$BWA star=$STAR
    
    rm $FQTMP1 $FQTMP2

else
    echo "please specify fastq files correctly"
fi

