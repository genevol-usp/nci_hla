#!/usr/bin/env bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1
#PBS -N StarMap
#PBS -j oe
#PBS -o star_map.log

cd $PBS_O_WORKDIR

STAR=/home/vitor/Libraries/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
SAMPLE=$( awk "FNR==$PBS_ARRAYID" samples.txt )
INDEX=starindex
FQ1=$( ls -v /raid/genevol/nci_rnaseq/phase1/fastq/${SAMPLE}*R1_001.fastq.gz | paste -s -d, - )
FQ2=$( ls -v /raid/genevol/nci_rnaseq/phase1/fastq/${SAMPLE}*R2_001.fastq.gz | paste -s -d, - )
OUT=/scratch/vitor/mapping/${SAMPLE}_

$STAR --runMode alignReads \
    --runThreadN $PBS_NUM_PPN \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMultimapNmax 20 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --outFileNamePrefix $OUT

