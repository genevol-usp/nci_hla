#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-50
#PBS -N coverage
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
BAM=/media/storage/genevol/vitor/bam/hlamapper/$SAMPLE/${SAMPLE}.adjusted.bam
MAPBAM=/scratch/vitor/${SAMPLE}.mappedReads.bam
SORTBAM=/scratch/vitor/${SAMPLE}.mappedReadsSorted.bam
BED=./hla.bed
OUT=./coverage/${SAMPLE}.cov

samtools view -b -@ $PBS_NUM_PPN -f 0x2 -F 0x100 $BAM > $MAPBAM

samtools sort -@ $PBS_NUM_PPN $MAPBAM -o $SORTBAM

samtools depth -a -m 0 -b $BED $SORTBAM > $OUT

rm $MAPBAM $SORTBAM
