#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-50
#PBS -N STAR-map
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

STAR=$HOME/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
INDEX=../../indices/STAR
FQ1=$HOME/simulation/fastq/${SAMPLE}/sample_01_1.fastq.gz
FQ2=$HOME/simulation/fastq/${SAMPLE}/sample_01_2.fastq.gz
OUT=$HOME/simulation/bam/${SAMPLE}

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
    --outFileNamePrefix ${OUT}_

samtools index ${OUT}_Aligned.sortedByCoord.out.bam
