#!/usr/bin/env bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N STAR-index
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

STAR=$HOME/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
GENOME=$HOME/gencode/GRCh38.primary_assembly.genome.fa
GTF=$HOME/gencode/gencode.v37.primary_assembly.annotation.gtf
OUT=STAR

mkdir -p $OUT

$STAR --runThreadN $PBS_NUM_PPN \
    --runMode genomeGenerate \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF
