#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N StarIndexer
#PBS -j oe
#PBS -o star_index.log

cd $PBS_O_WORKDIR

STAR=/home/ecastelli/STAR/STAR-2.7.1a/bin/Linux_x86_64_static/STAR
GENOME=/home/vitor/gencode_data/v36/GRCh38.primary_assembly.genome.fa
GTF=/home/vitor/gencode_data/v36/gencode.v36.annotation.gtf
OUT=starindex

mkdir -p $OUT

zcat $GENOME.gz > $GENOME
zcat $GTF.gz > $GTF

$STAR --runThreadN $PBS_NUM_PPN \
    --runMode genomeGenerate \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF

rm $GENOME $GTF
