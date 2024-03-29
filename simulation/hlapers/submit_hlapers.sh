#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N HLApers
#PBS -t 1-50
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
BAM=$HOME/simulation/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam
MHC=../../indices/hladb/mhc_coords.txt
OUT=./genotypes/${SAMPLE}

$HOME/Libraries/HLApers/hlapers bam2fq -b $BAM -m $MHC -o $OUT

INDEX=../../indices/HLAPERS
TRANSCRIPTS=../../indices/hladb/transcripts_MHC_HLAsupp.fa
FQ1=./genotypes/${SAMPLE}_mhc_unmap_1.fq
FQ2=./genotypes/${SAMPLE}_mhc_unmap_2.fq

$HOME/Libraries/HLApers/hlapers genotype -i $INDEX -t $TRANSCRIPTS \
    -1 $FQ1 -2 $FQ2 -p $PBS_NUM_PPN -o $OUT

rm -r $FQ1 $FQ2 ${OUT}_log

GENOTYPES=./genotypes/${SAMPLE}_genotypes.tsv
HLADB=../../indices/hladb
FQ1=$HOME/simulation/fastq/${SAMPLE}/sample_01_1.fastq.gz
FQ2=$HOME/simulation/fastq/${SAMPLE}/sample_01_2.fastq.gz
OUTPREFIX=./quant/$SAMPLE

$HOME/Libraries/HLApers/hlapers quant --salmonreads -t $HLADB -g $GENOTYPES \
    -1 $FQ1 -2 $FQ2 -o $OUTPREFIX -p $PBS_NUM_PPN
