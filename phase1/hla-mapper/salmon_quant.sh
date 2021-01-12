#!/usr/bin/env bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1
#PBS -N salmon
#PBS -j oe
#PBS -o $PBS_O_WORKDIR/log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" samples.txt )
BAM1=/scratch/vitor/mapping/${SAMPLE}_Aligned.toTranscriptome.out.bam
BAM2=/scratch/vitor/mapping/${SAMPLE}_Aligned.transcriptCoords.out.bam
TRANSCRIPTS=$HOME/gencode_data/v36/gencode.v36.primary_assembly.transcripts.fa
OUT1=./quant/STAR/${SAMPLE}
OUT2=./quant/xlate/${SAMPLE}

mkdir -p $OUT1
mkdir -p $OUT2

salmon quant -t $TRANSCRIPTS -l IU -a $BAM1 -p $PBS_NUM_PPN \
    -o $OUT1 \ 
    --seqBias --gcBias --posBias


salmon quant -t $TRANSCRIPTS -l IU -a $BAM2 -p $PBS_NUM_PPN \
    -o $OUT2 \
    --seqBias --gcBias --posBias
