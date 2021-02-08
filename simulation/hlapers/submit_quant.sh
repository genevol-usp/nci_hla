#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=12gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N HLApers
#PBS -t 1-50
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
GENOTYPES=./genotypes/${SAMPLE}.tsv 
FQ1=/scratch/vitor/simulation/${SAMPLE}/sample_01_1.fastq.gz
FQ2=/scratch/vitor/simulation/${SAMPLE}/sample_01_2.fastq.gz
OUTPREFIX=./quant/$SAMPLE
CPUS=$PBS_NUM_PPN 

./hlapers quant -t hladb -g $GENOTYPES -1 $FQ1 -2 $FQ2 -o $OUTPREFIX -p $CPUS 
