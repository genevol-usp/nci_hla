#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N Salmon-quant
#PBS -t 1-50
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
FQ1=/scratch/vitor/simulation/${SAMPLE}/sample_01_1.fastq.gz
FQ2=/scratch/vitor/simulation/${SAMPLE}/sample_01_2.fastq.gz
OUT=./quant_nocorrection/$SAMPLE
CPUS=$PBS_NUM_PPN 

salmon quant -i index -l A -1 $FQ1 -2 $FQ2 -p $CPUS -o $OUT