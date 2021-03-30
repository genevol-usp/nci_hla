#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N Salmon-quant
#PBS -t 1-107
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

INDEX=../../indices/SALMON
CPUS=$PBS_NUM_PPN 

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )

FQ1=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R1*.fastq.gz )
FQ2=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R2*.fastq.gz )

OUT=./quant/${SAMPLE}_t${PHASE}

salmon quant -i $INDEX -l A -1 <( zcat $FQ1 ) -2 <( zcat $FQ2 ) -p $CPUS -o $OUT \
    --seqBias --posBias --gcBias
