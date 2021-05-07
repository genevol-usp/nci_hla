#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=12gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-96
#PBS -N featureCounts
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )
OUT=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t${PHASE} 

featureCounts -p -B -M -O -f --primary \
    -T $PBS_NUM_PPN \
    -a ./gencode.v37.primary_assembly.annotation_HLAedit.gtf \
    -o ./quant/${SAMPLE}_t${PHASE}_exon.txt \
    $OUT/${SAMPLE}.adjusted.bam

