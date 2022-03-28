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
GTF=/raid/genevol/users/vitor/nci_hla/simulation/hla-mapper/gencode.v37.primary_assembly.annotation.EDIT.gtf
OUT=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t${PHASE} 

featureCounts -p -B -C -M --primary --fraction\
    -T $PBS_NUM_PPN\
    -g gene_id \
    -a $GTF \
    -o ./quant/${SAMPLE}_t${PHASE}_gene_rmPseudo.txt \
    $OUT/${SAMPLE}.adjusted.bam
