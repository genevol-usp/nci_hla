#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-96
#SBATCH --job-name=SamDepth
#SBATCH -o /home/vitor/log/samtools-depth.o%A.%a


SAMPLELIST=${SLURM_SUBMIT_DIR}/sample_ids_t1.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

BED=${SLURM_SUBMIT_DIR}/plot_data/hla.bed
BAM=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t1/${SAMPLE}.adjusted.bam
OUT=${SLURM_SUBMIT_DIR}/plot_data/coverage_${SAMPLE}.txt

samtools depth -a -m 1000000 -b $BED $BAM > $OUT
