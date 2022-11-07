#!/bin/bash

#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-96
#SBATCH --job-name=samtools
#SBATCH -o /scratch/vitor/log/getmapped-%a

cd $SLURM_SUBMIT_DIR

PHASE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ./samples.txt )
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' ./samples.txt )

BAM=/media/storage/genevol/vitor/nci/bam/${SAMPLE}_timepoint${PHASE}_Aligned.sortedByCoord.out.bam
READS=/scratch/vitor/results/${SAMPLE}_mapped.txt

samtools view -@ 2 -F 0x4 $BAM | cut -f 1 > $READS
Rscript count_reads.R $READS
