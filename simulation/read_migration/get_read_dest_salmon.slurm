#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-50
#SBATCH --job-name=ReadDest
#SBATCH -o /home/vitor/nci_hla/simulation/log/Read-dest-Salmon-%a


SAMPLELIST=/home/vitor/nci_hla/simulation/samples.txt

SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

cd $SLURM_SUBMIT_DIR

Rscript make_read_dest_salmon.R $SAMPLE
