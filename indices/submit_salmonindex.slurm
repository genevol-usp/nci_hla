#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --job-name=Salmon-index
#SBATCH -o /home/vitor/nci_hla/simulation/log/Salmon-index

cd $SLURM_SUBMIT_DIR

salmon index -t gencode.transcripts.fa -i SALMON
