#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=12gb
#SBATCH --time=48:00:00
#SBATCH -p long
#SBATCH --array=1-96
#SBATCH --job-name=SPLICE
#SBATCH -o /home/vitor/nci_hla/analysis/splicing/log/SPLICE-%a

# IO
SAMPLELIST=/home/vitor/nci_hla/analysis/samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SAMPLELIST )

FQ1=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t1/${SAMPLE}_HLA-A_R1.fastq
FQ2=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t1/${SAMPLE}_HLA-A_R2.fastq

INDEX=${SLURM_SUBMIT_DIR}/STARINDEX
OUT=/scratch/vitor/splicing/${SAMPLE}

# Mapping
#--twopassMode Basic \
#--readFilesCommand -  \

STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.99 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --outFileNamePrefix ${OUT}_

BAM=${OUT}_Aligned.sortedByCoord.out.bam
samtools index -@ $SLURM_CPUS_PER_TASK $BAM

JUNC=/scratch/vitor/splicing/${SAMPLE}.junc 

# Extract junctions
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $BAM -o $JUNC

