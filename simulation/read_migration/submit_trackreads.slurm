#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-50
#SBATCH --job-name=Track
#SBATCH -o /home/vitor/nci_hla/simulation/log/Track-reads-%a


SAMPLELIST=/home/vitor/nci_hla/simulation/samples.txt

SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
BED=/home/vitor/nci_hla/analysis/plot_data/hla.bed
BAM=/media/storage/genevol/vitor/simulation/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam
FQ=/media/storage/genevol/vitor/simulation/fastq/${SAMPLE}/sample_01_1.fastq.gz
TMP=/scratch/vitor/temp/${SAMPLE}_hla.sam
READS=/scratch/vitor/temp/${SAMPLE}_reads.txt
OUT=/scratch/vitor/temp/${SAMPLE}_hlareads.sam
FQOUT=/scratch/vitor/temp/${SAMPLE}_readids.txt

# simulated reads for HLA
samtools view $BAM |\
    awk '$1 ~ /read[0-9.]+_ENST[0-9.]+_[ABC]/' > $TMP

# other reads mapping to HLA
samtools view -L $BED $BAM |\
    awk '$1 !~ /read[0-9.]+_ENST[0-9.]+_[ABC]/' >> $TMP

awk '{ print $1 }' $TMP | sort | uniq > $READS

# union
samtools view $BAM | grep -F -f $READS > $OUT

# save their fastq header info
zcat $FQ | grep "^@" | grep -F -f $READS > $FQOUT

rm $TMP $READS
