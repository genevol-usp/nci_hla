#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-50
#SBATCH --job-name=ReadOrigins
#SBATCH -o /home/vitor/nci_hla/simulation/log/Read-origin-salmon-%a


SAMPLELIST=/home/vitor/nci_hla/simulation/samples.txt

SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
REFTX=/home/vitor/nci_hla/indices/hla_ref_transcripts.txt
BAM=/home/vitor/nci_hla/simulation/salmon/quant_salmon.v1.8/${SAMPLE}/mapping.bam
FQ=/media/storage/genevol/vitor/simulation/fastq/${SAMPLE}/sample_01_1.fastq.gz
TMP=/scratch/vitor/temp/${SAMPLE}_tmp_salmon.txt
READS=/scratch/vitor/temp/${SAMPLE}_reads_salmon.txt
OUT=/scratch/vitor/temp/${SAMPLE}_hlareads_salmon.sam
FQOUT=/scratch/vitor/temp/${SAMPLE}_readids_salmon.txt

# simulated reads for HLA
samtools view $BAM |\
    awk '$1 ~ /read[0-9.]+_ENST[0-9.]+_[ABC]/ { print $1 }' > $TMP

# other reads mapping to HLA
samtools view $BAM |\
    awk '$1 !~ /read[0-9.]+_ENST[0-9.]+_[ABC]/' |\
    grep -F -f $REFTX |\
    awk '{ print $1 }' >> $TMP

awk '{ print $1 }' $TMP | sort | uniq > $READS

# union
samtools view $BAM | grep -F -f $READS > $OUT

# save their fastq header info
zcat $FQ | grep "^@" | grep -F -f $READS > $FQOUT

rm $TMP $READS
