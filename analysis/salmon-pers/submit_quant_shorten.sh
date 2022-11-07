#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p long
#SBATCH --array=1-96
#SBATCH --job-name=Salmon
#SBATCH -o /home/vitor/nci_hla/analysis/salmon-pers/log/Salmon-short-%a

cd $SLURM_SUBMIT_DIR

SALMON=/home/vitor/Libraries/salmon-1.8.0_linux_x86_64/bin/salmon

PHASE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' ../samples.txt )

FQ1=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R1*.fastq.gz )
FQ2=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R2*.fastq.gz )

OUT=./quant-shorten/${SAMPLE}_t${PHASE}
CPUS=$SLURM_CPUS_PER_TASK

NOHLA=../../indices/transcripts_noABC.fa
HLA=$OUT/hla.fa
FASTA=$OUT/index.fa
INDEX=$OUT/INDEX

cat $NOHLA $HLA > $FASTA

$SALMON index -t $FASTA -i $INDEX

$SALMON quant -i $INDEX -l A -1 $FQ1 -2 $FQ2 -p $CPUS -o $OUT \
    --seqBias --posBias --gcBias

rm -r $HLA $FASTA $INDEX
