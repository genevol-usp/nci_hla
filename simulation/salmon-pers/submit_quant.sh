#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N Salmon-pers-quant
#PBS -t 1-50
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
FQ1=$HOME/simulation/fastq/${SAMPLE}/sample_01_1.fastq.gz
FQ2=$HOME/simulation/fastq/${SAMPLE}/sample_01_2.fastq.gz
OUT=./quant/$SAMPLE
CPUS=$PBS_NUM_PPN 

NOHLA=../../indices/transcripts_noABC.fa
HLA=$OUT/hla.fa
FASTA=$OUT/index.fa
INDEX=$OUT/INDEX

cat $NOHLA $HLA > $FASTA

salmon index -t $FASTA -i $INDEX

salmon quant -i $INDEX -l A -1 $FQ1 -2 $FQ2 -p $CPUS -o $OUT \
    --seqBias --posBias --gcBias

rm -r $HLA $FASTA $INDEX
