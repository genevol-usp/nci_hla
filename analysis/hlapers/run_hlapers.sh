#!/usr/bin/env bash

SAMPLE=$1
PHASE=1

HLADB=../../indices/hladb
BAM=/media/storage/genevol/vitor/nci/bam/${SAMPLE}_timepoint${PHASE}_Aligned.sortedByCoord.out.bam
OUT=./genotypes/${SAMPLE}_t${PHASE}

hlapers bam2fq -b $BAM -m $HLADB/mhc_coords.txt -o $OUT

#INDEX=../../indices/HLAPERS
#TRANSCRIPTS=$HLADB/transcripts_MHC_HLAsupp.fa
#FQ1=${OUT}_mhc_unmap_1.fq
#FQ2=${OUT}_mhc_unmap_2.fq
#
#$HOME/Libraries/HLApers/hlapers genotype -i $INDEX -t $TRANSCRIPTS \
#    -1 $FQ1 -2 $FQ2 -p $PBS_NUM_PPN -o $OUT
#
#rm -r $FQ1 $FQ2 ${OUT}_log
#
#FQS1=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R1*.fastq.gz )
#FQS2=$( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R2*.fastq.gz )
#
#TMPFQ1=/scratch/vitor/${SAMPLE}_t${PHASE}_R1.fq
#TMPFQ2=/scratch/vitor/${SAMPLE}_t${PHASE}_R2.fq
#
#zcat $FQS1 > $TMPFQ1 &
#zcat $FQS2 > $TMPFQ2
#wait
#
#GENOTYPES=${OUT}_genotypes.tsv
#OUTPREFIX=./quant/${SAMPLE}_t${PHASE}
#
#$HOME/Libraries/HLApers/hlapers quant --salmonreads -t $HLADB -g $GENOTYPES \
#    -1 $TMPFQ1 -2 $TMPFQ2 -o $OUTPREFIX -p $PBS_NUM_PPN
#
#rm $TMPFQ1 $TMPFQ2
