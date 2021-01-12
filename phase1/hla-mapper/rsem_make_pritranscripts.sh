#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -l walltime=01:00:00
#PBS -q short
#PBS -N RSEM-ref
#PBS -j oe
#PBS -o $PBS_O_WORKDIR/log/$PBS_JOBNAME

GENCODE=$HOME/gencode_data/v36
GENOME=$GENCODE/GRCh38.primary_assembly.genome.fa
GTF=$GENCODE/gencode.v36.primary_assembly.annotation.gtf
OUT=$GENCODE/gencode.v36.primary_assembly

zcat $GENOME.gz > $GENOME
zcat $GTF.gz > $GTF

rsem-prepare-reference --gtf $GTF $GENOME $OUT

rm $OUT.grp $OUT.ti $OUT.chrlist $OUT.seq $OUT.idx.fa $OUT.n2g.idx.fa 
rm $GENOME $GTF
