#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N HLApers-prepref
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

ANNOT=/home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf
TRANSCRIPTS=../salmon-reads/ensembl.transcripts.fa

hlapers prepare-ref -t $TRANSCRIPTS -a $ANNOT -i /home/vitor/IMGTHLA -o hladb 
