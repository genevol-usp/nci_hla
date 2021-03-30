#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=8gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-50
#PBS -N featureCounts
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
GTF=$HOME/gencode/gencode.v37.primary_assembly.annotation.gtf
BAM=$HOME/simulation/bam_hlamapper/${SAMPLE}/${SAMPLE}.adjusted.bam

featureCounts -p -B -C -M --primary --fraction\
    -T $PBS_NUM_PPN \
    -g gene_id \
    -a $GTF \
    -R SAM \
    -o ./quant/${SAMPLE}_counts.txt \
    $BAM

SAM=./quant/$( basename $BAM).featureCounts.sam
MAPREADS=./quant/${SAMPLE}_hlaMappedReads

grep -E read\[0-9\]+_ENST\[0-9.\]+_ $SAM > ./quant/${SAMPLE}_simulReads.sam
    
grep -F -f ./hla_geneids.txt $SAM |\
    awk '$16 ~ /Assigned/ { print $1 }' |\
    sort |\
    uniq > $MAPREADS

grep -F -f $MAPREADS $SAM |\
    awk '$16 ~ /Assigned/' > ./quant/${SAMPLE}_hlaMapped.sam

rm $SAM $MAPREADS
