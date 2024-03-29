#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-50
#PBS -N hla-mapper
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
GTF=/raid/genevol/users/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf
OUT=/raid/genevol/users/vitor/simulation/bam_hlamapper/${SAMPLE}

mkdir -p $OUT

/raid/genevol/users/vitor/hla-mapper_new/build/hla-mapper rna \
    bam=/raid/genevol/users/vitor/simulation/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    sample=$SAMPLE \
    threads=$PBS_NUM_PPN \
    db=/raid/genevol/users/vitor/hla-mapper_new/db/hla-mapper_db_004.1_HLA \
    output=$OUT \
    bwa=/raid/genevol/users/vitor/HLAMAPPER/bwa/bwa-0.7.17/bwa \
    samtools=/raid/genevol/users/vitor/HLAMAPPER/samtools-1.11/samtools \
    star=/raid/genevol/users/vitor/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

BAM=$OUT/${SAMPLE}.adjusted.bam

featureCounts -p -B -C -M --primary --fraction\
    -T $PBS_NUM_PPN \
    -g gene_id \
    -a $GTF \
    -R SAM \
    -o ./quant/${SAMPLE}_gene.txt \
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
