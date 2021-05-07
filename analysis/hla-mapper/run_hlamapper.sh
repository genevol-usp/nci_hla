#!/usr/bin/env bash

PBS_ARRAYID=$1

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )
BAM=/media/storage/genevol/vitor/nci/bam/${SAMPLE}_timepoint${PHASE}_Aligned.sortedByCoord.out.bam
OUT=/media/storage/genevol/vitor/nci/bam_hlamapper/${SAMPLE}_t${PHASE} 

mkdir -p $OUT

/raid/genevol/users/vitor/hla-mapper_new/build/hla-mapper rna \
    bam=$BAM \
    sample=$SAMPLE \
    threads=12 \
    db=/raid/genevol/users/vitor/hla-mapper_new/db/hla-mapper_db_004.1_HLA \
    output=$OUT \
    bwa=/raid/genevol/users/vitor/HLAMAPPER/bwa/bwa-0.7.17/bwa \
    samtools=/raid/genevol/users/vitor/HLAMAPPER/samtools-1.11/samtools \
    star=/raid/genevol/users/vitor/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

featureCounts -p -B -C -M --primary --fraction\
    -T 12 \
    -g gene_id \
    -a /home/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf \
    -o ./quant/${SAMPLE}_t${PHASE}_gene.txt \
    $OUT/${SAMPLE}.adjusted.bam

featureCounts -p -B -M -O -f --primary --fraction \
    -T 12 \
    -a ./gencode.v37.primary_assembly.annotation_HLAedit.gtf \
    -o ./quant/${SAMPLE}_t${PHASE}_exon.txt \
    $OUT/${SAMPLE}.adjusted.bam

