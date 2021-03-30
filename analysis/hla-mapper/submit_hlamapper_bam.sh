#!/usr/bin/env bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=120gb
#PBS -l walltime=48:00:00
#PBS -q bigmem
#PBS -t 1-107
#PBS -N hla-mapper
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )
BAM=/media/storage/genevol/vitor/nci/bam/${SAMPLE}_timepoint${PHASE}_Aligned.sortedByCoord.out.bam
OUT=/raid/genevol/users/vitor/nci_hla/analysis/hla-mapper/bam/${SAMPLE}_t${PHASE} 

mkdir -p $OUT

/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper rna \
    bam=$BAM \
    sample=$SAMPLE \
    threads=$PBS_NUM_PPN \
    db=/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper_db_004.1_HLA \
    output=$OUT \
    bwa=/raid/genevol/users/vitor/HLAMAPPER/bwa/bwa-0.7.17/bwa \
    samtools=/raid/genevol/users/vitor/HLAMAPPER/samtools-1.11/samtools \
    star=/raid/genevol/users/vitor/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

featureCounts -p -B -C -M -O --fraction -T $PBS_NUM_PPN -g transcript_id \
    -a /home/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf \
    -o ./quant/${SAMPLE}_t${PHASE}_counts.txt \
    $OUT/${SAMPLE}.adjusted.bam
