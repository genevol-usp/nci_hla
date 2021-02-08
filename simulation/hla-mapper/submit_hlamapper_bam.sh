#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-50
#PBS -N hla-mapper-bam
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk "FNR==$PBS_ARRAYID" ../samples.txt )
OUT=/media/storage/genevol/vitor/bam_hlamapper/${SAMPLE}

mkdir -p $OUT

/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper rna \
    bam=/media/storage/genevol/vitor/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    sample=$SAMPLE \
    threads=$PBS_NUM_PPN \
    db=/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper_db_004.1_HLA \
    output=$OUT \
    bwa=/raid/genevol/users/vitor/HLAMAPPER/bwa/bwa-0.7.17/bwa \
    samtools=/raid/genevol/users/vitor/HLAMAPPER/samtools-1.11/samtools \
    star=/raid/genevol/users/vitor/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

featureCounts -p -B -C -M -O --fraction --fracOverlap 0.2 -T $PBS_NUM_PPN -t exon -g transcript_id \
    -a /home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf \
    -o ./quant/${SAMPLE}_counts.txt \
    ${OUT}/${SAMPLE}.adjusted.bam
