#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=200gb
#PBS -l walltime=24:00:00
#PBS -q bigmem
#PBS -t 83
#PBS -N hla-mapper
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )
BAMFILE=${SAMPLE}_timepoint${PHASE}_Aligned.sortedByCoord.out.bam
BAMPATH=/media/storage/genevol/vitor/bam_nci \
OUTSCRATCH=/scratch/vitor/hlamapper/${SAMPLE}_t${PHASE} 
OUT=/media/storage/genevol/vitor/bam_nci/hlamapper

mkdir -p $OUTSCRATCH

cp $BAMPATH/$BAMFILE /scratch/vitor/hlamapper/
cp $BAMPATH/${BAMFILE}.bai /scratch/vitor/hlamapper/

/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper rna \
    bam=/scratch/vitor/hlamapper/$BAMFILE \
    sample=$SAMPLE \
    threads=$PBS_NUM_PPN \
    db=/raid/genevol/users/vitor/HLAMAPPER/hlamapper/hla-mapper_db_004.1_HLA \
    output=$OUTSCRATCH \
    bwa=/raid/genevol/users/vitor/HLAMAPPER/bwa/bwa-0.7.17/bwa \
    samtools=/raid/genevol/users/vitor/HLAMAPPER/samtools-1.11/samtools \
    star=/raid/genevol/users/vitor/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

featureCounts -p -B -C -M -O --fraction -d 0 -D 1000 -T $PBS_NUM_PPN -t exon -g transcript_id \
    -a /home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf \
    -o ./quant/${SAMPLE}_t${PHASE}_counts.txt \
    $OUTSCRATCH/${SAMPLE}.adjusted.bam

rm /scratch/vitor/hlamapper/$BAMFILE*
mv $OUTSCRATCH ${OUT}/
