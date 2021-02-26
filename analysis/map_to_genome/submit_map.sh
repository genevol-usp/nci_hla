#!/usr/bin/env bash

#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-107
#PBS -N STAR-map
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

PHASE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ../samples.txt )
SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $2 }' ../samples.txt )

FQS1=( $( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R1*.fastq.gz ) )
FQS2=( $( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R2*.fastq.gz ) )

if [[ "${#FQS1[@]}" -eq 1 ]]; then
    FQ1=$FQS1
    FQ2=$FQS2
elif [[ "${#FQS1[@]}" -eq 2 ]]; then
    FQ1=${FQS1[0]},${FQS1[1]}
    FQ2=${FQS2[0]},${FQS2[1]}
fi

STAR=$HOME/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
INDEX=../../simulation/map_to_genome/index
OUT=/media/storage/genevol/vitor/bam_nci/${SAMPLE}_timepoint${PHASE}

$STAR --runMode alignReads \
    --runThreadN $PBS_NUM_PPN \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMultimapNmax 20 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --quantMode TranscriptomeSAM \
    --quantTranscriptomeBan Singleend \
    --outFileNamePrefix ${OUT}_

samtools index ${OUT}_Aligned.sortedByCoord.out.bam
