#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=48:00:00
#SBATCH -p long
#SBATCH --array=1-96
#SBATCH --job-name=STARMAP
#SBATCH -o /home/vitor/nci_hla/analysis/map_to_genome/log/STARMAP-%a

# IO
SAMPLELIST=/home/vitor/nci_hla/analysis/samples.txt
PHASE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SAMPLELIST )

FQS1=( $( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R1*.fastq.gz ) )
FQS2=( $( ls -v /raid/genevol/nci_rnaseq/phase${PHASE}/fastq/${SAMPLE}*R2*.fastq.gz ) )

if [[ "${#FQS1[@]}" -eq 1 ]]; then
    FQ1=$FQS1
    FQ2=$FQS2
elif [[ "${#FQS1[@]}" -eq 2 ]]; then
    FQ1=${FQS1[0]},${FQS1[1]}
    FQ2=${FQS2[0]},${FQS2[1]}
fi

INDEX=/home/vitor/nci_hla/indices/STARLEAF
OUT=/scratch/vitor/bam_leaf/${SAMPLE}

# Mapping
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --outFileNamePrefix ${OUT}_

BAM=${OUT}_Aligned.sortedByCoord.out.bam
samtools index -@ $SLURM_CPUS_PER_TASK $BAM

## hla-mapper
#DB=/home/vitor/hla-mapper_new/db/hla-mapper_db_004.1_HLA
#HLAOUT=/scratch/vitor/bam_leaf/hlamapper/${SAMPLE}
#
#$HLAMAPPER rna \
#    bam=$BAM \
#    sample=$SAMPLE \
#    threads=$SLURM_CPUS_PER_TASK \
#    db=$DB \
#    output=$HLAOUT \
#    bwa=$BWA \
#    star=$STAR
#
#BAMHLA=${HLAOUT}/${SAMPLE}.adjusted.bam
#BAMREFCHR=${OUT}_refchrs.bam
#JUNC=/scratch/vitor/bam_leaf/hlamapper/junctions/${SAMPLE}.junc 
#
# Extract junctions
#samtools index -@ $SLURM_CPUS_PER_TASK $BAMHLA
#samtools view -@ $SLURM_CPUS_PER_TASK -b -h $BAM `(seq 1 22; echo X) | sed 's/^/chr/'` > $BAMREFCHR
#samtools index -@ $SLURM_CPUS_PER_TASK $BAMREFCHR
#regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $BAMREFCHR -o $JUNC
#
#rm $BAMREFCHR
#
