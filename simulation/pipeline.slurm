#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p short
#SBATCH --array=1-50
#SBATCH --job-name=simulpipe
#SBATCH -o /scratch/vitor/log/simulation.o-%A-%a

# Software
STAR=${HOME}/HLAMAPPER/STAR-2.7.3a/bin/Linux_x86_64_static/STAR
HLAMAPPER=/lab-share/IM-Gutierrez-e2/Public/vitor/nci/hlamapper/build/hla-mapper
SAMTOOLS=${HOME}/HLAMAPPER/samtools-1.11/samtools
BWA=${HOME}/HLAMAPPER/bwa/bwa-0.7.17/bwa
MUDSKIPPER=${HOME}/Libraries/mudskipper/target/release/mudskipper
SALMON=${HOME}/Libraries/salmon-1.8.0_linux_x86_64/bin/salmon

# References and annotations
INDEX=${HOME}/nci_hla/indices/STAR_v3.7.3a
ANNOT=/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf
TRANSCRIPTS=${HOME}/nci_hla/indices/gencode.transcripts.fa

# Sample input data
SAMPLELIST=${SLURM_SUBMIT_DIR}/samples.txt
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
FQ1=/media/storage/genevol/vitor/simulation/fastq/${SAMPLEID}/sample_01_1.fastq.gz
FQ2=/media/storage/genevol/vitor/simulation/fastq/${SAMPLEID}/sample_01_2.fastq.gz

# Output
OUT=/media/storage/genevol/vitor/simulation/bam/${SAMPLEID}
BAM=${OUT}_Aligned.sortedByCoord.out.bam
HLABAM=/media/storage/genevol/vitor/simulation/bam_hlamapper/${SAMPLEID}/${SAMPLEID}.adjusted.bam
HLAFIX=/scratch/vitor/${SAMPLEID}.fixed
HLASORTED=/scratch/vitor/${SAMPLEID}.sortedByName.bam 
BAMTRANSC=/media/storage/genevol/vitor/simulation/bam_transc/${SAMPLEID}.bam

# Main
$STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --outFileNamePrefix ${OUT}_

$SAMTOOLS index $BAM

$HLAMAPPER rna \
    bam=$BAM \
    sample=$SAMPLEID \
    threads=$SLURM_CPUS_PER_TASK \
    db=$HOME/hla-mapper_new/db/hla-mapper_db_004.1_HLA \
    output=/media/storage/genevol/vitor/simulation/bam_hlamapper/${SAMPLEID} \
    bwa=$BWA \
    samtools=$SAMTOOLS \
    star=$STAR

samtools view $HLABAM |\
    sort --parallel=$SLURM_CPUS_PER_TASK -u > ${HLAFIX}.sam 

samtools view -H $HLABAM |\
    cat - ${HLAFIX} |\
    samtools view -@ $SLURM_CPUS_PER_TASK -Sb - > ${HLAFIX}.bam

samtools sort -@ $SLURM_CPUS_PER_TASK -n -o $HLASORTED ${HLAFIX}.bam

$MUDSKIPPER bulk \
    --threads $SLURM_CPUS_PER_TASK \
    --index ${HOME}/nci_hla/indices/mudskipper_gencode_v37 \
    --alignment $HLASORTED \
    --out $BAMTRANSC

rm $HLASORTED ${HLAFIX}.sam ${HLAFIX}.bam

$SALMON quant \
    -p $SLURM_CPUS_PER_TASK \
    -t $TRANSCRIPTS \
    -l A \
    -a $BAMTRANSC \
    -o ${SLURM_SUBMIT_DIR}/pipeline_results/salmon/${SAMPLEID}
