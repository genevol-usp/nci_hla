#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N sam-xlate
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

#sort BAM by reference, then by read name
SORTSAM=$HOME/Libraries/jvarkit/dist/sortsamrefname.jar
BAMPREFIX=/scratch/vitor/mapping/66K00003_Aligned
BAM=${BAMPREFIX}.sortedByCoord.out.bam
BAMSORTED=${BAMPREFIX}.sortedByRefName.out.bam

#java -Xmx48g -jar $SORTSAM -o $BAMSORTED --outputsamformat BAM $BAM 2> /dev/null

#translate to transcript coordinates
UBUJAR=$HOME/Libraries/ubu-1.2-jar-with-dependencies.jar
BED=$PBS_O_WORKDIR/gencode_v36.bed
OUT=${BAMPREFIX}.transcriptCoords.out.bam

java -Xmx48g -jar $UBUJAR sam-xlate \
    --bed $BED \
    --in $BAMSORTED \
    --xgtags --reverse \
    --out $OUT
