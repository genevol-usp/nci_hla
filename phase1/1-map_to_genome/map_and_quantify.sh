#!/bin/bash

STAR=/home/vitor/Libraries/STAR
samtools=/home/vitor/Libraries/samtools-1.3.1/samtools
seqtk=/home/vitor/Libraries/seqtk/seqtk
salmon=/home/vitor/Libraries/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=8
indexDIR=/home/vitor/hlaexpression/index_genome/star/index
outMap=./mappings
runid=$(pwgen 5 1)
outPrefix=/scratch/genevol/users/vitor/${runid}_${sample}_

fqs1=(`ls -v ../data/fastq/${sample}*R1_001.fastq.gz`)
fqs2=(`ls -v ../data/fastq/${sample}*R2_001.fastq.gz`)

if [ "${#fqs1[@]}" == 1 ]; then
    fq1=$fqs1
    fq2=$fqs2
elif [ "${#fqs1[@]}" == 2 ]; then
    fq1=${fqs1[0]},${fqs1[1]} 
    fq2=${fqs2[0]},${fqs2[1]} 
else
    echo "wrong number of fastq files"
fi

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 20\
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within KeepPairs\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bamGenome=${outPrefix}Aligned.sortedByCoord.out.bam
readsalign=${outPrefix}readsAligned.txt
readsunmap=${outPrefix}readsUnmapped.txt
readids=${outPrefix}readids.txt

$samtools index $bamGenome $bamGenome.bai

$samtools view $bamGenome "chr6:29722000-33144000" |\
    cut -f1 |\
    sort |\
    uniq > $readsalign

$samtools view -F 0x2 $bamGenome |\
    cut -f1 |\
    sort |\
    uniq > $readsunmap

cat $readsalign $readsunmap |\
    sort |\
    uniq > $readids

fq12=./mappings/mhc_fqs/${sample}_1.fq
fq22=./mappings/mhc_fqs/${sample}_2.fq

if [ "${#fqs1[@]}" == 1 ]; then
    $seqtk subseq $fq1 $readids > $fq12 
    $seqtk subseq $fq2 $readids > $fq22 
elif [ "${#fqs1[@]}" == 2 ]; then
    $seqtk subseq ${fqs1[0]} $readids > $fq12 
    $seqtk subseq ${fqs1[1]} $readids >> $fq12 
    $seqtk subseq ${fqs2[0]} $readids > $fq22 
    $seqtk subseq ${fqs2[1]} $readids >> $fq22 
else
    echo "wrong number of fastq files for seqtk"
fi

bamTransc=${outPrefix}Aligned.toTranscriptome.out.bam
fasta=/home/vitor/hlaexpression/index_transcriptome/gencode.v25.PRI.transcripts.fa
out=./quantifications/$sample

if [ -d "$outTransc" ]; then
    rm -r $outTransc
fi

$salmon quant -t $fasta -l IU -a $bamTransc -o $out -p $CPUS --seqBias --gcBias

rm ${outPrefix}*
