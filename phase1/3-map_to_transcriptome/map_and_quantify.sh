#!/bin/bash

STAR=/home/vitor/Libraries/STAR
salmon=/home/vitor/Libraries/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=6
indexDIR=./sample_indices/$sample
fasta=/home/vitor/hlaexpression/imgt_index/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=../2-hla_typing/sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $indexDIR\
    --genomeFastaFiles $sample_fa\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${indexDIR}_

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
    --outFilterMultimapNmax 150\
    --winAnchorMultimapNmax 300\
    --alignIntronMax 0\
    --alignEndsType Local\
    --outSAMunmapped Within KeepPairs\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
outQuant=./quantifications
out=$outQuant/$sample

if [ -d "$out" ]; then
    rmdir $out
fi

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p $CPUS --seqBias --gcBias

rm -r $indexDIR ${indexDIR}_Log.out $sample_fa $outPrefix* 
