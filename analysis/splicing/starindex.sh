#!/usr/bin/env bash

GENOME=./chr6_masked.fasta
GTF=./hla.gtf
OUT=./STARINDEX

mkdir -p $OUT

STAR --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --sjdbOverhang 125 \
    --genomeSAindexNbases 12
