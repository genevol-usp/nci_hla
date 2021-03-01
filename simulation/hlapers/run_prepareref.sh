#!/usr/bin/env bash

ANNOT=/home/vitor/hisat2/grch38_snp_tran/Homo_sapiens.GRCh38.99.gtf
TRANSCRIPTS=../salmon-reads/ensembl.transcripts.fa
IMGT=$HOME/IMGTHLA

hlapers prepare-ref -t $TRANSCRIPTS -a $ANNOT -i $IMGT -o hladb 
