#!/usr/bin/env bash

ANNOT=$HOME/gencode/gencode.v37.primary_assembly.annotation.gtf
TRANSCRIPTS=./gencode.transcripts.fa
IMGT=$HOME/IMGTHLA

hlapers prepare-ref -t $TRANSCRIPTS -a $ANNOT -i $IMGT -o hladb 
