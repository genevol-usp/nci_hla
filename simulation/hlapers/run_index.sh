#!/usr/bin/env bash

TRANSCRIPTS=./hladb/transcripts_MHC_HLAsupp.fa

$HOME/Libraries/HLApers/hlapers index -t $TRANSCRIPTS -p 8 -o index
