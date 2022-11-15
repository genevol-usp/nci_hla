#!/usr/bin/env bash

ANNOT=./hla.gtf
EXON=./exon.txt.gz

bgzip -c $ANNOT > ${ANNOT}.gz

$HOME/Libraries/leafcutter/scripts/gtf_to_exons.R ${ANNOT}.gz $EXON
