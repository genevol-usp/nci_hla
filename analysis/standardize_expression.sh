#!/usr/bin/env bash

SAMPLES=./sample_ids.txt
BED=./salmon-pers/quants.bed
OUT=./salmon-pers/quants_std.bed

bgzip $BED && tabix -p bed $BED.gz

QTLtools correct --include-samples $SAMPLES --bed ${BED}.gz --normal --out $OUT

