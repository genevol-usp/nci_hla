#!/usr/bin/env bash

SAMPLES=$1
BED=$2
OUT=$( echo $BED | sed 's/\.bed/_std.bed/' )

bgzip $BED && tabix -p bed $BED.gz

QTLtools correct --include-samples $SAMPLES --bed ${BED}.gz --normal --out $OUT

