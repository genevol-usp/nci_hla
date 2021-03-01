#!/usr/bin/env bash

FASTADIR=$1

for i in {1..2}
do
    FA=${FASTADIR}/sample_01_${i}.fasta 
    FQ=${FASTADIR}/sample_01_${i}.fastq.gz

    cat $FA |\
	awk -v mate=$i 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {
	    print\
	    "@"gensub(/;/, " ", 1, gensub(/\//, "_", 1, $1))"/"mate"\n"\
	    $2"\n"\
	    "+\n"\
	    gensub(/./, "F", "g", $2)\
	}' | gzip -c > $FQ
done
