#!/usr/bin/env bash

# IO
JUNC=./juncfiles.txt
OUT=./results

mkdir -p $OUT

## Intron clustering
python $HOME/Libraries/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -m 10 \
    -l 500000 \
    -p 0.05 \
    -r ${OUT} \
    -o leafcutter

rm ${OUT}/*.leafcutter.sorted.gz

# Differential splicing analysis for each dataset and cell type
COUNT=${OUT}/leafcutter_perind_numers.counts.gz
GROUP=./group_file.txt
EXON=./exon.txt.gz

${HOME}/Libraries/leafcutter/scripts/leafcutter_ds.R \
    -i 5 -g 5 \
    -p 4 \
    -e $EXON \
    -o ${OUT}/leafcutter \
    $COUNT \
    $GROUP

# Prepare results for visualization with leafviz
${HOME}/Libraries/leafcutter/leafviz/prepare_results.R \
    --code leafcutter-a1 \
    --meta_data_file $GROUP \
    $COUNT \
    ${OUT}/leafcutter_cluster_significance.txt \
    ${OUT}/leafcutter_effect_sizes.txt \
    ./gencode/gencode \
    -o ${OUT}/a1.Rdata
