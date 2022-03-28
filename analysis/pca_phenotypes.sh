#!/usr/bin/env bash

SAMPLES=$1
BED=$2

QTLtools pca --include-samples $SAMPLES --bed BED --center --scale --out pca_salmon

