#!/usr/bin/env bash

mkdir -p gencode

$HOME/Libraries/leafcutter/leafviz/gtf2leafcutter.pl \
    -o gencode/gencode \
    ./hla.gtf
