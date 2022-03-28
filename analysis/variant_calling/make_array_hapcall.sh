#!/bin/bash

echo -n > array_hapcall.txt

for s in $(awk '{print $1}' bam_list.txt ) 
do 
    for j in $(seq 1 22; echo X)
    do 
	echo -e ${s}'\t'chr${j} >> array_hapcall.txt
    done
done
