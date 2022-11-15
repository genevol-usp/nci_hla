#!/usr/bin/bash

#BAM=/scratch/vitor/hlamapper/66K00003/66K00003.adjusted.bam

#samtools view $BAM chr6:29940000-29947000 > 66K00003.sam

#BAM=/scratch/vitor/bam_leaf/66K00003_Aligned.sortedByCoord.out.bam

#samtools view $BAM chr6:29940000-29947000 > 66K00003_Aligned.sam


BAM=/media/storage/genevol/vitor/nci/bam_hlamapper/66K00003_t1/66K00003.adjusted.bam

samtools view $BAM chr6:29940000-29947000 > 66K00003_originalmapper.sam

