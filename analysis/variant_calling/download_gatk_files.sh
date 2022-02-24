!/usr/bin/env bash


DBSNP_IN=/scratch/vitor/GCF_000001405.39.gz
DBSNP_OUT=./dbsnp_155.hg38.vcf.gz 

# Newer dbSNP from NCBI
wget -P /scratch/vitor https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
wget -P /scratch/vitor https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt 

tabix -p vcf $DBSNP_IN

# convert chr names to Gencode standard
Rscript prepare_gatk_files.R /scratch/vitor/GCF_000001405.39_GRCh38.p13_assembly_report.txt  

bcftools annotate \
    --rename-chrs /scratch/vitor/chr_dbsnpToGencode_names.txt \
    -O z -o $DBSNP_OUT \
    $DBSNP_IN

tabix -p vcf $DBSNP_OUT

# select chrs in GRCh38 Primary Assembly
DBSNP_PRI=./dbsnp_155.hg38pri.vcf.gz 
bcftools view -R ref_pri_chr.bed -O z -o $DBSNP_PRI $DBSNP_OUT 
tabix -p vcf $DBSNP_PRI

#rm $DBSNP_IN $DBSNP_OUT 
