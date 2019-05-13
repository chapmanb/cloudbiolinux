#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# merges per chr vcf files

export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/bin:$PATH
export TMPDIR=`pwd`

cd variation
bcftools concat gnomad_exome.chr1.vcf.gz gnomad_exome.chr2.vcf.gz gnomad_exome.chr3.vcf.gz gnomad_exome.chr4.vcf.gz gnomad_exome.chr5.vcf.gz gnomad_exome.chr6.vcf.gz \
gnomad_exome.chr7.vcf.gz gnomad_exome.chr8.vcf.gz gnomad_exome.chr9.vcf.gz gnomad_exome.chr10.vcf.gz gnomad_exome.chr11.vcf.gz gnomad_exome.chr12.vcf.gz \
gnomad_exome.chr13.vcf.gz gnomad_exome.chr14.vcf.gz gnomad_exome.chr15.vcf.gz gnomad_exome.chr16.vcf.gz gnomad_exome.chr17.vcf.gz gnomad_exome.chr18.vcf.gz \
gnomad_exome.chr19.vcf.gz gnomad_exome.chr20.vcf.gz gnomad_exome.chr21.vcf.gz gnomad_exome.chr22.vcf.gz gnomad_exome.chrX.vcf.gz gnomad_exome.chrY.vcf.gz -Ov | bgzip -c > gnomad_exome.vcf.gz

tabix -f -p vcf variation/gnomad_exome.vcf.gz
tabix -f -p vcf --csi variation/gnomad_exome.vcf.gz
cd ..
