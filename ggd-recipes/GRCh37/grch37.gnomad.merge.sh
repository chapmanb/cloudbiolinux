#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# to merge preprocessed per chr vcfs, resulting file is 
date
cd txtmp/variation
bcftools concat gnomad_genome.chr1.vcf.gz gnomad_genome.chr2.vcf.gz gnomad_genome.chr3.vcf.gz gnomad_genome.chr4.vcf.gz gnomad_genome.chr5.vcf.gz gnomad_genome.chr6.vcf.gz \
gnomad_genome.chr7.vcf.gz gnomad_genome.chr8.vcf.gz gnomad_genome.chr9.vcf.gz gnomad_genome.chr10.vcf.gz gnomad_genome.chr11.vcf.gz gnomad_genome.chr12.vcf.gz \
gnomad_genome.chr13.vcf.gz gnomad_genome.chr14.vcf.gz gnomad_genome.chr15.vcf.gz gnomad_genome.chr16.vcf.gz gnomad_genome.chr17.vcf.gz gnomad_genome.chr18.vcf.gz \
gnomad_genome.chr19.vcf.gz gnomad_genome.chr20.vcf.gz gnomad_genome.chr21.vcf.gz gnomad_genome.chr22.vcf.gz gnomad_genome.chrX.vcf.gz -Ov | bgzip -c > gnomad_genome.vcf.gz
tabix -f -p vcf gnomad_genome.vcf.gz
tabix -f -p vcf --csi gnomad_genome.vcf.gz
cd ../../
date