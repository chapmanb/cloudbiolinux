#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

date

set -eu -o pipefail
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/bin:$PATH
url_prefix=http://ftp.ensemblorg.ebi.ac.uk/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad/r2.1/exomes/
vcf_prefix=gnomad.exomes.r2.1.sites.chr
ref=../seq/GRCh37.fa
mkdir -p variation
export TMPDIR=`pwd`

vcf=${vcf_prefix}${chrom}_noVEP.vcf.gz
vcf_url=${url_prefix}${vcf}
#wget -c $vcf_url
#wget -c $vcf_url.tbi

fields_to_keep="INFO/"$(cat gnomad_fields_to_keep.txt | paste -s | sed s/"\t"/",INFO\/"/g)
gunzip -c $vcf | gsort -m 3000 /dev/stdin $ref.fai | bcftools view -f PASS -Ov | bcftools annotate -x "^$fields_to_keep" -Ov | grep -v "##contig="| bgzip -c >  variation/gnomad_exome.chr${chrom}.vcf.gz

tabix variation/gnomad_exome.chr${chrom}.vcf.gz

date