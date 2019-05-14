#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# process 1 chr : 3h per chr

date
set -eu -o pipefail
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/bin:$PATH
vcf_prefix=gnomad.genomes.r2.1.sites.chr
ref=../seq/GRCh37.fa
mkdir -p variation
export TMPDIR=`pwd`

# download before running as compute nodes have no internet
# gnomad_fields_to_keep_url=https://gist.githubusercontent.com/naumenko-sa/d20db928b915a87bba4012ba1b89d924/raw/cf343b105cb3347e966cc95d049e364528c86880/gnomad_fields_to_keep.txt
# no chr Y for gnomad genomes
vcf=${vcf_prefix}${chrom}_noVEP.vcf.gz

fields_to_keep="INFO/"$(cat gnomad_fields_to_keep.txt | paste -s | sed s/"\t"/",INFO\/"/g)
gunzip -c $vcf | gsort -m 3000 /dev/stdin $ref.fai | bcftools view -f PASS -Ov | bcftools annotate -x "^$fields_to_keep" -Ov | vt normalize -r $ref -n - | vt uniq - | grep -v "##contig="| bgzip -c >  variation/gnomad_genome.chr${chrom}.vcf.gz

tabix variation/gnomad_genome.chr${chrom}.vcf.gz
date