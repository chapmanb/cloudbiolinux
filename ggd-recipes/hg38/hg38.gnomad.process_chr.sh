#!/bin/bash

#PBS -l walltime=23:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# to process chromosomes in parallel
# processing time ?? for 1 chromosome

set -eu -o pipefail
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/bin:$PATH
vcf_prefix=gnomad.genomes.r2.1.sites.grch38.chr

remap_url=http://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2UCSC.txt
ref=../seq/hg38.fa
mkdir -p variation
# download before running - no internet on compute nodes
# wget --no-check-certificate -qO- $remap_url | awk '{ print length, $0 }' | sort -n -s -r | cut -d" " -f2- | awk '{if(!$2) print "/^"$1"/d"; else if($1!=$2) print "s/^"$1"/"$2"/g";}' > remap.sed

export TMPDIR=`pwd`

gnomad_fields_to_keep_url=https://gist.githubusercontent.com/naumenko-sa/d20db928b915a87bba4012ba1b89d924/raw/cf343b105cb3347e966cc95d049e364528c86880/gnomad_fields_to_keep.txt
# download before running
# wget --no-check-certificate -c $gnomad_fields_to_keep_url

# no chrY in gnomad genome in hg38
vcf=${vcf_prefix}${chrom}_noVEP.vcf.gz
#  wget -c $vcf_url
#  wget -c $vcf_url.tbi

fields_to_keep="INFO/"$(cat gnomad_fields_to_keep.txt | paste -s | sed s/"\t"/",INFO\/"/g)
# bcftools annotate is picky about vcf header and brakes if moved down the pipe after remap
gunzip -c $vcf | bcftools view -f PASS -Ov | bcftools annotate -x "^$fields_to_keep" -Ov | sed -f remap.sed | grep -v "##contig=" | gsort -m 3000 /dev/stdin $ref.fai | vt normalize -r $ref -n - | vt uniq - | bgzip -c >  variation/gnomad_genome.chr${chrom}.vcf.gz

tabix variation/gnomad_genome.chr${chrom}.vcf.gz
