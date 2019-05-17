#!/bin/bash

#PBS -l walltime=230:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

date
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/bin:$PATH
remap_url=https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_ensembl2UCSC.txt
ref=seq/hg19.fa
export TMPDIR=`pwd`
#wget --no-check-certificate -qO- $remap_url | awk '{if($1!=$2) print "s/^"$1"/"$2"/g"}' > remap.sed
gunzip -c gnomad_genome.grch37.vcf.gz | sed -f remap.sed | gsort -m 3000 /dev/stdin $ref.fai | vt normalize -r $ref -n - | vt uniq - | bgzip -c >  gnomad_genome.vcf.gz
tabix -f -p vcf gnomad_genome.vcf.gz
tabix -f -p vcf --csi gnomad_genome.vcf.gz
date
