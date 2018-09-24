#!/bin/bash
set -eu -o pipefail
version=20160215
civicv=2016-02-15
civicv38=2016-02-15

mkdir -p hg19/cancer
mkdir -p GRCh37/cancer
mkdir -p hg38/cancer

# prep work to do
# bcbio-prioritize create-civic -b GRCh37
# bcbio-prioritize create-civic -b GRCh38
# python az300_to_bed.py AZ300.txt
# python az300_to_bed.py AZ300_with_known.txt
# python az300_to_bed.py az-cancer-panel.txt

# CIViC
zcat civic-GRCh37-$civicv.bed.gz | sort -V | bgzip -c > GRCh37/cancer/civic-$civicv.bed.gz
tabix -f -p bed GRCh37/cancer/civic-$civicv.bed.gz
zcat civic-GRCh37-$civicv.bed.gz | sort -V | sed 's/^/chr/' | bgzip -c > hg19/cancer/civic-$civicv.bed.gz
tabix -f -p bed hg19/cancer/civic-$civicv.bed.gz
zcat civic-GRCh38-$civicv38.bed.gz | sort -V | sed 's/^/chr/' | bgzip -c > hg38/cancer/civic-$civicv38.bed.gz
tabix -f -p bed hg38/cancer/civic-$civicv38.bed.gz
tabix -f -p bed -m 12 --csi hg38/cancer/civic-$civicv38.bed.gz

# az300
cat AZ300-hg19.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg19/cancer/az300.bed.gz
tabix -f -p bed hg19/cancer/az300.bed.gz
cat AZ300-GRCh37.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > GRCh37/cancer/az300.bed.gz
tabix -f -p bed GRCh37/cancer/az300.bed.gz
cat AZ300-hg38.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg38/cancer/az300.bed.gz
tabix -f -p bed hg38/cancer/az300.bed.gz
tabix -f -p bed -m 12 --csi hg38/cancer/az300.bed.gz

# az300 with fusions
cat AZ300_with_known-hg19.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg19/cancer/az300-with-fusion.bed.gz
tabix -f -p bed hg19/cancer/az300-with-fusion.bed.gz
cat AZ300_with_known-GRCh37.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > GRCh37/cancer/az300-with-fusion.bed.gz
tabix -f -p bed GRCh37/cancer/az300-with-fusion.bed.gz
cat AZ300_with_known-hg38.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg38/cancer/az300-with-fusion.bed.gz
tabix -f -p bed hg38/cancer/az300-with-fusion.bed.gz
tabix -f -p bed -m 12 --csi hg38/cancer/az300-with-fusion.bed.gz

# az-cancer-panel
cat az-cancer-panel-hg19.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg19/cancer/az-cancer-panel.bed.gz
tabix -f -p bed hg19/cancer/az-cancer-panel.bed.gz
cat az-cancer-panel-GRCh37.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > GRCh37/cancer/az-cancer-panel.bed.gz
tabix -f -p bed GRCh37/cancer/az-cancer-panel.bed.gz
cat az-cancer-panel-hg38.bed | cut -f 1-4 | sort -V -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bgzip -c > hg38/cancer/az-cancer-panel.bed.gz
tabix -f -p bed hg38/cancer/az-cancer-panel.bed.gz
tabix -f -p bed -m 12 --csi hg38/cancer/az-cancer-panel.bed.gz

# tar up downloads
cd hg19
tar -czvpf ../prioritize-cancer-hg19-$version.tar.gz cancer
cd ../GRCh37
tar -czvpf ../prioritize-cancer-GRCh37-$version.tar.gz cancer
cd ../hg38
tar -czvpf ../prioritize-cancer-hg38-$version.tar.gz cancer
cd ..

aws s3 cp prioritize-cancer-hg19-$version.tar.gz s3://biodata/coverage/prioritize/ --acl public-read
aws s3 cp prioritize-cancer-GRCh37-$version.tar.gz s3://biodata/coverage/prioritize/ --acl public-read
aws s3 cp prioritize-cancer-hg38-$version.tar.gz s3://biodata/coverage/prioritize/ --acl public-read
