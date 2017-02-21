#!/usr/bin/env python
"""Prepare truth set of crowd sourced CNVs from GiaB samples.

http://biorxiv.org/content/early/2016/12/13/093526
"""
import requests

from bcbio.variation import vcfutils

url = "http://biorxiv.org/content/biorxiv/suppl/2016/12/13/093526.DC1/093526-3.txt"

out_base = "NA24385-crowd-dels-%s.bed"

r = requests.get(url, stream=True)
header = None
calls = []
for line in r.iter_lines():
    if not header:
        header = line.rstrip().split()
    else:
        cur = dict(zip(header, line.rstrip().split()))
        pos = (int(cur["chrom"].replace("chr", "")), int(cur["start"]))
        probs = sorted([(float(cur["%s_prob" % c]), c) for c in ["CN0", "CN1", "CN2"]], reverse=True)
        cn = probs[0][-1]
        if probs[0][0] > 0.8 and cn != "CN2":
            line = "\t".join([cur["chrom"], cur["start"], cur["end"], cur["site_id"], cn]) + "\n"
            calls.append((pos, line))
calls.sort()

hg19_bed = out_base % "hg19"
grch37_bed = out_base % "GRCh37"
with open(hg19_bed, "w") as hg19_handle:
    with open(grch37_bed, "w") as grch37_handle:
        for _, line in calls:
            hg19_handle.write(line)
            grch37_handle.write(line.replace("chr", ""))

vcfutils.bgzip_and_index(hg19_bed, {})
vcfutils.bgzip_and_index(grch37_bed, {})
