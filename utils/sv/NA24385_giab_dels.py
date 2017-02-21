#!/usr/bin/env python
"""Prepare SV truth set from NIST integrated calls for NA24385
"""
import os
import subprocess

from bcbio.variation import vcfutils

url = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_DraftIntegratedDeletionsgt19bp_v0.1.8/"
base = "GIAB_HG002_del_2tech_filt_%s_noTRs_v0.1.8.bed"
sizes = ["50to100bp", "100to1000bp", "1000to3000bp", "gt3000bp"]
fai_file = "/human/GRCh37/seq/GRCh37.fa.fai"
out_file = "NA24385-GIAB-2tech-dels-v0_1_8-GRCh37.bed"

for size in sizes:
    cur_url = url + base % size
    subprocess.check_call(["wget", "-c", cur_url])

with open("header.txt", "w") as out_handle:
    with open(os.path.basename(cur_url)) as in_handle:
        out_handle.write("# %s" % in_handle.readline())
subprocess.check_call(("cat GIAB_HG002* | grep -v ^X. | gsort /dev/stdin {fai_file} "
                       "| cat header.txt - > {out_file}".format(**locals())), shell=True)
out_file_hg19 = out_file.replace("GRCh37", "hg19")
with open(out_file) as in_handle:
    with open(out_file_hg19, "w") as out_handle:
        for line in in_handle:
            if not line.startswith("#"):
                line = "chr" + line
            out_handle.write(line)
vcfutils.bgzip_and_index(out_file, {})
vcfutils.bgzip_and_index(out_file_hg19, {})
