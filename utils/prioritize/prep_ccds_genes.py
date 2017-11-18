#!/usr/bin/env python
"""Prepare adjusted CCDS BED file.

Consensus CDS (CCDS https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi) regions with two adjustments:

- 2 bps added to internal introns to capture canonical splice acceptor/donor sites
- Multiple transcripts from a single gene are merged into a single all inclusive gene entry.
"""
import subprocess
import sys

def main(in_file, out_name, fai_file):
    out_file  = "%s.bed" % out_name
    with open(in_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                gene, chrom, coords, _ = line.split()
                chrom = get_chrom(chrom)
                for start, end in split_coords(coords):
                    out_handle.write("\t".join([chrom, start, end, gene]) + "\n")
    cmd = "gsort {out_file} {fai_file} | bgzip -c > {out_file}.gz"
    subprocess.check_call(cmd.format(**locals()), shell=True)

def get_chrom(chrom):
    assert chrom in set([str(x) for x in range(1, 23)] + ["X", "Y"]), chrom
    return "chr%s" % chrom

def split_coords(coords):
    assert coords[0] == "("
    assert coords[-1] == ")", coords
    out = []
    for c in coords[1:-1].split(","):
        start, end = c.split("..")
        out.append((start, end))
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
