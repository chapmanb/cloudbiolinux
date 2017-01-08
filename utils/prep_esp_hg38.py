#!/usr/bin/env python
"""Prepare hg38 compatible ESP from original download.

http://evs.gs.washington.edu/EVS/

The Download is a tarball of individual VCF files with hg38 coordinates in
an INFO key. Extract these into a full VCF, then sort, normalize and bgzip.
"""
import glob
import subprocess

from bcbio.bam import ref
from bcbio.variation import vcfutils
from bcbio.heterogeneity import chromhacks

def main():
    url = "http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz"
    ref_file = "../seq/hg38.fa"
    subprocess.check_call("wget -c -O esp-orig.tar.gz {url}".format(**locals()), shell=True)
    subprocess.check_call("tar -xzvpf esp-orig.tar.gz", shell=True)
    raw_file = "esp-raw.vcf"
    with open(raw_file, "w") as out_handle:
        for i, chrom in enumerate(range(1, 22) + ["X", "Y"]):
            fnames = glob.glob("*chr%s.snps_indels.vcf" % chrom)
            assert len(fnames) == 1, (chrom, fnames)
            with open(fnames[0]) as in_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        if i == 0:
                            if line.startswith("#CHROM"):
                                _add_contigs(out_handle, ref_file)
                            out_handle.write(line)
                    else:
                        parts = line.strip().split("\t")
                        key, val = parts[-1].split(";")[-1].split("=")
                        assert key == "GRCh38_POSITION"
                        if val != "-1":
                            new_chrom, new_pos = val.split(":")
                            if chromhacks.is_autosomal_or_sex(new_chrom):
                                parts[0] = "chr%s" % new_chrom
                                parts[1] = new_pos
                                out_handle.write("\t".join(parts) + "\n")
    out_file = "ESP6500SI-V2-hg38.vcf.gz"
    subprocess.check_call(("vt sort {raw_file} | vt decompose -s - | "
                           "vt normalize -n -r {ref_file} - | bgzip -c > {out_file}").format(**locals()),
                          shell=True)
    vcfutils.bgzip_and_index(out_file)

def _add_contigs(out_handle, ref_file):
    for contig in ref.file_contigs(ref_file):
        if chromhacks.is_autosomal_or_sex(contig.name):
            out_handle.write("##contig=<ID=%s,length=%s>\n" % (contig.name, contig.size))

if __name__ == "__main__":
    main()
