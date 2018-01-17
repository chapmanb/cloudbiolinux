#!/usr/bin/env python
"""Prepare combined VCF files of COSMIC resource for cancer variant calling.

This prepares the specified version and copies the updates into your bcbio installation
in the correct locations in the variation directory.

Because of licensing restrictions, download from COSMIC requires registration:

https://cancer.sanger.ac.uk/cosmic/download

Usage:

  export COSMIC_USER="your@registered.email.edu"
  export COSMIC_PASS="cosmic_password"
  bcbio_python prepare_cosmic.py <cosmic_version> </path/to/bcbio>

References:
http://gatkforums.broadinstitute.org/discussion/2226/cosmic-and-dbsnp-files-for-mutect
"""
import os
import subprocess
import sys

import paramiko

from bcbio import utils
from bcbio.variation import vcfutils

def main(cosmic_version, bcbio_genome_dir):
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cosmic-prep"))
    os.chdir(work_dir)

    for genome_build, bcbio_build, add_chr in [("grch37", "GRCh37", False), ("grch38", "hg38", True)]:
        bcbio_base = os.path.join(bcbio_genome_dir, "genomes", "Hsapiens", bcbio_build)
        if not os.path.exists(bcbio_base):
            continue
        bcbio_ref = os.path.join(bcbio_base, "seq", "%s.fa" % bcbio_build)
        sorted_inputs = []
        for fname in get_cosmic_files(genome_build, cosmic_version):
            sorted_inputs.append(sort_to_ref(fname, bcbio_ref, add_chr=add_chr))
        out_dir = utils.safe_makedir(os.path.join("v%s" % cosmic_version, "bcbio_ready", bcbio_build))
        out_file = os.path.join(out_dir, "cosmic.vcf.gz")
        ready_cosmic = combine_cosmic(sorted_inputs, bcbio_ref, out_file)
        utils.copy_plus(ready_cosmic, os.path.join(bcbio_base, "variation", os.path.basename(ready_cosmic)))
        if bcbio_build == "GRCh37":
            bcbio_base = os.path.join(bcbio_genome_dir, "genomes", "Hsapiens", "hg19")
            if not os.path.exists(bcbio_base):
                continue
            out_dir = utils.safe_makedir(os.path.join("v%s" % cosmic_version, "bcbio_ready", "hg19"))
            out_file = os.path.join(out_dir, "cosmic.vcf.gz")
            hg19_cosmic = map_coords_to_ucsc(ready_cosmic, bcbio_ref, out_file)
            utils.copy_plus(hg19_cosmic, os.path.join(bcbio_base, "variation", os.path.basename(ready_cosmic)))

def map_coords_to_ucsc(grc_cosmic, ref_file, out_file):
    hg19_ref_file = ref_file.replace("GRCh37", "hg19")
    if not os.path.exists(out_file):
        contig_cl = vcfutils.add_contig_to_header_cl(hg19_ref_file, out_file)
        cmd = ("zcat {grc_cosmic} | "
               r'sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g" '
               "| {contig_cl} "
               "| bgzip -c > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    if os.path.exists("%s-header.txt" % utils.splitext_plus(out_file)[0]):
        os.remove("%s-header.txt" % utils.splitext_plus(out_file)[0])
    return vcfutils.bgzip_and_index(out_file, {})

def _rename_to_ucsc(line):
    chrom, rest = line.split("\t", 1)
    if chrom == "MT":
        new_chrom = "chrM"
    else:
        new_chrom = "chr%s" % chrom
    return "%s\t%s" % (new_chrom, rest)

def combine_cosmic(fnames, ref_file, out_file):
    if not os.path.exists(out_file):
        cmd = ["picard", "MergeVcfs", "O=%s" % out_file, "D=%s" % ref_file.replace(".fa", ".dict")] + \
              ["I=%s" % x for x in fnames] + \
              ["USE_JDK_DEFLATER=true", "USE_JDK_INFLATER=true", "CREATE_INDEX=false"]
        subprocess.check_call(cmd)
    return vcfutils.bgzip_and_index(out_file, {})

def sort_to_ref(fname, ref_file, add_chr):
    """Match reference genome ordering.
    """
    out_file = "%s-prep.vcf.gz" % (fname.replace(".vcf.gz", ""))
    if not os.path.exists(out_file):
        if add_chr:
            fix_chrom = r'| sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g" '
        else:
            fix_chrom = ''
        contig_cl = vcfutils.add_contig_to_header_cl(ref_file, out_file)
        cmd = ("gunzip -c {fname} {fix_chrom} | "
               "gsort /dev/stdin {ref_file}.fai | {contig_cl} | "
               "bgzip -c > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    return vcfutils.bgzip_and_index(out_file, {})

def get_cosmic_files(genome_build, cosmic_version):
    """
    get /files/grch38/cosmic/v83/VCF/CosmicCodingMuts.vcf.gz
    get /files/grch38/cosmic/v83/VCF/CosmicNonCodingVariants.vcf.gz
    """
    out_dir = utils.safe_makedir(os.path.join("v%s" % cosmic_version, genome_build))
    fnames = []
    for ctype in ["CosmicCodingMuts", "CosmicNonCodingVariants"]:
        fname = os.path.join(out_dir, "%s.vcf.gz" % ctype)
        if not os.path.exists(fname):
            remote_path = "/files/%s/cosmic/v%s/VCF/%s.vcf.gz" % (genome_build, cosmic_version, ctype)
            t = paramiko.Transport(("sftp-cancer.sanger.ac.uk", 22))
            t.connect(username=os.environ["COSMIC_USER"], password=os.environ["COSMIC_PASS"])
            sftp = paramiko.SFTPClient.from_transport(t)
            print("Getting %s" % remote_path)
            sftp.get(remote_path, fname)
        fnames.append(fname)
    return fnames

if __name__ == "__main__":
    main(*sys.argv[1:])
