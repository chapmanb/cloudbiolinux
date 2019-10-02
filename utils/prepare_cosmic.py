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
import gzip
import logging
import os
import re
import requests
import subprocess
import sys
import tempfile
import shutil
from argparse import ArgumentParser

from bcbio import utils
from bcbio.variation import vcfutils

logging.basicConfig(format='%(asctime)s [%(levelname).1s] %(message)s', level=logging.INFO)


def main(cosmic_version, bcbio_genome_dir, overwrite=False, clean=False):
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cosmic-prep"))
    os.chdir(work_dir)

    for genome_build, bcbio_build, add_chr in [("GRCh37", "GRCh37", False), ("GRCh38", "hg38", True)]:
        bcbio_base = os.path.join(bcbio_genome_dir, "genomes", "Hsapiens", bcbio_build)
        installed_file = os.path.join(bcbio_base, "variation", f"cosmic-v{cosmic_version}.vcf.gz")
        installed_link = os.path.join(bcbio_base, "variation", "cosmic.vcf.gz")
        logging.info(f"Beginning COSMIC v{cosmic_version} prep for {genome_build}.")
        if not os.path.exists(bcbio_base):
            continue
        if os.path.exists(installed_file):
            if not overwrite:
                logging.info(f"{installed_file} exists, please use the --overwrite flag to overwrite the existing files if you want to reinstall.")
                continue
            else:
                logging.info(f"{installed_file} exists, removing.")
                remove_installed(installed_file, installed_link)
        bcbio_ref = os.path.join(bcbio_base, "seq", f"{bcbio_build}.fa")
        cosmic_vcf_files = get_cosmic_vcf_files(genome_build, cosmic_version, clean)
        sorted_inputs = []
        for fname in cosmic_vcf_files:
            sorted_inputs.append(sort_to_ref(fname, bcbio_ref, add_chr=add_chr))
        out_dir = utils.safe_makedir(os.path.join(f"v{cosmic_version}", "bcbio_ready", bcbio_build))
        out_file = os.path.join(out_dir, "cosmic.vcf.gz")
        ready_cosmic = combine_cosmic(sorted_inputs, bcbio_ref, out_file)
        variation_dir = utils.safe_makedir(os.path.join(bcbio_base, "variation"))
        utils.copy_plus(ready_cosmic, installed_file)
        logging.info(f"Created COSMIC v{cosmic_version} resource in {installed_file}.")
        logging.info(f"Linking {installed_file} as {installed_link}.")
        make_links(installed_file, installed_link)
        update_version_file(bcbio_base, cosmic_version)
        logging.info(f"Finished COSMIC v{cosmic_version} prep for {genome_build}.")
        # prepare hg19 from the GRCh37 file
        if bcbio_build == "GRCh37":
            genome_build = "hg19"
            logging.info(f"Prepping COSMIC v{cosmic_version} for {genome_build} from the GRCh37 preparation.")
            bcbio_base = os.path.join(bcbio_genome_dir, "genomes", "Hsapiens", genome_build)
            if not os.path.exists(bcbio_base):
                continue
            if os.path.exists(installed_file):
                installed_file = os.path.join(bcbio_base, "variation", f"cosmic-v{cosmic_version}.vcf.gz")
                installed_link = os.path.join(bcbio_base, "variation", "cosmic.vcf.gz")
                if not overwrite:
                    logging.info(f"{installed_file} exists, please use the --overwrite flag to overwrite the existing files if you want to reinstall.")
                    continue
                else:
                    logging.info(f"{installed_file} exists, removing.")
                    remove_installed(installed_file, installed_link)
            out_dir = utils.safe_makedir(os.path.join(f"v{cosmic_version}", "bcbio_ready", genome_build))
            out_file = os.path.join(out_dir, f"cosmic-v{cosmic_version}.vcf.gz")
            logging.info(f"Translating GRCh37 chromosome names to hg19 chromosome names.")
            hg19_cosmic = map_coords_to_ucsc(ready_cosmic, bcbio_ref, out_file)
            variation_dir = utils.safe_makedir(os.path.join(bcbio_base, "variation"))
            utils.copy_plus(hg19_cosmic, installed_file)
            logging.info(f"Created COSMIC v{cosmic_version} resource in {installed_file}.")
            logging.info(f"Linking {installed_file} as {installed_link}.")
            make_links(installed_file, installed_link)
            update_version_file(bcbio_base, cosmic_version)
            logging.info(f"Finished COSMIC v{cosmic_version} prep for {genome_build}.")


def remove_installed(installed_file, installed_link):
    logging.info(f"Removing {installed_file}.")
    if os.path.exists(installed_file):
        os.remove(installed_file)
    installed_index = installed_file + ".tbi"
    logging.info(f"Removing {installed_index}.")
    if os.path.exists(installed_index):
        os.remove(installed_index)
    logging.info(f"Removing {installed_link}.")
    if os.path.lexists(installed_link):
        os.remove(installed_link)
    installed_index = installed_link + ".tbi"
    logging.info(f"Removing {installed_index}.")
    if os.path.lexists(installed_index):
        os.remove(installed_index)


def make_links(installed_file, installed_link):
    if os.path.islink(installed_link):
        os.remove(installed_link)
        os.remove(installed_link + ".tbi")
    os.symlink(os.path.basename(installed_file), installed_link)
    os.symlink(os.path.basename(installed_file + ".tbi"), installed_link + ".tbi")


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
    logging.info(f"Combining COSMIC files to {out_file}.")
    if not os.path.exists(out_file):
        cmd = ["picard", "MergeVcfs", "O=%s" % out_file, "D=%s" % ref_file.replace(".fa", ".dict")] + \
              ["I=%s" % x for x in fnames] + \
              ["USE_JDK_DEFLATER=true", "USE_JDK_INFLATER=true", "CREATE_INDEX=false"]
        subprocess.check_call(cmd)
    return vcfutils.bgzip_and_index(out_file, {})


def sort_to_ref(fname, ref_file, add_chr):
    """Match reference genome ordering.
    """
    logging.info(f"Sorting {fname} to match the order of {ref_file}.")
    out_file = "%s-prep.vcf.gz" % (fname.replace(".vcf.gz", ""))
    if not os.path.exists(out_file):
        if add_chr:
            fix_chrom = r'| sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g" '
        else:
            fix_chrom = ''
        contig_cl = vcfutils.add_contig_to_header_cl(ref_file, out_file)
        cmd = ("gunzip -c {fname} {fix_chrom} | "
               "bcftools norm --check-ref s --do-not-normalize -f {ref_file} |"
               "bcftools view -e 'SNP=1' |"
               "gsort /dev/stdin {ref_file}.fai | {contig_cl} | "
               "bgzip -c > {out_file}")
        subprocess.check_call(cmd.format(**locals()), shell=True)
    logging.info(f"bgzipping and indexing {out_file}.")
    return vcfutils.bgzip_and_index(out_file, {})


def get_cosmic_vcf_files(genome_build, cosmic_version, clean):
    """Retrieve using new authentication based download approach.

    GRCh38/cosmic/v85/VCF/CosmicCodingMuts.vcf.gz
    GRCh38/cosmic/v85/VCF/CosmicNonCodingVariants.vcf.gz
    """
    vdir = os.path.join("v%s" % cosmic_version, genome_build)
    if os.path.exists(vdir):
        if not clean:
            logging.info(f"{vdir} files exist, please use the --clean flag to overwrite the existing files if you want to reinstall.")
        else:
            logging.info(f"{vdir} exists, removing.")
            remove_cosmic_directory(vdir)
    logging.info("Downloading COSMIC VCF files.")
    url = "https://cancer.sanger.ac.uk/cosmic/file_download/"
    out_dir = utils.safe_makedir(os.path.join("v%s" % cosmic_version, genome_build))
    fnames = []
    for ctype in ["CosmicCodingMuts", "CosmicNonCodingVariants"]:
        filename = os.path.join(out_dir, "%s.vcf.gz" % ctype)
        if not os.path.exists(filename):
            filepath = "%s/cosmic/v%s/VCF/%s.vcf.gz" % (genome_build, cosmic_version, ctype)
            logging.info("Downloading %s" % (url + filepath))
            try:
                r = requests.get(url + filepath, auth=(os.environ["COSMIC_USER"], os.environ["COSMIC_PASS"]))
            except KeyError as e:
                print("KeyError: {} not found. Be sure to export your COSMIC_USER and COSMIC_PASS before running in order to download the files".format(e))
                raise e
            download_url = r.json()["url"]
            r = requests.get(download_url)
            with open(filename, "wb") as f:
                f.write(r.content)
        fnames.append(filename)
    return fnames


def remove_cosmic_directory(installed_directory):
    logging.info(f"Removing {installed_directory}.")
    shutil.rmtree(installed_directory)

def update_version_file(bcbio_base, version):
    """
    update the version of cosmic used in the versions.csv file, adding it if it does not exist
    """
    versionfile = os.path.join(bcbio_base, "versions.csv")
    updatedfile = os.path.join(bcbio_base, "versions.csv-tmp")
    logging.info(f"Updating {versionfile}.")
    found = False
    with open(versionfile) as in_handle, open(updatedfile, "w") as out_handle:
        for line in in_handle:
            tokens = line.split(",")
            if tokens[0] != "cosmic":
                out_handle.write(line)
            else:
                # only write it once
                if not found:
                    out_handle.write(f"cosmic,{version}\n")
                found = True
        if not found:
            out_handle.write(f"cosmic,{version}\n")
    shutil.move(updatedfile, versionfile)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("cosmic_version", help="COSMIC version to install.", default="89")
    parser.add_argument("bcbio_directory", help="Path to bcbio installation. Should contain the 'genomes' directory.")
    parser.add_argument("--overwrite", action="store_true", default=False, help="Overwrite existing cosmic.vcf.gz files.")
    parser.add_argument("--clean", action="store_true", default=False, help="Remove pre-downloaded files, if available.")
    args = parser.parse_args()
    main(args.cosmic_version, args.bcbio_directory, args.overwrite, args.clean)
