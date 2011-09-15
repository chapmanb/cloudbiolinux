#!/usr/bin/env python
"""Prepare GFF transcript files for use as input to Cufflinks.

This supports running RNA-seq pipelines by providing reference annotation and
mask files needed by Cufflinks, ready for use with UCSC named genomes:

http://cufflinks.cbcb.umd.edu/manual.html

Usage:
  prepare_tx_gff.py <org_build>
"""
import os
import sys
import shutil
import platform
import collections
import subprocess

import MySQLdb
import rpy2.robjects as robjects

from bcbio.utils import chdir, safe_makedir

# ##  Version and retrieval details for Ensembl and UCSC
ensembl_release = "63"
base_ftp = "ftp://ftp.ensembl.org/pub/release-{release}/gtf"

Build = collections.namedtuple("Build", ["taxname", "fname", "biomart_name"])
build_info = {
    "hg19": Build("homo_sapiens", "Homo_sapiens.GRCh37.{release}.gtf.gz",
                  "hsapiens_gene_ensembl")}

ucsc_db= "genome-mysql.cse.ucsc.edu"
ucsc_user="genome"

# ## Main driver functions

def main(org_build):
    work_dir = os.path.join(os.getcwd(), "tmpcbl")
    safe_makedir(work_dir)
    with chdir(work_dir):
        build = build_info[org_build]
        tx_gff = prepare_tx_gff(build, org_build)
        mask_gff = prepare_mask_gff(tx_gff, build)
        upload_to_s3([tx_gff, mask_gff], org_build)

def upload_to_s3(fnames, org_build):
    final_dir = "{org}-rnaseq".format(org=org_build)
    final_tarball = "{0}.tar.xz".format(final_dir)
    if not os.path.exists(final_tarball):
        safe_makedir(final_dir)
        for fname in fnames:
            shutil.copy(fname, final_dir)
        subprocess.check_call("tar -cvpf - {dir} | xz -zc - > {tarball}".format(
            dir=final_dir, tarball=final_tarball), shell=True)
    python_exe = "python{0}.{1}".format(*platform.python_version_tuple()[:2])
    upload_script = os.path.join(os.path.dirname(__file__), "s3_multipart_upload.py")
    subprocess.check_call([python_exe, upload_script, final_tarball, "biodata",
                           os.path.join("annotation", os.path.basename(final_tarball)),
                           "--public"])

# ## Get transcripts to filter from Ensembl BioMart query

def prepare_mask_gff(base_gff, build):
    """Prepare GFF file with high abundance transcripts to mask.
    """
    out_file = "{0}-mask{1}".format(*os.path.splitext(base_gff))
    if not os.path.exists(out_file):
        tx_ids = high_abudance_transcripts(build)
        with open(base_gff) as in_handle, \
             open(out_file, "w") as out_handle:
            for line in in_handle:
                tx = line.split('transcript_id "')[1].split('"')[0]
                if tx in tx_ids:
                    out_handle.write(line)
    return out_file

def high_abudance_transcripts(build):
    """Get identifiers for high abundance transcriptions using R biomaRt.
    """
    robjects.r.assign("biomart.org", build.biomart_name)
    robjects.r('''
      library(biomaRt)
      mart <- useMart("ensembl", dataset=biomart.org)
      attrs <- c("ensembl_transcript_id")
      filters <- c("biotype")
      filter.types <- c("Mt_rRNA", "rRNA", "rRNA_pseudogene")
      result <- getBM(attributes=attrs, filters=filters, values=filter.types,
                      mart=mart)
      result <- unique(result)
    ''')
    # get first column of data frame: the transcript IDs
    return set(robjects.r["result"][0])

# ## Retrieve GFF file from Ensembl and map to UCSC coordinates

def prepare_tx_gff(build, org_name):
    """Prepare UCSC ready transcript file given build information.
    """
    ensembl_gff = _download_ensembl_gff(build)
    ucsc_name_map = _query_for_ucsc_ensembl_map(org_name)
    return _remap_gff(ensembl_gff, ucsc_name_map)

def _remap_gff(base_gff, name_map):
    """Remap chromosome names to UCSC instead of Ensembl
    """
    base, ext = os.path.splitext(base_gff)
    out_file = "{0}-ucsc{1}".format(base.replace(".", "-"), ext)
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle, \
             open(base_gff) as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                ucsc_name = name_map.get(parts[0], None)
                if ucsc_name:
                    out_handle.write("\t".join([ucsc_name] + parts[1:]))
    return out_file

def _query_for_ucsc_ensembl_map(org_name):
    """Retrieve UCSC to Ensembl name mappings from UCSC MySQL database.
    """
    db = MySQLdb.connect(host=ucsc_db, user=ucsc_user, db=org_name)
    cursor = db.cursor()
    cursor.execute("select * from ucscToEnsembl")
    ucsc_map = {}
    for ucsc, ensembl in cursor.fetchall():
        ucsc_map[ensembl] = ucsc
    return ucsc_map

def _download_ensembl_gff(build):
    """Given build details, download and extract the relevant ensembl GFF.
    """
    dl_url = "/".join([base_ftp, build.taxname, build.fname]).format(release=ensembl_release)
    out_file = os.path.splitext(os.path.basename(dl_url))[0]
    if not os.path.exists(out_file):
        subprocess.check_call(["wget", dl_url])
        subprocess.check_call(["gunzip", os.path.basename(dl_url)])
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
