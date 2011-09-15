"""Prepare supplemental files to work with RNA-seq transcriptome experiments.

Retrieves annotations in a format usable by Cufflinks, using repositories
provided by Illumina:

http://cufflinks.cbcb.umd.edu/igenomes.html
"""
import os

from fabric.api import *
from fabric.contrib.files import *

_orgname_map = {"Mmusculus": "Mus_musculus",
                "Hsapiens": "Homo_sapiens"}

def download_transcripts(genomes, env):
    folder_name = "rnaseq"
    base_url = "ftp://igenome:*password*@ftp.illumina.com/{orgname}/{data_source}" \
               "/{version}/{orgname}_{data_source}_{version}.tar.gz"
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("rnaseq", False)):
        cur_url = base_url.format(orgname=_orgname_map[orgname],
                                  version=gid, data_source=manager.data_source)
        tx_dir = os.path.join(genome_dir, orgname, gid, folder_name)
        if not exists(tx_dir):
            run("mkdir -p %s" % tx_dir)
        with cd(tx_dir):
            _download_igenome_bundle(cur_url)

def _download_igenome_bundle(url):
    """Download iGenome tarball and extract files of interest.
    """
    files_want = ["Genes/genes.gtf"]
    print url
