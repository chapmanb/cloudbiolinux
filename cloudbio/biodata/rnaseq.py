"""Prepare supplemental files to work with RNA-seq transcriptome experiments.

Retrieves annotations in a format usable by Cufflinks, using repositories
provided by Illumina:

http://cufflinks.cbcb.umd.edu/igenomes.html
"""
import os

from fabric.api import *
from fabric.contrib.files import *

def download_transcripts(genomes, env):
    folder_name = "rnaseq"
    base_url = "https://s3.amazonaws.com/biodata/annotation/{gid}-rnaseq.tar.xz"
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("rnaseq", False)):
        cur_url = base_url.format(orgname=_orgname_map[orgname],
                                  version=gid, data_source=manager.data_source)
        org_dir = os.path.join(genome_dir, orgname)
        tx_dir = os.path.join(org_dir, gid, folder_name)
        if not exists(tx_dir):
            with cd(org_dir):
                _download_annotation_bundle(url.format(gid=gid))

def _download_annotation_bundle(url):
    """Download bundle of RNA-seq data from S3 biodata/annotation
    """
    tarball = os.path.basename(url)
    if not exists(tarball):
        run("wget %s" % url)
    run("xz -dc %s | tar -xvpf -" % tarball)
    run("rm -f %s" % tarball)
