"""Prepare supplemental files to work with RNA-seq transcriptome experiments.

Retrieves annotations in a format usable by Cufflinks, using repositories
provided by Illumina:

http://cufflinks.cbcb.umd.edu/igenomes.html
"""
import os

from fabric.api import cd
from cloudbio.fabutils import warn_only


VERSIONS = {"GRCh37": "-2013-08-21"}

def download_transcripts(genomes, env):
    folder_name = "rnaseq"
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("rnaseq", False)):
        version = VERSIONS.get(gid, "")
        base_url = "https://s3.amazonaws.com/biodata/annotation/{gid}-rnaseq{version}.tar.xz"
        org_dir = os.path.join(genome_dir, orgname)
        tx_dir = os.path.join(org_dir, gid, folder_name)
        version_dir = tx_dir + version
        if not env.safe_exists(version_dir):
            with cd(org_dir):
                _download_annotation_bundle(env, base_url.format(gid=gid, version=version), gid)
                if version:
                    _symlink_version(env, tx_dir, version_dir)

def _symlink_version(env, tx_dir, version_dir):
    """Symlink the expected base output directory to our current version.
    """
    if env.safe_exists(tx_dir):
        env.safe_run("rm -rf %s" % tx_dir)
    env.safe_run("ln -s %s %s" % (version_dir, tx_dir))

def _download_annotation_bundle(env, url, gid):
    """Download bundle of RNA-seq data from S3 biodata/annotation
    """
    tarball = os.path.basename(url)
    if not env.safe_exists(tarball):
        with warn_only():
            env.safe_run("wget %s" % url)
    if env.safe_exists(tarball):
        env.safe_run("xz -dc %s | tar -xvpf -" % tarball)
        env.safe_run("rm -f %s" % tarball)
    else:
        env.logger.warn("RNA-seq transcripts not available for %s" % gid)
