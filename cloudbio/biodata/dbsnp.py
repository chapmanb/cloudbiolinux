"""Download variation data from dbSNP and install within directory structure.

Uses Broad's GATK resource bundles:

 http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle

Retrieves dbSNP plus training data for variant recalibration:
  - dbsnp_132.hg19.vcf.gz
  - hapmap_3.3.hg19.sites.vcf
  - 1000G_omni2.5.hg19.sites.vcf
  - Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
"""
import os

from fabric.api import *
from fabric.contrib.files import *

def download_dbsnp(genomes, bundle_version, dbsnp_version):
    """Download and install dbSNP variation data for supplied genomes.
    """
    folder_name = "variation"
    to_download = [("dbsnp_{ver}".format(ver=dbsnp_version), ""),
                   ("hapmap_3.3", ".sites"),
                   ("1000G_omni2.5", ".sites"),
                   ("Mills_and_1000G_gold_standard.indels", ".sites")]
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("dbsnp", False)):
        vrn_dir = os.path.join(genome_dir, orgname, gid, folder_name)
        if not exists(vrn_dir):
            run('mkdir -p %s' % vrn_dir)
        with cd(vrn_dir):
            for dl_name, dl_ext in to_download:
                _download_broad_bundle(manager.dl_name, bundle_version, dl_name, dl_ext)
            _download_background_vcf(gid)

def _download_broad_bundle(gid, bundle_version, name, ext):
    broad_fname = "{name}.{gid}{ext}.vcf".format(gid=gid, name=name, ext=ext)
    fname = broad_fname.replace(".{0}".format(gid), "").replace(".sites", "")
    base_url = "ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/bundle/" + \
               "{bundle}/{gid}/{fname}.gz".format(
                   bundle=bundle_version, fname=broad_fname, gid=gid)
    if not exists(fname):
        run("wget %s" % base_url)
        run("gunzip %s" % os.path.basename(base_url))
        run("mv %s %s" % (broad_fname, fname))
    return fname

def _download_background_vcf(gid):
    """Download background file of variant to use in calling.
    """
    base_url = "https://s3.amazonaws.com/biodata/variants/"
    base_name = "background-diversity-1000g.vcf"
    if gid in ["GRCh37"] and not exists("{0}.gz".format(base_name)):
        for ext in ["gz", "gz.tbi"]:
            run("wget {0}/{1}.{2}".format(base_url, base_name, ext))
