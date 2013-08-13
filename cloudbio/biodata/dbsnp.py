"""Download variation data from dbSNP and install within directory structure.

Uses Broad's GATK resource bundles:

 http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle

Retrieves dbSNP plus training data for variant recalibration:
  - dbsnp_132.hg19.vcf.gz
  - hapmap_3.3.hg19.sites.vcf
  - 1000G_omni2.5.hg19.sites.vcf
  - Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

For MuTect and cancer calling:
  - cosmic
"""
import os

from fabric.api import *
from fabric.contrib.files import *

def download_dbsnp(genomes, bundle_version, dbsnp_version):
    """Download and install dbSNP variation data for supplied genomes.
    """
    folder_name = "variation"
    to_download = ["dbsnp_{ver}".format(ver=dbsnp_version),
                   "hapmap_3.3",
                   "1000G_omni2.5",
                   "Mills_and_1000G_gold_standard.indels"]
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("dbsnp", False)):
        vrn_dir = os.path.join(genome_dir, orgname, gid, folder_name)
        if not env.safe_exists(vrn_dir):
            env.safe_run('mkdir -p %s' % vrn_dir)
        with cd(vrn_dir):
            for dl_name in to_download:
                for ext in ["", ".idx"]:
                    _download_broad_bundle(manager.dl_name, bundle_version, dl_name, ext)
            _download_cosmic(gid)
            # XXX Wait to get this by default until it is used more widely
            #_download_background_vcf(gid)

def _download_broad_bundle(gid, bundle_version, name, ext):
    broad_fname = "{name}.{gid}.vcf{ext}".format(gid=gid, name=name, ext=ext)
    fname = broad_fname.replace(".{0}".format(gid), "").replace(".sites", "")
    base_url = "ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/bundle/" + \
               "{bundle}/{gid}/{fname}.gz".format(
                   bundle=bundle_version, fname=broad_fname, gid=gid)
    if not env.safe_exists(fname):
        with warn_only():
            dl = env.safe_run("wget -c %s" % base_url)
        if dl.succeeded:
            env.safe_run("gunzip %s" % os.path.basename(base_url))
            env.safe_run("mv %s %s" % (broad_fname, fname))
        else:
            env.logger.warn("dbSNP resources not available for %s" % gid)
    return fname

def _download_cosmic(gid):
    base_url = "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/"
    base_name = "b37_cosmic_v54_120711.vcf"
    if gid in ["GRCh37"] and not env.safe_exists(base_name):
        env.safe_run("wget -c {0}/{1}".format(base_url, base_name))

def _download_background_vcf(gid):
    """Download background file of variant to use in calling.
    """
    base_url = "https://s3.amazonaws.com/biodata/variants"
    base_name = "background-diversity-1000g.vcf"
    if gid in ["GRCh37"] and not env.safe_exists("{0}.gz".format(base_name)):
        for ext in ["gz", "gz.tbi"]:
            env.safe_run("wget -c {0}/{1}.{2}".format(base_url, base_name, ext))
