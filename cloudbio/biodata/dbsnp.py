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

For structural variant calling and SNP/indel filtering
  - low complexity regions
  - centromere and telomere regions
"""
import os

from fabric.api import env
from fabric.contrib.files import cd

from cloudbio.custom import shared

def download_dbnsfp(genomes):
    """Back compatible download target for dbNSFP, to be moved to GGD recipes.
    """
    folder_name = "variation"
    genome_dir = os.path.join(env.data_files, "genomes")
    gids = set(["hg19", "GRCh37"])
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if g in gids and m.config.get("dbnsfp")):
        vrn_dir = os.path.join(genome_dir, orgname, gid, folder_name)
        if not env.safe_exists(vrn_dir):
            env.safe_run('mkdir -p %s' % vrn_dir)
        with cd(vrn_dir):
            _download_dbnsfp(env, gid, manager.config)

def _download_dbnsfp(env, gid, gconfig):
    """Download and prepare dbNSFP functional prediction resources if configured.

    Feeds into VEP for annotating VCF files:
    https://sites.google.com/site/jpopgen/dbNSFP
    https://github.com/ensembl-variation/VEP_plugins/blob/master/dbNSFP.pm
    """
    version = "3.0b2c"
    url = "ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv%s.zip" % version
    dl_file = "dbNSFPv%s.zip" % version
    if gconfig.get("dbnsfp"):
        outfile = "dbNSFP_v%s.gz" % (version)
        if gid == "GRCh37" or (gid == "hg19" and not env.safe_exists("../../GRCh37")):
            if not env.safe_exists(outfile):
                zipfile = shared._remote_fetch(env, url, out_file=dl_file, samedir=True)
                outdir = "dbNSFPv%s" % version
                env.safe_run("mkdir -p %s" % outdir)
                env.safe_run("7za x %s -y -o%s" % (zipfile, outdir))
                env.safe_run("cat %s/dbNSFP*_variant.chr* | bgzip -c > %s" % (outdir, outfile))
                env.safe_run("rm -f %s/* && rmdir %s" % (outdir, outdir))
                env.safe_run("rm -f %s" % (zipfile))
            if not env.safe_exists(outfile + ".tbi"):
                env.safe_run("tabix -s 1 -b 2 -e 2 -c '#' %s" % outfile)
        elif gid == "hg19":  # symlink to GRCh37 download
            if not env.safe_exists(outfile):
                env.safe_run("ln -sf ../../GRCh37/variation/%s %s" % (outfile, outfile))
            if not env.safe_exists(outfile + ".tbi"):
                env.safe_run("ln -sf ../../GRCh37/variation/%s.tbi %s.tbi" % (outfile, outfile))
