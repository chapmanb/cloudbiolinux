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

def download_dbsnp(genomes, bundle_version, dbsnp_version):
    """Download and install dbSNP variation data for supplied genomes.
    """
    folder_name = "variation"
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in ((o, g, m) for (o, g, m) in genomes
                                    if m.config.get("dbsnp", False)):
        vrn_dir = os.path.join(genome_dir, orgname, gid, folder_name)
        if not env.safe_exists(vrn_dir):
            env.safe_run('mkdir -p %s' % vrn_dir)
        with cd(vrn_dir):
            if gid in ["GRCh37", "hg19"]:
                _dbsnp_human(env, gid, manager, bundle_version, dbsnp_version)
            elif gid in ["mm10", "canFam3"]:
                _dbsnp_custom(env, gid)
                _download_lcrs_custom(env, gid)

def _dbsnp_custom(env, gid):
    """Retrieve resources for dbsnp builds from custom S3 biodata bucket.
    """
    remote_dir = "https://s3.amazonaws.com/biodata/variants/"
    files = {"mm10": ["mm10-dbSNP-2013-09-12.vcf.gz"],
             "canFam3": ["canFam3-dbSNP-2014-05-10.vcf.gz"]}
    for f in files[gid]:
        for ext in ["", ".tbi"]:
            fname = f + ext
            if not env.safe_exists(fname):
                shared._remote_fetch(env, "%s%s" % (remote_dir, fname))

def _dbsnp_human(env, gid, manager, bundle_version, dbsnp_version):
    """Retrieve resources for human variant analysis from Broad resource bundles.
    """
    to_download = ["dbsnp_{ver}".format(ver=dbsnp_version),
                   "hapmap_3.3",
                   "1000G_omni2.5",
                   "1000G_phase1.snps.high_confidence",
                   "Mills_and_1000G_gold_standard.indels"]
    for dl_name in to_download:
        for ext in [""]:
            _download_broad_bundle(manager.dl_name, bundle_version, dl_name, ext)
    _download_cosmic(gid)
    _download_repeats(gid)
    _download_dbnsfp(env, gid, manager.config)
    _download_ancestral(env, gid, manager.config)
    _download_qsignature(env, gid, manager.config)
    # XXX Wait to get this by default until it is used more widely
    #_download_background_vcf(gid)

def _download_broad_bundle(gid, bundle_version, name, ext):
    # Broad bundle directories have uneven use of ".sites" in VCF files
    # only present in hg19 for non-dbSNP resources
    sites = ".sites" if gid == "hg19" and not name.startswith("dbsnp") else ""
    broad_fname = "{name}.{gid}{sites}.vcf{ext}".format(gid=gid, name=name, sites=sites, ext=ext)
    fname = broad_fname.replace(".{0}".format(gid), "").replace(".sites", "") + ".gz"
    base_url = "ftp://gsapubftp-anonymous:@ftp.broadinstitute.org/bundle/" + \
               "{bundle}/{gid}/{fname}.gz".format(
                   bundle=bundle_version, fname=broad_fname, gid=gid)
    # compress and prepare existing uncompressed versions
    if env.safe_exists(fname.replace(".vcf.gz", ".vcf")):
        env.safe_run("bgzip %s" % fname.replace(".vcf.gz", ".vcf"))
        env.safe_run("tabix -f -p vcf %s" % fname)
    # otherwise, download and bgzip and tabix index
    if not env.safe_exists(fname):
        out_file = shared._remote_fetch(env, base_url)
        env.safe_run("gunzip -c %s | bgzip -c > %s" % (out_file, fname))
        env.safe_run("tabix -f -p vcf %s" % fname)
        env.safe_run("rm -f %s" % out_file)
    # clean up old files
    for ext in [".vcf", ".vcf.idx"]:
        if env.safe_exists(fname.replace(".vcf.gz", ext)):
            env.safe_run("rm -f %s" % (fname.replace(".vcf.gz", ext)))
    return fname

def _download_cosmic(gid):
    """Prepared versions of COSMIC, pre-sorted and indexed.
    utils/prepare_cosmic.py handles the work of creating the VCFs from standard
    COSMIC resources.
    """
    base_url = "https://s3.amazonaws.com/biodata/variants"
    version = "v68"
    supported = ["hg19", "GRCh37"]
    if gid in supported:
        url = "%s/cosmic-%s-%s.vcf.gz" % (base_url, version, gid)
        fname = os.path.basename(url)
        if not env.safe_exists(fname):
            shared._remote_fetch(env, url)
        if not env.safe_exists(fname + ".tbi"):
            shared._remote_fetch(env, url + ".tbi")

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

def _download_ancestral(env, gid, gconfig):
    """Download ancestral genome sequence for loss of function evaluation.

    Used by LOFTEE VEP plugin: https://github.com/konradjk/loftee
    """
    base_url = "https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
    if gid == "GRCh37" or (gid == "hg19" and not env.safe_exists("../../GRCh37")):
        for ext in ["", ".fai", ".gzi"]:
            outfile = os.path.basename(base_url) + ext
            if not env.safe_exists(outfile):
                shared._remote_fetch(env, base_url + ext, samedir=True)
    elif gid == "hg19":  # symlink to GRCh37 download
        for ext in ["", ".fai", ".gzi"]:
            outfile = os.path.basename(base_url) + ext
            if not env.safe_exists(outfile):
                env.safe_run("ln -sf ../../GRCh37/variation/%s %s" % (outfile, outfile))

def _download_qsignature(env, gid, gconfig):
    """Download qsignature position file to detect samples problems

    :param env
    :param gid: str genome id
    :param gconfig: 

    :returns: NULL
    """
    base_url = "http://downloads.sourceforge.net/project/adamajava/qsignature.tar.bz2"
    outfile = "qsignature.vcf"
    if gid == "GRCh37" or (gid == "hg19" and not env.safe_exists("../../GRCh37")):
        if not env.safe_exists(outfile):
            zipfile = shared._remote_fetch(env, base_url, samedir=True)
            outdir = "qsignature"
            env.safe_run("mkdir -p %s" % outdir)
            env.safe_run("tar -jxf %s -C %s" % (zipfile, outdir))
            env.safe_run("mv %s/qsignature_positions.txt %s" % (outdir, outfile))
            env.safe_run("rm -rf %s" % outdir)
            env.safe_run("rm -rf %s" % zipfile)
    elif gid == "hg19":  # symlink to GRCh37 download
        if not env.safe_exists(outfile):
            env.safe_run("ln -sf ../../GRCh37/variation/%s %s" % (outfile, outfile))

def _download_background_vcf(gid):
    """Download background file of variant to use in calling.
    """
    base_url = "https://s3.amazonaws.com/biodata/variants"
    base_name = "background-diversity-1000g.vcf"
    if gid in ["GRCh37"] and not env.safe_exists("{0}.gz".format(base_name)):
        for ext in ["gz", "gz.tbi"]:
            shared._remote_fetch(env, "{0}/{1}.{2}".format(base_url, base_name, ext))

def _download_repeats(gid):
    _download_sv_repeats(gid)
    _download_lcrs(gid)

def _download_sv_repeats(gid):
    """Retrieve telomere and centromere exclusion regions for structural variant calling.
    From Delly: https://github.com/tobiasrausch/delly
    """
    mere_url = "https://raw.githubusercontent.com/chapmanb/delly/master/human.hg19.excl.tsv"
    out_file = "sv_repeat_telomere_centromere.bed"
    if not env.safe_exists(out_file):
        def _select_by_gid(env, orig_file):
            if gid == "hg19":
                env.safe_run("grep ^chr %s > %s" % (orig_file, out_file))
            else:
                assert gid == "GRCh37"
                env.safe_run("grep -v ^chr %s > %s" % (orig_file, out_file))
            return out_file
        shared._remote_fetch(env, mere_url, fix_fn=_select_by_gid)

def _download_lcrs(gid):
    """Retrieve low complexity regions from Heng Li's variant analysis paper.
    """
    lcr_url = "https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz"
    out_file = "LCR.bed.gz"
    if not env.safe_exists(out_file):
        def _fix_chrom_names(env, orig_file):
            if gid == "hg19":
                convert_cmd = "| grep -v ^GL | grep -v ^NC | grep -v ^hs | sed 's/^/chr/'"
            else:
                assert gid == "GRCh37"
                convert_cmd = ""
            env.safe_run("zcat %s %s | bgzip -c > %s" % (orig_file, convert_cmd, out_file))
            return out_file
        shared._remote_fetch(env, lcr_url, fix_fn=_fix_chrom_names)
        env.safe_run("tabix -p vcf -f %s" % out_file)

def _download_lcrs_custom(env, gid):
    """Retrieve low complexity regions from other sources.

    mm10 from Brent Pedersen: http://figshare.com/articles/LCR_mm10_bed_gz/1180124
    """
    urls = {"mm10": "http://files.figshare.com/1688228/LCR_mm10.bed.gz"}
    out_file = "LCR.bed.gz"
    cur_url = urls.get(gid)
    if cur_url and not env.safe_exists(out_file):
        def _bgzip_file(env, orig_file):
            env.safe_run("zcat %s | bgzip -c > %s" % (orig_file, out_file))
            return out_file
        shared._remote_fetch(env, cur_url, fix_fn=_bgzip_file)
        env.safe_run("tabix -p vcf -f %s" % out_file)
