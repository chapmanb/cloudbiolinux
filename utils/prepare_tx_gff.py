#!/usr/bin/env python
"""Prepare GFF transcript files for use as input to RNA-seq pipelines

Usage, from within the main genome directory of your organism:
  prepare_tx_gff.py <org_build>

requires these python packages which may not be installed
---------------------------------------------------------
rnaseqlib
git@github.com:yarden/rnaseqlib.git


gffutils (the refactor branch):
https://github.com/daler/gffutils/tree/refactor

MISO:
git@github.com:yarden/MISO.git

mysql-python (via pip)


"""
import os
import sys
import shutil
import platform
import collections
import subprocess
import tempfile
import glob

import MySQLdb
import rnaseqlib.utils as utils
import rnaseqlib.events.defineEvents as def_events
import gffutils
import time


from bcbio.utils import chdir, safe_makedir, file_exists

def build_ucsc_map(ensembl_chrs, mt_chr):
    """Build mapping of Ensembl to UCSC names for standard chromosomes.
    """
    out = {}
    for c in ensembl_chrs:
        out[str(c)] = "chr{0}".format(c)
    out[mt_chr] = "chrM"
    return out

# ##  Version and retrieval details for Ensembl and UCSC
ensembl_release = "73"
base_ftp = "ftp://ftp.ensembl.org/pub/release-{release}/gtf"

Build = collections.namedtuple("Build", ["taxname", "fname", "biomart_name",
                                         "ucsc_map"])
build_info = {
    "hg19": Build("homo_sapiens", "Homo_sapiens.GRCh37.{release}.gtf.gz",
                  "hsapiens_gene_ensembl", None),
    "mm9": Build("mus_musculus", "Mus_musculus.NCBIM37.67.gtf.gz",
                 "mmusculus_gene_ensembl", None),
    "mm10": Build("mus_musculus", "Mus_musculus.GRCm38.{release}.gtf.gz",
                 "mmusculus_gene_ensembl", None)}

ucsc_db= "genome-mysql.cse.ucsc.edu"
ucsc_user="genome"

# ## Main driver functions

def main(org_build):
    work_dir = os.path.join(os.getcwd(), org_build, "tmpcbl")
    out_dir = os.path.join(os.getcwd(), org_build, "rnaseq")
    tophat_dir = os.path.join(out_dir, "tophat")
    safe_makedir(work_dir)
    with chdir(work_dir):
        build = build_info[org_build]
        tx_gff = prepare_tx_gff(build, org_build)
        gtf_to_refflat(tx_gff)
        mask_gff = prepare_mask_gtf(tx_gff)
        rrna_gtf = prepare_rrna_gtf(tx_gff)
        gtf_to_interval(rrna_gtf, org_build)
        make_miso_events(tx_gff, org_build)
        prepare_tophat_index(tx_gff, org_build)
        cleanup(work_dir, out_dir)
    tar_dirs = [out_dir]
    upload_to_s3(tar_dirs, org_build)

def cleanup(work_dir, out_dir):
    db_files = glob.glob(os.path.join(work_dir, "*.db"))
    map(os.remove, db_files)
    shutil.move(work_dir, out_dir)

def upload_to_s3(tar_dirs, org_build):
    tar_dirs = " ".join(tar_dirs)
    tarball = "{org}-rnaseq.tar.xz".format(org=org_build)
    if not os.path.exists(tarball):
        subprocess.check_call("tar -cvpf - {out_dir} | xz -zc - > {tarball}".format(
            out_dir=tar_dirs, tarball=tarball), shell=True)
    python_exe = "python{0}.{1}".format(*platform.python_version_tuple()[:2])
    upload_script = os.path.join(os.path.dirname(__file__), "s3_multipart_upload.py")
    subprocess.check_call([python_exe, upload_script, tarball, "biodata",
                           os.path.join("annotation", os.path.basename(tarball)),
                           "--public"])


def make_miso_annotation(tables_dir, output_dir, org_build):
    """
    Make GFF annotation. Takes GFF tables directory
    and an output directory.

    Adapted from
    https://github.com/yarden/rnaseqlib/
    """
    tables_dir = utils.pathify(tables_dir)
    output_dir = utils.pathify(output_dir)
    print "Making GFF alternative events annotation..."
    print " - UCSC tables read from: %s" % (tables_dir)
    print " - Output dir: %s" % (output_dir)
    t1 = time.time()
    table_fnames = def_events.load_ucsc_tables(tables_dir)
    num_tables = len(table_fnames)
    if num_tables == 0:
        raise Exception("No UCSC tables found in %s." % (tables_dir))
    print "Loaded %d UCSC tables." % (num_tables)
    def_events.defineAllSplicing(tables_dir, output_dir,
                                 flanking="commonshortest",
                                 multi_iso=False,
                                 sanitize=False,
                                 genome_label=org_build)
    t2 = time.time()
    print "Took %.2f minutes to make the annotation." \
        % ((t2 - t1)/60.)


def genepred_to_UCSC_table(genepred):
    header = ["#bin", "name", "chrom", "strand",
              "txStart", "txEnd", "cdsStart", "cdsEnd",
              "exonCount", "exonStarts", "exonEnds", "score",
              "name2",	"cdsStartStat",	"cdsEndStat",
              "exonFrames"]
    out_file = os.path.splitext(genepred)[0] + ".UCSCTable"
    if file_exists(out_file):
        return out_file
    with open(genepred) as in_handle, open(out_file, "w") as out_handle:
        counter = -1
        current_item = None
        out_handle.write("\t".join(header) + "\n")
        for l in in_handle:
            item = l.split("\t")[0]
            if current_item != item:
                current_item = item
                counter = counter + 1
            out_handle.write("\t".join([str(counter), l]))
    return out_file

def gtf_to_genepred(gtf):
    out_file = os.path.splitext(gtf)[0] + ".genePred"
    if file_exists(out_file):
        return out_file

    cmd = "gtfToGenePred -allErrors -genePredExt {gtf} {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

def gtf_to_refflat(gtf):
    out_file = os.path.splitext(gtf)[0] + ".refFlat"
    if file_exists(out_file):
        return out_file

    genepred = gtf_to_genepred(gtf)
    with open(genepred) as in_handle, open(out_file, "w") as out_handle:
        for l in in_handle:
            first = l.split("\t")[0]
            out_handle.write("\t".join([first, l]))

    return out_file

def make_miso_events(gtf, org_build):

    genepred = gtf_to_genepred(gtf)
    genepred = genepred_to_UCSC_table(genepred)
    pred_dir = tempfile.mkdtemp()
    miso_dir = os.path.join(os.path.dirname(gtf), "miso")
    tmp_pred = os.path.join(pred_dir, "ensGene.txt")
    os.symlink(os.path.abspath(genepred), tmp_pred)
    make_miso_annotation(pred_dir, miso_dir, org_build)

    gff_files = glob.glob(os.path.join(miso_dir, "commonshortest", "*.gff3"))

    cmd = "index_gff --index {f} {prefix}"

    for f in gff_files:
        prefix = f.split(".")[0] + "_indexed"
        print prefix
        print f
        print cmd.format(**locals())
        subprocess.check_call(cmd.format(**locals()), shell=True)

def prepare_tophat_index(gtf, org_build):
    tophat_dir = os.path.abspath(os.path.join(os.path.dirname(gtf), "tophat",
                                              org_build + "_transcriptome"))
    bowtie_dir = os.path.abspath(os.path.join(os.path.dirname(gtf),
                                              os.path.pardir, "bowtie2",
                                              org_build))
    out_dir = tempfile.mkdtemp()
    fastq = _create_dummy_fastq()
    cmd = ("tophat --transcriptome-index {tophat_dir} -G {gtf} "
           "-o {out_dir} {bowtie_dir} {fastq}")
    subprocess.check_call(cmd.format(**locals()), shell=True)
    shutil.rmtree(out_dir)
    os.remove(fastq)

def _create_dummy_fastq():
    read = ("@HWI-ST333_0178_FC:5:1101:1107:2112#ATCTCG/1\n"
            "GGNCTTTCCTGCTTCTATGTCTTGATCGCCTGTAGGCAGG\n"
            "+HWI-ST333_0178_FC:5:1101:1107:2112#ATCTCG/1\n"
            "[[BS\\a`ceeagfhhhhhaefhcdfhcf`efeg[cg_b__\n")
    fn = "dummy.fq"
    with open(fn, "w") as out_handle:
        out_handle.write(read)
    return fn

def gtf_to_interval(gtf, build):
    fa_dict = os.path.join(os.getcwd(), os.pardir, "seq", build + ".dict")
    if not file_exists(fa_dict):
        raise IOError("%s is not found, please make with "
                      "CreateSequenceDictionary.")
    db = _get_gtf_db(gtf)
    out_file = os.path.splitext(gtf)[0] + ".interval_list"
    if file_exists(out_file):
        return out_file

    with open(out_file, "w") as out_handle:
        with open(fa_dict) as in_handle:
            for l in in_handle:
                out_handle.write(l)

        for l in db.all_features():
            out_handle.write("\t".join([str(l.seqid), str(l.start),
                                        str(l.end), str(l.strand),
                                        str(l.attributes.get("transcript_id",
                                                             ["."])[0])]) + "\n")
    return out_file

def prepare_mask_gtf(gtf):
    """
    make a mask file of usually-masked RNA biotypes
    """

    mask_biotype = ["rRNA", "Mt_rRNA", "misc_RNA", "snRNA", "snoRNA",
                    "tRNA","Mt_tRNA"]
    mask_chrom = ["MT"]
    out_file = os.path.join(os.path.dirname(gtf), "ref-transcripts-mask.gtf")
    if file_exists(out_file):
        return out_file

    db = _get_gtf_db(gtf)
    with open(out_file, "w") as out_handle:
        for g in db.all_features():
            biotype = g.attributes.get("gene_biotype", None)
            if ((biotype and biotype[0] in mask_biotype) or
               g.chrom in mask_chrom):
                out_handle.write(str(g) + "\n")
    return out_file

def prepare_rrna_gtf(gtf):
    """
    extract out just the rRNA biotypes, for assessing rRNA contamination
    """
    mask_biotype = ["rRNA", "Mt_rRNA", "tRNA", "MT_tRNA"]

    out_file = os.path.join(os.path.dirname(gtf), "rRNA.gtf")
    if os.path.exists(out_file):
        return out_file

    db = _get_gtf_db(gtf)

    with open(out_file, "w") as out_handle:
        for g in db.all_features():
            biotype = g.attributes.get("gene_biotype", None)
            if biotype and biotype[0] in mask_biotype:
                out_handle.write(str(g) + "\n")

    return out_file

def gtf_to_genepred(gtf):
    out_file = os.path.splitext(gtf)[0] + ".genePred"
    if file_exists(out_file):
        return out_file

    cmd = "gtfToGenePred -allErrors -genePredExt {gtf} {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
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
      filter.types <- c("Mt_rRNA", "rRNA")
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
    ucsc_name_map = (_query_for_ucsc_ensembl_map(org_name)
                     if build.ucsc_map is None else build.ucsc_map)
    tx_gff = _remap_gff(ensembl_gff, ucsc_name_map)
    os.remove(ensembl_gff)
    return tx_gff

def _remap_gff(base_gff, name_map):
    """Remap chromosome names to UCSC instead of Ensembl
    """
    out_file = "ref-transcripts.gtf"
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

def _get_gtf_db(gtf):
    db_file = gtf + ".db"
    if not file_exists(db_file):
        gffutils.create_db(gtf, dbfn=db_file)

    return gffutils.FeatureDB(db_file)

if __name__ == "__main__":
    main(*sys.argv[1:])
