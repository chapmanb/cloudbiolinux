#!/usr/bin/env python
"""Prepare GFF transcript files for use as input to RNA-seq pipelines

Usage, from within the main genome directory of your organism:
  prepare_tx_gff.py <org_build>

requires these python packages which may not be installed
---------------------------------------------------------
mysql-python (via conda)
pandas (via conda)

"""
import os
import sys
import shutil
import collections
import datetime
import subprocess
import tempfile
import glob
from argparse import ArgumentParser

import gffutils

try:
    import MySQLdb
except:
    MySQLdb = None


from bcbio.utils import chdir, safe_makedir, file_exists


# ##  Version and retrieval details for Ensembl and UCSC
ensembl_release = "74"
base_ftp = "ftp://ftp.ensembl.org/pub/release-{release}/gtf"

# taxname:
# biomart_name: name of ensembl gene_id on biomart
# ucsc_map:
# fbase: the base filename for ensembl files using this genome

Build = collections.namedtuple("Build", ["taxname", "biomart_name",
                                         "ucsc_map", "fbase"])

def ucsc_ensembl_map_via_download(org_build):
    ensembl_dict_file = get_ensembl_dict(org_build)
    ucsc_dict_file = get_ucsc_dict(org_build)
    ensembl_dict = parse_sequence_dict(ensembl_dict_file)
    ucsc_dict = parse_sequence_dict(ucsc_dict_file)
    return ensembl_to_ucsc(ensembl_dict, ucsc_dict)

def ensembl_to_ucsc(ensembl_dict, ucsc_dict):
    name_map = {}
    for md5, name in ensembl_dict.items():
        name_map[name] = ucsc_dict.get(md5, None)
    return name_map

def ucsc_ensembl_map_via_query(org_build):
    """Retrieve UCSC to Ensembl name mappings from UCSC MySQL database.
    """
    # if MySQLdb is not installed, figure it out via download
    if not MySQLdb:
        return ucsc_ensembl_map_via_download(org_build)

    db = MySQLdb.connect(host=ucsc_db, user=ucsc_user, db=org_build)
    cursor = db.cursor()
    cursor.execute("select * from ucscToEnsembl")
    ucsc_map = {}
    for ucsc, ensembl in cursor.fetchall():
        # workaround for GRCh37/hg19 additional haplotype contigs.
        # Coordinates differ between builds so do not include these regions.
        if org_build == "hg19" and "hap" in ucsc:
            continue
        else:
            ucsc_map[ensembl] = ucsc
    return ucsc_map


build_info = {
    "hg19": Build("homo_sapiens", "hsapiens_gene_ensembl",
                  ucsc_ensembl_map_via_query,
                  "Homo_sapiens.GRCh37." + ensembl_release),
    "mm9": Build("mus_musculus", "mmusculus_gene_ensembl",
                 ucsc_ensembl_map_via_query,
                 "Mus_musculus.NCBIM37.67"),
    "mm10": Build("mus_musculus", "mmusculus_gene_ensembl",
                  ucsc_ensembl_map_via_query,
                  "Mus_musculus.GRCm38." + ensembl_release),
    "rn5": Build("rattus_norvegicus", None,
                 ucsc_ensembl_map_via_download,
                 "Rattus_norvegicus.Rnor_5.0." + ensembl_release),
    "GRCh37": Build("homo_sapiens", "hsapiens_gene_ensembl",
                    None,
                    "Homo_sapiens.GRCh37." + ensembl_release)}

ucsc_db = "genome-mysql.cse.ucsc.edu"
ucsc_user = "genome"


def parse_sequence_dict(fasta_dict):
    def _tuples_from_line(line):
        name = line.split("\t")[1].split(":")[1]
        md5 = line.split("\t")[4].split(":")[1]
        return md5, name
    with open(fasta_dict) as dict_handle:
        tuples = [_tuples_from_line(x) for x in dict_handle if "@SQ" in x]
        md5_dict = {x[0]: x[1] for x in tuples}
    return md5_dict

class SequenceDictParser(object):

    def __init__(self, fname):
        self.fname = fname

    def _get_sequences_in_genome_dict(self):
        with open(self.fname) as genome_handle:
            sequences = [self._sequence_from_line(x) for x in genome_handle if "@SQ" in x]
        return sequences

    def _sequence_from_line(self, line):
        name = line.split("\t")[1].split(":")[1]
        md5 = line.split("\t")[4].split(":")[1]
        return md5, name


def get_ensembl_dict(org_build):
    genome_dict = org_build + ".dict"
    if not os.path.exists(genome_dict):
        genome = _download_ensembl_genome(org_build)
        org_fa = org_build + ".fa"
        shutil.move(genome, org_fa)
        genome_dict = make_fasta_dict(org_fa)
    return genome_dict

def get_ucsc_dict(org_build):
    fa_dict = os.path.join(os.getcwd(), os.pardir, "seq", org_build + ".dict")
    if not file_exists(fa_dict):
        fa_file = os.path.splitext(fa_dict)[0] + ".fa"
        fa_dict = make_fasta_dict(fa_file)
    return fa_dict


def make_fasta_dict(fasta_file):
    dict_file = os.path.splitext(fasta_file)[0] + ".dict"
    if not os.path.exists(dict_file):
        picard_jar = os.path.join(PICARD_DIR, "CreateSequenceDictionary.jar")
        subprocess.check_call("java -jar {picard_jar} R={fasta_file} "
                              "O={dict_file}".format(**locals()), shell=True)
    return dict_file


def _download_ensembl_genome(org_build):
    build = build_info[org_build]
    fname = build.fbase + ".dna_sm.toplevel.fa.gz"
    dl_url = ("ftp://ftp.ensembl.org/pub/release-{release}/"
                   "fasta/{taxname}/dna/{fname}").format(release=ensembl_release,
                                                         taxname=build.taxname,
                                                         fname=fname)
    out_file = os.path.splitext(os.path.basename(dl_url))[0]
    if not os.path.exists(out_file):
        subprocess.check_call(["wget", dl_url])
        subprocess.check_call(["gunzip", os.path.basename(dl_url)])
    return out_file

def prepare_gff_db(gff_file):
    """
    make a database of a GTF file with gffutils
    """
    dbfn = gff_file + ".db"
    if not os.path.exists(dbfn):
        db = gffutils.create_db(gff_file, dbfn=dbfn, keep_order=False,
                                merge_strategy='merge', force=False,
                                infer_gene_extent=False)
    return dbfn

# ## Main driver functions

def main(org_build):
    work_dir = os.path.join(os.getcwd(), org_build, "tmpcbl")
    out_dir = os.path.join(os.getcwd(), org_build,
                           "rnaseq-%s" % datetime.datetime.now().strftime("%Y-%m-%d"))
    tophat_dir = os.path.join(out_dir, "tophat")
    safe_makedir(work_dir)
    with chdir(work_dir):
        build = build_info[org_build]
        tx_gff = prepare_tx_gff(build, org_build)
        db = prepare_gff_db(tx_gff)
        gtf_to_refflat(tx_gff)
        mask_gff = prepare_mask_gtf(tx_gff)
        rrna_gtf = prepare_rrna_gtf(tx_gff)
        gtf_to_interval(rrna_gtf, org_build)
        prepare_tophat_index(tx_gff, org_build)
        cleanup(work_dir, out_dir, org_build)
    tar_dirs = [out_dir]
    upload_to_s3(tar_dirs, org_build)

def cleanup(work_dir, out_dir, org_build):
    try:
        os.remove(os.path.join(work_dir, org_build + ".dict"))
        os.remove(os.path.join(work_dir, org_build + ".fa"))
    except:
        pass
    shutil.move(work_dir, out_dir)

def upload_to_s3(tar_dirs, org_build):
    str_tar_dirs = " ".join(os.path.relpath(d) for d in tar_dirs)
    tarball = "{org}-{dir}.tar.xz".format(org=org_build, dir=os.path.basename(tar_dirs[0]))
    if not os.path.exists(tarball):
        subprocess.check_call("tar -cvpf - {out_dir} | xz -zc - > {tarball}".format(
            out_dir=str_tar_dirs, tarball=tarball), shell=True)
    upload_script = os.path.join(os.path.dirname(__file__), "s3_multipart_upload.py")
    subprocess.check_call([sys.executable, upload_script, tarball, "biodata",
                           os.path.join("annotation", os.path.basename(tarball)),
                           "--public"])

def genepred_to_UCSC_table(genepred):
    header = ["#bin", "name", "chrom", "strand",
              "txStart", "txEnd", "cdsStart", "cdsEnd",
              "exonCount", "exonStarts", "exonEnds", "score",
              "name2", "cdsStartStat", "cdsEndStat",
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
        if not file_exists(prefix):
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
    make_large_exons_gtf(gtf)
    shutil.rmtree(out_dir)
    os.remove(fastq)


def make_large_exons_gtf(gtf_file):
    """
    Save all exons > 1000 bases to a separate file for estimating the
    insert size distribution
    """
    out_dir = os.path.abspath(os.path.join(os.path.dirname(gtf_file), "tophat"))
    out_file = os.path.join(out_dir, "large_exons.gtf")

    if file_exists(out_file):
        return out_file

    dbfn = gtf_file + ".db"
    if not file_exists(dbfn):
        db = gffutils.create_db(gtf_file, dbfn=dbfn, keep_order=True,
                                merge_strategy='merge', force=False,
                                infer_gene_extent=False)
    else:
        db = gffutils.FeatureDB(dbfn)
    processed_count = 0
    kept_exons = []
    for exon in db.features_of_type('exon'):
        processed_count += 1
        if processed_count % 10000 == 0:
            print("Processed %d exons." % processed_count)
        if exon.end - exon.start > 1000:
            kept_exons.append(exon)

    with open(out_file, "w") as out_handle:
        print("Writing %d large exons to %s." % (processed_count,
                                                 out_file))
        for exon in kept_exons:
            out_handle.write(str(exon) + "\n")
    return out_file


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
    fa_dict = get_ucsc_dict(build)
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
                    "tRNA", "Mt_tRNA"]
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

def prepare_tx_gff(build, org_name):
    """Prepare UCSC ready transcript file given build information.
    """
    ensembl_gff = _download_ensembl_gff(build)
    # if we need to do the name remapping
    if build.ucsc_map:
        ucsc_name_map = build.ucsc_map(org_name)
        tx_gff = _remap_gff(ensembl_gff, ucsc_name_map)
        os.remove(ensembl_gff)
    else:
        tx_gff = "ref-transcripts.gtf"
        os.rename(ensembl_gff, tx_gff)
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

def _download_ensembl_gff(build):
    """Given build details, download and extract the relevant ensembl GFF.
    """
    fname = build.fbase + ".gtf.gz"
    dl_url = "/".join([base_ftp, build.taxname, fname]).format(release=ensembl_release)
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
    parser = ArgumentParser(description="Prepare the transcriptome files for an "
                            "organism.")
    parser.add_argument("picard",
                        help="Path to Picard")
    parser.add_argument("org_build", help="Build of organism to run.",
                        choices=build_info.keys())
    args = parser.parse_args()
    global PICARD_DIR
    PICARD_DIR = args.picard
    main(args.org_build)
