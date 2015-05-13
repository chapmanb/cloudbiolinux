#!/usr/bin/env python
"""Prepare GFF transcript files for use as input to RNA-seq pipelines

Usage, from within the main genome directory of your organism:
  prepare_tx_gff.py <org_build>

requires these python packages which may not be installed
---------------------------------------------------------
mysql-python
gffutils
requests

requires a picard installation with a `picard` shell script, installed
by homebrew-science and bcbio
"""
import csv
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
import requests

try:
    import MySQLdb
except:
    MySQLdb = None


from bcbio.utils import chdir, safe_makedir, file_exists


# ##  Version and retrieval details for Ensembl and UCSC
ensembl_release = "79"
base_ftp = "ftp://ftp.ensembl.org/pub/release-{release}/gtf"
supported_oldbuilds = {"GRCh37": "75", "hg19": "75"}
build_subsets = {"hg38-noalt": "hg38"}

ucsc_db = "genome-mysql.cse.ucsc.edu"
ucsc_user = "genome"

# Chromosome name remappings thanks to Devon Ryan
# https://github.com/dpryan79/ChromosomeMappings
manual_remaps = {"hg38":
                 "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2UCSC.txt"}

def manual_ucsc_ensembl_map(org_build):
    org_build = build_subsets.get(org_build, org_build)
    requests.packages.urllib3.disable_warnings()
    r = requests.get(manual_remaps[org_build], verify=False)
    out = {}
    for line in r.text.split("\n"):
        try:
            ensembl, ucsc = line.split()
            out[ensembl] = ucsc
        except ValueError:
            pass
    return out

def ucsc_ensembl_map_via_download(org_build):
    """Compare .dict files by md5, then length to compare two builds.
    """
    ensembl_dict_file = get_ensembl_dict(org_build)
    ucsc_dict_file = get_ucsc_dict(org_build)
    ensembl_dict = parse_sequence_dict(ensembl_dict_file)
    ucsc_dict = parse_sequence_dict(ucsc_dict_file)
    return ensembl_to_ucsc(ensembl_dict, ucsc_dict, org_build)

def ensembl_to_ucsc(ensembl_dict, ucsc_dict, org_build):
    name_map = {}
    for md5, name in ensembl_dict.items():
        if ucsc_dict.get("md5"):
            name_map[name] = ucsc_dict["md5"]
    map_file = "%s-map.csv" % (org_build)
    with open(map_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["ensembl", "ucsc"])
        for md5, name in ensembl_dict.items():
            ucsc = ucsc_dict.get(md5)
            if ucsc is not None:
                writer.writerow([name, ucsc])
    return name_map

def ucsc_ensembl_map_via_query(org_build):
    """Retrieve UCSC to Ensembl name mappings from UCSC MySQL database.
    """
    org_build = build_subsets.get(org_build, org_build)
    # if MySQLdb is not installed, figure it out via download
    if not MySQLdb:
        return ucsc_ensembl_map_via_download(org_build)

    db = MySQLdb.connect(host=ucsc_db, user=ucsc_user, db=org_build)
    cursor = db.cursor()
    cursor.execute("select * from ucscToEnsembl")
    ucsc_map = {}
    for fields in cursor.fetchall():
        ucsc = fields[0]
        ensembl = fields[-1]
        # workaround for GRCh37/hg19 additional haplotype contigs.
        # Coordinates differ between builds so do not include these regions.
        if org_build == "hg19" and "hap" in ucsc:
            continue
        else:
            ucsc_map[ensembl] = ucsc
    return ucsc_map

# taxname:
# biomart_name: name of ensembl gene_id on biomart
# ucsc_map:
# fbase: the base filename for ensembl files using this genome

Build = collections.namedtuple("Build", ["taxname", "biomart_name",
                                         "ucsc_map", "fbase"])

build_info = {
    "hg19": Build("homo_sapiens", "hsapiens_gene_ensembl",
                  ucsc_ensembl_map_via_query,
                  "Homo_sapiens.GRCh37." + supported_oldbuilds["GRCh37"]),
    "GRCh37": Build("homo_sapiens", "hsapiens_gene_ensembl",
                    None,
                    "Homo_sapiens.GRCh37." + supported_oldbuilds["hg19"]),
    "mm9": Build("mus_musculus", "mmusculus_gene_ensembl",
                 ucsc_ensembl_map_via_query,
                 "Mus_musculus.NCBIM37.67"),
    "mm10": Build("mus_musculus", "mmusculus_gene_ensembl",
                  ucsc_ensembl_map_via_query,
                  "Mus_musculus.GRCm38." + ensembl_release),
    "rn5": Build("rattus_norvegicus", None,
                 ucsc_ensembl_map_via_download,
                 "Rattus_norvegicus.Rnor_5.0." + ensembl_release),
    "hg38": Build("homo_sapiens", "hsapiens_gene_ensembl",
                  manual_ucsc_ensembl_map,
                  "Homo_sapiens.GRCh38." + ensembl_release),
    "hg38-noalt": Build("homo_sapiens", "hsapiens_gene_ensembl",
                        manual_ucsc_ensembl_map,
                        "Homo_sapiens.GRCh38." + ensembl_release),
    "canFam3": Build("canis_familiaris", None,
                     ucsc_ensembl_map_via_download,
                     "Canis_familiaris.CanFam3.1." + ensembl_release),
    "sacCer3": Build("saccharomyces_cerevisiae", None,
                     ucsc_ensembl_map_via_download,
                     "Saccharomyces_cerevisiae.R64-1-1." + ensembl_release),
    "WBcel235": Build("caenorhabditis_elegans", None,
                      ucsc_ensembl_map_via_download,
                      "Caenorhabditis_elegans.WBcel235." + ensembl_release),
    "dm3": Build("drosophila_melanogaster", None,
                 ucsc_ensembl_map_via_download,
                 "Drosophila_melanogaster.BDGP5." + ensembl_release),
    "Zv9": Build("danio_rerio", None,
                 ucsc_ensembl_map_via_download,
                 "Danio_rerio.Zv9." + ensembl_release),
    "xenTro3": Build("xenopus_tropicalis", None,
                     ucsc_ensembl_map_via_download,
                     "Xenopus_tropicalis.JGI_4.2." + ensembl_release),
}


def parse_sequence_dict(fasta_dict):
    def _tuples_from_line(line):
        attrs = {}
        for tag, val in [x.split(":", 1) for x in line.strip().split("\t")[1:]]:
            attrs[tag] = val
        return attrs["SN"], attrs["LN"], attrs["M5"]
    out = {}
    with open(fasta_dict) as dict_handle:
        for name, length, md5 in [_tuples_from_line(x) for x in dict_handle if x.startswith("@SQ")]:
            out[md5] = name
    return out

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
        org_fa = org_build + ".fa.gz"
        if not os.path.exists(org_fa):
            genome = _download_ensembl_genome(org_build)
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
    dict_file = os.path.splitext(fasta_file.replace(".fa.gz", ".fa"))[0] + ".dict"
    if not os.path.exists(dict_file):
        subprocess.check_call("picard CreateSequenceDictionary R={fasta_file} "
                              "O={dict_file}".format(**locals()), shell=True)
    return dict_file


def _download_ensembl_genome(org_build):
    build = build_info[org_build]
    # reference files do not use the ensembl_release version so split it off
    fname = os.path.splitext(build.fbase)[0] + ".dna_sm.toplevel.fa.gz"
    dl_url = ("ftp://ftp.ensembl.org/pub/release-{release}/"
              "fasta/{taxname}/dna/{fname}").format(release=ensembl_release,
                                                    taxname=build.taxname,
                                                    fname=fname)
    out_file = os.path.basename(dl_url)
    if not os.path.exists(out_file):
        subprocess.check_call(["wget", "-c", dl_url])
    return out_file

# ## Main driver functions

def main(org_build, gtf_file=None):
    build_dir = os.path.abspath(os.path.join(os.curdir, org_build))
    work_dir = os.path.join(build_dir, "tmpcbl")
    out_dir = os.path.join(build_dir,
                           "rnaseq-%s" % datetime.datetime.now().strftime("%Y-%m-%d"))
    tophat_dir = os.path.join(out_dir, "tophat")
    safe_makedir(work_dir)
    with chdir(work_dir):
        if not gtf_file:
            build = build_info[org_build]
            gtf_file = prepare_tx_gff(build, org_build)
        else:
            work_gtf = os.path.join(work_dir, "ref-transcripts.gtf")
            if not os.path.exists(work_gtf):
                shutil.copy(gtf_file, work_gtf)
            gtf_file = work_gtf

        gtf_file = clean_gtf(gtf_file, org_build)
        db = _get_gtf_db(gtf_file)
        gtf_to_refflat(gtf_file)
        gtf_to_bed(gtf_file)
        prepare_dexseq(gtf_file)
        mask_gff = prepare_mask_gtf(gtf_file)
        rrna_gtf = prepare_rrna_gtf(gtf_file)
        if rrna_gtf:
            gtf_to_interval(rrna_gtf, org_build)
        prepare_tophat_index(gtf_file, org_build)
        cleanup(work_dir, out_dir, org_build)
        rnaseq_dir = os.path.join(build_dir, "rnaseq")
        if os.path.exists(rnaseq_dir):
            if os.path.islink(rnaseq_dir):
                os.unlink(rnaseq_dir)
            else:
                shutil.rmtree(rnaseq_dir)
        os.symlink(out_dir, rnaseq_dir)

    tar_dirs = [out_dir]
    tarball = create_tarball(tar_dirs, org_build)


def clean_gtf(gtf_file, org_build):
    """
    remove transcripts that have the following properties
    1) don't have a corresponding ID in the reference
    2) are bugged in the gencode release (Selenocysteine)
    3) are not associated with a gene (no gene_id field)
    """
    temp_gtf = tempfile.NamedTemporaryFile(suffix=".gtf").name
    fa_names = get_fasta_names(org_build)
    with open(gtf_file) as in_gtf, open(temp_gtf, "w") as out_gtf:
        for line in in_gtf:
            if line.startswith("#"):
                continue
            # these are bugged in the gencode release
            if "Selenocysteine" in line:
                continue
            if line.split()[0].strip() not in fa_names:
                continue
            if 'gene_id' not in line:
                continue
            out_gtf.write(line)
    shutil.move(temp_gtf, gtf_file)
    return gtf_file

def get_fasta_names(org_build):
    fa_dict = os.path.abspath(os.path.join(os.curdir, os.pardir, "seq", org_build + ".fa.fai"))
    with open(fa_dict) as in_handle:
        return [line.split("\t")[0] for line in in_handle]
    return seqs

def cleanup(work_dir, out_dir, org_build):
    try:
        os.remove(os.path.join(work_dir, org_build + ".dict"))
        os.remove(os.path.join(work_dir, org_build + ".fa"))
    except:
        pass
    shutil.move(work_dir, out_dir)

def create_tarball(tar_dirs, org_build):
    str_tar_dirs = " ".join(os.path.relpath(d) for d in tar_dirs)
    tarball = "{org}-{dir}.tar.xz".format(org=org_build, dir=os.path.basename(tar_dirs[0]))
    if not os.path.exists(tarball):
        subprocess.check_call("tar -cvpf - {out_dir} | xz -zc - > {tarball}".format(
            out_dir=str_tar_dirs, tarball=tarball), shell=True)
    return tarball

def upload_to_s3(tarball):
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

    cmd = "gtfToGenePred -allErrors -ignoreGroupsWithoutExons -genePredExt {gtf} {out_file}"
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

def gtf_to_bed(gtf):
    db = _get_gtf_db(gtf)
    out_file = os.path.splitext(gtf)[0] + ".bed"
    if file_exists(out_file):
        return out_file
    with open(out_file, "w") as out_handle:
        for feature in db.features_of_type('transcript'):
            chrom = feature.chrom
            start = feature.start
            end = feature.end
            attributes = feature.attributes.keys()
            strand = feature.strand
            name = (feature['gene_name'][0] if 'gene_name' in attributes else
                    feature['gene_id'][0])
            line = "\t".join(map(str, [chrom, start, end, name, ".", strand]))
            out_handle.write(line + "\n")
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
    biotype_lookup = _biotype_lookup_fn(gtf)
    db = _get_gtf_db(gtf)
    with open(out_file, "w") as out_handle:
        for g in db.all_features():
            biotype = biotype_lookup(g)
            if (biotype in mask_biotype) or (g.chrom in mask_chrom):
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
    biotype_lookup = _biotype_lookup_fn(gtf)
    # if we can't find a biotype column, skip this
    if not biotype_lookup:
        return None
    with open(out_file, "w") as out_handle:
        for feature in db.all_features():
            biotype = biotype_lookup(feature)
            if biotype in mask_biotype:
                out_handle.write(str(feature) + "\n")
    return out_file

def _biotype_lookup_fn(gtf):
    """
    return a function that will look up the biotype of a feature
    this checks for either gene_biotype or biotype being set or for the source
    column to have biotype information
    """
    db = _get_gtf_db(gtf)
    sources = set([feature.source for feature in db.all_features()])
    gene_biotypes = set([feature.attributes.get("gene_biotype", [None])[0]
                         for feature in db.all_features()])
    biotypes = set([feature.attributes.get("biotype", [None])[0]
                    for feature in db.all_features()])
    if "protein_coding" in sources:
        return lambda feature: feature.source
    elif "protein_coding" in biotypes:
        return lambda feature: feature.attributes.get("biotype", [None])[0]
    elif "protein_coding" in gene_biotypes:
        return lambda feature: feature.attributes.get("gene_biotype", [None])[0]
    else:
        return None

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
    ensembl_gff = _download_ensembl_gff(build, org_name)
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
    wrote_missing = set([])
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle, \
             open(base_gff) as in_handle:
            for line in in_handle:
                parts = line.split("\t")
                ucsc_name = name_map.get(parts[0], None)
                if ucsc_name:
                    out_handle.write("\t".join([ucsc_name] + parts[1:]))
                elif parts[0] not in wrote_missing and not line.startswith("#"):
                    print "Missing", parts[0]
                    wrote_missing.add(parts[0])
    return out_file

def _download_ensembl_gff(build, org_name):
    """Given build details, download and extract the relevant ensembl GFF.
    """
    fname = build.fbase + ".gtf.gz"
    dl_url = "/".join([base_ftp, build.taxname, fname]).format(
        release=supported_oldbuilds.get(org_name, ensembl_release))
    out_file = os.path.splitext(os.path.basename(dl_url))[0]
    if not os.path.exists(out_file):
        subprocess.check_call(["wget", dl_url])
        subprocess.check_call(["gunzip", os.path.basename(dl_url)])
    return out_file

def guess_infer_extent(gtf_file):
    """
    guess if we need to use the gene extent option when making a gffutils
    database by making a tiny database of 1000 lines from the original
    GTF and looking for all of the features
    """
    _, ext = os.path.splitext(gtf_file)
    tmp_out = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False).name
    with open(tmp_out, "w") as out_handle:
        count = 0
        in_handle = open(gtf_file) if ext != ".gz" else gzip.open(gtf_file)
        for line in in_handle:
            if count > 1000:
                break
            out_handle.write(line)
            count += 1
        in_handle.close()
    db = gffutils.create_db(tmp_out, dbfn=":memory:", infer_gene_extent=False)
    os.remove(tmp_out)
    features = [x for x in db.featuretypes()]
    if "gene" in features and "transcript" in features:
        return False
    else:
        return True

def _get_gtf_db(gtf):
    db_file = gtf + ".db"
    if not file_exists(db_file):
        print "Creating gffutils database for %s." % (gtf)
        infer_extent = guess_infer_extent(gtf)
        if infer_extent:
            print ("'transcript' and 'gene' entries not found, so inferring"
                   "their extent. This can be very slow.")
        gffutils.create_db(gtf, dbfn=db_file, infer_gene_extent=infer_extent)
    return gffutils.FeatureDB(db_file)

def _dexseq_preparation_path():
    PREP_FILE = "python_scripts/dexseq_prepare_annotation.py"
    try:
        cmd = "Rscript -e 'find.package(\"DEXSeq\")'"
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        return None
    for line in output.split("\n"):
        if line.startswith("["):
            dirname = line.split("[1]")[1].replace("\"", "").strip()
            path = os.path.join(dirname, PREP_FILE)
            if os.path.exists(path):
                return path
    return None

def prepare_dexseq(gtf):
    out_file = os.path.splitext(gtf)[0] + ".dexseq.gff3"
    if file_exists(out_file):
        return out_file

    dexseq_path = _dexseq_preparation_path()
    if not dexseq_path:
        return None
    executable = sys.executable
    cmd = "{executable} {dexseq_path} {gtf} {out_file}"
    subprocess.check_call(cmd.format(**locals()), shell=True)
    return out_file

if __name__ == "__main__":
    parser = ArgumentParser(description="Prepare the transcriptome files for an "
                            "organism.")
    parser.add_argument("--gtf",
                        help="Optional GTF file (instead of downloading from Ensembl.",
                        default=None),
    parser.add_argument("org_build", help="Build of organism to run.")
    args = parser.parse_args()
    main(args.org_build, args.gtf)
