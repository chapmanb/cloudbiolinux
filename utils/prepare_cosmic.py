#!/usr/bin/env python
"""Prepare combined VCF files of COSMIC resource for cancer variant calling.

http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/
ftp://ngs.sanger.ac.uk/production/cosmic/
http://gatkforums.broadinstitute.org/discussion/2226/cosmic-and-dbsnp-files-for-mutect
"""
import os
import subprocess
import sys

FTP_DIR = "ftp://ngs.sanger.ac.uk/production/cosmic/"
VERSION = "v68"
BCBIO_NEXTGEN_BASE = "/usr/local"

def main():
    work_dir = "tmp-cosmic-GRCh37"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    ref_file = BCBIO_NEXTGEN_BASE + "/share/bcbio_nextgen/genomes/Hsapiens/GRCh37/seq/GRCh37.fa"
    fnames = [reorder_reference(x, ref_file) for x in get_cosmic_files()]
    grc_cosmic = combine_cosmic(fnames, ref_file)
    hg_cosmic = map_coords_to_ucsc(grc_cosmic, ref_file)
    for ready_file in [bgzip_vcf(x) for x in [grc_cosmic, hg_cosmic]]:
        upload_to_s3(ready_file)
        upload_to_s3(ready_file.replace(".gz", ".idx"))
        upload_to_s3(ready_file + ".tbi")

def upload_to_s3(fname):
    upload_script = os.path.join(os.path.dirname(__file__), "s3_multipart_upload.py")
    subprocess.check_call([sys.executable, upload_script, fname, "biodata",
                           "variants/%s" % os.path.basename(fname), "--public"])

def map_coords_to_ucsc(grc_cosmic, ref_file):
    hg19_ref_file = ref_file.replace("GRCh37", "hg19")
    out_file = grc_cosmic.replace("GRCh37.vcf", "hg19.vcf")
    if not os.path.exists(out_file):
        tmp_file = "%s-raw%s" % os.path.splitext(out_file)
        with open(tmp_file, "w") as out_handle:
            # header
            with open(grc_cosmic) as in_handle:
                for line in in_handle:
                    if line.startswith("#") and not line.startswith("##contig"):
                        out_handle.write(line)
            # chromsome M
            with open(grc_cosmic) as in_handle:
                for line in in_handle:
                    if line.startswith("MT"):
                        line = _rename_to_ucsc(line)
                        out_handle.write(line)
            # rest
            with open(grc_cosmic) as in_handle:
                for line in in_handle:
                    if not line.startswith(("MT", "#")):
                        line = _rename_to_ucsc(line)
                        out_handle.write(line)
        # Create clean VCF and index for upload
        subprocess.check_call(["gatk-framework", "-R", hg19_ref_file, "-T", "SelectVariants",
                               "--variant", tmp_file, "--out", out_file])
    return out_file

def _rename_to_ucsc(line):
    chrom, rest = line.split("\t", 1)
    if chrom == "MT":
        new_chrom = "chrM"
    else:
        new_chrom = "chr%s" % chrom
    return "%s\t%s" % (new_chrom, rest)

def combine_cosmic(fnames, ref_file):
    out_file = "cosmic-%s-GRCh37.vcf" % VERSION
    if not os.path.exists(out_file):
        cmd = ["gatk-framework", "-T", "CombineVariants", "-R", ref_file, "--out", out_file,
               "--suppressCommandLineHeader", "--setKey", "null"]
        for v in fnames:
            cmd += ["--variant", v]
        subprocess.check_call(cmd)
    return out_file

def reorder_reference(fname, ref_file):
    """Move mitochondrial calls to end to match GATK reference ordering.
    """
    out_file = "%s-prep%s" % os.path.splitext(fname)
    bcbiov_jar = BCBIO_NEXTGEN_BASE + "/share/java/bcbio_variation/bcbio.variation-0.1.6-SNAPSHOT-standalone.jar"
    if not os.path.exists(out_file):
        cmd = ["java", "-jar", bcbiov_jar, "variant-utils", "sort-vcf", fname, ref_file, ""]
        subprocess.check_call(cmd)
    return out_file

def get_cosmic_files():
    fnames = []
    for ctype in ["CodingMuts", "NonCodingVariants"]:
        to_get = "%sCosmic%s_%s.vcf.gz" % (FTP_DIR, ctype, VERSION)
        fname_gz = os.path.basename(to_get)
        fname = os.path.splitext(fname_gz)[0]
        if not os.path.exists(fname):
            if not os.path.exists(fname_gz):
                subprocess.check_call(["wget", to_get])
            subprocess.check_call(["gunzip", fname_gz])
        fnames.append(fname)
    return fnames

def bgzip_vcf(in_file):
    out_file = in_file + ".gz"
    if not os.path.exists(out_file):
        subprocess.check_call("bgzip -c %s > %s" % (in_file, out_file), shell=True)
    tabix_file = out_file + ".tbi"
    if not os.path.exists(tabix_file):
        subprocess.check_call(["tabix", "-p", "vcf", out_file])
    return out_file

if __name__ == "__main__":
    main()
