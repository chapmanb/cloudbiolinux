"""Prepare sorted and consolidated dbSNP resources for mouse mm10/GRCh38 in VCF format.
"""
import datetime
import ftplib
import gzip
import os
import subprocess
from argparse import ArgumentParser
import re
import shutil

FTP = "ftp.ncbi.nih.gov"

REMOTES = {"mm10": "snp/organisms/mouse_10090/VCF",
           "canFam3": "snp/organisms/dog_9615/VCF/"}

def main(org):
    work_dir = "tmp-dbsnp-%s" % org
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    conn = ftplib.FTP(FTP, "anonymous", "me@example.com")
    conn.cwd(REMOTES[org])

    os.chdir(work_dir)
    files = []
    def add_files(x):
        if x.endswith("vcf.gz"):
            files.append(get_file(x, REMOTES[org], conn))
    conn.retrlines("NLST", add_files)
    out_file = "%s-dbSNP-%s.vcf" % (org, datetime.datetime.now().strftime("%Y-%m-%d"))
    with open(out_file, "w") as out_handle:
        for i, f in enumerate(karyotype_sort(files)):
            with gzip.open(f) as in_handle:
                for line in in_handle:
                    if line.startswith("#"):
                        if i == 0:
                            out_handle.write(line)
                    else:
                        out_handle.write("\t".join(fix_info(fix_chrom(line.rstrip().split("\t")))) + "\n")
    subprocess.check_call(["bgzip", out_file])
    shutil.move(out_file + ".gz", os.path.join(os.pardir, out_file + ".gz"))
    os.chdir(os.pardir)
    subprocess.check_call(["tabix", "-p", "vcf", out_file + ".gz"])
    shutil.rmtree(work_dir)

multi_whitespace = re.compile(r"\s+")

def fix_info(parts):
    """Fix the INFO file to remove whitespace.
    """
    parts[7] = multi_whitespace.sub("_", parts[7])
    return parts

def fix_chrom(parts):
    MAX_CHROMOSOMES = 50
    if parts[0] in [str(x) for x in range(1, MAX_CHROMOSOMES)] + ["X", "Y"]:
        new_chrom = "chr%s" % parts[0]
    elif parts[0] == "MT":
        new_chrom = "chrM"
    else:
        raise NotImplementedError(parts)
    parts[0] = new_chrom
    return parts

def get_file(x, ftp_dir, conn):
    if not os.path.exists(x):
        print "Retrieving %s" % x
        with open(x, "wb") as out_handle:
            conn = ftplib.FTP(FTP, "anonymous", "me@example.com")
            conn.cwd(ftp_dir)
            conn.retrbinary("RETR %s" % x, out_handle.write)
    return x

def karyotype_sort(xs):
    """Sort in karyotypic order to work with GATK's defaults.
    """
    def karyotype_keyfn(x):
        for suffix in [".gz"]:
            if x.endswith(suffix):
                x = x[:-len(suffix)]
        base = os.path.splitext(os.path.basename(x))[0]
        for prefix in ["chr", "vcf_chr_"]:
            if base.startswith(prefix):
                base = base[len(prefix):]
        parts = base.split("_")
        try:
            parts[0] =  int(parts[0])
        except ValueError:
            pass
        # unplaced at the very end
        if isinstance(parts[0], basestring) and parts[0].startswith(("Un", "Alt", "Multi", "NotOn")):
            parts.insert(0, "z")
        # mitochondrial special case -- after X/Y
        elif parts[0] in ["M", "MT"]:
            parts.insert(0, "x")
        # sort random and extra chromosomes after M
        elif len(parts) > 1:
            parts.insert(0, "y")
        return parts
    return sorted(xs, key=karyotype_keyfn)

if __name__ == "__main__":
    parser = ArgumentParser(description="Prepare a dbSNP file from NCBI.")
    parser.add_argument("org_build", choices=REMOTES.keys(),
                        help="genome build")
    args = parser.parse_args()
    main(args.org_build)
