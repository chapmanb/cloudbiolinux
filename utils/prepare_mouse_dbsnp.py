"""Prepare sorted and consolidated dbSNP resources for mouse mm10/GRCh38 in VCF format.
"""
import datetime
import ftplib
import gzip
import os
import subprocess

REMOTES = {"ftp": "ftp.ncbi.nih.gov",
           "ftp_dir": "snp/organisms/mouse_10090/VCF/",
           "org": "mm10"}

def main():
    work_dir = "tmp-dbsnp-%s" % REMOTES["org"]
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    conn = ftplib.FTP(REMOTES["ftp"], "anonymous", "me@example.com")
    conn.cwd(REMOTES["ftp_dir"])

    os.chdir(work_dir)
    files = []
    def add_files(x):
        if x.endswith("vcf.gz"):
            files.append(get_file(x, REMOTES, conn))
    conn.retrlines("NLST", add_files)
    out_file = "%s-dbSNP-%s.vcf" % (REMOTES["org"], datetime.datetime.now().strftime("%Y-%m-%d"))
    with open(out_file, "w") as out_handle:
        for i, f in enumerate(karyotype_sort(files)):
            with gzip.open(f) as in_handle:
                test = 0
                for line in in_handle:
                    if line.startswith("#"):
                        if i == 0:
                            out_handle.write(line)
                    else:
                        out_handle.write(fix_chrom(line))
    subprocess.check_call(["bgzip", out_file])

def fix_chrom(line):
    parts = line.split("\t")
    if parts[0] in [str(x) for x in range(1, 23)] + ["X", "Y"]:
        new_chrom = "chr%s" % parts[0]
    elif parts[0] == "MT":
        new_chrom = "chrM"
    else:
        raise NotImplementedError(line)
    parts[0] = new_chrom
    return "\t".join(parts)

def get_file(x, remotes, conn):
    if not os.path.exists(x):
        print "Retrieving %s" % x
        with open(x, "wb") as out_handle:
            conn = ftplib.FTP(remotes["ftp"], "anonymous", "me@example.com")
            conn.cwd(remotes["ftp_dir"])
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
    main()
