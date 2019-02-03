"""Infrastructure for RNA-seq supporting files.
"""
import os
import subprocess

from cloudbio.custom import shared

def finalize(genomes, data_filedir):
    """Provide symlinks back to reference genomes so tophat avoids generating FASTA genomes.
    """
    genome_dir = os.path.join(data_filedir, "genomes")
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname)
        for aligner in ["bowtie", "bowtie2"]:
            aligner_dir = os.path.join(org_dir, gid, aligner)
            if os.path.exists(aligner_dir):
                with shared.chdir(aligner_dir):
                    for ext in ["", ".fai"]:
                        orig_seq = os.path.join(os.pardir, "seq", "%s.fa%s" % (gid, ext))
                        if os.path.exists(orig_seq) and not os.path.exists(os.path.basename(orig_seq)):
                            subprocess.check_call("ln -sf %s" % orig_seq, shell=True)

def cleanup(genomes, data_filedir):
    """Cleanup for GGD recipe installation, removing old rnaseq symlinks.
    """
    folder_name = "rnaseq"
    genome_dir = os.path.join(data_filedir, "genomes")
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname)
        target_dir = os.path.join(org_dir, gid, folder_name)
        if os.path.lexists(target_dir) and os.path.islink(target_dir):
            os.remove(target_dir)
