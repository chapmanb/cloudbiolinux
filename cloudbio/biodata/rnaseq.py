"""Infrastructure for RNA-seq supporting files.
"""
import os

from fabric.api import cd

def finalize(genomes, env):
    """Provide symlinks back to reference genomes so tophat avoids generating FASTA genomes.
    """
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname)
        for aligner in ["bowtie", "bowtie2"]:
            aligner_dir = os.path.join(org_dir, gid, aligner)
            if env.safe_exists(aligner_dir):
                with cd(aligner_dir):
                    for ext in ["", ".fai"]:
                        orig_seq = os.path.join(os.pardir, "seq", "%s.fa%s" % (gid, ext))
                        if env.safe_exists(orig_seq) and not env.safe_exists(os.path.basename(orig_seq)):
                            env.safe_run("ln -sf %s" % orig_seq)

def cleanup(genomes, env):
    """Cleanup for GGD recipe installation, removing old rnaseq symlinks.
    """
    folder_name = "rnaseq"
    genome_dir = os.path.join(env.data_files, "genomes")
    for (orgname, gid, manager) in genomes:
        org_dir = os.path.join(genome_dir, orgname)
        target_dir = os.path.join(org_dir, gid, folder_name)
        if env.hosts == ["localhost"]:
            if os.path.lexists(target_dir) and os.path.islink(target_dir):
                os.remove(target_dir)
