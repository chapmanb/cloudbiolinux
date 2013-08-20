"""Install instructions for non-packaged phyologeny programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from cloudbio.custom.shared import _if_not_installed, _make_tmp_dir

def install_tracer(env):
    """A program for analysing results from Bayesian MCMC programs such as BEAST & MrBayes.
    http://tree.bio.ed.ac.uk/software/tracer/
    """
    version = "1.5"
    install_dir = os.path.join(env.system_install, "bioinf")
    final_exe = os.path.join(env.system_install, "bin", "tracer")
    if env.safe_exists(final_exe):
        return
    if not exists(final_exe):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget -O Tracer_v{0}.tgz 'http://tree.bio.ed.ac.uk/download.php?id=80&num=3'".format(
                        version))
                run("tar xvzf Tracer_v{0}.tgz".format(version))
                run("chmod a+x Tracer_v{0}/bin/tracer".format(version))
                env.safe_sudo("mkdir -p %s" % install_dir)
                env.safe_sudo("rm -rvf %s/tracer" % install_dir)
                env.safe_sudo("mv -f Tracer_v%s %s/tracer" % (version, install_dir))
                env.safe_sudo("ln -sf %s/tracer/bin/tracer %s" % (install_dir, final_exe))

@_if_not_installed("beast -help")
def install_beast(env):
    """BEAST: Bayesian MCMC analysis of molecular sequences.
    http://beast.bio.ed.ac.uk
    """
    version = "1.7.4"
    install_dir = os.path.join(env.system_install, "bioinf")
    final_exe = os.path.join(env.system_install, "bin", "beast")
    if not exists(final_exe):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget http://beast-mcmc.googlecode.com/files/BEASTv%s.tgz" % version)
                run("tar xvzf BEASTv%s.tgz" % version)
                env.safe_sudo("mkdir -p %s" % install_dir)
                env.safe_sudo("rm -rvf %s/beast" % install_dir)
                env.safe_sudo("mv -f BEASTv%s %s/beast" % (version, install_dir))
                for l in ["beast","beauti","loganalyser","logcombiner","treeannotator","treestat"]:
                    env.safe_sudo("ln -sf %s/beast/bin/%s %s/bin/%s" % (install_dir, l,
                                                                        env.system_install, l))

