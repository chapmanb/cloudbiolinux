"""Install instructions for non-packaged phyologeny programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from cloudbio.custom.shared import _if_not_installed, _make_tmp_dir

@_if_not_installed("tracer")
def install_tracer(env):
    """Tracer tool for phylogeny
    """
    install_dir = os.path.join(env.system_install, "bioinf")
    final_exe = os.path.join(env.system_install, "bin", "tracer")
    if not exists(final_exe):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget http://bio4.dnsalias.net/download/biolinux/packages/Tracer_v1.5.tgz")
                run("tar xvzf Tracer_v1.5.tgz")
                run("chmod a+x Tracer_v1.5/bin/tracer")
                env.safe_sudo("mkdir -p %s" % install_dir)
                env.safe_sudo("rm -rvf %s/tracer" % install_dir)
                env.safe_sudo("mv -f Tracer_v1.5 %s/tracer" % install_dir)
                env.safe_sudo("ln -sf %s/tracer/bin/tracer %s" % (install_dir, final_exe))

@_if_not_installed("beast")
def install_beast(env):
    """BEAST for phylogeny
    """
    install_dir = os.path.join(env.system_install, "bioinf")
    final_exe = os.path.join(env.system_install, "bin", "beast")
    if not exists(final_exe):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                run("wget http://beast-mcmc.googlecode.com/files/BEASTv1.6.2.tgz")
                run("tar xvzf BEASTv1.6.2.tgz")
                env.safe_sudo("mkdir -p %s" % install_dir)
                env.safe_sudo("rm -rvf %s/beast" % install_dir)
                env.safe_sudo("mv -f BEASTv1.6.2 %s/beast" % install_dir)
                for l in ["beast","beauti","loganalyser","logcombiner","treeannotator","treestat"]:
                    env.safe_sudo("ln -sf %s/beast/bin/%s %s/bin/%s" % (install_dir, l,
                                                                        env.system_install, l))

