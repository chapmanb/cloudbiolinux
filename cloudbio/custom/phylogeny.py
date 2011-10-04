"""Install instructions for non-packaged phyologeny programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed

@_if_not_installed("tracer")
def install_tracer(env):
    """Tracer tool for phylogeny
    """
    if not exists("/usr/local/bin/tracer"):
        run("wget http://bio4.dnsalias.net/download/biolinux/packages/Tracer_v1.5.tgz")
        run("tar xvzf Tracer_v1.5.tgz")
        run("chmod a+x Tracer_v1.5/bin/tracer")
        env.safe_sudo("mkdir -p /usr/local/bioinf")
        env.safe_sudo("rm -rvf /usr/local/bioinf/tracer")
        env.safe_sudo("mv -f Tracer_v1.5 /usr/local/bioinf/tracer")
        env.safe_sudo("ln -sf /usr/local/bioinf/tracer/bin/tracer /usr/local/bin/tracer")

