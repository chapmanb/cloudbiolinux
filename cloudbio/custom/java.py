"""Install instructions for non-packaged java programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed

@_if_not_installed("cljr")
def install_cljr(env):
    """Install the clojure package manager cljr

    http://github.com/liebke/cljr
    """
    run("wget http://incanter.org/downloads/cljr-installer.jar")
    run("java -jar cljr-installer.jar")
    env.safe_sudo("ln -s .cljr/bin/cljr /usr/bin")
    run("rm cljr-installer.jar")

@_if_not_installed("lein")
def install_leinengin(env):
    """Standard clojure build tool: http://github.com/technomancy/leiningen
    """
    run("wget --no-check-certificate https://github.com/technomancy/leiningen/raw/stable/bin/lein")
    run("chmod a+rwx lein")
    env.safe_sudo("mv lein %s" % os.path.join(env.system_install, "bin"))
    run("lein self-install")

@_if_not_installed("tracer")
def install_tracer(env):
    """Tracer tool for phylogeny
    """
    run("wget http://bio4.dnsalias.net/download/biolinux/packages/Tracer_v1.5.tgz")
    run("tar xvzf Tracer_v1.5.tgz")
    run("chmod a+x Tracer_v1.5/bin/tracer")
    env.safe_sudo("ln -s `pwd`/Tracer_v1.5/bin/tracer /usr/local/bin")

