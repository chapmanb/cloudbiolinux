"""Install instructions for non-packaged java programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_installed, _make_tmp_dir)

@_if_not_installed("cljr")
def install_cljr(env):
    """Clojure package manager, cljr.

    http://github.com/liebke/cljr
    """
    run("wget http://incanter.org/downloads/cljr-installer.jar")
    run("java -jar cljr-installer.jar")
    env.safe_sudo("ln -s .cljr/bin/cljr /usr/bin")
    run("rm cljr-installer.jar")

@_if_not_installed("lein")
def install_leiningen(env):
    """Clojure tool for project configuration and automation.
    http://github.com/technomancy/leiningen
    """
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            run("wget --no-check-certificate https://raw.github.com/technomancy/leiningen/stable/bin/lein")
            run("chmod a+rwx lein")
            env.safe_sudo("mv lein %s" % os.path.join(env.system_install, "bin"))
            run("lein")
