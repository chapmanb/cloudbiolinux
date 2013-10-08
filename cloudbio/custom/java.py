"""Install instructions for non-packaged java programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_installed, _make_tmp_dir)
from cloudbio.custom import shared

@_if_not_installed("lein")
def install_leiningen(env):
    """Clojure tool for project configuration and automation.
    http://github.com/technomancy/leiningen
    """
    bin_dir = os.path.join(env.system_install, "bin")
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            shared._remote_fetch(env,"https://raw.github.com/technomancy/leiningen/stable/bin/lein") 
            env.safe_run("chmod a+rwx lein")
            env.safe_sudo("mv lein %s" % bin_dir)
            env.safe_run("%s/lein" % bin_dir)
