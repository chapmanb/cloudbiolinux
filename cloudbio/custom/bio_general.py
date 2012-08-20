"""Custom installs for biological packages.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_installed, _get_install, _configure_make, _java_install)

@_if_not_installed("embossversion")
def install_emboss(env):
    """EMBOSS: A high-quality package of free, Open Source software for molecular biology.
    http://emboss.sourceforge.net/
    Emboss target for platforms without packages (CentOS -- rpm systems).
    """
    default_version = "6.3.1"
    version = env.get("tool_version", default_version)
    url = "ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-%s.tar.gz" % version
    _get_install(url, env, _configure_make)

@_if_not_installed("PGDSpider2.sh")
def install_pgdspider(env):
    """PGDSpider format conversion for population genetics programs.
    http://www.cmpg.unibe.ch/software/PGDSpider/
    """
    version = "2.0.1.2"
    url = "http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_{v}.zip".format(
        v=version)
    def _install_fn(env, install_dir):
        env.safe_sudo("mv *.jar %s" % install_dir)
        bin_dir = os.path.join(env.system_install, "bin")
        exe_file = "PGDSpider2.sh"
        jar = "PGDSpider2.jar"
        sed(exe_file, jar, "{dir}/{jar}".format(dir=install_dir, jar=jar))
        run("chmod a+x {0}".format(exe_file))
        env.safe_sudo("mv {exe} {bin}".format(exe=exe_file, bin=bin_dir))
    _java_install("PGDSpider", version, url, env, install_fn=_install_fn)

def install_bio4j(env):
    """Bio4j graph based database built on Neo4j with UniProt, GO, RefSeq and more.
    http://www.bio4j.com/
    """
    version = "0.7"
    url = "https://s3-eu-west-1.amazonaws.com/bio4j-public/releases/" \
          "{v}/bio4j-{v}.zip".format(v=version)
    def _install_fn(env, install_dir):
        targets = ["conf", "doc", "jars", "lib", "README"]
        for x in targets:
            env.safe_sudo("mv {0} {1}".format(x, install_dir))
    _java_install("bio4j", version, url, env, install_fn=_install_fn)
