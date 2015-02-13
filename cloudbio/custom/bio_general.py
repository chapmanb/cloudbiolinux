"""Custom installs for biological packages.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from cloudbio.custom import shared
from shared import (_if_not_installed, _get_install, _configure_make, _java_install,
                    _make_tmp_dir)

def install_anaconda(env):
    """Pre-packaged Anaconda Python installed from Continuum.
    http://docs.continuum.io/anaconda/index.html
    """
    version = "2.0.0"
    outdir = os.path.join(env.system_install, "anaconda")
    if env.distribution in ["ubuntu", "centos", "scientificlinux", "debian", "arch", "suse"]:
        platform = "Linux"
    elif env.distribution in ["macosx"]:
        platform = "MacOSX"
    else:
        raise ValueError("Unexpected distribution: %s" % env.distribution)
    url = "http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/" \
          "Anaconda-%s-%s-x86_64.sh" % (version, platform)
    if not env.safe_exists(outdir):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                installer = shared._remote_fetch(env, url)
                env.safe_sed(os.path.basename(url), "more <<EOF", "cat  <<EOF")
                env.safe_sudo("echo -e '\nyes\n%s\nyes\n' | bash %s" % (outdir, installer))
                env.safe_sudo("chown -R %s %s" % (env.user, outdir))
                comment_line = "# added by Ananconda %s installer" % version
                if not env.safe_contains(env.shell_config, comment_line):
                    env.safe_append(env.shell_config, comment_line)
                    env.safe_append(env.shell_config, "export PATH=%s/bin:$PATH" % outdir)
                # remove curl library with broken certificates
                env.safe_run("%s/bin/conda remove --yes curl" % outdir)
                env.safe_run("%s/bin/conda install --yes pip" % outdir)

@_if_not_installed("embossversion")
def install_emboss(env):
    """EMBOSS: A high-quality package of free, Open Source software for molecular biology.
    http://emboss.sourceforge.net/
    Emboss target for platforms without packages (CentOS -- rpm systems).
    """
    default_version = "6.6.0"
    version = env.get("tool_version", default_version)
    url = "https://science-annex.org/pub/emboss/EMBOSS-%s.tar.gz" % version
    #url = "ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-%s.tar.gz" % version
    _get_install(url, env, _configure_make)

def install_pgdspider(env):
    """PGDSpider format conversion for population genetics programs.
    http://www.cmpg.unibe.ch/software/PGDSpider/
    """
    if os.path.exists(os.path.join(shared._get_bin_dir(env), "PGDSpider2.sh")):
        return
    version = "2.0.2.0"
    url = "http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_{v}.zip".format(
        v=version)
    def _install_fn(env, install_dir):
        env.safe_sudo("mv *.jar %s" % install_dir)
        bin_dir = shared._get_bin_dir(env)
        exe_file = "PGDSpider2.sh"
        jar = "PGDSpider2.jar"
        env.safe_sed(exe_file, jar, "{dir}/{jar}".format(dir=install_dir, jar=jar))
        env.safe_run("chmod a+x {0}".format(exe_file))
        env.safe_sudo("mv {exe} {bin}".format(exe=exe_file, bin=bin_dir))
    _java_install("PGDSpider", version, url, env, install_fn=_install_fn)

def install_bio4j(env):
    """Bio4j graph based database built on Neo4j with UniProt, GO, RefSeq and more.
    http://www.bio4j.com/
    """
    version = "0.8"
    url = "https://s3-eu-west-1.amazonaws.com/bio4j-public/releases/" \
          "{v}/bio4j-{v}.zip".format(v=version)
    def _install_fn(env, install_dir):
        targets = ["conf", "doc", "jars", "lib", "README"]
        for x in targets:
            env.safe_sudo("mv {0} {1}".format(x, install_dir))
    _java_install("bio4j", version, url, env, install_fn=_install_fn)
