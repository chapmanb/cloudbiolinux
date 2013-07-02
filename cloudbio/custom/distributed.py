"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_installed, _if_not_python_lib,
                    _pip_cmd, _get_install, _configure_make)

@_if_not_installed("parallel -h")
def install_gnu_parallel(env):
    """GNU parallel: build and execute command lines from standard input in parallel.
    https://savannah.gnu.org/projects/parallel/
    """
    version = "20130122"
    url = "ftp://ftp.gnu.org/gnu/parallel/parallel-%s.tar.bz2" % version
    _get_install(url, env, _configure_make)

@_if_not_python_lib("pydoop")
def install_pydoop(env):
    """pydoop; provides Hadoop access for Python.
    http://pydoop.sourceforge.net/docs/
    """
    java_home = env.java_home if env.has_key("java_home") else os.environ["JAVA_HOME"]
    export_str = "export JAVA_HOME=%s" % (java_home)
    env.safe_sudo("%s && %s install pydoop" % (export_str, _pip_cmd(env)))

@_if_not_python_lib("bl.mr.seq.seqal")
def install_seal(env):
    """Install seal: process high-throughput sequencing with Hadoop.

    http://biodoop-seal.sf.net/
    """
    install_pydoop(env)

    java_home = env.java_home if env.has_key("java_home") else os.environ["JAVA_HOME"]
    export_str = "export JAVA_HOME=%s" % (java_home)
    env.safe_sudo("%s && %s install seal" % (export_str, _pip_cmd(env)))
