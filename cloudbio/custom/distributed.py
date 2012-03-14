"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir, _if_not_python_lib, _fetch_and_unpack


def _pip_cmd(env):
    if env.has_key("pip_cmd") and env.pip_cmd:
        return env.pip_cmd
    else:
        return "pip"

@_if_not_python_lib("pydoop")
def install_pydoop(env):
    """Install pydoop; provides Hadoop access for Python.
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
