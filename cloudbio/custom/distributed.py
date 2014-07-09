"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import (_if_not_python_lib, _pip_cmd, _is_anaconda)

@_if_not_python_lib("pydoop")
def install_pydoop(env):
    """pydoop; provides Hadoop access for Python.
    http://pydoop.sourceforge.net/docs/
    """
    java_home = env.java_home if "java_home" in env else os.environ["JAVA_HOME"]
    export_str = "export JAVA_HOME=%s" % (java_home)
    cmd = env.safe_run if _is_anaconda(env) else env.safe_sudo
    cmd("%s && %s install pydoop" % (export_str, _pip_cmd(env)))

@_if_not_python_lib("bl.mr.seq.seqal")
def install_seal(env):
    """Install seal: process high-throughput sequencing with Hadoop.

    http://biodoop-seal.sf.net/
    """
    install_pydoop(env)

    java_home = env.java_home if "java_home" in env else os.environ["JAVA_HOME"]
    export_str = "export JAVA_HOME=%s" % (java_home)
    cmd = env.safe_run if _is_anaconda(env) else env.safe_sudo
    cmd("%s && %s install --pre seal" % (export_str, _pip_cmd(env)))
