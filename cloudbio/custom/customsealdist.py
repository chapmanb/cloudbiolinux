"""Install instructions for Python components of seal flavour.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir, _if_not_python_lib, _fetch_and_unpack

@_if_not_python_lib("pydoop")
def install_pydoop(env):
    """Install pydoop; provides Hadoop access for Python.
    """
    java_home = env.java_home if env.has_key("java_home") else os.environ["JAVA_HOME"]
    export_str = "export JAVA_HOME=%s" % (java_home)
    env.safe_sudo("%s && pip-python install pydoop" % export_str)

@_if_not_python_lib("bl.mr.seq.seqal")
def install_seal(env):
    """Install seal: process high-throughput sequencing with Hadoop.

    http://biodoop-seal.sf.net/
    """
    java_home = env.java_home if env.has_key("java_home") else os.environ["JAVA_HOME"]
    hadoop_home = env.hadoop_home if env.has_key("hadoop_home") else os.environ["HADOOP_HOME"]

    try:
      import pydoop
    except ImportError:
      env.logger.fatal("Cannot load pydoop module!")
      raise ImportError("Error installing Seal: missing Pydoop.  Please install Pydoop first")

    export_str = "export JAVA_HOME='%s'" % java_home
    env.safe_sudo("%s && pip-python install seal" % export_str)
