"""Install instructions for python libraries not ready for easy_install.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_python_lib, _get_install

def _python_make(env):
    run("python%s setup.py build" % env.python_version_ext)
    sudo("python%s setup.py install --skip-build" % env.python_version_ext)
    sudo("rm -rf dist")
    sudo("rm -rf build")
    sudo("rm -rf lib/*.egg-info")

@_if_not_python_lib("bx")
def install_bx_python(env):
    url = "hg clone http://bitbucket.org/james_taylor/bx-python"
    _get_install(url, env, _python_make)

@_if_not_python_lib("matplotlib")
def install_matplotlib(env):
    version = "1.0.1"
    url = "http://downloads.sourceforge.net/project/matplotlib/matplotlib/" \
          "matplotlib-%s/matplotlib-%s.tar.gz" % (version, version)
    _get_install(url, env, _python_make)

@_if_not_python_lib("rpy")
def install_rpy(env):
    version = "1.0.3"
    ext = "a"
    url = "http://downloads.sourceforge.net/project/rpy/rpy/" \
          "%s/rpy-%s%s.zip" % (version, version, ext)
    _get_install(url, env, _python_make)
