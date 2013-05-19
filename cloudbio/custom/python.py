"""Install instructions for python libraries not ready for easy_install.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_python_lib, _get_install, _python_make, _pip_cmd

@_if_not_python_lib("bx")
def install_bx_python(env):
    """Tools for manipulating biological data, particularly multiple sequence alignments
    https://bitbucket.org/james_taylor/bx-python/wiki/Home
    """
    version = "bitbucket"
    url = "https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2"
    env.safe_sudo("%s install --upgrade %s" % (_pip_cmd(env), url))

@_if_not_python_lib("rpy")
def install_rpy(env):
    """RPy is a very simple, yet robust, Python interface to the R Programming Language.
    http://rpy.sourceforge.net/
    """
    version = "1.0.3"
    ext = "a"
    url = "http://downloads.sourceforge.net/project/rpy/rpy/" \
          "%s/rpy-%s%s.zip" % (version, version, ext)
    def _fix_libraries(env):
        run("""sed -i.bak -r -e "s/,'Rlapack'//g" setup.py""")
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):
        result = run("R --version")
        if result.failed:
            return
    _get_install(url, env, _python_make, post_unpack_fn=_fix_libraries)

@_if_not_python_lib("netsa")
def install_netsa_python(env):
    """A suite of open source tools for monitoring large-scale networks using flow data.
    http://tools.netsa.cert.org/index.html
    """
    version = "1.3"
    url = "http://tools.netsa.cert.org/releases/netsa-python-%s.tar.gz" % version
    env.safe_sudo("%s install %s" % (_pip_cmd(env), url))
