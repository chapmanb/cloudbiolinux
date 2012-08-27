"""Install instructions for python libraries not ready for easy_install.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_python_lib, _get_install, _python_make

@_if_not_python_lib("bx")
def install_bx_python(env):
    """Tools for manipulating biological data, particularly multiple sequence alignments
    https://bitbucket.org/james_taylor/bx-python/wiki/Home
    """
    version = "bitbucket"
    url = "hg clone http://bitbucket.org/james_taylor/bx-python"
    _get_install(url, env, _python_make)

@_if_not_python_lib("matplotlib")
def install_matplotlib(env):
    """matplotlib is a python 2D plotting library which produces publication quality figures
    http://matplotlib.sourceforge.net/
    """
    version = "1.1.1"
    url = "http://downloads.sourceforge.net/project/matplotlib/matplotlib/" \
          "matplotlib-%s/matplotlib-%s.tar.gz" % (version, version)
    _get_install(url, env, _python_make)

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
