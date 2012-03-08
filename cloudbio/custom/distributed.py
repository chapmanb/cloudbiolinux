"""Install instructions for distributed MapReduce style programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed, _make_tmp_dir, _if_not_python_lib, _fetch_and_unpack

def install_mahout(env):
    # ToDo setup mahout, must be checked out from repo ATM:
    # https://cwiki.apache.org/MAHOUT/mahoutec2.html
    pass
