"""
"""

from fabric.api import env

from cloudbio.galaxy.tools import _install_tools


def purge_tools():
    env.safe_sudo("rm -rf %s" % env.install_dir)


def install_tools(tools_conf):
    """Deploy a Galaxy server along with some tools.
    """
    _install_tools(env, tools_conf)
