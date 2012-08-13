import yaml

from cloudbio.galaxy.applications import *
from cloudbio.galaxy.r import _install_r_packages
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean

from fabric.api import sudo
from fabric.contrib.files import exists


def _install_tools(env, tools_conf=None):
    """
    Install tools needed for Galaxy along with tool configuration
    directories needed by Galaxy.
    """

    if not tools_conf:
        tools_conf = _load_tools_conf(env)

    if _read_boolean(env, "galaxy_install_dependencies", False):
       # Need to ensure the install dir exists and is owned by env.galaxy_user
        _setup_install_dir(env)
        _install_applications(env, tools_conf)

    if _read_boolean(env, "galaxy_install_r_packages", False):
        _install_r_packages(tools_conf)


def _load_tools_conf(env):
    tools_conf_path = env.get("galaxy_tools_conf", "contrib/cloudman/tools.yaml")
    with open(tools_conf_path) as in_handle:
        full_data = yaml.load(in_handle)
    return full_data


def _setup_install_dir(env):
    """Sets up install dir and ensures its owned by Galaxy"""
    if not exists(env.galaxy_tools_dir):
        sudo("mkdir -p %s" % env.galaxy_tools_dir)
        _chown_galaxy(env, env.galaxy_tools_dir)
    if not exists(env.galaxy_jars_dir):
        sudo("mkdir -p %s" % env.galaxy_jars_dir)
        _chown_galaxy(env, env.galaxy_jars_dir)


def _install_applications(env, tools_conf):
    """Install external tools (galaxy tool dependencies).
    """
    applications = tools_conf["applications"] or {}
    for (name, version) in applications.iteritems():
        eval("install_%s" % name)(env, version)
