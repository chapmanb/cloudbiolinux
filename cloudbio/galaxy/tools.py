import os
import yaml

from cloudbio.custom.bio_general import *
from cloudbio.custom.bio_nextgen import *
from cloudbio.custom.shared import _set_default_config
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
        _chown_galaxy(env, env.galaxy_tools_dir)
        _chown_galaxy(env, env.galaxy_jars_dir)

    if _read_boolean(env, "galaxy_install_r_packages", False):
        _install_r_packages(tools_conf)


def _load_tools_conf(env):
    tools_conf_path = env.get("galaxy_tools_conf",
                              os.path.join(env.config_dir, os.pardir,
                                           "contrib", "flavor", "cloudman", "tools.yaml"))
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
    for (name, versions) in applications.iteritems():
        if type(versions) is str:
            versions = [versions]
        for version_info in versions:
            if type(version_info) is str:
                _install_tool(env, name, version_info)
            else:
                version = version_info["version"]
                bin_dirs = version_info.get("bin_dirs", ["bin"])
                tool_env = _install_tool(env, name, version, bin_dirs)
                symlink_versions = version_info.get("symlink_versions", [])
                if type(symlink_versions) is str:
                    symlink_versions = [symlink_versions]
                for symlink_version in symlink_versions:
                    _set_default_config(tool_env, tool_env["system_install"], symlink_version)


def _install_tool(env, name, version, bin_dirs=["bin"]):
    tool_env = _build_tool_env(env, name, version)
    eval("install_%s" % name)(tool_env)
    _install_galaxy_config(tool_env, bin_dirs)
    return tool_env


def _build_tool_env(env, name, version):
    """ Build new env to have tool installed for Galaxy instead of into /usr. """
    tool_env = {"tool_version": version,
                "galaxy_tool_install": True}
    for key, value in env.iteritems():
        tool_env[key] = value
    tool_env["system_install"] = os.path.join(env.galaxy_tools_dir, name, version)
    return AttributeDict(tool_env)


class AttributeDict(dict):
    """
    Dictionary that allows attribute access to values.

    This is needed because cloudbio.custom.* accesses env extensively via
    attributes (e.g. env.system_install).

    http://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute-in-python
    """
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


def _install_galaxy_config(tool_env, bin_dirs):
    """
    Setup galaxy tool config files (env.sh-es) and default version
    symbolic links.
    """
    install_dir = tool_env["system_install"]
    env_path = os.path.join(install_dir, "env.sh")
    bin_paths = [os.path.join(install_dir, bin_dir) for bin_dir in bin_dirs]
    path_pieces = [bin_path for bin_path in bin_paths if exists(bin_path)]
    if len(path_pieces) > 0 and not exists(env_path):
        path_addtion = ":".join(path_pieces)
        # Standard bin install, just add it to path
        sudo("echo 'PATH=%s:$PATH' > %s/env.sh" % (path_addtion, install_dir))
        sudo("chmod +x %s/env.sh" % install_dir)
    _set_default_config(tool_env, install_dir)
