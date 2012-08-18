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
    for (name, versions) in applications.iteritems():
        if type(versions) is str:
            versions = [versions]
        for version in versions:
            tool_env = _build_tool_env(env, name, version)
            eval("install_%s" % name)(tool_env)
            _install_galaxy_config(tool_env)


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


def _install_galaxy_config(tool_env):
    """
    Setup galaxy tool config files (env.sh-es) and default version
    symbolic links.
    """
    install_dir = tool_env["system_install"]
    bin_dir = os.path.join(install_dir, "bin")
    env_path = os.path.join(install_dir, "env.sh")
    if exists(bin_dir) and not exists(env_path):
        # Standard bin install, just add it to path
        sudo("echo 'PATH=%s/bin:$PATH' > %s/env.sh" % (install_dir, install_dir))
        sudo("chmod +x %s/env.sh" % install_dir)
    _set_default_config(tool_env, install_dir)
