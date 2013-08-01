import os
from string import Template

import yaml

from cloudbio.custom.bio_general import *
from cloudbio.custom.bio_nextgen import *
from cloudbio.custom.bio_proteomics import *
from cloudbio.custom.shared import _set_default_config, _add_to_profiles
from cloudbio.galaxy.applications import *
from cloudbio.galaxy.r import _install_r_packages
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean

from fabric.api import sudo
from fabric.contrib.files import exists

FAILED_INSTALL_MESSAGE = \
    "Failed to install application %s as a Galaxy application. This may be a transient problem (e.g. mirror for download is currently unavailable) or misconfiguration. The contents of the CloudBioLinux temporary working directory may need to be deleted."


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
        _install_configured_applications(env, tools_conf)
        _chown_galaxy(env, env.galaxy_tools_dir)
        _chown_galaxy(env, env.galaxy_jars_dir)

    if _read_boolean(env, "galaxy_install_r_packages", False):
        _install_r_packages(tools_conf)


def _tools_conf_path(env):
    """
    Load path to galaxy_tools_conf file from env, allowing expansion of $__contrib_dir__.
    Default to $__contrib_dir__/flavor/cloudman/tools.yaml.
    """
    contrib_dir = os.path.join(env.config_dir, os.pardir, "contrib")
    default_tools_conf_path = os.path.join(contrib_dir, "flavor", "cloudman", "tools.yaml")
    tools_conf_path = env.get("galaxy_tools_conf", default_tools_conf_path)
    ## Allow expansion of __config_dir__ in galaxy_tools_conf property.
    return Template(tools_conf_path).safe_substitute({"__contrib_dir__": contrib_dir})


def _load_tools_conf(env):
    with open(_tools_conf_path(env)) as in_handle:
        full_data = yaml.load(in_handle)
    return full_data


def _setup_install_dir(env):
    """Sets up install dir and ensures its owned by Galaxy"""
    if not exists(env.galaxy_tools_dir):
        sudo("mkdir -p %s" % env.galaxy_tools_dir)
        _chown_galaxy(env, env.galaxy_tools_dir)
    # Create a general-purpose ``bin`` directory under the galaxy_tools_dir
    # and put it on the PATH so users can more easily add custom tools
    bin_dir = os.path.join(env.galaxy_tools_dir, 'bin')
    if not exists(bin_dir):
        sudo("mkdir -p %s" % bin_dir)
        _chown_galaxy(env, bin_dir)
        line = "export PATH={0}:$PATH".format(bin_dir)
        _add_to_profiles(line)
    if not exists(env.galaxy_jars_dir):
        sudo("mkdir -p %s" % env.galaxy_jars_dir)
        _chown_galaxy(env, env.galaxy_jars_dir)


def _install_configured_applications(env, tools_conf):
    """
    Install external tools defined by YAML or dictionary data structure.  Instead of
    installing in system_install (e.g. /usr), these custom tools will be installed as
    Galaxy dependency applications.
    """
    applications = tools_conf["applications"] or {}
    for (name, tool_conf) in applications.iteritems():
        try:
            _install_application(name, tool_conf)
        except:
            env.logger.warn(FAILED_INSTALL_MESSAGE % name)
            raise


def _install_application(name, versions):
    """
    Install single custom tool as Galaxy dependency application.

    TODO: Rename versions and document options.
    """
    if type(versions) is str:
        versions = [versions]
    for version_info in versions:
        if type(version_info) is str:
            _install_tool(env, name, version=version_info, requirement_name=name)
        else:
            version = version_info["version"]
            bin_dirs = version_info.get("bin_dirs", ["bin"])
            env_vars = version_info.get("env_vars", {})
            provides = version_info.get("provides", [])
            if isinstance(provides, (str, unicode, basestring)):
                provides = [provides]

            # Some requirements (e.g. blast+) maybe not have valid python
            # identifiers as name. Use install_blast to setup but override
            # requirement directory name with requirement_name field.
            requirement_name = version_info.get("requirement_name", name)
            tool_env = _install_tool(env, name, version, bin_dirs=bin_dirs, env_vars=env_vars, requirement_name=requirement_name)
            symlink_versions = version_info.get("symlink_versions", [])
            if type(symlink_versions) is str:
                symlink_versions = [symlink_versions]
            for symlink_version in symlink_versions:
                _set_default_config(tool_env, tool_env["system_install"], symlink_version)

            if provides:
                install_dir = tool_env["system_install"]
                ## Create additional symlinked packages from this one.
                tool_dir = "%s/.." % install_dir
                tools_dir = "%s/.." % tool_dir
                for package in provides:
                    link_dir = "%s/%s" % (tools_dir, package)
                    env.safe_sudo("ln -f -s '%s' '%s'" % (requirement_name, link_dir))


def _install_tool(env, name, version, requirement_name, bin_dirs=["bin"], env_vars={}):
    tool_env = _build_tool_env(env, requirement_name, version)
    env.logger.debug("Installing a Galaxy tool via 'install_%s'" % name)
    eval("install_%s" % name)(tool_env)
    _install_galaxy_config(tool_env, bin_dirs, env_vars=env_vars)
    return tool_env


def _build_tool_env(env, name, version):
    """ Build new env to have tool installed for Galaxy instead of into /usr. """
    tool_env = {"tool_version": version,
                "galaxy_tool_install": True}
    for key, value in env.iteritems():
        tool_env[key] = value
    tool_env["system_install"] = os.path.join(env.galaxy_tools_dir, name, version)
    tool_env["local_install"] = os.path.join(env.galaxy_tools_dir, name, version)
    tool_env["venv_directory"] = "%s/%s" % (tool_env["system_install"], "venv")
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


def _install_galaxy_config(tool_env, bin_dirs, env_vars):
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
        sudo("echo 'PATH=%s:$PATH' > %s" % (path_addtion, env_path))
        venv_path = "%s/%s" % (install_dir, "venv")
        if exists(venv_path):
            #  Have env.sh activate virtualdirectory
            sudo("echo '. %s/bin/activate' >> %s" % (venv_path, env_path))
        sudo("chmod +x %s" % env_path)
        for env_var, env_var_value in env_vars.iteritems():
            env_var_template = Template(env_var_value)
            expanded_env_var_value = env_var_template.substitute(tool_env)
            sudo("echo 'export %s=%s' >> %s" % (env_var, expanded_env_var_value, env_path))
        env.logger.debug("Added Galaxy env.sh file: %s" % env_path)

    _set_default_config(tool_env, install_dir)
    if _read_boolean(tool_env, "autoload_galaxy_tools", False) and exists(env_path):
        # In this case, the web user (e.g. ubuntu) should auto-load all of
        # galaxy's default env.sh files so they are available for direct use
        # as well.
        _add_to_profiles(". %s" % env_path, profiles=["~/.bashrc"])
