import os
from string import Template

import six
import yaml

from cloudbio.custom.bio_general import *
from cloudbio.custom.bio_nextgen import *
from cloudbio.custom.bio_proteomics import *
from cloudbio.custom.shared import _set_default_config, _add_to_profiles
from cloudbio.galaxy.applications import *
from cloudbio.galaxy.r import _install_r_packages
from cloudbio.galaxy.utils import _chown_galaxy, _read_boolean

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
        full_data = yaml.safe_load(in_handle)
    return full_data


def _setup_install_dir(env):
    """Sets up install dir and ensures its owned by Galaxy"""
    if not env.safe_exists(env.galaxy_tools_dir):
        env.safe_sudo("mkdir -p %s" % env.galaxy_tools_dir)
        _chown_galaxy(env, env.galaxy_tools_dir)
    # Create a general-purpose ``bin`` directory under the galaxy_tools_dir
    # and put it on the PATH so users can more easily add custom tools
    bin_dir = os.path.join(env.galaxy_tools_dir, 'bin')
    if not env.safe_exists(bin_dir):
        env.safe_sudo("mkdir -p %s" % bin_dir)
        _chown_galaxy(env, bin_dir)
        line = "export PATH={0}:$PATH".format(bin_dir)
        _add_to_profiles(line)
    if not env.safe_exists(env.galaxy_jars_dir):
        env.safe_sudo("mkdir -p %s" % env.galaxy_jars_dir)
        _chown_galaxy(env, env.galaxy_jars_dir)


def _install_configured_applications(env, tools_conf):
    """
    Install external tools defined by YAML or dictionary data structure.  Instead of
    installing in system_install (e.g. /usr), these custom tools will be installed as
    Galaxy dependency applications.
    """
    applications = tools_conf["applications"] or {}
    # Changing the default behavior here to install all tools and
    # just record exceptions as they occur, but wait until the end
    # raise an exception out of this block. Disable this behavior
    # by setting galaxay_tool_defer_errors to False.
    defer_errors = env.get("galaxy_tool_defer_errors", True)
    exceptions = {}
    for (name, tool_conf) in applications.iteritems():
        if not __check_conditional(tool_conf):
            continue

        try:
            _install_application(name, tool_conf)
        except BaseException as e:
            exceptions[name] = e
            if not defer_errors:
                break

    if exceptions:
        for name, exception in exceptions.iteritems():
            env.logger.warn(FAILED_INSTALL_MESSAGE % name)
        first_exception = list(exceptions.values())[0]
        raise first_exception


def _install_application(name, versions, tool_install_dir=None):
    """
    Install single custom tool as Galaxy dependency application.

    TODO: Rename versions and document options.
    """
    if type(versions) is str:
        versions = [versions]
    for version_info in versions:
        if type(version_info) is str:
            _install_tool(env, name, version=version_info, requirement_name=name, tool_install_dir=tool_install_dir)
        else:
            version = version_info["version"]
            bin_dirs = version_info.get("bin_dirs", ["bin"])
            env_vars = version_info.get("env_vars", {})
            provides = version_info.get("provides", [])
            if isinstance(provides, (str, unicode, six.string_types)):
                provides = [provides]
            for provide_conf in provides[:]:
                if isinstance(provide_conf, dict):
                    provides.remove(provide_conf)
                    if __check_conditional(provide_conf):
                        provies.append(provide_conf["name"])

            # Some requirements (e.g. blast+) maybe not have valid python
            # identifiers as name. Use install_blast to setup but override
            # requirement directory name with requirement_name field.
            requirement_name = version_info.get("requirement_name", name)
            tool_env = _install_tool(env, name, version, bin_dirs=bin_dirs, env_vars=env_vars, requirement_name=requirement_name, tool_install_dir=tool_install_dir)
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


def _install_tool(env, name, version, requirement_name, bin_dirs=["bin"], env_vars={}, tool_install_dir=None):
    tool_env = _build_tool_env(env, requirement_name, version, tool_install_dir)
    env.logger.debug("Installing a Galaxy tool via 'install_%s'" % name)
    eval("install_%s" % name)(tool_env)
    _install_galaxy_config(tool_env, bin_dirs, env_vars=env_vars)
    return tool_env


def _build_tool_env(env, name, version, tool_install_dir):
    """ Build new env to have tool installed for Galaxy instead of into /usr. """
    tool_env = {"tool_version": version,
                "galaxy_tool_install": True}
    for key, value in env.iteritems():
        tool_env[key] = value
    if not tool_install_dir:
        tool_install_dir = os.path.join(env.galaxy_tools_dir, name, version)
    tool_env["system_install"] = tool_install_dir
    tool_env["local_install"] = tool_install_dir
    tool_env["venv_directory"] = "%s/%s" % (tool_env["system_install"], "venv")
    return AttributeDict(tool_env)


def __check_conditional(conf_dict):
    passes = True
    try:
        if "if" in conf_dict:
            value = conf_dict["if"]
            passes = _read_boolean(env, value, False)
        elif "unless" in conf_dict:
            value = conf_dict["unless"]
            passes = not _read_boolean(env, value, False)
    except TypeError:
        # configuration is not a dictionary, default to True
        pass
    return passes


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
    path_pieces = [bin_path for bin_path in bin_paths if env.safe_exists(bin_path)]
    if len(path_pieces) > 0 and not env.safe_exists(env_path):
        path_addtion = ":".join(path_pieces)
        # Standard bin install, just add it to path
        env.safe_sudo("echo 'PATH=%s:$PATH' > %s" % (path_addtion, env_path))
        venv_path = "%s/%s" % (install_dir, "venv")
        if env.safe_exists(venv_path):
            #  Have env.sh activate virtualdirectory
            env.safe_sudo("echo '. %s/bin/activate' >> %s" % (venv_path, env_path))
        env.safe_sudo("chmod +x %s" % env_path)
        for env_var, env_var_value in env_vars.iteritems():
            env_var_template = Template(env_var_value)
            expanded_env_var_value = env_var_template.substitute(tool_env)
            env.safe_sudo("echo 'export %s=%s' >> %s" % (env_var, expanded_env_var_value, env_path))
        env.logger.debug("Added Galaxy env.sh file: %s" % env_path)

    # TODO: If a direct install (i.e. tool_install_dir specified instead of galaxy_tools_dir)
    # default is still setup. This is not really desired.
    _set_default_config(tool_env, install_dir)
    if _read_boolean(tool_env, "autoload_galaxy_tools", False) and env.safe_exists(env_path):
        # In this case, the web user (e.g. ubuntu) should auto-load all of
        # galaxy's default env.sh files so they are available for direct use
        # as well.
        _add_to_profiles(". %s" % env_path, profiles=["~/.bashrc"])
