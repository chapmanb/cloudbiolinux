"""Main Fabric deployment file for CloudBioLinux distribution.

This installs a standard set of useful biological applications on a remote
server. It is designed for bootstrapping a machine from scratch, as with new
Amazon EC2 instances.

Usage:

    fab -H hostname -i private_key_file install_biolinux

which will call into the 'install_biolinux' method below. See the README for
more examples. hostname can be a named host in ~/.ssh/config

Requires:
    Fabric http://docs.fabfile.org
    PyYAML http://pyyaml.org/wiki/PyYAMLDocumentation
"""
import os
import sys
from datetime import datetime

from fabric.api import *
from fabric.contrib.files import *
import yaml

# use local cloudbio directory
for to_remove in [p for p in sys.path if p.find("cloudbiolinux-") > 0]:
    sys.path.remove(to_remove)
sys.path.append(os.path.dirname(__file__))
import cloudbio

from cloudbio import libraries
from cloudbio.utils import _setup_logging, _configure_fabric_environment
from cloudbio.cloudman import _cleanup_ec2, _configure_cloudman
from cloudbio.cloudbiolinux import _cleanup_space, _freenx_scripts
from cloudbio.custom import shared
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.package import brew, conda
from cloudbio.package import (_configure_and_install_native_packages,
                              _connect_native_packages, _print_shell_exports)
from cloudbio.package.nix import _setup_nix_sources, _nix_packages
from cloudbio.flavor.config import get_config_file
from cloudbio.config_management.puppet import _puppet_provision
from cloudbio.config_management.chef import _chef_provision, chef, _configure_chef

# ### Shared installation targets for all platforms

def install_biolinux(target=None, flavor=None):
    """Main entry point for installing BioLinux on a remote server.

    `flavor` allows customization of CloudBioLinux behavior. It can either
    be a flavor name that maps to a corresponding directory in contrib/flavor
    or the path to a custom directory. This can contain:

      - alternative package lists (main.yaml, packages.yaml, custom.yaml)
      - custom python code (nameflavor.py) that hooks into the build machinery

    `target` allows running only particular parts of the build process. Valid choices are:

      - packages     Install distro packages
      - custom       Install custom packages
      - chef_recipes Provision chef recipes
      - libraries    Install programming language libraries
      - post_install Setup CloudMan, FreeNX and other system services
      - cleanup      Remove downloaded files and prepare images for AMI builds
    """
    _setup_logging(env)
    time_start = _print_time_stats("Config", "start")
    _check_fabric_version()
    if env.ssh_config_path and os.path.isfile(os.path.expanduser(env.ssh_config_path)):
        env.use_ssh_config = True
    _configure_fabric_environment(env, flavor,
                                  ignore_distcheck=(target is not None
                                                    and target in ["libraries", "custom"]))
    env.logger.debug("Target is '%s'" % target)
    env.logger.debug("Flavor is '%s'" % flavor)
    _perform_install(target, flavor)
    _print_time_stats("Config", "end", time_start)
    if hasattr(env, "keep_isolated") and env.keep_isolated:
        _print_shell_exports(env)

def _perform_install(target=None, flavor=None, more_custom_add=None):
    """
    Once CBL/fabric environment is setup, this method actually
    runs the required installation procedures.

    See `install_biolinux` for full details on arguments
    `target` and `flavor`.
    """
    pkg_install, lib_install, custom_ignore, custom_add = _read_main_config()
    if more_custom_add:
        if custom_add is None:
            custom_add = {}
        for k, vs in more_custom_add.iteritems():
            if k in custom_add:
                custom_add[k].extend(vs)
            else:
                custom_add[k] = vs
    if target is None or target == "packages":
        env.keep_isolated = getattr(env, "keep_isolated", "false").lower() in ["true", "yes"]
        # Only touch system information if we're not an isolated installation
        if not env.keep_isolated:
            # can only install native packages if we have sudo access or are root
            if env.use_sudo or env.safe_run_output("whoami").strip() == "root":
                _configure_and_install_native_packages(env, pkg_install)
            else:
                _connect_native_packages(env, pkg_install, lib_install)
        if env.nixpkgs:  # ./doc/nixpkgs.md
            _setup_nix_sources()
            _nix_packages(pkg_install)
    if target is None or target == "custom":
        _custom_installs(pkg_install, custom_ignore, custom_add)
    if target is None or target == "chef_recipes":
        _provision_chef_recipes(pkg_install, custom_ignore)
    if target is None or target == "puppet_classes":
        _provision_puppet_classes(pkg_install, custom_ignore)
    if target is None or target == "brew":
        install_brew(flavor=flavor, automated=True)
    if target is None or target == "conda":
        install_conda(flavor=flavor, automated=True)
    if target is None or target == "libraries":
        _do_library_installs(lib_install)
    if target is None or target == "post_install":
        env.flavor.post_install()
        if "is_ec2_image" in env and env.is_ec2_image.upper() in ["TRUE", "YES"]:
            _freenx_scripts(self.env)
            if pkg_install is not None and 'cloudman' in pkg_install:
                _configure_cloudman(self.env)
    if target is None or target == "cleanup":
        _cleanup_space(env)
        if "is_ec2_image" in env and env.is_ec2_image.upper() in ["TRUE", "YES"]:
            _cleanup_ec2(env)

def _print_time_stats(action, event, prev_time=None):
    """ A convenience method for displaying time event during configuration.

    :type action: string
    :param action: Indicates type of action (eg, Config, Lib install, Pkg install)

    :type event: string
    :param event: The monitoring event (eg, start, stop)

    :type prev_time: datetime
    :param prev_time: A timeststamp of a previous event. If provided, duration between
                      the time the method is called and the time stamp is included in
                      the printout

    :rtype: datetime
    :return: A datetime timestamp of when the method was called
    """
    time = datetime.utcnow()
    s = "{0} {1} time: {2}".format(action, event, time)
    if prev_time: s += "; duration: {0}".format(str(time-prev_time))
    env.logger.info(s)
    return time

def _check_fabric_version():
    """Checks for fabric version installed
    """
    version = env.version
    if int(version.split(".")[0]) < 1:
        raise NotImplementedError("Please install fabric version 1 or higher")

def _custom_installs(to_install, ignore=None, add=None):
    if not env.safe_exists(env.local_install) and env.local_install:
        env.safe_run("mkdir -p %s" % env.local_install)
    pkg_config = get_config_file(env, "custom.yaml").base
    packages, pkg_to_group = _yaml_to_packages(pkg_config, to_install)
    packages = [p for p in packages if ignore is None or p not in ignore]
    if add is not None:
        for key, vals in add.iteritems():
            for v in vals:
                pkg_to_group[v] = key
                packages.append(v)
    for p in env.flavor.rewrite_config_items("custom", packages):
        install_custom(p, True, pkg_to_group)


def _provision_chef_recipes(to_install, ignore=None):
    """
    Much like _custom_installs, read config file, determine what to install,
    and install it.
    """
    pkg_config = get_config_file(env, "chef_recipes.yaml").base
    packages, _ = _yaml_to_packages(pkg_config, to_install)
    packages = [p for p in packages if ignore is None or p not in ignore]
    recipes = [recipe for recipe in env.flavor.rewrite_config_items("chef_recipes", packages)]
    if recipes:  # Don't bother running chef if nothing to configure
        install_chef_recipe(recipes, True)


def _provision_puppet_classes(to_install, ignore=None):
    """
    Much like _custom_installs, read config file, determine what to install,
    and install it.
    """
    pkg_config = get_config_file(env, "puppet_classes.yaml").base
    packages, _ = _yaml_to_packages(pkg_config, to_install)
    packages = [p for p in packages if ignore is None or p not in ignore]
    classes = [recipe for recipe in env.flavor.rewrite_config_items("puppet_classes", packages)]
    if classes:  # Don't bother running chef if nothing to configure
        install_puppet_class(classes, True)


def install_chef_recipe(recipe, automated=False, flavor=None):
    """Install one or more chef recipes by name.

    Usage: fab [-i key] [-u user] -H host install_chef_recipe:recipe

    :type recipe:  string or list
    :param recipe: TODO

    :type automated:  bool
    :param automated: If set to True, the environment is not loaded.
    """
    _setup_logging(env)
    if not automated:
        _configure_fabric_environment(env, flavor)

    time_start = _print_time_stats("Chef provision for recipe(s) '{0}'".format(recipe), "start")
    _configure_chef(env, chef)
    recipes = recipe if isinstance(recipe, list) else [recipe]
    for recipe_to_add in recipes:
        chef.add_recipe(recipe_to_add)
    _chef_provision(env, recipes)
    _print_time_stats("Chef provision for recipe(s) '%s'" % recipe, "end", time_start)


def install_puppet_class(classes, automated=False, flavor=None):
    """Install one or more puppet classes by name.

    Usage: fab [-i key] [-u user] -H host install_puppet_class:class

    :type classes:  string or list
    :param classes: TODO

    :type automated:  bool
    :param automated: If set to True, the environment is not loaded.
    """
    _setup_logging(env)
    if not automated:
        _configure_fabric_environment(env, flavor)

    time_start = _print_time_stats("Puppet provision for class(es) '{0}'".format(classes), "start")
    classes = classes if isinstance(classes, list) else [classes]
    _puppet_provision(env, classes)
    _print_time_stats("Puppet provision for classes(s) '%s'" % classes, "end", time_start)


def install_custom(p, automated=False, pkg_to_group=None, flavor=None):
    """
    Install a single custom program or package by name.

    This method fetches program name from ``config/custom.yaml`` and delegates
    to a method in ``custom/*name*.py`` to proceed with the installation.
    Alternatively, if a program install method is defined in the appropriate
    package, it will be called directly (see param ``p``).

    Usage: fab [-i key] [-u user] -H host install_custom:program_name

    :type p:  string
    :param p: A name of the custom program to install. This has to be either a name
              that is listed in ``custom.yaml`` as a subordinate to a group name or a
              program name whose install method is defined in either ``cloudbio`` or
              ``custom`` packages
              (e.g., ``cloudbio/custom/cloudman.py -> install_cloudman``).

    :type automated:  bool
    :param automated: If set to True, the environment is not loaded and reading of
                      the ``custom.yaml`` is skipped.
    """
    p = p.lower() # All packages listed in custom.yaml are in lower case
    if not automated:
        _setup_logging(env)
        _configure_fabric_environment(env, flavor, ignore_distcheck=True)
        pkg_config = get_config_file(env, "custom.yaml").base
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)
    time_start = _print_time_stats("Custom install for '{0}'".format(p), "start")
    fn = _custom_install_function(env, p, pkg_to_group)
    fn(env)
    ## TODO: Replace the previous 4 lines with the following one, barring
    ## objections. Slightly different behavior because pkg_to_group will be
    ## loaded regardless of automated if it is None, but IMO this shouldn't
    ## matter because the following steps look like they would fail if
    ## automated is True and pkg_to_group is None.
    # _install_custom(p, pkg_to_group)
    _print_time_stats("Custom install for '%s'" % p, "end", time_start)


def _install_custom(p, pkg_to_group=None):
    if pkg_to_group is None:
        pkg_config = get_config_file(env, "custom.yaml").base
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)
    fn = _custom_install_function(env, p, pkg_to_group)
    fn(env)

def install_brew(p=None, version=None, flavor=None, automated=False):
    """Top level access to homebrew/linuxbrew packages.
    p is a package name to install, or all configured packages if not specified.
    """
    if not automated:
        _setup_logging(env)
        _configure_fabric_environment(env, flavor, ignore_distcheck=True)
    if p is not None:
        if version:
            p = "%s==%s" % (p, version)
        brew.install_packages(env, packages=[p])
    else:
        pkg_install = _read_main_config()[0]
        brew.install_packages(env, to_install=pkg_install)

def install_conda(p=None, flavor=None, automated=False):
    if not automated:
        _setup_logging(env)
        _configure_fabric_environment(env, flavor, ignore_distcheck=True)
    if p is not None:
        conda.install_packages(env, packages=[p])
    else:
        pkg_install = _read_main_config()[0]
        conda.install_packages(env, to_install=pkg_install)

def _custom_install_function(env, p, pkg_to_group):
    """
    Find custom install function to execute based on package name to
    pkg_to_group dict.
    """
    try:
        # Allow direct calling of a program install method, even if the program
        # is not listed in the custom list (ie, not contained as a key value in
        # pkg_to_group). For an example, see 'install_cloudman' or use p=cloudman.
        mod_name = pkg_to_group[p] if p in pkg_to_group else p
        env.logger.debug("Importing module cloudbio.custom.%s" % mod_name)
        mod = __import__("cloudbio.custom.%s" % mod_name,
                         fromlist=["cloudbio", "custom"])
    except ImportError:
        raise ImportError("Need to write module cloudbio.custom.%s" %
                pkg_to_group[p])
    replace_chars = ["-"]
    try:
        for to_replace in replace_chars:
            p = p.replace(to_replace, "_")
        env.logger.debug("Looking for custom install function %s.install_%s"
            % (mod.__name__, p))
        fn = getattr(mod, "install_%s" % p)
    except AttributeError:
        raise ImportError("Need to write a install_%s function in custom.%s"
                % (p, pkg_to_group[p]))
    return fn


def _read_main_config():
    """Pull a list of groups to install based on our main configuration YAML.

    Reads 'main.yaml' and returns packages and libraries
    """
    yaml_file = get_config_file(env, "main.yaml").base
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    packages = full_data.get('packages', [])
    packages = env.flavor.rewrite_config_items("main_packages", packages)
    libraries = full_data.get('libraries', [])
    custom_ignore = full_data.get('custom_ignore', [])
    custom_add = full_data.get("custom_additional")
    if packages is None: packages = []
    if libraries is None: libraries = []
    if custom_ignore is None: custom_ignore = []
    env.logger.info("Meta-package information from {2}\n- Packages: {0}\n- Libraries: "
            "{1}".format(",".join(packages), ",".join(libraries), yaml_file))
    return packages, sorted(libraries), custom_ignore, custom_add

# ### Library specific installation code

def _python_library_installer(config):
    """Install python specific libraries using pip, conda and easy_install.
    Handles using isolated anaconda environments.
    """
    if shared._is_anaconda(env):
        conda_bin = shared._conda_cmd(env)
        for pname in env.flavor.rewrite_config_items("python", config.get("conda", [])):
            env.safe_run("{0} install --yes {1}".format(conda_bin, pname))
        cmd = env.safe_run
        with settings(warn_only=True):
            cmd("%s -U distribute" % os.path.join(os.path.dirname(conda_bin), "easy_install"))
    else:
        pip_bin = shared._pip_cmd(env)
        ei_bin = pip_bin.replace("pip", "easy_install")
        env.safe_sudo("%s -U pip" % ei_bin)
        with settings(warn_only=True):
            env.safe_sudo("%s -U distribute" % ei_bin)
        cmd = env.safe_sudo
    for pname in env.flavor.rewrite_config_items("python", config['pypi']):
        cmd("{0} install --upgrade {1} --allow-unverified {1} --allow-external {1}".format(shared._pip_cmd(env), pname)) # fixes problem with packages not being in pypi

def _ruby_library_installer(config):
    """Install ruby specific gems.
    """
    gem_ext = getattr(env, "ruby_version_ext", "")
    def _cur_gems():
        with settings(
                hide('warnings', 'running', 'stdout', 'stderr')):
            gem_info = env.safe_run_output("gem%s list --no-versions" % gem_ext)
        return [l.rstrip("\r") for l in gem_info.split("\n") if l.rstrip("\r")]
    installed = _cur_gems()
    for gem in env.flavor.rewrite_config_items("ruby", config['gems']):
        # update current gems only to check for new installs
        if gem not in installed:
            installed = _cur_gems()
        if gem in installed:
            env.safe_sudo("gem%s update %s" % (gem_ext, gem))
        else:
            env.safe_sudo("gem%s install %s" % (gem_ext, gem))

def _perl_library_installer(config):
    """Install perl libraries from CPAN with cpanminus.
    """
    with shared._make_tmp_dir() as tmp_dir:
        with cd(tmp_dir):
            env.safe_run("wget --no-check-certificate -O cpanm "
                         "https://raw.github.com/miyagawa/cpanminus/master/cpanm")
            env.safe_run("chmod a+rwx cpanm")
            env.safe_sudo("mv cpanm %s/bin" % env.system_install)
    sudo_str = "--sudo" if env.use_sudo else ""
    for lib in env.flavor.rewrite_config_items("perl", config['cpan']):
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        env.safe_run("cpanm %s --skip-installed --notest %s < /dev/null" % (sudo_str, lib))

def _haskell_library_installer(config):
    """Install haskell libraries using cabal.
    """
    run("cabal update")
    for lib in config["cabal"]:
        sudo_str = "--root-cmd=sudo" if env.use_sudo else ""
        env.safe_run("cabal install %s --global %s" % (sudo_str, lib))

lib_installers = {
    "r-libs" : libraries.r_library_installer,
    "python-libs" : _python_library_installer,
    "ruby-libs" : _ruby_library_installer,
    "perl-libs" : _perl_library_installer,
    "haskell-libs": _haskell_library_installer,
    }

def install_libraries(language):
    """High level target to install libraries for a specific language.
    """
    _setup_logging(env)
    _check_fabric_version()
    _configure_fabric_environment(env, ignore_distcheck=True)
    _do_library_installs(["%s-libs" % language])

def _do_library_installs(to_install):
    for iname in to_install:
        yaml_file = get_config_file(env, "%s.yaml" % iname).base
        with open(yaml_file) as in_handle:
            config = yaml.load(in_handle)
        lib_installers[iname](config)
