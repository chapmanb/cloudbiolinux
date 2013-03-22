"""Main Fabric deployment file for CloudBioLinux distribution.

This installs a standard set of useful biological applications on a remote
server. It is designed for bootstrapping a machine from scratch, as with new
Amazon EC2 instances.

Usage:

    fab -H hostname -i private_key_file install_biolinux

which will call into the 'install_biolinux' method below. See the README for
more examples.

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
from cloudbio.utils import _setup_logging, _update_biolinux_log, _configure_fabric_environment
from cloudbio.cloudman import _cleanup_ec2
from cloudbio.cloudbiolinux import _cleanup_space
from cloudbio.custom.shared import _make_tmp_dir, _pip_cmd
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.package import (_configure_and_install_native_packages,
                              _connect_native_packages)
from cloudbio.package.nix import _setup_nix_sources, _nix_packages
from cloudbio.flavor.config import get_config_file

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
      - libraries    Install programming language libraries
      - post_install Setup CloudMan, FreeNX and other system services
      - cleanup      Remove downloaded files and prepare images for AMI builds
    """
    _setup_logging(env)
    time_start = _print_time_stats("Config", "start")
    _check_fabric_version()
    _configure_fabric_environment(env, flavor)
    env.logger.debug("Target is '%s'" % target)
    _perform_install(target, flavor)
    _print_time_stats("Config", "end", time_start)

def _perform_install(target=None, flavor=None):
    """
    Once CBL/fabric environment is setup, this method actually
    runs the required installation procedures.

    See `install_biolinux` for full details on arguments
    `target` and `flavor`.
    """
    pkg_install, lib_install, custom_ignore = _read_main_config()
    if target is None or target == "packages":
        # can only install native packages if we have sudo access.
        if env.use_sudo:
            _configure_and_install_native_packages(env, pkg_install)
        else:
            _connect_native_packages(env, pkg_install)
        if env.nixpkgs:  # ./doc/nixpkgs.md
            _setup_nix_sources()
            _nix_packages(pkg_install)
        _update_biolinux_log(env, target, flavor)
    if target is None or target == "custom":
        _custom_installs(pkg_install, custom_ignore)
    if target is None or target == "libraries":
        _do_library_installs(lib_install)
    if target is None or target == "post_install":
        env.edition.post_install(pkg_install=pkg_install)
        env.flavor.post_install()
    if target is None or target == "cleanup":
        _cleanup_space(env)
        if env.has_key("is_ec2_image") and env.is_ec2_image.upper() in ["TRUE", "YES"]:
            if env.distribution in ["ubuntu"]:
                # For the time being (Dec 2012), must install development version
                # of cloud-init because of a boto & cloud-init bug:
                # https://bugs.launchpad.net/cloud-init/+bug/1068801
                sudo('wget --output-document=cloud-init_0.7.1-0ubuntu4_all.deb ' +
                    'https://launchpad.net/ubuntu/+archive/primary/+files/cloud-init_0.7.1-0ubuntu4_all.deb')
                sudo("dpkg -i cloud-init_0.7.1-0ubuntu4_all.deb")
                sudo("rm -f cloud-init_0.7.1-0ubuntu4_all.deb")
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

def _custom_installs(to_install, ignore=None):
    if not exists(env.local_install):
        run("mkdir -p %s" % env.local_install)
    pkg_config = get_config_file(env, "custom.yaml").base
    packages, pkg_to_group = _yaml_to_packages(pkg_config, to_install)
    packages = [p for p in packages if ignore is None or p not in ignore]
    for p in env.flavor.rewrite_config_items("custom", packages):
        install_custom(p, True, pkg_to_group)

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
    _setup_logging(env)
    p = p.lower() # All packages listed in custom.yaml are in lower case
    time_start = _print_time_stats("Custom install for '{0}'".format(p), "start")
    if not automated:
        _configure_fabric_environment(env, flavor)
        pkg_config = get_config_file(env, "custom.yaml").base
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)

    try:
        env.logger.debug("Import %s" % p)
        # Allow direct calling of a program install method, even if the program
        # is not listed in the custom list (ie, not contained as a key value in
        # pkg_to_group). For an example, see 'install_cloudman' or use p=cloudman.
        mod_name = pkg_to_group[p] if p in pkg_to_group else p
        mod = __import__("cloudbio.custom.%s" % mod_name,
                         fromlist=["cloudbio", "custom"])
    except ImportError:
        raise ImportError("Need to write a %s module in custom." %
                pkg_to_group[p])
    replace_chars = ["-"]
    try:
        for to_replace in replace_chars:
            p = p.replace(to_replace, "_")
        fn = getattr(mod, "install_%s" % p)
    except AttributeError:
        raise ImportError("Need to write a install_%s function in custom.%s"
                % (p, pkg_to_group[p]))
    fn(env)
    _print_time_stats("Custom install for '%s'" % p, "end", time_start)

def _read_main_config():
    """Pull a list of groups to install based on our main configuration YAML.

    Reads 'main.yaml' and returns packages and libraries
    """
    yaml_file = get_config_file(env, "main.yaml").base
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    packages = full_data.get('packages', [])
    libraries = full_data.get('libraries', [])
    custom_ignore = full_data.get('custom_ignore', [])
    if packages is None: packages = []
    if libraries is None: libraries = []
    if custom_ignore is None: custom_ignore = []
    env.logger.info("Meta-package information from {2}\n- Packages: {0}\n- Libraries: "
            "{1}".format(",".join(packages), ",".join(libraries), yaml_file))
    return packages, sorted(libraries), custom_ignore

# ### Library specific installation code

def _python_library_installer(config):
    """Install python specific libraries using easy_install.
    """
    version_ext = "-%s" % env.python_version_ext if env.python_version_ext else ""
    env.safe_sudo("easy_install%s -U pip" % version_ext)
    for pname in env.flavor.rewrite_config_items("python", config['pypi']):
        env.safe_sudo("{0} install --upgrade {1}".format(_pip_cmd(env), pname))

def _ruby_library_installer(config):
    """Install ruby specific gems.
    """
    gem_ext = getattr(env, "ruby_version_ext", "")
    def _cur_gems():
        with settings(
                hide('warnings', 'running', 'stdout', 'stderr')):
            gem_info = run("gem%s list --no-versions" % gem_ext)
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
    with _make_tmp_dir() as tmp_dir:
        with cd(tmp_dir):
            run("wget --no-check-certificate -O cpanm "
                "https://raw.github.com/miyagawa/cpanminus/master/cpanm")
            run("chmod a+rwx cpanm")
            env.safe_sudo("mv cpanm %s/bin" % env.system_install)
    sudo_str = "--sudo" if env.use_sudo else ""
    for lib in env.flavor.rewrite_config_items("perl", config['cpan']):
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        run("cpanm %s --skip-installed --notest %s < /dev/null" % (sudo_str, lib))

def _haskell_library_installer(config):
    """Install haskell libraries using cabal.
    """
    run("cabal update")
    for lib in config["cabal"]:
        sudo_str = "--root-cmd=sudo" if env.use_sudo else ""
        run("cabal install %s --global %s" % (sudo_str, lib))

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
    _configure_fabric_environment(env)
    _do_library_installs(["%s-libs" % language])

def _do_library_installs(to_install):
    for iname in to_install:
        yaml_file = get_config_file(env, "%s.yaml" % iname).base
        with open(yaml_file) as in_handle:
            config = yaml.load(in_handle)
        lib_installers[iname](config)
