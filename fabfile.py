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

from fabric.main import load_settings
from fabric.api import *
from fabric.contrib.files import *
import yaml

# use local cloudbio directory
for to_remove in [p for p in sys.path if p.find("cloudbiolinux-") > 0]:
    sys.path.remove(to_remove)
sys.path.append(os.path.dirname(__file__))
import cloudbio

from cloudbio.edition import _setup_edition
from cloudbio.distribution import _setup_distribution_environment
from cloudbio.utils import _setup_logging, _update_biolinux_log
from cloudbio.cloudman import _cleanup_ec2
from cloudbio.cloudbiolinux import _cleanup_space
from cloudbio.custom.shared import _make_tmp_dir
from cloudbio.package.shared import _yaml_to_packages
from cloudbio.package.deb import (_apt_packages, _add_apt_gpg_keys,
                                  _setup_apt_automation, _setup_apt_sources)
from cloudbio.package.rpm import (_yum_packages, _setup_yum_bashrc,
                                  _setup_yum_sources)
from cloudbio.package.nix import _setup_nix_sources, _nix_packages

# ## Utility functions for establishing our build environment

def _parse_fabricrc():
    """Defaults from fabricrc.txt file; loaded if not specified at commandline.
    """
    # ## General setup
    env.config_dir = os.path.join(os.path.dirname(__file__), "config")

    if not env.has_key("distribution") and not env.has_key("system_install"):
        env.logger.info("Reading default fabricrc.txt")
        config_file = os.path.join(env.config_dir, "fabricrc.txt")
        if os.path.exists(config_file):
            env.update(load_settings(config_file))
    else:
        env.logger.warn("Skipping fabricrc.txt as distribution is already defined")

def _create_local_paths():
    """Expand any paths defined in terms of shell shortcuts (like ~).
    """
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):
        # This is the first point we call into a remote host - make sure
        # it does not fail silently by calling a dummy run
        env.logger.info("Now, testing connection to host...")
        test = run("pwd")
        # If there is a connection failure, the rest of the code is (sometimes) not
        # reached - for example with Vagrant the program just stops after above run
        # command.
        if test != None:
          env.logger.info("Connection to host appears to work!")
        else:
          raise NotImplementedError("Connection to host failed")
        env.logger.debug("Expand paths")
        if env.has_key("local_install"):
            if not exists(env.local_install):
                run("mkdir -p %s" % env.local_install)
            with cd(env.local_install):
                result = run("pwd")
                env.local_install = result

def _setup_flavor(flavor, environment=None):
    """Setup flavor
    """
    if not flavor:
        flavor = env.get("flavor", None)
    if not environment:
        environment = env.get("environment", None)
    if environment:
        env.environment = environment
        env.logger.info("Environment %s" % env.environment)
    if flavor:
        # import a flavor defined through parameters flavor and flavor_path
        flavor_path = env.get("flavor_path", None)
        if flavor_path == None:
          raise ImportError("You need to define the flavor_path for %s!" % flavor)
        # Add path for flavors
        sys.path.append(os.path.join(os.path.dirname(__file__), "contrib", "flavor"))
        env.logger.info("Flavor %s loaded from %s" % (flavor, flavor_path))
        try:
            mod = __import__(flavor_path, fromlist=[flavor])
        except ImportError:
            raise ImportError("Failed to import %s" % flavor)
    else:
        # import default Flavor
        from cloudbio.flavor import Flavor
    env.logger.info("This is a %s" % env.flavor.name)

# ### Shared installation targets for all platforms

def install_biolinux(target=None, packagelist=None, flavor=None, environment=None):
    """Main entry point for installing Biolinux on a remote server.

    This allows a different main package list (the main YAML file is passed in),
    and/or use of Flavor. So you can say:

      install_biolinux:packagelist=contrib/mylist/main.yaml,flavor=specialflavor

    Both packagelist and flavor, as well as the Edition, can also be passed in
    through the fabricrc file.

    target can also be supplied on the fab CLI. Special targets are:

      - packages     Install distro packages
      - custom       Install custom packages
      - libraries    Install programming language libraries
      - post_install Setup CloudMan, FreeNX and other system services
      - cleanup      Remove downloaded files and prepare images for AMI builds

    environment allows adding additional information on the command line -
    usually for defining environments, for example environment=testing, or
    environment=production, will set the deployment environment and tune
    post-installation settings.
    """
    _setup_logging(env)
    _check_fabric_version()
    _parse_fabricrc()
    _setup_edition(env)
    _setup_flavor(flavor, environment)
    _setup_distribution_environment() # get parameters for distro, packages etc.
    _create_local_paths()
    env.logger.info("packagelist=%s" % packagelist)
    pkg_install, lib_install = _read_main_config(packagelist)  # read yaml
    env.logger.info("Target=%s" % target)
    if target is None or target == "packages":
        if env.distribution in ["debian", "ubuntu"]:
            _setup_apt_sources()
            _setup_apt_automation()
            _add_apt_gpg_keys()
            _apt_packages(pkg_install)
        elif env.distibution in ["centos"]:
            _setup_yum_sources()
            _yum_packages(pkg_install)
            _setup_yum_bashrc()
        else:
            raise NotImplementedError("Unknown target distribution")
        if env.nixpkgs: # ./doc/nixpkgs.md
            _setup_nix_sources()
            _nix_packages(pkg_install)
        _update_biolinux_log(env, target, flavor)
    if target is None or target == "custom":
        _custom_installs(pkg_install)
    if target is None or target == "libraries":
        _do_library_installs(lib_install)
    if target is None or target == "post_install":
        env.edition.post_install()
        env.flavor.post_install()
    if target is None or target == "cleanup":
        _cleanup_space(env)
        if env.has_key("is_ec2_image") and env.is_ec2_image.upper() in ["TRUE", "YES"]:
            _cleanup_ec2(env)

def _check_fabric_version():
    """Checks for fabric version installed
    """
    version = env.version
    if int(version.split(".")[0]) < 1:
        raise NotImplementedError("Please install fabric version 1 or higher")

def _custom_installs(to_install):
    if not exists(env.local_install):
        run("mkdir -p %s" % env.local_install)
    pkg_config = os.path.join(env.config_dir, "custom.yaml")
    packages, pkg_to_group = _yaml_to_packages(pkg_config, to_install)
    for p in env.flavor.rewrite_config_items("custom", packages):
        install_custom(p, True, pkg_to_group)

def install_custom(p, automated=False, pkg_to_group=None):
    """Install a single custom package by name.

    This method fetches names from custom.yaml that delegate to a method
    in the custom/name.py program.

    fab install_custom_package:package_name
    """
    _setup_logging(env)
    env.logger.info("Install custom software packages")
    if not automated:
        _parse_fabricrc()
        _setup_edition(env)
        _setup_distribution_environment()
        _create_local_paths()
        pkg_config = os.path.join(env.config_dir, "custom.yaml")
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)
    try:
        env.logger.debug("Import %s" % p)
        mod = __import__("cloudbio.custom.%s" % pkg_to_group[p],
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

def _read_main_config(yaml_file=None):
    """Pull a list of groups to install based on our main configuration YAML.

    Reads 'main.yaml' and returns packages and libraries
    """
    if yaml_file is None:
        yaml_file = os.path.join(env.config_dir, "main.yaml")
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    packages = full_data['packages']
    packages = packages if packages else []
    libraries = full_data['libraries']
    libraries = libraries if libraries else []
    env.logger.info("Meta-package information")
    env.logger.info(",".join(packages))
    env.logger.info(",".join(libraries))
    return packages, sorted(libraries)

# ### Library specific installation code

def _r_library_installer(config):
    """Install R libraries using CRAN and Bioconductor.
    """
    # Create an Rscript file with install details.
    out_file = "install_packages.R"
    if exists(out_file):
        run("rm -f %s" % out_file)
    run("touch %s" % out_file)
    repo_info = """
    cran.repos <- getOption("repos")
    cran.repos["CRAN" ] <- "%s"
    options(repos=cran.repos)
    source("%s")
    """ % (config["cranrepo"], config["biocrepo"])
    append(out_file, repo_info)
    install_fn = """
    repo.installer <- function(repos, install.fn) {
      update.or.install <- function(pname) {
        if (pname %in% installed.packages())
          update.packages(lib.loc=c(pname), repos=repos, ask=FALSE)
        else
          install.fn(pname)
      }
    }
    """
    append(out_file, install_fn)
    std_install = """
    std.pkgs <- c(%s)
    std.installer = repo.installer(cran.repos, install.packages)
    lapply(std.pkgs, std.installer)
    """ % (", ".join('"%s"' % p for p in config['cran']))
    append(out_file, std_install)
    bioc_install = """
    bioc.pkgs <- c(%s)
    bioc.installer = repo.installer(biocinstallRepos(), biocLite)
    lapply(bioc.pkgs, bioc.installer)
    """ % (", ".join('"%s"' % p for p in config['bioc']))
    append(out_file, bioc_install)
    final_update = """
    update.packages(repos=biocinstallRepos(), ask=FALSE)
    update.packages(ask=FALSE)
    """
    append(out_file, final_update)
    # run the script and then get rid of it
    env.safe_sudo("Rscript %s" % out_file)
    run("rm -f %s" % out_file)

def _python_library_installer(config):
    """Install python specific libraries using easy_install.
    """
    version_ext = "-%s" % env.python_version_ext if env.python_version_ext else ""
    env.safe_sudo("easy_install%s -U pip" % version_ext)
    for pname in env.flavor.rewrite_config_items("python", config['pypi']):
        env.safe_sudo("easy_install%s -U %s" % (version_ext, pname))
        # Use pip when it doesn't re-download even if latest package installed
        # https://bitbucket.org/ianb/pip/issue/13/upgrade-always-downloads-most-recent
        #sudo("pip%s install -U %s" % (version_ext,  pname))

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

    # No need to prevent TOCTTOU, nothing critical is going to be touched
    if not os.path.isfile("%s/bin/cpanm" % env.system_install):
        with _make_tmp_dir() as tmp_dir:
            with cd(tmp_dir):
                cpanm_header = ''
                while cpanm_header.find('perl') == -1:
                    run("wget --no-check-certificate http://xrl.us/cpanm -O cpanm")
                    cpanm_file = open("%s/bin/cpanm" % env.system_install)
                    cpanm_header = cpanm_file.readline()

                run("chmod a+rwx cpanm")
                env.safe_sudo("mv cpanm %s/bin" % env.system_install)

    # TODO: Check cpanm file, making sure it is a legitimate perl file, not an HTML page.
    #    cpanm_file = open("%s/bin/cpanm" % env.system_install, 'r')
    #    cpanm_header = cpanm_file.readline()
    #    if cpanm_header.find('perl') == -1:
        # Retry to download the file

    sudo_str = "--sudo" if env.use_sudo else ""
    for lib in env.flavor.rewrite_config_items("perl", config['cpan']):
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        run("cpanm %s --skip-installed --notest %s < /dev/null" % (sudo_str, lib))

def _clojure_library_installer(config):
    """Install clojure libraries using cljr.
    """
    for lib in config['cljr']:
        run("cljr install %s" % lib)

def _haskell_library_installer(config):
    """Install haskell libraries using cabal.
    """
    run("cabal update")
    for lib in config["cabal"]:
        sudo_str = "--root-cmd=sudo" if env.use_sudo else ""
        run("cabal install %s --global %s" % (sudo_str, lib))

lib_installers = {
    "r-libs" : _r_library_installer,
    "python-libs" : _python_library_installer,
    "ruby-libs" : _ruby_library_installer,
    "perl-libs" : _perl_library_installer,
    "clojure-libs": _clojure_library_installer,
    "haskell-libs": _haskell_library_installer,
    }

def install_libraries(language):
    """High level target to install libraries for a specific language.
    """
    _setup_logging(env)
    _check_fabric_version()
    _parse_fabricrc()
    _setup_edition(env)
    _setup_flavor(None)
    _setup_distribution_environment()
    _create_local_paths()
    _do_library_installs(["%s-libs" % language])

def _do_library_installs(to_install):
    for iname in to_install:
        yaml_file = os.path.join(env.config_dir, "%s.yaml" % iname)
        with open(yaml_file) as in_handle:
            config = yaml.load(in_handle)
        lib_installers[iname](config)
