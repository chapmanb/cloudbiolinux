"""Main Fabric deployment file for BioLinux distributionrepo=.

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
from cloudbio.utils import _setup_logging
from cloudbio.cloudman import (_configure_ec2_autorun, _cleanup_ec2)

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

def _setup_flavor(flavor):
    """Setup flavor
    """
    if flavor == None:
        flavor = env.get("flavor", None)
        flavor_path = env.get("flavor_path", None)
    if flavor != None:
        # Add path for flavors
        sys.path.append(os.path.join(os.path.dirname(__file__), "contrib", "flavor"))
        env.logger.info("Flavor %s loaded from %s" % (flavor, flavor_path))
        try:
            mod = __import__(flavor_path, fromlist=[flavor])
        except ImportError:
            raise ImportError("Failed to import %s" % flavor)
    else:
        from cloudbio.flavor import Flavor
    env.logger.info("This is a %s" % env.flavor.name)

# ### Shared installation targets for all platforms

def install_biolinux(target=None, packagelist=None, flavor=None):
    """Main entry point for installing Biolinux on a remote server.

    This allows a different main package list (the main YAML file is passed in),
    and/or use of Flavor. So you can say:

      install_bare:packagelist=contrib/mylist/main.yaml,flavor=specialflavor

    Both packagelist and flavor, as well as the Edition, can also be passed in
    through the fabricrc file.

    target can also be supplied on the fab CLI. Special targets are:

      - packages     Install distro packages
      - custom       Install custom packages
      - libraries    Install programming language libraries
      - finalize     Setup freenx and clean-up environment
    """
    _setup_logging(env)
    _check_fabric_version()
    _parse_fabricrc()
    _setup_edition(env)
    _setup_flavor(flavor)
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
    if target is None or target == "custom":
        _custom_installs(pkg_install)
    if target is None or target == "libraries":
        _do_library_installs(lib_install)
    if target is None or target == "finalize":
        env.edition.post_install()
        env.flavor.post_install()
        _cleanup_space()
        if env.has_key("is_ec2_image") and env.is_ec2_image.upper() in ["TRUE", "YES"]:
            _freenx_scripts()
            _configure_ec2_autorun(env)
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
    for p in env.flavor.rewrite_custom_list(packages):
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
        print env
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

def _yaml_to_packages(yaml_file, to_install, subs_yaml_file = None):
    """Read a list of packages from a nested YAML configuration file.
    """
    # allow us to check for packages only available on 64bit machines
    machine = run("uname -m")
    is_64bit = machine.find("_64") > 0
    env.logger.info("Reading %s" % yaml_file)
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    if subs_yaml_file is not None:
        with open(subs_yaml_file) as in_handle:
            subs = yaml.load(in_handle)
    else:
        subs = {}
    # filter the data based on what we have configured to install
    data = [(k, v) for (k,v) in full_data.iteritems()
            if to_install is None or k in to_install]
    data.sort()
    packages = []
    pkg_to_group = dict()
    while len(data) > 0:
        cur_key, cur_info = data.pop(0)
        if cur_info:
            if isinstance(cur_info, (list, tuple)):
                packages.extend(_filter_subs_packages(cur_info, subs))
                for p in cur_info:
                    pkg_to_group[p] = cur_key
            elif isinstance(cur_info, dict):
                for key, val in cur_info.iteritems():
                    # if we are okay, propagate with the top level key
                    if key != 'needs_64bit' or is_64bit:
                        data.append((cur_key, val))
            else:
                raise ValueError(cur_info)
    env.logger.debug("Packages:")
    env.logger.debug(",".join(packages))
    return packages, pkg_to_group

def _filter_subs_packages(initial, subs):
    """Rename and filter package list with subsitutions; for similar systems.
    """
    final = []
    for p in initial:
        try:
            new_p = subs[p]
        except KeyError:
            new_p = p
        if new_p:
            final.append(new_p)
    return sorted(final)

def _read_main_config(yaml_file=None):
    """Pull a list of groups to install based on our main configuration YAML.

    Reads 'main.yaml' and returns packages and libraries
    """
    if yaml_file==None:
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
    sudo("Rscript %s" % out_file)
    run("rm -f %s" % out_file)

def _python_library_installer(config):
    """Install python specific libraries using easy_install.
    """
    version_ext = "-%s" % env.python_version_ext if env.python_version_ext else ""
    sudo("easy_install%s -U pip" % version_ext)
    for pname in env.flavor.rewrite_python_egg_list(config['pypi']):
        sudo("easy_install%s -U %s" % (version_ext, pname))
        # Use pip when it doesn't re-download even if latest package installed
        # https://bitbucket.org/ianb/pip/issue/13/upgrade-always-downloads-most-recent
        #sudo("pip%s install -U %s" % (version_ext,  pname))

def _ruby_library_installer(config):
    """Install ruby specific gems.
    """
    def _cur_gems():
        with settings(
                hide('warnings', 'running', 'stdout', 'stderr')):
            gem_info = run("gem list --no-versions")
        return [l.rstrip("\r") for l in gem_info.split("\n") if l.rstrip("\r")]
    installed = _cur_gems()
    for gem in env.flavor.rewrite_ruby_gem_list(config['gems']):
        # update current gems only to check for new installs
        if gem not in installed:
            installed = _cur_gems()
        if gem in installed:
            sudo("gem update %s" % gem)
        else:
            sudo("gem install %s" % gem)

def _perl_library_installer(config):
    """Install perl libraries from CPAN with cpanminus.
    """
    run("wget --no-check-certificate http://xrl.us/cpanm")
    run("chmod a+rwx cpanm")
    sudo("mv cpanm %s/bin" % env.system_install)
    for lib in env.flavor.rewrite_perl_cpan_list(config['cpan']):
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        run("cpanm --sudo --skip-installed --notest %s < /dev/null" % (lib))

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
        run("cabal install --root-cmd=sudo --global %s" % lib)

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

# ### Automated installation on apt systems

def _apt_packages(to_install):
    """Install packages available via apt-get.
    """
    env.logger.info("Update and install all packages")
    pkg_config_file = os.path.join(env.config_dir, "packages.yaml")
    subs_pkg_config_file = os.path.join(env.config_dir, "packages-%s.yaml" %
                                   env.distribution)
    if not os.path.exists(subs_pkg_config_file): subs_pkg_config_file = None
    sudo("apt-get update") # Always update
    env.edition.apt_upgrade_system()
    # Retrieve final package names
    (packages, _) = _yaml_to_packages(pkg_config_file, to_install,
                                      subs_pkg_config_file)
    # At this point allow the Flavor to rewrite the package list
    packages = env.flavor.rewrite_packages_list(packages)

    # A single line install is much faster - note that there is a max
    # for the command line size, so we do 30 at a time
    group_size = 30
    i = 0
    env.logger.info("Updating %i packages" % len(packages))
    while i < len(packages):
        sudo("apt-get -y --force-yes install %s" % " ".join(packages[i:i+group_size]))
        i += group_size
    sudo("apt-get clean")

def _add_apt_gpg_keys():
    """Adds GPG keys from all repositories
    """
    env.logger.info("Update GPG keys for repositories")
    standalone = [
        "http://archive.cloudera.com/debian/archive.key",
        'http://download.virtualbox.org/virtualbox/debian/oracle_vbox.asc'
    ]
    keyserver = [
            ("keyserver.ubuntu.com", "7F0CEB10"),
            ("keyserver.ubuntu.com", "E084DAB9"),
            ("keyserver.ubuntu.com", "D67FC6EAE2A11821"),
        ]
    standalone, keyserver = env.edition.rewrite_apt_keys(standalone, keyserver)
    for key in standalone:
        sudo("wget -q -O- %s | apt-key add -" % key)
    for url, key in keyserver:
        sudo("apt-key adv --keyserver %s --recv %s" % (url, key))

def _setup_apt_automation():
    """Setup the environment to be fully automated for tricky installs.

    Sun Java license acceptance:
    http://www.davidpashley.com/blog/debian/java-license

    MySQL root password questions; install with empty root password:
    http://snowulf.com/archives/540-Truly-non-interactive-unattended-apt-get-install.html

    Postfix, setup for no configuration. See more on issues here:
    http://www.uluga.ubuntuforums.org/showthread.php?p=9120196
    """
    interactive_cmd = "export DEBIAN_FRONTEND=noninteractive"
    if not contains(env.shell_config, interactive_cmd):
        append(env.shell_config, interactive_cmd)
    package_info = [
            "postfix postfix/main_mailer_type select No configuration",
            "postfix postfix/mailname string notusedexample.org",
            "mysql-server-5.1 mysql-server/root_password string '(password omitted)'",
            "mysql-server-5.1 mysql-server/root_password_again string '(password omitted)'",
            "sun-java6-jdk shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-jre shared/accepted-sun-dlj-v1-1 select true",
            "sun-java6-bin shared/accepted-sun-dlj-v1-1 select true",
            "grub-pc grub2/linux_cmdline string ''",
            "grub-pc grub-pc/install_devices_empty boolean true",
            "acroread acroread/default-viewer boolean false",
            ]
    package_info = env.edition.rewrite_apt_automation(package_info)
    cmd = ""
    for l in package_info:
        cmd += "echo %s | /usr/bin/debconf-set-selections ; " % l
    sudo(cmd)

def _setup_apt_sources():
    """Add sources for retrieving library packages.
       Using add-apt-repository allows processing PPAs (on Ubuntu)

       This method modifies the apt sources file.

       Uses python-software-properties, which provides an abstraction of apt repositories
    """
    env.edition.check_packages_source()

    comment = "# This file was modified for "+ env.edition.name
    if not exists(env.sources_file):
        sudo("touch %s" % env.sources_file)
    if not contains(env.sources_file, comment):
        append(env.sources_file, comment, use_sudo=True)
    for source in env.edition.rewrite_apt_sources_list(env.std_sources):
        env.logger.debug("Source %s" % source)
        if source.startswith("ppa:"):
            sudo("apt-get install -y --force-yes python-software-properties")
            sudo("add-apt-repository '%s'" % source)
        elif not contains(env.sources_file, source):
            append(env.sources_file, source, use_sudo=True)

# ### Automated installation on systems with yum package manager

def _yum_packages(to_install):
    """Install rpm packages available via yum.
    """
    pkg_config = os.path.join(env.config_dir, "packages-yum.yaml")
    with settings(warn_only=True):
        sudo("yum check-update")
    sudo("yum -y upgrade")
    # Retrieve packages to get and install each of them
    (packages, _) = _yaml_to_packages(pkg_config, to_install)
    # At this point allow the Flavor to rewrite the package list
    packages = env.flavor.rewrite_packages_list(packages)
    for package in packages:
        sudo("yum -y install %s" % package)

def _setup_yum_bashrc():
    """Fix the user bashrc to update compilers.
    """
    to_include = ["export CC=gcc44", "export CXX=g++44", "export FC=gfortran44",
                  "export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/usr/lib/pkgconfig"]
    fname = run("ls %s" % env.shell_config)
    for line in to_include:
        if not contains(fname, line.split("=")[0]):
            append(fname, line)

def _setup_yum_sources():
    """Add additional useful yum repositories.
    """
    repos = ["http://download.fedora.redhat.com/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm"]
    for repo in repos:
        with settings(warn_only=True):
            sudo("rpm -Uvh %s" % repo)

# ### CloudBioLinux specific scripts

def _freenx_scripts():
    """Provide graphical access to clients via FreeNX.
    """
    setup_script = "setupnx.sh"
    remote_setup = "%s/bin/%s" % (env.system_install, setup_script)
    install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
    if not exists(remote_setup):
        put(os.path.join(install_file_dir, setup_script), setup_script,
                mode=0777)
        sudo("mv %s %s" % (setup_script, remote_setup))
    remote_login = "configure_freenx.sh"
    if not exists(remote_login):
        put(os.path.join(install_file_dir, 'bash_login'), remote_login,
                mode=0777)

def _cleanup_space():
    """Cleanup to recover space from builds and packages.
    """
    env.logger.info("Cleaning up space from package builds")
    sudo("rm -rf .cpanm")
    sudo("rm -f /var/crash/*")

