"""Main Fabric deployment file for BioLinux distribution.

This installs a standard set of useful biological applications on a remote
server. It is designed for bootstrapping a machine from scratch, as with new
Amazon EC2 instances.

Usage:
    fab -H hostname -i private_key_file install_biolinux

which will call into the 'install_biolinux' method below.

Requires:
    Fabric http://docs.fabfile.org
    PyYAML http://pyyaml.org/wiki/PyYAMLDocumentation
"""
import os
import sys
import subprocess

from fabric.main import load_settings
from fabric.api import *
from fabric.contrib.files import *
import yaml
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
# add formatter to ch
ch.setFormatter(formatter)

logger.addHandler(ch)

env.config_dir = os.path.join(os.path.dirname(__file__), "config")

# ### Configuration details for different server types

def _setup_distribution_environment():
    """Setup distribution environment
    """
    _parse_fabricrc()
    logger.info("distribution=%s" % env.distribution)
    if env.hosts == ["vagrant"]:
        _setup_vagrant_environment()
    elif env.hosts == ["localhost"]:
        _setup_local_environment()
    if env.distribution == "ubuntu":
        _setup_ubuntu()
    elif env.distribution == "centos":
        _setup_centos()
    elif env.distribution == "debian":
        _setup_debian()
    else:
        raise ValueError("Unexpected distribution %s" % env.distribution)
    _expand_shell_paths()

def _validate_target_distribution():
    """Check target matches environment setting (for sanity)

    Throws exception on error
    """
    logger.debug(env.distribution)
    if env.distribution == "debian":
        tag = run("cat /proc/version")
        if tag.find('ebian') == -1:
           raise ValueError("Debian does not match target, are you using correct fabconfig?")
    if env.distribution == "ubuntu":
        tag = run("cat /proc/version")
        if tag.find('buntu') == -1:
           raise ValueError("Ubuntu does not match target, are you using correct fabconfig?")
    else:
        logger.debug("Unknown target distro")

def _setup_ubuntu():
    logger.info("Ubuntu setup")
    shared_sources = _setup_deb_general()
    # package information. This is ubuntu/debian based and could be generalized.
    version = env.dist_name
    sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb-src http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen", # mongodb
      "deb http://cran.stat.ucla.edu/bin/linux/ubuntu %s/", # lastest R versions
      "deb http://archive.cloudera.com/debian %s-cdh3 contrib", # Hadoop
      "ppa:sun-java-community-team/sun-java6", # sun-java
    ] + shared_sources
    env.std_sources = _add_source_versions(version, sources)

def _setup_debian():
    logger.info("Debian setup")
    shared_sources = _setup_deb_general()
    version = env.dist_name
    if not env.get('debian_repository'):
      main_repository = 'http://ftp.us.debian.org/debian/'
    else:
      main_repository = env.debian_repository

    sources = [
      "deb {repo} %s main contrib non-free".format(repo=main_repository),
      "deb {repo} %s-updates main contrib non-free".format(repo=main_repository),
      "deb http://downloads-distro.mongodb.org/repo/debian-sysvinit dist 10gen", # mongodb
      "deb http://cran.stat.ucla.edu/bin/linux/debian %s-cran/", # lastest R versions
      "deb http://archive.cloudera.com/debian lenny-cdh3 contrib", # Hadoop
    ] + shared_sources
    env.std_sources = _add_source_versions(version, sources)

def _setup_deb_general():
    """Shared settings for different debian based/derived distributions.
    """
    logger.debug("Debian-shared setup")
    env.sources_file = "/etc/apt/sources.list"
    env.python_version_ext = ""
    if not env.has_key("java_home"):
        # XXX look for a way to find JAVA_HOME automatically
        env.java_home = "/usr/lib/jvm/java-6-openjdk"
    shared_sources = [
      "deb http://nebc.nox.ac.uk/bio-linux/ unstable bio-linux", # Bio-Linux
      "deb http://ppa.launchpad.net/freenx-team/ppa/ubuntu lucid main", # FreeNX PPA
      "deb http://download.virtualbox.org/virtualbox/debian %s contrib", # virtualbox
    ]
    return shared_sources

def _setup_centos():
    logger.info("CentOS setup")
    env.python_version_ext = "2.6"
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"

def _parse_fabricrc():
    """Defaults from fabricrc.txt file; loaded if not specified at commandline.
    """
    if not env.has_key("distribution"):
        logger.info("Reading default fabricrc.txt")
        config_file = os.path.join(env.config_dir, "fabricrc.txt")
        if os.path.exists(config_file):
            env.update(load_settings(config_file))
    else:
        logger.warn("Skipping fabricrc.txt as distribution is already defined")

def _expand_shell_paths():
    """Expand any paths defined in terms of shell shortcuts (like ~).
    """
    logger.debug("Expand paths")
    if env.has_key("local_install"):
        if not exists(env.local_install):
            run("mkdir -p %s" % env.local_install)
        with cd(env.local_install):
            with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                          warn_only=True):
                result = run("pwd")
                env.local_install = result

def _setup_local_environment():
    """Setup a localhost environment based on system variables.
    """
    logger.info("Get local environment")
    if not env.has_key("user"):
        env.user = os.environ["USER"]
    if not env.has_key("java_home"):
        env.java_home = os.environ.get("JAVA_HOME", "/usr/lib/jvm/java-6-openjdk")

def _setup_vagrant_environment():
    """Use vagrant commands to get connection information.
    https://gist.github.com/1d4f7c3e98efdf860b7e
    """
    logger.info("Get vagrant environment")
    raw_ssh_config = subprocess.Popen(["vagrant", "ssh-config"],
                                      stdout=subprocess.PIPE).communicate()[0]
    ssh_config = dict([l.strip().split() for l in raw_ssh_config.split("\n") if l])
    env.user = ssh_config["User"]
    env.hosts = [ssh_config["HostName"]]
    env.port = ssh_config["Port"]
    env.host_string = "%s@%s:%s" % (env.user, env.hosts[0], env.port)
    env.key_filename = ssh_config["IdentityFile"]
    logger.debug("ssh %s" % env.host_string)

def _add_source_versions(version, sources):
    """Patch package source strings for version, e.g. Debian 'stable'
    """
    name = version
    logger.debug("Set source=%s" % name)
    final = []
    for s in sources:
        if s.find("%s") > 0:
            s = s % name
        final.append(s)
    return final

# ### Shared installation targets for all platforms

def install_biolinux(target=None):
    """Main entry point for installing Biolinux on a remote server.

    target is supplied on the fab CLI. Special targets are:

      - packages     Install distro packages (default)
      - custom
      - libraries
      - finalize     Setup freenx
    """
    _check_fabric_version()
    _setup_distribution_environment() # get parameters for distro, packages etc.
    pkg_install, lib_install = _read_main_config()  # read main.yaml
    _validate_target_distribution()
    logger.info("Target=%s" % target)
    if target is None or target == "packages":
        if env.distribution in ["ubuntu", "debian"]:
            _setup_apt_sources()
            _setup_apt_automation()
            _add_apt_gpg_keys()
            _apt_packages(pkg_install)
        elif env.distribution in ["centos"]:
            _setup_yum_sources()
            _yum_packages(pkg_install)
            _setup_yum_bashrc()
    if target is None or target == "custom":
        _custom_installs(pkg_install)
    if target is None or target == "libraries":
        _do_library_installs(lib_install)
    if target is None or target == "finalize":
        _freenx_scripts()
        _cleanup()

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
    sys.path.append(os.path.split(__file__)[0])
    for p in packages:
        install_custom(p, True, pkg_to_group)

def install_custom(p, automated=False, pkg_to_group=None):
    """Install a single custom package by name.

    This method fetches names from custom.yaml that delegate to a method
    in the custom/name.py program.

    fab install_custom_package:package_name
    """
    logger.info("Install custom software packages")
    if not automated:
        if not env.has_key("system_install"):
            _parse_fabricrc()
        pkg_config = os.path.join(env.config_dir, "custom.yaml")
        packages, pkg_to_group = _yaml_to_packages(pkg_config, None)
        sys.path.append(os.path.split(__file__)[0])
    try:
        logger.debug("Import %s" % p)
        mod = __import__("custom.%s" % pkg_to_group[p], fromlist=["custom"])
    except ImportError:
        raise ImportError("Need to write a %s module in custom." %
                pkg_to_group[p])
    try:
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
    logger.info("Reading %s" % yaml_file)
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

def _read_main_config():
    """Pull a list of groups to install based on our main configuration YAML.

    Reads 'main.yaml' and returns packages and libraries
    """
    yaml_file = os.path.join(env.config_dir, "main.yaml")
    with open(yaml_file) as in_handle:
        full_data = yaml.load(in_handle)
    packages = full_data['packages']
    packages = packages if packages else []
    libraries = full_data['libraries']
    libraries = libraries if libraries else []
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
    for pname in config['pypi']:
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
    for gem in config['gems']:
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
    for lib in config['cpan']:
        # Need to hack stdin because of some problem with cpanminus script that
        # causes fabric to hang
        # http://agiletesting.blogspot.com/2010/03/getting-past-hung-remote-processes-in.html
        run("cpanm --sudo --skip-installed --notest %s < /dev/null" % (lib))

def _clojure_library_installer(config):
    """Install clojure libraries using cljr.
    """
    for lib in config['cljr']:
        run("cljr install %s" % lib)

lib_installers = {
        "r-libs" : _r_library_installer,
        "python-libs" : _python_library_installer,
        "ruby-libs" : _ruby_library_installer,
        "perl-libs" : _perl_library_installer,
        "clojure-libs": _clojure_library_installer,
        }

def install_libraries(language):
    """High level target to install libraries for a specific language.
    """
    _check_fabric_version()
    _setup_distribution_environment()
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
    logger.info("Update and install all packages")
    pkg_config = os.path.join(env.config_dir, "packages.yaml")
    subs_pkg_config = os.path.join(env.config_dir, "packages-%s.yaml" %
                                   env.distribution)
    if not os.path.exists(subs_pkg_config): subs_pkg_config = None
    sudo("apt-get update")
    sudo("apt-get -y --force-yes upgrade")
    # Retrieve packages to get and install each of them
    (packages, _) = _yaml_to_packages(pkg_config, to_install,
                                      subs_pkg_config)
    # for package in packages:
    #     sudo("apt-get -y --force-yes install %s" % package)
    # A single line install is much faster - note that there is a max
    # for the command line size, so we do 30 at a time
    i = 0
    while i < len(packages):
      list = packages[i:i+30]
      sudo("apt-get -y --force-yes install %s" % " ".join(list))
      i += 30
    sudo("apt-get clean")

def _add_apt_gpg_keys():
    """Adds GPG keys from all repositories
    """
    logger.info("Update GPG keys for repositories")
    standalone = [
        "http://archive.cloudera.com/debian/archive.key",
        "http://download.virtualbox.org/virtualbox/debian/oracle_vbox.asc"]
    keyserver = [
        ("keyserver.ubuntu.com", "7F0CEB10"),
        ("keyserver.ubuntu.com", "E084DAB9"),
        ("keyserver.ubuntu.com", "D67FC6EAE2A11821"),
    ]
    for url, key in keyserver:
        sudo("apt-key adv --keyserver %s --recv %s" % (url, key))
    for key in standalone:
        sudo("wget -q -O- %s | apt-key add -" % key)

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
    cmd = ""
    for l in package_info:
        #     sudo("echo %s | /usr/bin/debconf-set-selections" % l)
        cmd += "echo %s | /usr/bin/debconf-set-selections ; " % l

    sudo(cmd)

def _setup_apt_sources():
    """Add sources for retrieving library packages.
       Using add-apt-repository allows processing PPAs (on Ubuntu)

       This method modifies the apt sources file.

       Uses python-software-properties, which provides an abstraction of apt repositories
    """
    sudo("apt-get install -y --force-yes python-software-properties")
    comment = "# This file was modified for BioLinux"
    if not contains(env.sources_file, comment):
        append(env.sources_file, comment, use_sudo=True)
    for source in env.std_sources:
        logger.debug("Source %s" % source)
        if source.startswith("ppa:"):
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

def _cleanup():
    """Clean up any extra files after building.
    """
    logger.info("Cleaning up")
    run("rm -f .bash_history")
    sudo("rm -f /var/crash/*")
    sudo("rm -f /var/log/firstboot.done")
    sudo("rm -f .nx_setup_done")
    sudo("rm -rf .cpanm")
    # RabbitMQ fails to start if its database is embedded into the image
    # because it saves the current IP address or host name so delete it now.
    # When starting up, RabbitMQ will recreate that directory.
    sudo('/etc/init.d/rabbitmq-server stop')
    for db_location in ['/var/lib/rabbitmq/mnesia', '/mnesia']:
        if exists(db_location):
            sudo('rm -rf %s' % db_location)
