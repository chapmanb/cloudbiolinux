"""Configuration details for specific server types.

This module contains functions that help with initializing a Fabric environment
for standard server types.
"""
import os
import subprocess

from fabric.api import env

from cloudbio.fabutils import quiet
from cloudbio.fabutils import configure_runsudo
from cloudbio.custom import system

def _setup_distribution_environment(ignore_distcheck=False):
    """Setup distribution environment.

    In low-level terms, this method attempts to populate various values in the fabric
    env data structure for use other places in CloudBioLinux.
    """
    if "distribution" not in env:
        env.distribution = "__auto__"
    if "dist_name" not in env:
        env.dist_name = "__auto__"
    env.logger.info("Distribution %s" % env.distribution)

    if env.hosts == ["vagrant"]:
        _setup_vagrant_environment()
    elif env.hosts == ["localhost"]:
        _setup_local_environment()
    configure_runsudo(env)
    if env.distribution == "__auto__":
        env.distribution = _determine_distribution(env)
    if env.distribution == "ubuntu":
        ## TODO: Determine if dist_name check works with debian.
        if env.dist_name == "__auto__":
            env.dist_name = _ubuntu_dist_name(env)
        _setup_ubuntu()
    elif env.distribution == "centos":
        _setup_centos()
    elif env.distribution == "scientificlinux":
        _setup_scientificlinux()
    elif env.distribution == "debian":
        if env.dist_name == "__auto__":
            env.dist_name = _debian_dist_name(env)
        _setup_debian()
    elif env.distribution == "arch":
        pass  # No package support for Arch yet
    elif env.distribution == "macosx":
        _setup_macosx(env)
        ignore_distcheck = True
    else:
        raise ValueError("Unexpected distribution %s" % env.distribution)
    if not ignore_distcheck:
        _validate_target_distribution(env.distribution, env.get('dist_name', None))
    _cloudman_compatibility(env)
    _setup_nixpkgs()
    _setup_fullpaths(env)
    # allow us to check for packages only available on 64bit machines
    machine = env.safe_run_output("uname -m")
    env.is_64bit = machine.find("_64") > 0


def _setup_fullpaths(env):
    home_dir = env.safe_run_output("echo $HOME")
    for attr in ["data_files", "galaxy_home", "local_install"]:
        if hasattr(env, attr):
            x = getattr(env, attr)
            if x.startswith("~"):
                x = x.replace("~", home_dir)
                setattr(env, attr, x)


def _cloudman_compatibility(env):
    """Environmental variable naming for compatibility with CloudMan.
    """
    env.install_dir = env.system_install


def _validate_target_distribution(dist, dist_name=None):
    """Check target matches environment setting (for sanity)

    Throws exception on error
    """
    env.logger.debug("Checking target distribution " + env.distribution)
    if dist in ["debian", "ubuntu"]:
        tag = env.safe_run_output("cat /proc/version")
        if tag.lower().find(dist) == -1:
           # hmmm, test issue file
            tag2 = env.safe_run_output("cat /etc/issue")
            if tag2.lower().find(dist) == -1:
                raise ValueError("Distribution does not match machine; are you using correct fabconfig for " + dist)
        if env.edition.short_name in ["minimal"]:
            # "minimal editions don't actually change any of the apt
            # source except adding biolinux, so won't cause this
            # problem and don't need to match dist_name"
            return
        if not dist_name:
            raise ValueError("Must specify a dist_name property when working with distribution %s" % dist)
        # Does this new method work with CentOS, do we need this.
        actual_dist_name = _ubuntu_dist_name(env)
        if actual_dist_name != dist_name:
            raise ValueError("Distribution does not match machine; are you using correct fabconfig for " + dist)
    else:
        env.logger.debug("Unknown target distro")


def _setup_ubuntu():
    env.logger.info("Ubuntu setup")
    shared_sources = _setup_deb_general()
    # package information. This is ubuntu/debian based and could be generalized.
    sources = [
      "deb http://us.archive.ubuntu.com/ubuntu/ %s universe",  # unsupported repos
      "deb http://us.archive.ubuntu.com/ubuntu/ %s multiverse",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates universe",
      "deb http://us.archive.ubuntu.com/ubuntu/ %s-updates multiverse",
      "deb http://archive.canonical.com/ubuntu %s partner",  # partner repositories
      "deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen",  # mongodb
      "deb http://cran.fhcrc.org/bin/linux/ubuntu %s/",  # lastest R versions
      "deb http://archive.cloudera.com/debian maverick-cdh3 contrib",  # Hadoop
      "deb http://archive.canonical.com/ubuntu %s partner",  # sun-java
      "deb http://ppa.launchpad.net/freenx-team/ppa/ubuntu precise main",  # Free-NX
      "deb http://ppa.launchpad.net/nebc/bio-linux/ubuntu precise main",  # Free-NX
      "deb [arch=amd64 trusted=yes] http://research.cs.wisc.edu/htcondor/debian/stable/ squeeze contrib"  # HTCondor
    ] + shared_sources
    env.std_sources = _add_source_versions(env.dist_name, sources)


def _setup_debian():
    env.logger.info("Debian setup")
    unstable_remap = {"sid": "squeeze"}
    shared_sources = _setup_deb_general()
    sources = [
        "deb http://downloads-distro.mongodb.org/repo/debian-sysvinit dist 10gen",  # mongodb
        "deb http://cran.fhcrc.org/bin/linux/debian %s-cran/",  # lastest R versions
        "deb http://archive.cloudera.com/debian lenny-cdh3 contrib"  # Hadoop
        ] + shared_sources
    # fill in %s
    dist_name = unstable_remap.get(env.dist_name, env.dist_name)
    env.std_sources = _add_source_versions(dist_name, sources)


def _setup_deb_general():
    """Shared settings for different debian based/derived distributions.
    """
    env.logger.debug("Debian-shared setup")
    env.sources_file = "/etc/apt/sources.list.d/cloudbiolinux.list"
    env.global_sources_file = "/etc/apt/sources.list"
    env.apt_preferences_file = "/etc/apt/preferences"
    if not hasattr(env, "python_version_ext"):
        env.python_version_ext = ""
    if not hasattr(env, "ruby_version_ext"):
        env.ruby_version_ext = "1.9.1"
    if not env.has_key("java_home"):
        # Try to determine java location from update-alternatives
        java_home = "/usr/lib/jvm/java-7-openjdk-amd64"
        with quiet():
            java_info = env.safe_run_output("update-alternatives --display java")
        for line in java_info.split("\n"):
            if line.strip().startswith("link currently points to"):
                java_home = line.split()[-1].strip()
                java_home = java_home.replace("/jre/bin/java", "")
        env.java_home = java_home
    shared_sources = [
        "deb http://nebc.nerc.ac.uk/bio-linux/ unstable bio-linux",  # Bio-Linux
        "deb http://download.virtualbox.org/virtualbox/debian %s contrib",  # virtualbox
    ]
    return shared_sources


def _setup_centos():
    env.logger.info("CentOS setup")
    if not hasattr(env, "python_version_ext"):
        # use installed anaconda version instead of package 2.6
        #env.python_version_ext = "2.6"
        env.python_version_ext = ""
    #env.pip_cmd = "pip-python"
    if not hasattr(env, "ruby_version_ext"):
        env.ruby_version_ext = ""
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"


def _setup_scientificlinux():
    env.logger.info("ScientificLinux setup")
    if not hasattr(env, "python_version_ext"):
        env.python_version_ext = ""
    env.pip_cmd = "pip-python"
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"

def _setup_macosx(env):
    # XXX Only framework in place; needs testing
    env.logger.info("MacOSX setup")
    # XXX Ensure XCode is installed and provide useful directions if not
    system.install_homebrew(env)
    # XXX find java correctly
    env.java_home = ""

def _setup_nixpkgs():
    # for now, Nix packages are only supported in Debian - it can
    # easily be done for others - just get Nix installed from the .rpm
    nixpkgs = False
    if env.has_key("nixpkgs"):
        if env.distribution in ["debian", "ubuntu"]:
            if env.nixpkgs == "True":
                nixpkgs = True
            else:
                nixpkgs = False
        else:
            env.logger.warn("NixPkgs are currently not supported for " + env.distribution)
    if nixpkgs:
        env.logger.info("NixPkgs: supported")
    else:
        env.logger.debug("NixPkgs: Ignored")
    env.nixpkgs = nixpkgs


def _setup_local_environment():
    """Setup a localhost environment based on system variables.
    """
    env.logger.info("Get local environment")
    if not env.has_key("user"):
        env.user = os.environ["USER"]


def _setup_vagrant_environment():
    """Use vagrant commands to get connection information.
    https://gist.github.com/1d4f7c3e98efdf860b7e
    """
    env.logger.info("Get vagrant environment")
    raw_ssh_config = subprocess.Popen(["vagrant", "ssh-config"],
                                      stdout=subprocess.PIPE).communicate()[0]
    env.logger.info(raw_ssh_config)
    ssh_config = dict([l.strip().split() for l in raw_ssh_config.split("\n") if l])
    env.user = ssh_config["User"]
    env.hosts = [ssh_config["HostName"]]
    env.port = ssh_config["Port"]
    env.host_string = "%s@%s:%s" % (env.user, env.hosts[0], env.port)
    env.key_filename = ssh_config["IdentityFile"].replace('"', '')
    env.logger.debug("ssh %s" % env.host_string)


def _add_source_versions(version, sources):
    """Patch package source strings for version, e.g. Debian 'stable'
    """
    name = version
    env.logger.debug("Source=%s" % name)
    final = []
    for s in sources:
        if s.find("%s") > 0:
            s = s % name
        final.append(s)
    return final


def _ubuntu_dist_name(env):
    """
    Determine Ubuntu dist name (e.g. precise or quantal).
    """
    return env.safe_run_output("cat /etc/*release | grep DISTRIB_CODENAME | cut -f 2 -d =")


def _debian_dist_name(env):
    """
    Determine Debian dist name (e.g. squeeze).
    """
    return env.safe_run_output("lsb_release -a | grep Codename | cut -f 2")


def _determine_distribution(env):
    """
    Attempt to automatically determine the distribution of the target machine.

    Currently works for Ubuntu, CentOS, Debian, Scientific Linux and Mac OS X.
    """
    with quiet():
        output = env.safe_run_output("cat /etc/*release").lower()
    if output.find("distrib_id=ubuntu") >= 0:
        return "ubuntu"
    elif output.find("centos release") >= 0:
        return "centos"
    elif output.find("red hat enterprise linux server release") >= 0:
        return "centos"
    elif output.find("fedora release") >= 0:
        return "centos"
    elif output.find("scientific linux release") >= 0:
        return "scientificlinux"
    elif env.safe_exists("/etc/debian_version"):
        return "debian"
    elif output.find("id=arch"):
        return "arch"
    # check for file used by Python's platform.mac_ver
    elif env.safe_exists("/System/Library/CoreServices/SystemVersion.plist"):
        return "macosx"
    else:
        raise Exception("Attempt to automatically determine Linux distribution of target machine failed, please manually specify distribution in fabricrc.txt")
