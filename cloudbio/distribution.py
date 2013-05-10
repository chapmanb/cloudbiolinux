"""Configuration details for specific server types.

This module contains functions that help with initializing a Fabric environment
for standard server types.
"""
import os
import subprocess

from fabric.api import env

from cloudbio.fabutils import configure_runsudo

def _setup_distribution_environment(ignore_distcheck=False):
    """Setup distribution environment
    """
    env.logger.info("Distribution %s" % env.distribution)

    if env.hosts == ["vagrant"]:
        _setup_vagrant_environment()
    elif env.hosts == ["localhost"]:
        _setup_local_environment()
    if env.distribution == "ubuntu":
        _setup_ubuntu()
    elif env.distribution == "centos":
        _setup_centos()
    elif env.distribution == "scientificlinux":
        _setup_scientificlinux()
    elif env.distribution == "debian":
        _setup_debian()
    else:
        raise ValueError("Unexpected distribution %s" % env.distribution)
    configure_runsudo(env)
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
        actual_dist_name = env.safe_run_output("cat /etc/*release | grep DISTRIB_CODENAME | cut -f 2 -d =")
        if actual_dist_name.lower().find(dist_name) == -1:
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
      "deb http://watson.nci.nih.gov/cran_mirror/bin/linux/ubuntu %s/",  # lastest R versions
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
        "deb http://downloads-distro.mongodb.org/repo/debian-sysvinit dist 10gen", # mongodb
        "deb http://watson.nci.nih.gov/cran_mirror/bin/linux/debian %s-cran/", # lastest R versions
        "deb http://archive.cloudera.com/debian lenny-cdh3 contrib" # Hadoop
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
    env.python_version_ext = ""
    env.ruby_version_ext = "1.9.1"
    if not env.has_key("java_home"):
        # XXX look for a way to find JAVA_HOME automatically
        env.java_home = "/usr/lib/jvm/java-6-openjdk"
    shared_sources = [
        "deb http://nebc.nerc.ac.uk/bio-linux/ unstable bio-linux", # Bio-Linux
        "deb http://download.virtualbox.org/virtualbox/debian %s contrib", # virtualbox
    ]
    return shared_sources

def _setup_centos():
    env.logger.info("CentOS setup")
    env.python_version_ext = "2.6"
    env.ruby_version_ext = ""
    env.pip_cmd = "pip-python"
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"

def _setup_scientificlinux():
    env.logger.info("ScientificLinux setup")
    env.python_version_ext = ""
    env.pip_cmd = "pip-python"
    if not env.has_key("java_home"):
        env.java_home = "/etc/alternatives/java_sdk"

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
    if not env.has_key("java_home"):
        env.java_home = os.environ.get("JAVA_HOME", "/usr/lib/jvm/java-6-openjdk")

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
    env.key_filename = ssh_config["IdentityFile"]
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
