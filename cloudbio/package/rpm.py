"""Automated installation on RPM systems with the yum package manager.
"""
from fabric.api import *
from fabric.contrib.files import *

from cloudbio.package.shared import _yaml_to_packages
from cloudbio.flavor.config import get_config_file

def _yum_packages(to_install):
    """Install rpm packages available via yum.
    """
    if env.distribution == "scientificlinux":
        package_file = "packages-scientificlinux.yaml"
    else:
        package_file = "packages-yum.yaml"
    pkg_config = get_config_file(env, package_file).base
    with settings(warn_only=True):
        env.safe_sudo("yum check-update")
    env.safe_sudo("yum -y upgrade")
    # Retrieve packages to get and install each of them
    (packages, _) = _yaml_to_packages(pkg_config, to_install)
    # At this point allow the Flavor to rewrite the package list
    packages = env.flavor.rewrite_config_items("packages", packages)
    for package in packages:
        env.safe_sudo("yum -y install %s" % package)

def _setup_yum_bashrc():
    """Fix the user bashrc to update compilers.
    """
    if env.distribution in ["centos"]:
        to_include = ["export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/usr/lib/pkgconfig"]
        # gcc fixes no longer necessary on recent CentOS versions
        #"export CC=gcc44", "export CXX=g++44", "export FC=gfortran44",
        fname = env.safe_run_output("ls %s" % env.shell_config)
        for line in to_include:
            if not env.safe_contains(fname, line.split("=")[0]):
                env.safe_append(fname, line)

def _setup_yum_sources():
    """Add additional useful yum repositories.
    """
    repos = [
      "http://dl.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm",
      "http://archive.cloudera.com/redhat/6/x86_64/cdh/cdh3-repository-1.0-1.noarch.rpm"
    ]
    for repo in repos:
        with settings(warn_only=True):
            env.safe_sudo("rpm -Uvh %s" % repo)
