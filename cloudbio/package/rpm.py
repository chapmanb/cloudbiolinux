"""Automated installation on RPM systems with the yum package manager.
"""
from fabric.api import *
from fabric.contrib.files import *

from cloudbio.package.shared import _yaml_to_packages

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
    packages = env.flavor.rewrite_config_items("packages", packages)
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
