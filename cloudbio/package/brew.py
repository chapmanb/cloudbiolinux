"""Install packages via the MacOSX Homebrew and Linux Linuxbrew package manager.
https://github.com/mxcl/homebrew
https://github.com/Homebrew/linuxbrew
"""
import os

from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages

from fabric.api import quiet

def install_packages(env, to_install=None, packages=None):
    """Install packages using the home brew package manager.

    Handles upgrading brew, tapping required repositories and installing or upgrading
    packages as appropriate.

    `to_install` is a CloudBioLinux compatible set of top level items to add,
    alternatively `packages` is a list of raw package names.
    """
    config_file = get_config_file(env, "packages-homebrew.yaml")
    if to_install:
        (packages, _) = _yaml_to_packages(config_file.base, to_install, config_file.dist)
    brew_cmd = _brew_cmd(env)
    formula_repos = ["homebrew/science"]
    env.safe_run("%s update" % brew_cmd)
    current_taps = set([x.strip() for x in env.safe_run_output("%s tap" % brew_cmd).split()])
    for repo in formula_repos:
        if repo not in current_taps:
            env.safe_run("%s tap %s" % (brew_cmd, repo))
    current_pkgs = set([x.strip() for x in env.safe_run_output("%s list" % brew_cmd).split()])
    outdated_pkgs = set([x.strip() for x in env.safe_run_output("%s outdated" % brew_cmd).split()])
    for pkg in packages:
        if pkg in outdated_pkgs:
            brew_subcmd = "upgrade"
        elif pkg in current_pkgs:
            brew_subcmd = None
        else:
            brew_subcmd = "install"
        if brew_subcmd:
            env.safe_run("%s %s %s" % (brew_cmd, brew_subcmd, pkg))

def _brew_cmd(env):
    """Retrieve brew command for installing homebrew packages.
    """
    local_brew = os.path.join(env.local_install, "bin", "brew")
    for cmd in [local_brew, "brew"]:
        with quiet():
            test_version = env.safe_run("%s --version" % cmd)
        if test_version.succeeded:
            return cmd
    raise ValueError("Did not find installation of Homebrew")
