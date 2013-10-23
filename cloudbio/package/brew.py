"""Install packages via the MacOSX Homebrew and Linux Linuxbrew package manager.
https://github.com/mxcl/homebrew
https://github.com/Homebrew/linuxbrew
"""
import os

from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages

from fabric.api import quiet, cd

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
    ipkgs = {"outdated": set([x.strip() for x in env.safe_run_output("%s outdated" % brew_cmd).split()]),
             "current" : _get_current_pkgs(env, brew_cmd)}
    _install_brew_baseline(env, brew_cmd, ipkgs)
    for pkg_str in packages:
        _install_pkg(env, pkg_str, brew_cmd, ipkgs)

def _get_current_pkgs(env, brew_cmd):
    out = {}
    with quiet():
        which_out = env.safe_run_output("{brew_cmd} which".format(**locals()))
    for line in which_out.split("\n"):
        if line:
            pkg, version = line.rstrip().split()
            if pkg.endswith(":"):
                pkg = pkg[:-1]
            out[pkg] = version
    return out

def _install_pkg(env, pkg_str, brew_cmd, ipkgs):
    """Install a specific brew package, handling versioning and existing packages.
    """
    pkg, version = _get_pkg_and_version(pkg_str)
    if version:
        _install_pkg_version(env, pkg, version, brew_cmd, ipkgs)
    else:
        _install_pkg_latest(env, pkg, brew_cmd, ipkgs)

def _install_pkg_version(env, pkg, version, brew_cmd, ipkgs):
    """Install a specific version of a package by retrieving from git history.
    https://gist.github.com/gcatlin/1847248
    Handles both global packages and those installed via specific taps.
    """
    if ipkgs["current"].get(pkg) == version:
        return
    git_cmd = _git_cmd_for_pkg_version(env, brew_cmd, pkg, version)
    git_fname = git_cmd.split()[-1]
    brew_prefix = env.safe_run_output("{brew_cmd} --prefix".format(**locals()))
    if git_fname.startswith("{brew_prefix}/Library/Taps/".format(**locals())):
        brew_prefix = os.path.dirname(git_fname)
    print brew_prefix, git_fname
    with cd(brew_prefix):
        env.safe_run(git_cmd)
    if pkg in ipkgs["current"]:
        env.safe_run("{brew_cmd} unlink {pkg}".format(**locals()))
    env.safe_run("{brew_cmd} install {pkg}".format(**locals()))
    env.safe_run("{brew_cmd} switch {pkg} {version}".format(**locals()))
    env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))
    # reset Git back to latest
    with cd(brew_prefix):
        cmd_parts = git_cmd.split()
        cmd_parts[2] = "--"
        env.safe_run(" ".join(cmd_parts))

def _git_cmd_for_pkg_version(env, brew_cmd, pkg, version):
    """Retrieve git command to check out a specific version from homebrew.
    """
    git_cmd = None
    for git_line in env.safe_run_output("{brew_cmd} versions {pkg}".format(**locals())).split("\n"):
        if git_line.startswith(version):
            git_cmd = " ".join(git_line.rstrip().split()[1:])
            break
    if git_cmd is None:
        raise ValueError("Did not find version %s for %s" % (version, pkg))
    return git_cmd

def _latest_pkg_version(env, brew_cmd, pkg):
    """Retrieve the latest available version of a package.
    """
    for git_line in env.safe_run_output("{brew_cmd} versions {pkg}".format(**locals())).split("\n"):
        if git_line.strip():
            return git_line.split()[0].strip()

def _install_pkg_latest(env, pkg, brew_cmd, ipkgs):
    """Install the latest version of the given package.
    """
    if pkg in ipkgs["outdated"]:
        brew_subcmd = "upgrade"
    elif pkg in ipkgs["current"]:
        brew_subcmd = None
        pkg_version =  _latest_pkg_version(env, brew_cmd, pkg)
        if ipkgs["current"][pkg] != pkg_version:
            _install_pkg_version(env, pkg, pkg_version, brew_cmd, ipkgs)
    else:
        brew_subcmd = "install"
    if brew_subcmd:
        perl_setup = "export PERL5LIB=%s/lib/perl5:${PERL5LIB}" % env.system_install
        env.safe_run("%s && %s %s %s" % (perl_setup, brew_cmd, brew_subcmd, pkg))
        env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))

def _get_pkg_and_version(pkg_str):
    """Uses Python style package==0.1 version specifications.
    """
    parts = pkg_str.split("==")
    if len(parts) == 1:
        return parts[0], None
    else:
        assert len(parts) == 2
        return parts

def _install_brew_baseline(env, brew_cmd, ipkgs):
    """Install baseline brew components not handled by dependency system.

    Handles installation of required Perl libraries.
    """
    _install_pkg_latest(env, "cpanminus", brew_cmd, ipkgs)
    cpanm_cmd = os.path.join(os.path.dirname(brew_cmd), "cpanm")
    for perl_lib in ["Statistics::Descriptive"]:
        env.safe_run("%s -i --notest --local-lib=%s '%s'" % (cpanm_cmd, env.system_install, perl_lib))

def _brew_cmd(env):
    """Retrieve brew command for installing homebrew packages.
    """
    local_brew = os.path.join(env.system_install, "bin", "brew")
    for cmd in [local_brew, "brew"]:
        with quiet():
            test_version = env.safe_run("%s --version" % cmd)
        if test_version.succeeded:
            return cmd
    raise ValueError("Did not find installation of Homebrew")
