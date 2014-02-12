"""Install packages via the MacOSX Homebrew and Linux Linuxbrew package manager.
https://github.com/mxcl/homebrew
https://github.com/Homebrew/linuxbrew
"""
import contextlib
from distutils.version import LooseVersion
import os

from cloudbio.flavor.config import get_config_file
from cloudbio.fabutils import quiet, find_cmd
from cloudbio.package.shared import _yaml_to_packages

from fabric.api import cd, settings

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
    formula_repos = ["homebrew/science", "chapmanb/cbl"]
    current_taps = set([x.strip() for x in env.safe_run_output("%s tap" % brew_cmd).split()])
    _safe_update(env, brew_cmd, formula_repos, current_taps)
    for repo in formula_repos:
        if repo not in current_taps:
            env.safe_run("%s tap %s" % (brew_cmd, repo))
    env.safe_run("%s tap --repair" % brew_cmd)
    ipkgs = {"outdated": set([x.strip() for x in env.safe_run_output("%s outdated" % brew_cmd).split()]),
             "current": _get_current_pkgs(env, brew_cmd)}
    _install_brew_baseline(env, brew_cmd, ipkgs, packages)
    for pkg_str in packages:
        _install_pkg(env, pkg_str, brew_cmd, ipkgs)

def _safe_update(env, brew_cmd, formula_repos, cur_taps):
    """Revert any taps if we fail to update due to local changes.
    """
    with settings(warn_only=True):
        out = env.safe_run("%s update" % brew_cmd)
    if out.failed:
        for repo in formula_repos:
            if repo in cur_taps:
                env.safe_run("%s untap %s" % (brew_cmd, repo))
        env.safe_run("%s update" % brew_cmd)

def _get_current_pkgs(env, brew_cmd):
    out = {}
    with quiet():
        which_out = env.safe_run_output("{brew_cmd} which".format(**locals()))
    for line in which_out.split("\n"):
        if line:
            try:
                pkg, version = line.rstrip().split()
                if pkg.endswith(":"):
                    pkg = pkg[:-1]
                    out[pkg] = version
            except:
                print(line)
                continue
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
    if ipkgs["current"].get(pkg.split("/")[-1]) == version:
        return
    with _git_pkg_version(env, brew_cmd, pkg, version):
        if pkg.split("/")[-1] in ipkgs["current"]:
            with settings(warn_only=True):
                env.safe_run("{brew_cmd} unlink {pkg}".format(
                    brew_cmd=brew_cmd, pkg=pkg.split("/")[-1]))
        # if we have a more recent version, uninstall that first
        cur_version_parts = env.safe_run_output("{brew_cmd} list --versions {pkg}".format(
            brew_cmd=brew_cmd, pkg=pkg.split("/")[-1])).strip().split()
        if len(cur_version_parts) > 1 and LooseVersion(cur_version_parts[1]) > LooseVersion(version):
            with settings(warn_only=True):
                env.safe_run("{brew_cmd} uninstall {pkg}".format(**locals()))
        if version == "HEAD":
            env.safe_run("{brew_cmd} install --HEAD {pkg}".format(**locals()))
        else:
            env.safe_run("{brew_cmd} install {pkg}".format(**locals()))
            with settings(warn_only=True):
                env.safe_run("{brew_cmd} switch {pkg} {version}".format(**locals()))
        env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))

@contextlib.contextmanager
def _git_pkg_version(env, brew_cmd, pkg, version):
    """Convert homebrew Git to previous revision to install a specific package version.
    """
    git_cmd = _git_cmd_for_pkg_version(env, brew_cmd, pkg, version)
    git_fname = git_cmd.split()[-1]
    brew_prefix = env.safe_run_output("{brew_cmd} --prefix".format(**locals()))
    if git_fname.startswith("{brew_prefix}/Library/Taps/".format(**locals())):
        brew_prefix = os.path.dirname(git_fname)
    try:
        with cd(brew_prefix):
            if version != "HEAD":
                env.safe_run(git_cmd)
        yield
    finally:
        # reset Git back to latest
        with cd(brew_prefix):
            if version != "HEAD":
                cmd_parts = git_cmd.split()
                env.safe_run("%s reset HEAD %s" % (cmd_parts[0], cmd_parts[-1]))
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
        pkg_version = _latest_pkg_version(env, brew_cmd, pkg)
        if ipkgs["current"][pkg] != pkg_version:
            _install_pkg_version(env, pkg, pkg_version, brew_cmd, ipkgs)
    else:
        brew_subcmd = "install"
    if brew_subcmd:
        perl_setup = "export PERL5LIB=%s/lib/perl5:${PERL5LIB}" % env.system_install
        compiler_setup = "export CC=${CC:-`which gcc`} && export CXX=${CXX:-`which g++`}"
        env.safe_run("%s && %s && %s %s --env=inherit %s" % (compiler_setup, perl_setup, brew_cmd, brew_subcmd, pkg))
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

def _install_brew_baseline(env, brew_cmd, ipkgs, packages):
    """Install baseline brew components not handled by dependency system.

    Handles installation of required Perl libraries.
    """
    for dep in ["cpanminus", "expat"]:
        _install_pkg_latest(env, dep, brew_cmd, ipkgs)
    # if installing samtools, avoid conflicts with cbl and homebrew-science versions
    if len([x for x in packages if x.find("samtools") >= 0]):
        with settings(warn_only=True):
            has_bcftools = int(env.safe_run_output("{brew_cmd} list samtools | grep -c bcftools".format(
                brew_cmd=brew_cmd)))
            if has_bcftools:
                env.safe_run("{brew_cmd} uninstall {pkg}".format(brew_cmd=brew_cmd, pkg="samtools"))
    cpanm_cmd = os.path.join(os.path.dirname(brew_cmd), "cpanm")
    for perl_lib in ["Statistics::Descriptive"]:
        env.safe_run("%s -i --notest --local-lib=%s '%s'" % (cpanm_cmd, env.system_install, perl_lib))
    # Ensure paths we may have missed on install are accessible to regular user
    if env.use_sudo:
        paths = ["share", "share/java"]
        for path in paths:
            with quiet():
                test_access = env.safe_run("test -d %s/%s && test -O %s/%s" % (env.system_install, path,
                                                                               env.system_install, path))
            if test_access.failed and env.safe_exists("%s/%s" % (env.system_install, path)):
                env.safe_sudo("chown %s %s/%s" % (env.user, env.system_install, path))

def _brew_cmd(env):
    """Retrieve brew command for installing homebrew packages.
    """
    cmd = find_cmd(env, "brew", "--version")
    if cmd is None:
        raise ValueError("Did not find working installation of Linuxbrew/Homebrew")
    else:
        return cmd
