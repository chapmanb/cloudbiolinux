"""Install packages via the MacOSX Homebrew and Linux Linuxbrew package manager.
https://github.com/mxcl/homebrew
https://github.com/Homebrew/linuxbrew
"""
import contextlib
from distutils.version import LooseVersion
import os
import sys

from cloudbio.custom import system, shared
from cloudbio.flavor.config import get_config_file
from cloudbio.fabutils import quiet, find_cmd
from cloudbio.package import cpan
from cloudbio.package.shared import _yaml_to_packages

from fabric.api import cd, settings

BOTTLE_URL = "https://s3.amazonaws.com/cloudbiolinux/brew_bottles/{pkg}-{version}.x86_64-linux.bottle.tar.gz"
BOTTLE_SUPPORTED = set(["isaac-aligner", "isaac-variant-caller", "cmake"])

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
    # if we have no packages to install, do not try to install or update brew
    if len(packages) == 0:
        return
    system.install_homebrew(env)
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
    ipkgs = {"outdated": set([x.strip() for x in env.safe_run_output("%s outdated" % brew_cmd).split()]),
             "current": _get_current_pkgs(env, brew_cmd)}
    for pkg_str in packages:
        _install_pkg(env, pkg_str, brew_cmd, ipkgs)
    for pkg_str in ["pkg-config", "openssl", "cmake"]:
        _safe_unlink_pkg(env, pkg_str, brew_cmd)
    for pkg_str in ["curl"]:
        _safe_uninstall_pkg(env, pkg_str, brew_cmd)


def _safe_update(env, brew_cmd, formula_repos, cur_taps):
    """Revert any taps if we fail to update due to local changes.
    """
    with quiet():
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
        which_out = env.safe_run_output("{brew_cmd} list --versions".format(**locals()))
    for line in which_out.split("\n"):
        if line:
            parts = line.rstrip().split()
            if len(parts) == 2:
                pkg, version = line.rstrip().split()
                if pkg.endswith(":"):
                    pkg = pkg[:-1]
                out[pkg] = version
    return out

def _safe_unlink_pkg(env, pkg_str, brew_cmd):
    """Unlink packages which can cause issues with a Linux system.
    """
    with settings(warn_only=True):
        with quiet():
            env.safe_run("{brew_cmd} unlink {pkg_str}".format(**locals()))

def _safe_uninstall_pkg(env, pkg_str, brew_cmd):
    """Uninstall packages which get pulled in even when unlinked by brew.
    """
    with settings(warn_only=True):
        with quiet():
            env.safe_run("{brew_cmd} uninstall {pkg_str}".format(**locals()))

def _install_pkg(env, pkg_str, brew_cmd, ipkgs):
    """Install a specific brew package, handling versioning and existing packages.
    """
    pkg, version, args = _get_pkg_version_args(pkg_str)
    installed = False
    if version:
        _install_pkg_version(env, pkg, args, version, brew_cmd, ipkgs)
        installed = True
    if pkg in BOTTLE_SUPPORTED and not env.distribution == "macosx":
        installed = _install_bottle(env, brew_cmd, pkg, ipkgs)
    if not installed:
        _install_pkg_latest(env, pkg, args, brew_cmd, ipkgs)

def _install_pkg_version(env, pkg, args, version, brew_cmd, ipkgs):
    """Install a specific version of a package by retrieving from git history.
    https://gist.github.com/gcatlin/1847248
    Handles both global packages and those installed via specific taps.
    """
    if ipkgs["current"].get(pkg.split("/")[-1]) == version:
        return
    if version == "HEAD":
        args = " ".join(args)
        brew_install = _get_brew_install_cmd(brew_cmd, env, pkg)
        env.safe_run("{brew_install} {args} --HEAD {pkg}".format(**locals()))
    else:
        raise ValueError("Cannot currently handle installing brew packages by version.")
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

def _latest_pkg_version(env, brew_cmd, pkg, devel=False):
    """Retrieve the latest available version of a package and if it is linked.
    """
    i = 0
    version, is_linked = None, False
    with settings(warn_only=True):
	info_str = env.safe_run_output("{brew_cmd} info {pkg}".format(**locals()))
    for i, git_line in enumerate(info_str.split("\n")):
        if git_line.strip():
            if i == 0:
                _, version_str = git_line.split(":")
                versions = version_str.split(",")
                if devel:
                    dev_strs = [x for x in versions if x.strip().startswith("devel")]
                    version = dev_strs[0].split()[-1].strip()
                else:
                    version = versions[0].replace("(bottled)", "").split()[-1].strip()
            elif git_line.find("Cellar/%s" % pkg) > 0 and git_line.find(" files,") > 0:
                is_linked = git_line.strip().split()[-1] == "*"
    return version, is_linked

def _get_brew_install_cmd(brew_cmd, env, pkg):
    perl_setup = "export PERL5LIB=%s/lib/perl5:${PERL5LIB}" % env.system_install
    compiler_setup = "export CC=${CC:-`which gcc`} && export CXX=${CXX:-`which g++`}"
    extra_args = ""
    if pkg in ["cmake"]:
        extra_args += " --without-docs"
    if pkg in ["lumpy-sv", "bamtools", "freebayes"]:
        extra_args += " --ignore-dependencies"
    return "%s && %s && %s install --env=inherit %s" % (compiler_setup, perl_setup, brew_cmd, extra_args)

def _install_pkg_latest(env, pkg, args, brew_cmd, ipkgs):
    """Install the latest version of the given package.
    """
    short_pkg = pkg.split("/")[-1]
    do_install = True
    is_linked = True
    remove_old = False
    if pkg in ipkgs["outdated"] or short_pkg in ipkgs["outdated"]:
        remove_old = True
    elif pkg in ipkgs["current"] or short_pkg in ipkgs["current"]:
        do_install = False
        pkg_version, is_linked = _latest_pkg_version(env, brew_cmd, pkg, devel="--devel" in args)
        cur_version = ipkgs["current"].get(pkg, ipkgs["current"][short_pkg])
        if cur_version != pkg_version and cur_version.split("_")[0] != pkg_version:
            remove_old = True
            do_install = True
    if do_install:
        if remove_old:
            env.safe_run("{brew_cmd} remove --force {short_pkg}".format(**locals()))
        flags = " ".join(args)
        with settings(warn_only=True):
            cmd = "%s %s %s" % (_get_brew_install_cmd(brew_cmd, env, pkg), flags, pkg)
            result = env.safe_run_output(cmd)
            if result.failed and not result.find("Could not symlink") > 0:
                sys.tracebacklimit = 1
                raise ValueError("Failed to install brew formula: %s\n" % pkg +
                                 "To debug, please try re-running the install command with verbose output:\n" +
                                 cmd.replace("brew install", "brew install -v"))
        env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))
    # installed but not linked
    elif not is_linked:
        env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))

def _get_pkg_version_args(pkg_str):
    """Uses Python style package==0.1 version specifications and args separated with ';'
    """
    arg_parts = pkg_str.split(";")
    if len(arg_parts) == 1:
        args = []
    else:
        pkg_str = arg_parts[0]
        args = arg_parts[1:]
    parts = pkg_str.split("==")
    if len(parts) == 1:
        return parts[0], None, args
    else:
        assert len(parts) == 2
        name, version = parts
        return name, version, args

def _install_bottle(env, brew_cmd, pkg, ipkgs):
    """Install Linux bottles for brew packages that can be tricky to build.
    """
    if env.distribution == "macosx":  # Only Linux bottles, build away on Mac
        return False
    pkg_version, is_linked = _latest_pkg_version(env, brew_cmd, pkg)
    install_version = ipkgs["current"].get(pkg)
    if pkg_version == install_version:  # Up to date
        if not is_linked:
            env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))
        return True
    elif install_version or pkg in ipkgs["outdated"]:
        env.safe_run("{brew_cmd} remove --force {pkg}".format(**locals()))
    url = BOTTLE_URL.format(pkg=pkg, version=pkg_version)
    brew_cachedir = env.safe_run_output("%s --cache" % brew_cmd)
    brew_cellar = os.path.join(env.safe_run_output("%s --prefix" % brew_cmd), "Cellar")
    with quiet():
        env.safe_run("mkdir -p %s" % brew_cellar)
    out_file = os.path.join(brew_cachedir, os.path.basename(url))
    if env.safe_exists(out_file):
        env.safe_run("rm -f %s" % out_file)
    bottle_file = shared._remote_fetch(env, url, out_file=out_file,
                                       allow_fail=True, samedir=True)
    if bottle_file:
        with cd(brew_cellar):
            env.safe_run("tar -xf %s" % bottle_file)
        env.safe_run("%s link --overwrite %s" % (brew_cmd, pkg))
        return True
    else:
        return False

def _install_brew_baseline(env, brew_cmd, ipkgs, packages):
    """Install baseline brew components not handled by dependency system.

    - Installation of required Perl libraries.
    - Ensures installed samtools does not overlap with bcftools
    - Upgrades any package dependencies
    """
    for dep in ["expat", "cmake", "pkg-config"]:
        _install_pkg(env, dep, brew_cmd, ipkgs)
    for dep in ["sambamba"]:  # Avoid conflict with homebrew-science sambamba
        env.safe_run("{brew_cmd} remove --force {dep}".format(**locals()))
    # if installing samtools, avoid bcftools conflicts
    if len([x for x in packages if x.find("samtools") >= 0]):
        with settings(warn_only=True):
            def _has_prog(prog):
                try:
                    return int(env.safe_run_output("{brew_cmd} list samtools | grep -c {prog} | cat".format(
                        brew_cmd=brew_cmd, prog=prog)))
                except ValueError:
                    return 0
            if any(_has_prog(p) for p in ["bctools", "vcfutils.pl"]):
                env.safe_run("{brew_cmd} uninstall {pkg}".format(brew_cmd=brew_cmd, pkg="samtools"))
                ipkgs["current"].pop("samtools", None)
        _install_pkg_latest(env, "samtools", ["--without-curses"], brew_cmd, ipkgs)
    for dependency in ["htslib"]:
        if dependency in packages:
            if (dependency in ipkgs["outdated"] or "chapmanb/cbl/%s" % dependency in ipkgs["outdated"]
                  or dependency not in ipkgs["current"]):
                _install_pkg_latest(env, dependency, [], brew_cmd, ipkgs)
    if "cpanminus" in packages:
        _install_pkg_latest(env, "cpanminus", [], brew_cmd, ipkgs)
        _install_pkg_latest(env, "samtools-library-0.1", [], brew_cmd, ipkgs)
        cpan.install_packages(env)
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
        raise ValueError("Did not find working installation of Linuxbrew/Homebrew. "
                         "Please check if you have ruby available.")
    else:
        return cmd
