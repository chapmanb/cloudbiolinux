"""Install software and configure package managers.
"""
import os

from fabric.api import env, cd

from cloudbio.custom.shared import _make_tmp_dir
from cloudbio.package import brew
from cloudbio.package.deb import (_apt_packages, _add_apt_gpg_keys,
                                  _setup_apt_automation, _setup_apt_sources)
from cloudbio.package.rpm import (_yum_packages, _setup_yum_bashrc,
                                  _setup_yum_sources)

def _configure_and_install_native_packages(env, pkg_install):
    """
    Setups up native package repositories, determines list
    of native packages to install, and installs them.
    """
    home_dir = env.safe_run("echo $HOME")
    if home_dir:
        if env.shell_config.startswith("~"):
            nonhome = env.shell_config.split("~/", 1)[-1]
            env.shell_config = os.path.join(home_dir, nonhome)
    if env.distribution in ["debian", "ubuntu"]:
        _setup_apt_sources()
        _setup_apt_automation()
        _add_apt_gpg_keys()
        _apt_packages(pkg_install)
    elif env.distribution in ["centos", "scientificlinux"]:
        _setup_yum_sources()
        _yum_packages(pkg_install)
        _setup_yum_bashrc()
    elif env.distribution in ["arch", "suse"]:
        pass  # No package support for Arch, SUSE yet
    elif env.distribution == "macosx":
        brew.install_packages(env, pkg_install)
    else:
        raise NotImplementedError("Unknown target distribution")

def _connect_native_packages(env, pkg_install, lib_install):
    """Connect native installed packages to local versions.

    This helps setup a non-sudo environment to handle software
    that needs a local version in our non-root directory tree.
    """
    bin_dir = os.path.join(env.system_install, "bin")
    exports = _get_shell_exports(env)
    path = env.safe_run_output("echo $PATH")
    comment_line = "# CloudBioLinux PATH updates"
    if not env.safe_contains(env.shell_config, comment_line):
        env.safe_append(env.shell_config, "\n" + comment_line)
    if bin_dir not in path and env.safe_exists(env.shell_config):
        if not env.safe_contains(env.shell_config, exports["path"]):
            env.safe_append(env.shell_config, exports["path"])
    if not env.safe_contains(env.shell_config, exports["ld_library"]):
        env.safe_append(env.shell_config, exports["ld_library"])
    if not env.safe_contains(env.shell_config, exports["perl"]):
        env.safe_append(env.shell_config, exports["perl"])
    if "python" in pkg_install and "python" in lib_install:
        _create_local_virtualenv(env.system_install)

def _get_shell_exports(env):
    bin_dir = os.path.join(env.system_install, "bin")
    ldlib_path = os.path.join(env.system_install, "lib")
    return {"path": "export PATH=%s:$PATH" % bin_dir,
            "ld_library": "export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH" % ldlib_path,
            "perl": "export PERL5LIB=%s/lib/perl5:${PERL5LIB}" % env.system_install}

def _print_shell_exports(env):
    """Print a set of exports to add to shell in isolated installations.
    """
    exports = _get_shell_exports(env)
    print "\nIsolated installation: add the following to your shell to include installed files:"
    for e in ["path", "ld_library", "perl"]:
        print exports[e]

def _create_local_virtualenv(target_dir):
    """Create virtualenv in target directory for non-sudo installs.
    """
    url = "https://raw.github.com/pypa/virtualenv/master/virtualenv.py"
    if not os.path.exists(os.path.join(target_dir, "bin", "python")):
        with _make_tmp_dir() as work_dir:
            with cd(work_dir):
                env.safe_run("wget --no-check-certificate %s" % url)
                env.safe_run("python virtualenv.py %s" % target_dir)
