"""Install packages via the Conda package manager: http://conda.pydata.org/
"""
import json
import os
import yaml

from cloudbio.custom import shared
from cloudbio.fabutils import quiet
from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages

from fabric.api import settings

def install_packages(env, to_install=None, packages=None):
    if shared._is_anaconda(env):
        conda_bin = shared._conda_cmd(env)
        config_file = get_config_file(env, "packages-conda.yaml")
        if config_file.base is None and packages is None:
            packages = []
        else:
            if to_install:
                (packages, _) = _yaml_to_packages(config_file.base, to_install, config_file.dist)
            with open(config_file.base) as in_handle:
                channels = " ".join(["-c %s" % x for x in yaml.safe_load(in_handle).get("channels", [])])
        conda_info = json.loads(env.safe_run_output("{conda_bin} info --json".format(**locals())))
        if len(packages) > 0:
            pkgs_str = " ".join(packages)
            env.safe_run("{conda_bin} install -y {channels} {pkgs_str}".format(**locals()))
            for package in packages:
                _link_bin(package, env, conda_info, conda_bin)
        for pkg in ["python", "conda", "pip"]:
            _link_bin(pkg, env, conda_info, conda_bin, [pkg], "bcbio_")
        # remove packages we want the system to supply
        # curl https://github.com/ContinuumIO/anaconda-issues/issues/72
        system_packages = ["curl"]
        pkgs_str = " ".join(system_packages)
        with settings(warn_only=True):
            env.safe_run("{conda_bin} uninstall -y {pkgs_str}".format(**locals()))

def _link_bin(package, env, conda_info, conda_bin, files=None, prefix=""):
    """Link files installed in the bin directory into the install directory.

    This is imperfect but we're trying not to require injecting everything in the anaconda
    directory into a user's path.
    """
    package = package.split("=")[0]
    final_bindir = os.path.join(env.system_install, "bin")
    base_bindir = os.path.dirname(conda_bin)
    # resolve any symlinks in the final and base heirarchies
    with quiet():
        final_bindir = env.safe_run_output("cd %s && pwd -P" % final_bindir)
        base_bindir = env.safe_run_output("cd %s && pwd -P" % base_bindir)
    for pkg_subdir in json.loads(env.safe_run_output("{conda_bin} list --json -f {package}".format(**locals()))):
        for pkg_dir in conda_info["pkgs_dirs"]:
            pkg_bindir = os.path.join(pkg_dir, pkg_subdir, "bin")
            if env.safe_exists(pkg_bindir):
                if not files:
                    with quiet():
                        files = env.safe_run_output("ls -1 {pkg_bindir}".format(**locals())).split()
                for fname in files:
                    # symlink to the original file in the /anaconda/bin directory
                    # this could be a hard or soft link
                    base_fname = os.path.join(base_bindir, fname)
                    if os.path.exists(base_fname) and os.path.lexists(base_fname):
                        _do_link(base_fname,
                                 os.path.join(final_bindir, "%s%s" % (prefix, fname)))

def _do_link(orig_file, final_file):
    """Perform a soft link of the original file into the final location.

    We need the symlink to point to /anaconda/bin directory, not the real location
    in the pkgs directory so conda can resolve LD_LIBRARY_PATH and the interpreters.
    """
    needs_link = True
    # working symlink, check if already in the right place or remove it
    if os.path.exists(final_file):
        if (os.path.realpath(final_file) == os.path.realpath(orig_file) and
              orig_file == os.path.normpath(os.path.join(os.path.dirname(final_file), os.readlink(final_file)))):
            needs_link = False
        else:
            os.remove(final_file)
    # broken symlink
    elif os.path.lexists(final_file):
        os.unlink(final_file)
    if needs_link:
        os.symlink(os.path.relpath(orig_file, os.path.dirname(final_file)), final_file)
