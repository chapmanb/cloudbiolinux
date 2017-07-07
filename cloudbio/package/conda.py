"""Install packages via the Conda package manager: http://conda.pydata.org/
"""
import collections
import json
import os
import yaml

from cloudbio.custom import shared
from cloudbio.fabutils import quiet
from cloudbio.flavor.config import get_config_file
from cloudbio.package.shared import _yaml_to_packages

def install_packages(env, to_install=None, packages=None):
    if shared._is_anaconda(env):
        conda_bin = shared._conda_cmd(env)
        if hasattr(env, "conda_yaml"):
            Config = collections.namedtuple("Config", "base dist")
            config_file = Config(base=env.conda_yaml, dist=None)
        else:
            config_file = get_config_file(env, "packages-conda.yaml")
        if config_file.base is None and packages is None:
            packages = []
        else:
            if to_install:
                (packages, _) = _yaml_to_packages(config_file.base, to_install, config_file.dist)
            with open(config_file.base) as in_handle:
                channels = " ".join(["-c %s" % x for x in yaml.safe_load(in_handle).get("channels", [])])
        conda_info = json.loads(env.safe_run_output("{conda_bin} info --json".format(**locals())))
        # Uninstall old R packages that clash with updated versions
        # Temporary fix to allow upgrades from older versions that have migrated
        # r-tximport is now bioconductor-tximport
        for problem in ["r-tximport"]:
            cur_packages = [x["name"] for x in
                            json.loads(env.safe_run_output("{conda_bin} list --json {problem}".format(**locals())))]
            if problem in cur_packages:
                env.safe_run("{conda_bin} remove --force -y {problem}".format(**locals()))
        # install our customized packages
        if len(packages) > 0:
            pkgs_str = " ".join(packages)
            env.safe_run("{conda_bin} install --quiet -y {channels} {pkgs_str}".format(**locals()))
            for package in packages:
                _link_bin(package, env, conda_info, conda_bin)
        for pkg in ["python", "conda", "pip"]:
            _link_bin(pkg, env, conda_info, conda_bin, [pkg], "bcbio_")

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
        # New style conda info (4.3+) is a dictionary, old style is a string. Handle both cases
        if isinstance(pkg_subdir, dict):
            pkg_subdir = pkg_subdir["dist_name"]
        pkg_subdir = pkg_subdir.split("::")[-1]
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
