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
        conda_envs = _create_environments(env, conda_bin)
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
            for env_name, env_packages in _split_by_condaenv(packages):
                if env_name:
                    assert env_name in conda_envs, (env_name, conda_envs)
                    env_str = "-n %s" % env_name
                else:
                    env_str = ""
                pkgs_str = " ".join(env_packages)
                env.safe_run("{conda_bin} install --quiet -y {env_str} {channels} {pkgs_str}".format(**locals()))
                conda_pkg_list = json.loads(env.safe_run_output(
                    "{conda_bin} list --json {env_str}".format(**locals())))
                for package in env_packages:
                    _link_bin(package, env, conda_info, conda_bin, conda_pkg_list,
                              conda_envdir=conda_envs.get(env_name))
        conda_pkg_list = json.loads(env.safe_run_output("{conda_bin} list --json".format(**locals())))
        for pkg in ["python", "conda", "pip"]:
            _link_bin(pkg, env, conda_info, conda_bin, conda_pkg_list, files=[pkg], prefix="bcbio_")

def _link_bin(package, env, conda_info, conda_bin, conda_pkg_list, files=None, prefix="", conda_env=None,
              conda_envdir=None):
    """Link files installed in the bin directory into the install directory.

    This is imperfect but we're trying not to require injecting everything in the anaconda
    directory into a user's path.
    """
    package = package.split("=")[0]
    final_bindir = os.path.join(env.system_install, "bin")
    if conda_envdir:
        base_bindir = os.path.join(conda_envdir, "bin")
    else:
        base_bindir = os.path.dirname(conda_bin)
    # resolve any symlinks in the final and base heirarchies
    with quiet():
        final_bindir = env.safe_run_output("cd %s && pwd -P" % final_bindir)
        base_bindir = env.safe_run_output("cd %s && pwd -P" % base_bindir)
    for pkg_subdir in [x for x in conda_pkg_list if x["name"] == package]:
        pkg_subdir = pkg_subdir["dist_name"].split("::")[-1]
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

def _split_by_condaenv(packages):
    """Split packages into those requiring special conda environments.
    """
    out = collections.defaultdict(list)
    for p in packages:
        parts = p.split(";")
        name = parts[0]
        metadata = parts[1:]
        condaenv = None
        for k, v in [x.split("=") for x in metadata]:
            if k == "env":
                condaenv = v
        out[condaenv].append(name)
    return dict(out).items()

def _create_environments(env, conda_bin):
    """Create a custom local build environment for tools. Handles python2/python3 divide.

    This is an initial step towards transitioning to more python3 tool support.
    """
    out = {}
    conda_envs = json.loads(env.safe_run_output("{conda_bin} info --envs --json".format(**locals())))["envs"]
    if not any(x.endswith("/python3") for x in conda_envs):
        env.safe_run("{conda_bin} create -y --name python3 python=3".format(**locals()))
    out["python3"] = [x for x in conda_envs if x.endswith("/python3")][0]
    return out
