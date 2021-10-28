"""Install packages via the Conda package manager: http://conda.pydata.org/
"""
from __future__ import print_function
import collections
import json
import os
import shutil
import subprocess
import yaml

from cloudbio.package.shared import _yaml_to_packages

ENV_PY_VERSIONS = collections.defaultdict(lambda: "python=3")
ENV_PY_VERSIONS[None] = "python=3"
ENV_PY_VERSIONS["python2"] = "python=2"
ENV_PY_VERSIONS["python3.6"] = "python=3.6"
ENV_PY_VERSIONS["dv"] = "python=3.6"
ENV_PY_VERSIONS["samtools0"] = "python=3"
ENV_PY_VERSIONS["r35"] = "python=3"
ENV_PY_VERSIONS["rbcbiornaseq"] = "python=3"
ENV_PY_VERSIONS["htslib1.9"] = "python=3"
ENV_PY_VERSIONS["htslib1.10"] = "python=3"
ENV_PY_VERSIONS["htslib1.11"] = "python=3"
ENV_PY_VERSIONS["htslib1.12"] = "python=3"
ENV_PY_VERSIONS["htslib1.12_py3.9"] = "python=3.9"
ENV_PY_VERSIONS["java"] = "python=3"
ENV_PY_VERSIONS["bwakit"] = "python=3"

def install_packages(env, to_install=None, packages=None):
    """Old installation, based on pre-configured fabric inputs.
    """
    from cloudbio.flavor.config import get_config_file
    from cloudbio.custom import shared

    if shared._is_anaconda(env):
        conda_bin = shared._conda_cmd(env)
        if hasattr(env, "conda_yaml"):
            Config = collections.namedtuple("Config", "base dist")
            config_file = Config(base=env.conda_yaml, dist=None)
        else:
            config_file = get_config_file(env, "packages-conda.yaml")
        install_in(conda_bin, env.system_install, config_file.base, packages)

def _install_env_pkgs(env_name, env_packages, conda_bin, conda_envs, channels):
    """Install packages into the given environment.

    Uses mamba for initial install for speed, following by conda for completeness.

    TODO: currently duplicates mamba base code in _initial_base_install to make it
    easy to remove or roll back general mamba usage. We can refactor _initial_base_install
    in favor of this after further testing.

    conda_bin could refer to mamba
    """
    mamba_bin = os.path.join(os.path.dirname(conda_bin), "mamba")
    conda_bin = os.path.join(os.path.dirname(mamba_bin), "conda")
    if env_name:
        assert env_name in conda_envs, (env_name, conda_envs)
        env_str = "-n %s" % env_name
    else:
        env_str = ""
    pkgs_str = " ".join(["'%s'" % x for x in sorted(env_packages)])
    py_version = ENV_PY_VERSIONS[env_name]
    if "deepvariant" in env_packages:
        # Ignore /etc/boto.cfg which creates conflicts with conda gsutils
        # https://github.com/GoogleCloudPlatform/gsutil/issues/516
        exports = "export BOTO_CONFIG=/ignoreglobal && "
    else:
        exports = ""
    if os.path.exists(mamba_bin):
        try:
            subprocess.check_call("{mamba_bin} install -q -y {env_str} {channels} "
                                  "{py_version} {pkgs_str}".format(**locals()), shell=True)
        except subprocess.CalledProcessError:
            # Fall back to standard conda install when we have system specific issues
            # https://github.com/bcbio/bcbio-nextgen/issues/2871
            subprocess.check_call("{exports}{conda_bin} install -q -y {env_str} {channels} "
                                "{py_version} {pkgs_str}".format(**locals()), shell=True)
    conda_pkg_list = json.loads(subprocess.check_output(
        "{conda_bin} list --json {env_str} -q".format(**locals()), shell=True))
    return conda_pkg_list

def install_in(conda_bin, system_installdir, config_file=None, packages=None):
    """Install packages inside a given anaconda directory.

    New approach, local only and not dependent on fabric.
    conda_bin could refer to mamba
    """
    if config_file is None and packages is None:
        packages = []
        check_channels = []
    else:
        (packages, _) = _yaml_to_packages(config_file)
        with open(config_file) as in_handle:
            check_channels = yaml.safe_load(in_handle).get("channels", [])
    channels = " ".join(["-c %s" % x for x in check_channels])
    conda_envs = _create_environments(conda_bin, packages)
    for env_dir in conda_envs.values():
        _clean_environment(env_dir)
    conda_info = json.loads(subprocess.check_output("{conda_bin} info --json -q".format(**locals()), shell=True))
    # Uninstall old R packages that clash with updated versions
    # Temporary fix to allow upgrades from older versions that have migrated
    # r-tximport is now bioconductor-tximport
    # py2cairo is incompatible with r 3.4.1+
    problems = ["r-tximport", "py2cairo"]
    # Add packages migrated into separate environments, like python2
    for env_name, env_packages in _split_by_condaenv(packages):
        if env_name:
            problems += env_packages
    if problems:
        print("Checking for problematic or migrated packages in default environment")
        cur_packages = [x["name"] for x in
                        json.loads(subprocess.check_output("%s list --json -q" % (conda_bin), shell=True))
                        if x["name"] in problems and x["channel"] in check_channels]
        if cur_packages:
            print("Found packages that moved from default environment: %s" % ", ".join(cur_packages))
            problems = " ".join(cur_packages)
            subprocess.check_call("{conda_bin} remove {channels} -y {problems}".format(**locals()), shell=True)
    _initial_base_install(conda_bin, [ps for (n, ps) in _split_by_condaenv(packages) if n is None][0],
                          check_channels)
    # install our customized packages
    if len(packages) > 0:
        for env_name, env_packages in _split_by_condaenv(packages):
            print("# Installing into conda environment %s: %s" % (env_name or "default", ", ".join(env_packages)))
            conda_pkg_list = _install_env_pkgs(env_name, env_packages, conda_bin, conda_envs, channels)
            for package in env_packages:
                _link_bin(package, system_installdir, conda_info, conda_bin, conda_pkg_list,
                            conda_envdir=conda_envs.get(env_name))
    conda_pkg_list = json.loads(subprocess.check_output("{conda_bin} list --json -q".format(**locals()), shell=True))
    for pkg in ["python", "conda", "pip"]:
        _link_bin(pkg, system_installdir, conda_info, conda_bin, conda_pkg_list, files=[pkg], prefix="bcbio_")

def _initial_base_install(conda_bin, env_packages, check_channels):
    """Provide a faster initial installation of base packages, avoiding dependency issues.

    Uses mamba (https://github.com/QuantStack/mamba) to provide quicker package resolution
    and avoid dependency conflicts with base install environment. Bootstraps the initial
    installation of all tools when key inputs that cause conflicts are missing.

    TODO: we could remove mamba package running code here in favor of _install_env_pkgs general
    mamba usage once that is further tested.
    """
    initial_package_targets = {None: ["r-base"]}
    env_name = None
    env_str = ""
    channels = " ".join(["-c %s" % x for x in check_channels])
    cur_ps = [x["name"] for x in
              json.loads(subprocess.check_output("{conda_bin} list --json {env_str} -q".format(**locals()), shell=True))
              if x["channel"] in check_channels]
    have_package_targets = env_name in initial_package_targets and any([p for p in cur_ps
                                                                        if p in initial_package_targets[env_name]])
    if not have_package_targets:
        print("Initalling initial set of packages for %s environment with mamba" % (env_name or "default"))
        py_version = ENV_PY_VERSIONS[env_name]
        pkgs_str = " ".join(["'%s'" % x for x in sorted(env_packages)])
        if "mamba" not in cur_ps:
            subprocess.check_call("{conda_bin} install -y {env_str} {channels} "
                                  "{py_version} mamba".format(**locals()), shell=True)
        mamba_bin = os.path.join(os.path.dirname(conda_bin), "mamba")
        pkgs_str = " ".join(["'%s'" % x for x in sorted(env_packages)])
        # Skip in favor of _install_env_pkgs, can be removed later after testing
        if False:
            try:
                subprocess.check_call("{mamba_bin} install -y {env_str} {channels} "
                                      "{py_version} {pkgs_str}".format(**locals()), shell=True)
            except subprocess.CalledProcessError:
                # Fall back to standard conda install when we have system specific issues
                # https://github.com/bcbio/bcbio-nextgen/issues/2871
                pass

def _link_bin(package, system_installdir, conda_info, conda_bin, conda_pkg_list, files=None,
              prefix="", conda_env=None, conda_envdir=None):
    """Link files installed in the bin directory into the install directory.

    This is imperfect but we're trying not to require injecting everything in the anaconda
    directory into a user's path.
    """
    package = package.split("=")[0].split(">")[0]
    final_bindir = os.path.join(system_installdir, "bin")
    if conda_envdir:
        base_bindir = os.path.join(conda_envdir, "bin")
    else:
        base_bindir = os.path.dirname(conda_bin)
    # resolve any symlinks in the final and base heirarchies
    final_bindir = subprocess.check_output("cd %s && pwd -P" % final_bindir, shell=True).decode().strip()
    base_bindir = subprocess.check_output("cd %s && pwd -P" % base_bindir, shell=True).decode().strip()
    for pkg_subdir in [x for x in conda_pkg_list if x["name"] == package]:
        pkg_subdir = pkg_subdir["dist_name"].split("::")[-1]
        for pkg_dir in conda_info["pkgs_dirs"]:
            pkg_bindir = os.path.join(os.path.realpath(pkg_dir), pkg_subdir, "bin")
            python_bindir = os.path.join(os.path.dirname(pkg_bindir), "python-scripts")
            if (os.path.commonprefix([pkg_bindir, base_bindir]).find("anaconda") > 0 and
                    (os.path.exists(python_bindir) or os.path.exists(pkg_bindir))):
                if not files:
                    if not os.path.exists(python_bindir):
                        python_bindir = ""
                    if not os.path.exists(pkg_bindir):
                        pkg_bindir = ""
                    files = subprocess.check_output("ls -1 {pkg_bindir} {python_bindir}"
                                                    .format(**locals()), shell=True).decode().split()
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
    envs = set([])
    for p in packages:
        parts = p.split(";")
        name = parts[0]
        metadata = parts[1:]
        condaenv = None
        for k, v in [x.split("=") for x in metadata]:
            if k == "env":
                condaenv = v
        envs.add(condaenv)
        out[condaenv].append(name)
    envs = [None] + sorted([x for x in list(envs) if x])
    return [(e, out[e]) for e in envs]

def _get_conda_envs(conda_bin):
    info = json.loads(subprocess.check_output("{conda_bin} info --envs --json -q".format(**locals()), shell=True))
    return [e for e in info["envs"] if e.startswith(info["conda_prefix"])]

def _create_environments(conda_bin, packages):
    """Creates custom local environments that conflict with global dependencies.

    Available environments:

    - python2 -- tools that require python 2 and are not compatible with python3.
      The goal is to move all other installs into a default python 3 base environment.
    - python3 -- support tools that require python 3. This will get deprecated but for
      and removed as we move to an all python3 install. For now it collects tools that
      require 3 or some other specific requirements.
    - samtools0 -- For tools that require older samtools 0.1.19
    - dv -- DeepVariant, which requires a specific version of numpy and tensorflow
    - r36 -- R3.6 for PureCN
    - htslib1.10 -- htslib 1.10 for mosdepth and other packages that require it until bioconda
      switches off of 1.9.
    - htslib1.11 -- htslib 1.11 for scramble until bioconda switched off of 1.9
    - bwakit -- requires a specific (old) version of samblaster
    """
    env_names = set([e for e, ps in _split_by_condaenv(packages) if e])
    out = {}
    conda_envs = _get_conda_envs(conda_bin)
    for addenv in ["python3.6", "samtools0", "dv", "python2", "r35", 
                   "htslib1.9", "htslib1.10", "htslib1.11", "htslib1.12", "htslib1.12_py3.9", 
                   "bwakit", "java", "rbcbiornaseq"]:
        if addenv in env_names:
            if not any(x.endswith("/%s" % addenv) for x in conda_envs):
                print("Creating conda environment: %s" % addenv)
                py_version = ENV_PY_VERSIONS[addenv]
                subprocess.check_call("{conda_bin} create -q --no-default-packages -y --name {addenv} {py_version} nomkl"
                                      .format(**locals()), shell=True)
                conda_envs = _get_conda_envs(conda_bin)
            out[addenv] = [x for x in conda_envs if x.endswith("/%s" % addenv)][0]
    return out

def _clean_environment(env_dir):
    """Remove problem elements in environmental directories.

    - Get rid of old history comment lines that cause parsing failures:
      https://github.com/bcbio/bcbio-nextgen/issues/2431
    """
    history_file = os.path.join(env_dir, "conda-meta", "history")
    if os.path.exists(history_file):
        has_problem = False
        cleaned_lines = []
        with open(history_file) as in_handle:
            for line in in_handle:
                # Remove lines like `# create specs:` which have no information after colon
                if line.startswith("#") and len([x for x in line.strip().split(":") if x]) == 1:
                    has_problem = True
                else:
                    cleaned_lines.append(line)
        if has_problem:
            shutil.copy(history_file, history_file + ".orig")
            with open(history_file, "w") as out_handle:
                for line in cleaned_lines:
                    out_handle.write(line)
