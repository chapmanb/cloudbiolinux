"""Process GGD (Get Genomics Data) configurations for installation in biodata directories.

Builds off work done by Aaron Quinlan to define and install genomic data:

https://github.com/arq5x/ggd
"""
import collections
import contextlib
from distutils.version import LooseVersion
import os
import shutil
import subprocess

from fabric.api import env
import yaml

def install_recipe(base_dir, recipe_file):
    """Install data in a biodata directory given instructions from GGD YAML recipe.
    """
    assert env.hosts == ["localhost"], "GGD recipes only work for local runs"
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    recipe = _read_recipe(recipe_file)
    if not version_uptodate(base_dir, recipe):
        if _has_required_programs(recipe["recipe"]["full"].get("required", [])):
            with tx_tmpdir(base_dir) as tmpdir:
                with chdir(tmpdir):
                    print("Running GGD recipe: %s" % recipe["attributes"]["name"])
                    _run_recipe(tmpdir, recipe["recipe"]["full"]["recipe_cmds"],
                                recipe["recipe"]["full"]["recipe_type"])
                _move_files(tmpdir, base_dir, recipe["recipe"]["full"]["recipe_outfiles"])
            add_version(base_dir, recipe)

def _has_required_programs(programs):
    """Ensure the provided programs exist somewhere in the current PATH.

    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    for p in programs:
        found = False
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, p)
            if is_exe(exe_file):
                found = True
                break
        if not found:
            return False
    return True

def _run_recipe(work_dir, recipe_cmds, recipe_type):
    """Create a bash script and run the recipe to download data.
    """
    assert recipe_type == "bash", "Can only currently run bash recipes"
    run_file = os.path.join(work_dir, "ggd-run.sh")
    with open(run_file, "w") as out_handle:
        out_handle.write("#!/bin/bash\nset -eu -o pipefail\n")
        out_handle.write("\n".join(recipe_cmds))
    subprocess.check_output(["bash", run_file])

def _move_files(tmp_dir, final_dir, out_files):
    for out_file in out_files:
        orig = os.path.join(tmp_dir, out_file)
        final = os.path.join(final_dir, out_file)
        assert os.path.exists(orig), ("Did not find expected output file %s in %s" %
                                      (out_file, tmp_dir))
        cur_dir = os.path.dirname(final)
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)
        os.rename(orig, final)

def _read_recipe(in_file):
    in_file = os.path.abspath(os.path.expanduser(in_file))
    with open(in_file) as in_handle:
        recipe = yaml.safe_load(in_handle)
    return recipe

# ## Versioning

def version_uptodate(base_dir, recipe):
    """Check if we have an up to date GGD installation in this directory.
    """
    versions = _get_versions(base_dir)
    return (recipe["attributes"]["name"] in versions and
            LooseVersion(versions[recipe["attributes"]["name"]]) >=
            LooseVersion(str(recipe["attributes"]["version"])))

def add_version(base_dir, recipe):
    versions = _get_versions(base_dir)
    versions[recipe["attributes"]["name"]] = recipe["attributes"]["version"]
    with open(_get_version_file(base_dir), "w") as out_handle:
        for n, v in versions.items():
            out_handle.write("%s,%s\n" % (n, v))

def _get_versions(base_dir):
    version_file = _get_version_file(base_dir)
    versions = collections.OrderedDict()
    if os.path.exists(version_file):
        with open(version_file) as in_handle:
            for line in in_handle:
                name, version = line.strip().split(",")
                versions[name] = version
    return versions

def _get_version_file(base_dir):
    return os.path.join(base_dir, "versions.csv")

# ## Transactional utilities

@contextlib.contextmanager
def tx_tmpdir(base_dir):
    """Context manager to create and remove a transactional temporary directory.
    """
    tmp_dir = os.path.join(base_dir, "txtmp")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    yield tmp_dir
    shutil.rmtree(tmp_dir, ignore_errors=True)

@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)
