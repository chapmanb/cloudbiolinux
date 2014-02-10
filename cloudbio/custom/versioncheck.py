"""Tool specific version checking to identify out of date dependencies.

This provides infrastructure to check version strings against installed
tools, enabling re-installation if a version doesn't match. This is a
lightweight way to avoid out of date dependencies.
"""
from distutils.version import LooseVersion

from cloudbio.custom import shared
from cloudbio.fabutils import quiet

def _parse_from_stdoutflag(out, flag, stdout_index=-1):
    """Extract version information from a flag in verbose stdout.

    flag -- text information to identify the line we should split for a version
    stdout_index -- Position of the version information in the split line. Defaults
    to the last item.
    """
    for line in out.split("\n") + out.stderr.split("\n"):
        if line.find(flag) >= 0:
            parts = line.split()
            return parts[stdout_index].strip()
    raise IOError("Did not find version information with flag %s from: \n %s"
                  % (flag, out))

def _clean_version(x):
    if x.startswith("upstream/"):
        x = x.replace("upstream/", "")
    if x.startswith("("):
        x = x[1:].strip()
    if x.endswith(")"):
        x = x[:-1].strip()
    if x.startswith("v"):
        x = x[1:].strip()
    return x

def up_to_date(env, cmd, version, args=None, stdout_flag=None,
               stdout_index=-1):
    iversion = get_installed_version(env, cmd, version, args, stdout_flag,
                                     stdout_index)
    if not iversion:
        return False
    else:
        return LooseVersion(iversion) >= LooseVersion(version)

def is_version(env, cmd, version, args=None, stdout_flag=None,
               stdout_index=-1):
    iversion = get_installed_version(env, cmd, version, args, stdout_flag,
                                     stdout_index)
    if not iversion:
        return False
    else:
        return LooseVersion(iversion) == LooseVersion(version)

def get_installed_version(env, cmd, version, args=None, stdout_flag=None,
                          stdout_index=-1):
    """Check if the given command is up to date with the provided version.
    """
    if shared._executable_not_on_path(cmd):
        return False
    if args:
        cmd = cmd + " " + " ".join(args)
    with quiet():
        path_safe = ("export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:{s}/lib/pkgconfig && "
                     "export PATH=$PATH:{s}/bin && "
                     "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{s}/lib && ".format(s=env.system_install))
        out = env.safe_run_output(path_safe + cmd)
    if stdout_flag:
        iversion = _parse_from_stdoutflag(out, stdout_flag, stdout_index)
    else:
        iversion = out.strip()
    iversion = _clean_version(iversion)
    if " not found in the pkg-config search path" in iversion:
        return False
    return iversion
