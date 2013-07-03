"""Tool specific version checking to identify out of date dependencies.

This provides infrastructure to check version strings against installed
tools, enabling re-installation if a version doesn't match. This is a
lightweight way to avoid out of date dependencies.
"""
from distutils.version import LooseVersion

from fabric.api import quiet

from cloudbio.custom import shared

def _parse_from_stdoutflag(stdout, flag):
    """Extract version information from a flag in verbose stdout.
    """
    for line in stdout.split("\n"):
        if line.find(flag) >= 0:
            parts = [x for x in line.split() if not x.startswith(flag)]
            return parts[0]
    return ""

def up_to_date(env, cmd, version, args=None, stdout_flag=None):
    """Check if the given command is up to date with the provided version.
    """
    if shared._executable_not_on_path(cmd):
        return False
    if args:
        cmd = cmd + " " + " ".join(args)
    with quiet():
        out = env.safe_run_output(cmd)
    if stdout_flag:
        iversion = _parse_from_stdoutflag(out, stdout_flag)
    else:
        iversion = out.strip()
    return LooseVersion(iversion) >= LooseVersion(version)
