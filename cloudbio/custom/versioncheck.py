"""Tool specific version checking to identify out of date dependencies.

This provides infrastructure to check version strings against installed
tools, enabling re-installation if a version doesn't match. This is a
lightweight way to avoid out of date dependencies.
"""
from distutils.version import LooseVersion

from cloudbio.custom import shared
from cloudbio.fabutils import quiet

def _parse_from_stdoutflag(out, flag):
    """Extract version information from a flag in verbose stdout.
    """
    for line in out.split("\n") + out.stderr.split("\n"):
        if line.find(flag) >= 0:
            parts = [x for x in line.split() if not x.startswith(flag)]
            return parts[0]
    raise IOError("Did not find version information with flag %s from: \n %s"
                  % (flag, out))

def up_to_date(env, cmd, version, args=None, stdout_flag=None):
    """Check if the given command is up to date with the provided version.
    """
    if shared._executable_not_on_path(cmd):
        return False
    if args:
        cmd = cmd + " " + " ".join(args)
    with quiet():
        path_safe = "export PATH=$PATH:%s/bin && "
        out = env.safe_run_output(path_safe + cmd)
    if stdout_flag:
        iversion = _parse_from_stdoutflag(out, stdout_flag)
    else:
        iversion = out.strip()
    return LooseVersion(iversion) >= LooseVersion(version)
