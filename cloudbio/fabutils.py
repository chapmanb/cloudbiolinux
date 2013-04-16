"""Utilities to generalize usage of fabric for local and remote builds.

Handles:
  - Providing a local equivalent of standard functions that avoid
    the need to ssh to a local machine.
    Adds generalized targets to the global `env` object which cleanly
    handle local and remote execution:
    - safe_run: Run a command
    - safe_run_output: Run a command, capturing the output
    - safe_sudo: Run a command as sudo user
    - safe_exists: Check for existence of a file.
    - safe_sed: Run sed command.
"""
import hashlib

from fabric.api import env, run, sudo, local, settings, hide
from fabric.contrib.files import exists

# ## Local non-ssh access

def local_exists(path):
    cmd = 'test -e "$(echo %s)"' % path
    with settings(hide('everything'), warn_only=True):
        env.lcwd = env.cwd
        return not local(cmd).failed

def run_local(use_sudo=False, capture=False):
    def _run(command, *args, **kwags):
        if use_sudo:
            command = "sudo " + command
        env.lcwd = env.cwd
        return local(command, capture=capture)
    return _run

def local_sed(filename, before, after, limit='', use_sudo=False, backup='.bak',
              flags='', shell=False):
    """ Run a search-and-replace on ``filename`` with given regex patterns.

    From main fabric contrib, modified to handle local.
    """
    func = env.safe_sudo if use_sudo else env.safe_run
    # Characters to be escaped in both
    for char in "/'":
        before = before.replace(char, r'\%s' % char)
        after = after.replace(char, r'\%s' % char)
    # Characters to be escaped in replacement only (they're useful in regexen
    # in the 'before' part)
    for char in "()":
        after = after.replace(char, r'\%s' % char)
    if limit:
        limit = r'/%s/ ' % limit
    context = {
        'script': r"'%ss/%s/%s/%sg'" % (limit, before, after, flags),
        'filename': '"$(echo %s)"' % filename,
        'backup': backup
    }
    # Test the OS because of differences between sed versions

    with hide('running', 'stdout'):
        platform = env.safe_run("uname")
    if platform in ('NetBSD', 'OpenBSD', 'QNX'):
        # Attempt to protect against failures/collisions
        hasher = hashlib.sha1()
        hasher.update(env.host_string)
        hasher.update(filename)
        context['tmp'] = "/tmp/%s" % hasher.hexdigest()
        # Use temp file to work around lack of -i
        expr = r"""cp -p %(filename)s %(tmp)s \
&& sed -r -e %(script)s %(filename)s > %(tmp)s \
&& cp -p %(filename)s %(filename)s%(backup)s \
&& mv %(tmp)s %(filename)s"""
    else:
        context['extended_regex'] = '-E' if platform == 'Darwin' else '-r'
        expr = r"sed -i%(backup)s %(extended_regex)s -e %(script)s %(filename)s"
    command = expr % context
    return func(command, shell=shell)

def configure_runsudo(env):
    """Setup env variable with safe_sudo and safe_run,
    supporting non-privileged users and local execution.
    """
    env.is_local = env.hosts == ["localhost"]
    env.safe_sed = local_sed
    if env.is_local:
        env.safe_exists = local_exists
        env.safe_run = run_local()
        env.safe_run_output = run_local(capture=True)
    else:
        env.safe_exists = exists
        env.safe_run = run
        env.safe_run_output = run
    if getattr(env, "use_sudo", "true").lower() in ["true", "yes"]:
        env.use_sudo = True
        if env.is_local:
            env.safe_sudo = run_local(True)
        else:
            env.safe_sudo = sudo
    else:
        env.use_sudo = False
        if env.is_local:
            env.safe_sudo = run_local()
        else:
            env.safe_sudo = run
