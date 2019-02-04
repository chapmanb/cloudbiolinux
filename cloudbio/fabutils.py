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
import os
import re
import shutil

from fabric.api import env, run, sudo, local, settings, hide, put
from fabric.contrib.files import exists, sed, contains, append, comment

import six


SUDO_ENV_KEEPS = []  # Environment variables passed through to sudo environment when using local sudo.
SUDO_ENV_KEEPS += ["http_proxy", "https_proxy"]  # Required for local sudo to work behind a proxy.


# ## Local non-ssh access
def local_exists(path, use_sudo=False):
    func = env.safe_sudo if use_sudo else env.safe_run
    cmd = 'test -e "$(echo %s)"' % path
    cmd_symbolic = 'test -h "$(echo %s)"' % path
    with settings(hide('everything'), warn_only=True):
        env.lcwd = env.cwd
        # We do not use cmd_symbolic so we avoid rescuing broken symlinks
        return not func(cmd).failed

def run_local(use_sudo=False, capture=False):
    def _run(command, *args, **kwags):
        if use_sudo:
            sudo_env = " ".join(["%s=$%s" % (keep, keep) for keep in SUDO_ENV_KEEPS])
            sudo_to = ""
            if "user" in kwags:
                sudo_to = "su - %s" % kwags["user"]
            sudo_prefix = "sudo %s %s bash -c " % (sudo_env, sudo_to)

            command = sudo_prefix + '"%s"' % command.replace('"', '\\"')
        env.lcwd = env.cwd
        return local(command, capture=capture)
    return _run

def local_put(orig_file, new_file):
    shutil.copyfile(orig_file, new_file)

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

def local_comment(filename, regex, use_sudo=False, char='#', backup='.bak', shell=False):
    carot, dollar = '', ''
    if regex.startswith('^'):
        carot = '^'
        regex = regex[1:]
    if regex.endswith('$'):
        dollar = '$'
        regex = regex[:-1]
    regex = "%s(%s)%s" % (carot, regex, dollar)
    return local_sed(
        filename,
        before=regex,
        after=r'%s\1' % char,
        use_sudo=use_sudo,
        backup=backup,
        shell=shell
    )

def _escape_for_regex(text):
    """Escape ``text`` to allow literal matching using egrep"""
    regex = re.escape(text)
    # Seems like double escaping is needed for \
    regex = regex.replace('\\\\', '\\\\\\')
    # Triple-escaping seems to be required for $ signs
    regex = regex.replace(r'\$', r'\\\$')
    # Whereas single quotes should not be escaped
    regex = regex.replace(r"\'", "'")
    return regex

def _expand_path(path):
    return '"$(echo %s)"' % path

def local_contains(filename, text, exact=False, use_sudo=False, escape=True,
    shell=False):
    func = use_sudo and env.safe_sudo or env.safe_run
    if escape:
        text = _escape_for_regex(text)
        if exact:
            text = "^%s$" % text
    with settings(hide('everything'), warn_only=True):
        egrep_cmd = 'egrep "%s" %s' % (text, _expand_path(filename))
        return func(egrep_cmd, shell=shell).succeeded

def local_append(filename, text, use_sudo=False, partial=False, escape=True, shell=False):
    func = use_sudo and env.safe_sudo or env.safe_run
    # Normalize non-list input to be a list
    if isinstance(text, six.string_types):
        text = [text]
    for line in text:
        regex = '^' + _escape_for_regex(line)  + ('' if partial else '$')
        if (env.safe_exists(filename, use_sudo=use_sudo) and line
            and env.safe_contains(filename, regex, use_sudo=use_sudo, escape=False,
                                  shell=shell)):
            continue
        line = line.replace("'", r"'\\''") if escape else line
        func("echo '%s' >> %s" % (line, _expand_path(filename)))

def run_output(*args, **kwargs):
    if not 'shell' in kwargs:
        kwargs['shell'] = False
    return run(*args, **kwargs)

def configure_runsudo(env):
    """Setup env variable with safe_sudo and safe_run,
    supporting non-privileged users and local execution.
    """
    env.is_local = env.hosts == ["localhost"]
    if env.is_local:
        env.safe_put = local_put
        env.safe_sed = local_sed
        env.safe_comment = local_comment
        env.safe_contains = local_contains
        env.safe_append = local_append
        env.safe_exists = local_exists
        env.safe_run = run_local()
        env.safe_run_output = run_local(capture=True)
    else:
        env.safe_put = put
        env.safe_sed = sed
        env.safe_comment = comment
        env.safe_contains = contains
        env.safe_append = append
        env.safe_exists = exists
        env.safe_run = run
        env.safe_run_output = run_output
    if isinstance(getattr(env, "use_sudo", "true"), six.string_types):
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

def find_cmd(env, cmd, args):
    """Retrieve location of a command, checking in installation directory.
    """
    local_cmd = os.path.join(env.system_install, "bin", cmd)
    for cmd in [local_cmd, cmd]:
        with quiet():
            test_version = env.safe_run("%s %s" % (cmd, args))
        if test_version.succeeded:
            return cmd
    return None

try:
    from fabric.api import quiet
except ImportError:
    def quiet():
        return settings(hide('warnings', 'running', 'stdout', 'stderr'), warn_only=True)

try:
    from fabric.api import warn_only
except ImportError:
    def warn_only():
        return settings(warn_only=True)
